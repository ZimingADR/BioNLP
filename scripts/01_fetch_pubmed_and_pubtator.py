# -*- coding: utf-8 -*-
"""
Step 1: 文献获取与 PubTator 标注批量下载 (支持断点续传 + 批量拼接请求 + 严格遵守课件 Page 31 的 sleep 5.8s)
检索式: ("Oryza sativa"[Title/Abstract] OR "Oryza sativa"[MeSH Terms] OR "rice"[Title/Abstract]) AND "drought"[Title/Abstract]
"""

import argparse
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

QUERY = '("Oryza sativa"[Title/Abstract] OR "Oryza sativa"[MeSH Terms] OR "rice"[Title/Abstract]) AND "drought"[Title/Abstract]'
XML_FILE = DATA_DIR / "pmid_raw.xml"
ALL_PMID_FILE = DATA_DIR / "pmid_list_all.txt"
SUBSET_PMID_FILE = DATA_DIR / "pmid_list_100.txt"
PUBTATOR_100_FILE = DATA_DIR / "result_Pubtator_100.txt"
PUBTATOR_ALL_FILE = DATA_DIR / "result_Pubtator_all.txt"
CACHE_DIR = DATA_DIR / "pubtator_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)


def fetch_pmid_list(limit: int = 9999, reuse_existing: bool = True):
    if reuse_existing and ALL_PMID_FILE.exists() and ALL_PMID_FILE.stat().st_size > 100:
        pmids = [line.strip() for line in ALL_PMID_FILE.read_text(encoding="utf-8").splitlines() if line.strip()]
        total_found = len(pmids)
        subset = pmids[:limit] if limit > 0 else pmids
        print(f"[1/2] 已复用本地检索列表 {ALL_PMID_FILE.name} (共 {total_found:,} 篇，本次目标: {len(subset):,} 篇)")
        return subset, total_found

    params = urllib.parse.urlencode({
        "db": "pubmed",
        "term": QUERY,
        "retmax": 9999,
        "usehistory": "y",
    })
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{params}"
    print(f"[1/2] 正在向 NCBI E-utilities 请求检索式:\n      {QUERY}")
    req = urllib.request.Request(url, headers={"User-Agent": "HZAU-BioNLP-Course/1.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        xml_text = resp.read().decode("utf-8", errors="replace")

    XML_FILE.write_text(xml_text, encoding="utf-8")

    count_match = re.search(r"<Count>(\d+)</Count>", xml_text)
    total_found = int(count_match.group(1)) if count_match else 0
    pmids = re.findall(r"<Id>(\d+)</Id>", xml_text)

    ALL_PMID_FILE.write_text("\n".join(pmids) + "\n", encoding="utf-8")
    subset = pmids[:limit] if limit > 0 else pmids
    if limit == 100:
        SUBSET_PMID_FILE.write_text("\n".join(subset) + "\n", encoding="utf-8")

    print(f"      PubMed 共命中匹配文献总数: {total_found:,} 篇")
    print(f"      已提取前 {len(pmids):,} 个 PMID 存入: {ALL_PMID_FILE.name}")
    return subset, total_found


def split_pubtator_blocks(raw_text: str):
    """将批量返回的 PubTator 文本按文章块切分，并提取每篇的 PMID"""
    records_by_pmid = {}
    blocks = re.split(r"\n\s*\n", raw_text.strip())
    for block in blocks:
        block = block.strip()
        if not block:
            continue
        m = re.match(r"^(\d+)\|t\|", block)
        if m:
            pmid = m.group(1)
            records_by_pmid[pmid] = block
    return records_by_pmid


def fetch_pubtator_batch_with_sleep(pmids, batch_size: int = 80, sleep_sec: float = 5.8, out_file: Path = PUBTATOR_ALL_FILE):
    total = len(pmids)
    # 先扫描本地缓存，筛出真正需要联网请求的 PMID
    uncached_pmids = []
    for pmid in pmids:
        cache_file = CACHE_DIR / f"{pmid}.pubtator.txt"
        unindexed_file = CACHE_DIR / f"{pmid}.unindexed"
        if (cache_file.exists() and cache_file.stat().st_size > 20) or unindexed_file.exists():
            continue
        uncached_pmids.append(pmid)

    num_batches = (len(uncached_pmids) + batch_size - 1) // batch_size
    print(
        f"\n[2/2] Let us curate the results! 目标文献: {total:,} 篇 | 已缓存: {total - len(uncached_pmids):,} 篇 | "
        f"待请求: {len(uncached_pmids):,} 篇 (每批 {batch_size} 篇，共 {num_batches} 批，批次间隔严格 sleep {sleep_sec}s)...",
        flush=True,
    )

    for b_idx in range(num_batches):
        batch = uncached_pmids[b_idx * batch_size : (b_idx + 1) * batch_size]
        pmid_str = ",".join(batch)
        url = f"https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/pubtator?pmids={pmid_str}"

        text = ""
        success = False
        for attempt in range(4):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "HZAU-BioNLP-Course/1.0"})
                with urllib.request.urlopen(req, timeout=45) as resp:
                    text = resp.read().decode("utf-8", errors="replace").strip()
                success = True
                break
            except Exception as e:
                print(f"  [警告] 第 {b_idx+1}/{num_batches} 批请求异常 (尝试 {attempt+1}/4): {e}", file=sys.stderr, flush=True)
                time.sleep(3.0)

        if success and text:
            records_map = split_pubtator_blocks(text)
            for pmid in batch:
                if pmid in records_map:
                    (CACHE_DIR / f"{pmid}.pubtator.txt").write_text(records_map[pmid] + "\n\n", encoding="utf-8")
                else:
                    (CACHE_DIR / f"{pmid}.unindexed").write_text("unindexed\n", encoding="utf-8")
            print(
                f"Batch [{b_idx+1}/{num_batches}] completed: requested {len(batch)} PMIDs -> fetched {len(records_map)} valid PubTator records.",
                flush=True,
            )
        else:
            print(f"Batch [{b_idx+1}/{num_batches}] FAILED after retries.", file=sys.stderr, flush=True)

        if b_idx + 1 < num_batches:
            time.sleep(sleep_sec)

    # 汇总所有已缓存的记录写入最终文件
    combined_records = []
    for pmid in pmids:
        cache_file = CACHE_DIR / f"{pmid}.pubtator.txt"
        if cache_file.exists() and cache_file.stat().st_size > 20:
            combined_records.append(cache_file.read_text(encoding="utf-8", errors="replace").strip())

    out_file.write_text("\n\n".join(combined_records) + "\n\n", encoding="utf-8")
    print(f"\n全部汇总完成! 有效 PubTator 文献: {len(combined_records):,}/{total:,} 篇 -> 已写入: {out_file}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PubMed & PubTator 批量获取脚本 (遵守 5.8s 限流)")
    parser.add_argument("--limit", type=int, default=9999, help="下载文献数量上限 (默认 9999 即全量 9,586 篇)")
    parser.add_argument("--batch-size", type=int, default=80, help="单次请求拼接的 PMID 数量 (默认 80，设为 1 即单篇请求)")
    parser.add_argument("--sleep", type=float, default=5.8, help="每次请求之间的强制休眠秒数 (默认 5.8s)")
    args = parser.parse_args()

    subset_pmids, _ = fetch_pmid_list(limit=args.limit, reuse_existing=True)
    target_out = PUBTATOR_100_FILE if args.limit == 100 else PUBTATOR_ALL_FILE
    fetch_pubtator_batch_with_sleep(subset_pmids, batch_size=args.batch_size, sleep_sec=args.sleep, out_file=target_out)
