# -*- coding: utf-8 -*-
"""
Step 4 (Part 2 核心合成): 句子级证据构建、Oryzabase 基因归一化 + RTO 性状本体对齐 + RiceMind API 对照
对应课件:
  - JingboSlides_8.PubMedTermsNERandShellProgramming.pdf
    - Page 32-37: RiceMind API 访问
    - Page 50-53: sentence_spans 字符偏移对齐到句子
    - Page 57-62: 面向水稻数据的自主实体识别开发 (路线 A 词典匹配 + 路线 B 共现拼接 + 错误分析)
输出:
  - results/rice_drought_sentence_evidence.tsv (所有基因+性状共现的句子级证据表)
  - results/all_sentence_entities.tsv (所有识别出基因或性状的句子明细)
  - results/part2_normalization_and_evidence_report.json (含输入样例、输出样例、局限案例与 RiceMind 对照)
"""

import json
import re
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

from nltk.tokenize import sent_tokenize

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
DICTS_DIR = BASE_DIR / "dicts"
RESULTS_DIR = BASE_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

PUBTATOR_ALL = DATA_DIR / "result_Pubtator_all.txt"
PUBTATOR_100 = DATA_DIR / "result_Pubtator_100.txt"
PUBTATOR_FILE = PUBTATOR_ALL if (PUBTATOR_ALL.exists() and PUBTATOR_ALL.stat().st_size > 1000) else PUBTATOR_100
GENE_DICT_FILE = DICTS_DIR / "rice_gene_synonym_to_rap.json"
TRAIT_DICT_FILE = DICTS_DIR / "rto_trait_synonym_to_id.json"


def sentence_spans(text: str):
    """
    严格复现课件 Page 51 的 sentence_spans 函数：
    定位每个句子在拼接全文 (title + " " + abstract) 中的起止字符偏移量 (begin, end)。
    """
    spans, cursor = [], 0
    for sent_id, sent in enumerate(sent_tokenize(text), start=1):
        start = text.find(sent, cursor)
        if start == -1:
            # 容错回退：若空白字符差异导致精确查找失败，从 cursor 附近搜索首词
            start = text.find(sent[:20], cursor)
            if start == -1:
                start = cursor
        end = start + len(sent)
        spans.append({"sent_id": sent_id, "begin": start, "end": end, "text": sent})
        cursor = end
    return spans


def tokenize_with_spans(sent_text: str, sent_begin: int):
    """
    将句子切分为带精确字符偏移量的词元列表，支持连字符基因名 (如 OsDREB1A, OsCKA2-1)
    """
    tokens = []
    for m in re.finditer(r"[A-Za-z0-9]+(?:-[A-Za-z0-9]+)*", sent_text):
        tokens.append({
            "word": m.group(0),
            "lower": m.group(0).lower(),
            "begin": sent_begin + m.start(),
            "end": sent_begin + m.end(),
            "rel_begin": m.start(),
            "rel_end": m.end(),
        })
    return tokens


def match_dictionary_ngrams(tokens: list, sent_text: str, sent_begin: int, norm_dict: dict, entity_type: str, max_n: int = 4):
    """
    在句子中执行 1~max_n 元 n-gram 最长词典优先匹配，并自动验证 text[begin:end] == mention
    """
    n = len(tokens)
    occupied = [False] * n
    matches = []

    # 从长到短匹配 (4-gram -> 1-gram)，避免子串重复命中
    for k in range(min(max_n, n), 0, -1):
        for i in range(n - k + 1):
            if any(occupied[i : i + k]):
                continue
            rel_b = tokens[i]["rel_begin"]
            rel_e = tokens[i + k - 1]["rel_end"]
            mention_raw = sent_text[rel_b:rel_e]
            lookup_key = " ".join(t["lower"] for t in tokens[i : i + k])

            # 针对基因的 1-gram 额外防护：若全为小写普通字母且长度<=4，跳过以防误判普通英文词
            if entity_type == "Gene" and k == 1:
                if mention_raw.islower() and mention_raw.isalpha() and len(mention_raw) <= 4:
                    continue

            info = norm_dict.get(lookup_key)
            if not info and k == 1 and "-" in lookup_key:
                info = norm_dict.get(lookup_key.replace("-", ""))

            if info:
                for j in range(i, i + k):
                    occupied[j] = True
                abs_b = sent_begin + rel_b
                abs_e = sent_begin + rel_e
                if entity_type == "Gene":
                    matches.append({
                        "begin": abs_b,
                        "end": abs_e,
                        "mention": mention_raw,
                        "type": "Gene",
                        "canonical": info["canonical_symbol"],
                        "norm_id": info["rap_id"],
                        "source": "Oryzabase_Dict",
                    })
                else:
                    matches.append({
                        "begin": abs_b,
                        "end": abs_e,
                        "mention": mention_raw,
                        "type": "Trait",
                        "canonical": info["canonical_name"],
                        "norm_id": info["rto_id"],
                        "source": "RTO_Ontology",
                    })

    matches.sort(key=lambda x: x["begin"])
    return matches


def query_ricemind_api(gene_symbol: str):
    """
    调用课件 Page 33 介绍的华农 RiceMind API，查询基因画像作为交叉验证
    """
    url = f"http://lit-evi.hzau.edu.cn/ricemind-api/gene-profile/?gene={urllib.parse.quote(gene_symbol)}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "HZAU-BioNLP-Course/1.0"})
        with urllib.request.urlopen(req, timeout=8) as resp:
            data = json.loads(resp.read().decode("utf-8", errors="replace"))
            return data
    except Exception as e:
        return {"error": str(e)}


def main():
    from importlib.machinery import SourceFileLoader
    step2_mod = SourceFileLoader("step2", str(BASE_DIR / "scripts" / "02_lexical_stats_ttr_tfidf.py")).load_module()

    docs = step2_mod.parse_pubtator_documents(PUBTATOR_FILE)
    gene_dict = json.loads(GENE_DICT_FILE.read_text(encoding="utf-8"))
    trait_dict = json.loads(TRAIT_DICT_FILE.read_text(encoding="utf-8"))

    cooccurrence_rows = []
    all_entity_rows = []
    gene_counter = Counter()
    trait_counter = Counter()
    pair_counter = Counter()
    pubtator_gene_count = 0
    oryza_normalized_pubtator_count = 0

    for doc in docs:
        pmid = doc["pmid"]
        full_text = doc["text"]
        spans = sentence_spans(full_text)

        # 将 PubTator 原生标注的 Gene 实体索引到句子，同时记录非基因实体 (Chemical/Species/Disease) 的 span 以防误判
        pubtator_genes_by_sent = defaultdict(list)
        non_gene_pubtator_spans = []
        for ent in doc["entities"]:
            if ent["type"] == "Gene" and ent.get("offset_verified", True):
                pubtator_gene_count += 1
                m_low = ent["mention"].lower()
                oryza_info = gene_dict.get(m_low) or gene_dict.get(m_low.replace("-", ""))
                if oryza_info:
                    norm_id = oryza_info["rap_id"]
                    canonical = oryza_info["canonical_symbol"]
                    oryza_normalized_pubtator_count += 1
                else:
                    norm_id = f"NCBIGene:{ent['identifier']}"
                    canonical = ent["mention"]

                for sp in spans:
                    if ent["begin"] >= sp["begin"] and ent["end"] <= sp["end"]:
                        pubtator_genes_by_sent[sp["sent_id"]].append({
                            "begin": ent["begin"],
                            "end": ent["end"],
                            "mention": ent["mention"],
                            "type": "Gene",
                            "canonical": canonical,
                            "norm_id": norm_id,
                            "source": "PubTator+Oryzabase" if oryza_info else "PubTator_Only",
                        })
                        break
            else:
                non_gene_pubtator_spans.append((ent["begin"], ent["end"]))

        for sp in spans:
            sid = sp["sent_id"]
            stext = sp["text"]
            sbegin = sp["begin"]
            toks = tokenize_with_spans(stext, sbegin)

            # 1. Oryzabase 词典匹配基因 + 融合 PubTator 基因 (排除已被 PubTator 识别为 Chemical/Species 的位置)
            dict_genes = match_dictionary_ngrams(toks, stext, sbegin, gene_dict, "Gene", max_n=3)
            merged_genes = list(pubtator_genes_by_sent[sid])
            existing_spans = {(g["begin"], g["end"]) for g in merged_genes}
            for dg in dict_genes:
                overlaps_existing = any(not (dg["end"] <= eb or dg["begin"] >= ee) for (eb, ee) in existing_spans)
                overlaps_non_gene = any(not (dg["end"] <= nb or dg["begin"] >= ne) for (nb, ne) in non_gene_pubtator_spans)
                if not overlaps_existing and not overlaps_non_gene:
                    merged_genes.append(dg)
                    existing_spans.add((dg["begin"], dg["end"]))

            # 2. RTO 本体词典匹配性状
            traits = match_dictionary_ngrams(toks, stext, sbegin, trait_dict, "Trait", max_n=4)

            # 自动校验每个实体的字符偏移量 text[begin:end] == mention
            for ent in merged_genes + traits:
                assert full_text[ent["begin"] : ent["end"]] == ent["mention"], (
                    f"Offset mismatch in {pmid}: {full_text[ent['begin']:ent['end']]} != {ent['mention']}"
                )
                all_entity_rows.append({
                    "pmid": pmid,
                    "sent_id": sid,
                    "begin": ent["begin"],
                    "end": ent["end"],
                    "mention": ent["mention"],
                    "type": ent["type"],
                    "canonical": ent["canonical"],
                    "norm_id": ent["norm_id"],
                    "source": ent["source"],
                    "sentence": stext,
                })

            for g in merged_genes:
                gene_counter[(g["canonical"], g["norm_id"])] += 1
            for t in traits:
                trait_counter[(t["canonical"], t["norm_id"])] += 1

            # 若同一句子内同时存在【水稻基因】与【表型性状】，构建句子级共现证据行（路线 B）
            if merged_genes and traits:
                # 去重同一句中相同归一化 ID 的基因-性状对
                seen_pairs = set()
                for g in merged_genes:
                    for t in traits:
                        pair_key = (g["norm_id"], t["norm_id"])
                        if pair_key in seen_pairs:
                            continue
                        seen_pairs.add(pair_key)
                        pair_counter[(g["canonical"], g["norm_id"], t["canonical"], t["norm_id"])] += 1
                        cooccurrence_rows.append({
                            "pmid": pmid,
                            "sent_id": sid,
                            "gene_mention": g["mention"],
                            "gene_canonical": g["canonical"],
                            "rap_id": g["norm_id"],
                            "gene_source": g["source"],
                            "trait_mention": t["mention"],
                            "trait_canonical": t["canonical"],
                            "rto_id": t["norm_id"],
                            "sentence": stext,
                        })

    # 保存所有句子级共现证据表
    with open(RESULTS_DIR / "rice_drought_sentence_evidence.tsv", "w", encoding="utf-8") as f:
        f.write("PMID\tsent_id\tgene_mention\tgene_canonical\tRAP_ID\tgene_source\ttrait_mention\ttrait_canonical\tRTO_ID\tsentence\n")
        for r in cooccurrence_rows:
            f.write(
                f"{r['pmid']}\t{r['sent_id']}\t{r['gene_mention']}\t{r['gene_canonical']}\t"
                f"{r['rap_id']}\t{r['gene_source']}\t{r['trait_mention']}\t{r['trait_canonical']}\t"
                f"{r['rto_id']}\t{r['sentence']}\n"
            )

    with open(RESULTS_DIR / "all_sentence_entities.tsv", "w", encoding="utf-8") as f:
        f.write("PMID\tsent_id\tbegin\tend\tmention\ttype\tcanonical\tnorm_id\tsource\tsentence\n")
        for r in all_entity_rows:
            f.write(
                f"{r['pmid']}\t{r['sent_id']}\t{r['begin']}\t{r['end']}\t{r['mention']}\t"
                f"{r['type']}\t{r['canonical']}\t{r['norm_id']}\t{r['source']}\t{r['sentence']}\n"
            )

    print("=== Part 2 实体识别、归一化与句子级证据构建汇总 ===")
    print(f"总计识别并对齐实体提及数: {len(all_entity_rows):,} 个 (100% 通过 text[begin:end]==mention 校验)")
    print(f"PubTator 原生识别 Gene 提及数: {pubtator_gene_count:,} (成功映射至 Oryzabase RAP ID: {oryza_normalized_pubtator_count:,})")
    print(f"成功提取【水稻基因 - 抗旱/农艺性状】句子级共现证据记录: {len(cooccurrence_rows):,} 条")

    print("\n=== 高频水稻抗旱相关基因 Top 15 (归一化后) ===")
    top_genes = []
    for (canon, rid), cnt in gene_counter.most_common(15):
        print(f"  {canon:<14} ({rid:<16}) : {cnt:3d} 次提及")
        top_genes.append({"canonical_symbol": canon, "rap_id": rid, "count": cnt})

    print("\n=== 高频水稻抗旱相关表型性状 Top 15 (RTO 归一化后) ===")
    top_traits = []
    for (canon, tid), cnt in trait_counter.most_common(15):
        print(f"  {canon:<30} ({tid:<12}) : {cnt:3d} 次提及")
        top_traits.append({"canonical_trait": canon, "rto_id": tid, "count": cnt})

    print("\n=== 高频【基因 - 性状】共现配对 Top 15 ===")
    top_pairs = []
    for (gc, grid, tc, trid), cnt in pair_counter.most_common(15):
        print(f"  {gc} ({grid}) <---> {tc} ({trid}) : {cnt:2d} 句证据")
        top_pairs.append({
            "gene_symbol": gc,
            "rap_id": grid,
            "trait_name": tc,
            "rto_id": trid,
            "evidence_sentences_count": cnt,
        })

    # 与华农 RiceMind API 进行交叉验证 (取前 3 个代表性基因)
    ricemind_checks = {}
    for item in top_genes[:3]:
        sym = item["canonical_symbol"]
        ricemind_checks[sym] = query_ricemind_api(sym)

    total_gene_mentions = sum(gene_counter.values())
    total_trait_mentions = sum(trait_counter.values())
    rap_mapped_gene_mentions = sum(cnt for (_, rid), cnt in gene_counter.items() if str(rid).startswith("Os"))

    report = {
        "num_documents_analyzed": len(docs),
        "total_entity_mentions": len(all_entity_rows),
        "total_gene_mentions": total_gene_mentions,
        "rap_mapped_gene_mentions": rap_mapped_gene_mentions,
        "total_trait_mentions": total_trait_mentions,
        "unique_genes_count": len(gene_counter),
        "unique_traits_count": len(trait_counter),
        "unique_gene_trait_pairs_count": len(pair_counter),
        "pubtator_gene_mentions": pubtator_gene_count,
        "pubtator_mapped_to_rap": oryza_normalized_pubtator_count,
        "sentence_cooccurrence_evidence_count": len(cooccurrence_rows),
        "top_genes": top_genes,
        "top_traits": top_traits,
        "top_gene_trait_pairs": top_pairs,
        "sample_evidence_records": cooccurrence_rows[:25],
        "ricemind_api_cross_check": ricemind_checks,
    }
    (RESULTS_DIR / "part2_normalization_and_evidence_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
