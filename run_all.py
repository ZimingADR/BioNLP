# -*- coding: utf-8 -*-
"""
BioNLP workspace 一键复现脚本（阶段一至阶段四 + 结果包导出）
================================================================
默认复用 data/ 下已缓存的 PubTator 语料，不重新下载（下载 3,306 篇需约 5.3 小时）。
若需要从零下载，先运行 scripts/01_fetch_pubmed_and_pubtator.py。

执行顺序:
  02 词汇统计/TTR/TF-IDF
  03 Oryzabase + RTO 词典
  04 句子级实体与共现证据
  06 阶段一/二图表
  07 SDP 因果抽取（水稻-only）
  08 SDP 审计评测（默认只评测；审计样本与标签已冻结）
  09 TO/RTO 富集
  10 RiceMind 交叉验证 + 网络 + 注释缺口候选
  11 阶段三/四图表
  12 导出表格/key_numbers/数据清单

用法:
  python run_all.py
  python run_all.py --skip-figures
  python run_all.py --with-audit-sample   # 会重新抽 100 条并覆盖旧审计表，需重新标注
"""

import argparse
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
SCRIPTS = BASE_DIR / "scripts"
RESULTS = BASE_DIR / "results"
LOG = RESULTS / "run_all.log"


def run_step(name, args, log_handle):
    print(f"\n{'=' * 70}\n[{datetime.now().strftime('%H:%M:%S')}] {name}\n{'=' * 70}")
    log_handle.write(f"\n{'=' * 70}\n[{datetime.now().isoformat(timespec='seconds')}] {name}\n{'=' * 70}\n")
    log_handle.flush()
    t0 = time.time()
    proc = subprocess.run([sys.executable, str(SCRIPTS / args[0]), *args[1:]],
                          cwd=str(BASE_DIR), text=True, encoding="utf-8", errors="replace")
    dt = time.time() - t0
    log_handle.write(f"[exit={proc.returncode}] elapsed={dt:.1f}s\n")
    log_handle.flush()
    if proc.returncode != 0:
        raise SystemExit(f"步骤失败: {name} (exit {proc.returncode})")
    print(f"[done] {name}  ({dt:.1f}s)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-figures", action="store_true")
    parser.add_argument("--with-audit-sample", action="store_true")
    args = parser.parse_args()

    RESULTS.mkdir(exist_ok=True)
    if not (BASE_DIR / "data" / "result_Pubtator_all.txt").exists():
        raise SystemExit("缺少 data/result_Pubtator_all.txt；请先运行 scripts/01_fetch_pubmed_and_pubtator.py")
    if not (BASE_DIR / "dicts" / "OryzabaseGeneListEn.txt").exists() or not (BASE_DIR / "dicts" / "RTO-1.0.obo").exists():
        raise SystemExit("缺少 dicts/ 下的 Oryzabase 或 RTO 文件")

    steps = [
        ("02_lexical_stats_ttr_tfidf.py", []),
        ("03_build_oryza_and_rto_dicts.py", []),
        ("04_extract_sentence_evidence.py", []),
    ]
    if not args.skip_figures:
        steps.append(("06_generate_bibm_nature_figures.py", []))
    steps.append(("07_sdp_relation_extraction.py", ["--rice-only"]))
    if args.with_audit_sample:
        steps.append(("08_sdp_audit_eval.py", ["--sample"]))
    steps.append(("08_sdp_audit_eval.py", ["--evaluate"]))
    steps.append(("09_enrichment_validation.py", []))
    steps.append(("10_ricemind_novelty.py", []))
    if not args.skip_figures:
        steps.append(("11_generate_part34_figures.py", []))
    steps.append(("12_export_results_package.py", []))

    print(f"共 {len(steps)} 步；日志: {LOG}")
    with open(LOG, "a", encoding="utf-8") as log:
        log.write(f"\n\n########## RUN {datetime.now().isoformat(timespec='seconds')} ##########\n")
        for name, argv in steps:
            run_step(name, [name, *argv], log)
    print("\n全部完成。结果见 results/，表格见 results/tables/。")


if __name__ == "__main__":
    main()