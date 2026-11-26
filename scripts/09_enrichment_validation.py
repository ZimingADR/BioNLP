# -*- coding: utf-8 -*-
"""
Step 9 (阶段四): 基于 SDP 高置信基因集的 TO/RTO 超几何富集分析
==========================================================================
输入:
  - results/sdp_relations.tsv                       阶段三 SDP 关系
  - lab03/task3_to_enrichment/results/RiceMind_TO_RTO_background.gmt
      （已由 lab03 任务③构建并验证：N=2,218 背景基因，272 个 TO/RTO 条目）

方法（严格复现白板“坛子抽球”四步）:
  Step 1  固定一个 TO/RTO 条目 T_j
  Step 2  取 SDP 因果调控基因集作为待观测小规模 GeneSet (k)
  Step 3  在背景基因库 (N) 中统计与 T_j 关联的黑球数 m
  Step 4  x = GeneSet 中命中 T_j 的基因数; p = hypergeom.sf(x-1, N, m, k)
          并对全部测试条目做 Benjamini-Hochberg FDR

输出:
  - results/part4_enrichment_all_terms.csv
  - results/part4_enrichment_significant.csv
  - results/part4_enrichment_summary.json
"""

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import hypergeom

BASE_DIR = Path(__file__).resolve().parent.parent
LAB03_DIR = BASE_DIR.parent / "lab03"
RESULTS_DIR = BASE_DIR / "results"
SDP_FILE = RESULTS_DIR / "sdp_relations.tsv"
GMT_FILE = LAB03_DIR / "task3_to_enrichment" / "results" / "RiceMind_TO_RTO_background.gmt"

CAUSAL = {"positive_regulation", "negative_regulation", "regulation"}


def benjamini_hochberg(p_values):
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1)
    adjusted = p * n / ranks
    adjusted = np.minimum.accumulate(adjusted[order][::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adjusted, 0.0, 1.0)
    return out


def load_gmt(path):
    term_to_genes = {}
    term_meta = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        term = parts[0].strip()
        source = parts[1].strip()
        genes = {g.strip() for g in parts[2:] if g.strip()}
        term_to_genes[term] = genes
        term_meta[term] = {"source": source}
    return term_to_genes, term_meta


def load_sdp_gene_stats():
    rows = list(csv.DictReader(open(SDP_FILE, encoding="utf-8"), delimiter="\t"))
    causal = [r for r in rows if r["relation"] in CAUSAL]
    gene_stats = defaultdict(lambda: {"relations": 0, "pairs": set(), "traits": Counter(),
                                      "positive": 0, "negative": 0, "neutral": 0,
                                      "high_conf": 0, "pmids": set()})
    for r in causal:
        g = r["RAP_ID"]
        s = gene_stats[g]
        s["relations"] += 1
        s["pairs"].add((r["RAP_ID"], r["RTO_ID"]))
        s["traits"][r["trait_canonical"]] += 1
        s["pmids"].add(r["PMID"])
        if r["relation"] == "positive_regulation":
            s["positive"] += 1
        elif r["relation"] == "negative_regulation":
            s["negative"] += 1
        else:
            s["neutral"] += 1
        if r["confidence"] == "high":
            s["high_conf"] += 1
    return rows, causal, gene_stats


def run_ora(query_genes, background_genes, term_to_genes):
    query = {g for g in query_genes if g in background_genes}
    N = len(background_genes)
    k = len(query)
    records = []
    for term, genes in term_to_genes.items():
        term_in_bg = genes & background_genes
        m = len(term_in_bg)
        if m == 0:
            continue
        x = len(query & term_in_bg)
        if x == 0:
            p = 1.0
        else:
            p = float(hypergeom.sf(x - 1, N, m, k))
        records.append({
            "Term": term,
            "N_background": N,
            "m_term_in_background": m,
            "k_query_mapped": k,
            "x_overlap": x,
            "expected": (m * k / N) if N else 0.0,
            "fold_enrichment": (x / (m * k / N)) if m and k else 0.0,
            "P-value": p,
            "Genes": ";".join(sorted(query & term_in_bg)),
        })
    if not records:
        return []
    pvals = [r["P-value"] for r in records]
    qvals = benjamini_hochberg(pvals)
    for r, q in zip(records, qvals):
        r["Adjusted P-value"] = float(q)
    records.sort(key=lambda r: (r["P-value"], -r["x_overlap"]))
    return records


def main():
    print("[1/4] 读取 SDP 关系与 TO/RTO 背景 GMT ...")
    rows, causal, gene_stats = load_sdp_gene_stats()
    term_to_genes, term_meta = load_gmt(GMT_FILE)
    background = set().union(*term_to_genes.values())
    print(f"    SDP 因果关系 {len(causal):,} 条；基因 {len(gene_stats):,} 个")
    print(f"    TO/RTO 条目 {len(term_to_genes):,} 个；背景基因 {len(background):,} 个")

    # 排序：关系数优先，高置信数次之
    ranked = sorted(
        gene_stats.items(),
        key=lambda kv: (len(kv[1]["pairs"]), kv[1]["high_conf"], kv[1]["relations"]),
        reverse=True,
    )
    top100 = [g for g, s in ranked if g in background][:100]
    top50 = top100[:50]
    top20 = top100[:20]
    core14 = ["Os09g0522200", "Os07g0154100", "Os02g0766700", "Os03g0815100",
              "Os01g0884300", "Os06g0211200", "Os09g0522000", "Os07g0448800",
              "Os05g0542500", "Os05g0455500", "Os01g0165000", "Os01g0859300",
              "Os03g0786400", "Os09g0439800"]

    print("[2/4] 对 top100 / top50 / top20 / core14 分别执行超几何富集 ...")
    sets = {"top100": top100, "top50": top50, "top20": top20, "core14": core14}
    all_results = {}
    for name, genes in sets.items():
        records = run_ora(genes, background, term_to_genes)
        all_results[name] = records
        sig = [r for r in records if r["Adjusted P-value"] < 0.05 and r["x_overlap"] > 0]
        print(f"    {name}: query={len(genes)}, mapped={len(set(genes) & background)}, "
              f"hit={len([r for r in records if r['x_overlap'] > 0])}, sig(q<0.05)={len(sig)}")

    print("[3/4] 写出 CSV 与汇总 ...")
    out_fields = ["Term", "N_background", "m_term_in_background", "k_query_mapped",
                  "x_overlap", "expected", "fold_enrichment", "P-value",
                  "Adjusted P-value", "Genes"]
    all_csv = RESULTS_DIR / "part4_enrichment_all_terms.csv"
    with open(all_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene_set"] + out_fields)
        w.writeheader()
        for name, records in all_results.items():
            for r in records:
                w.writerow({"gene_set": name, **r})

    sig_csv = RESULTS_DIR / "part4_enrichment_significant.csv"
    with open(sig_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene_set"] + out_fields)
        w.writeheader()
        for name, records in all_results.items():
            for r in records:
                if r["Adjusted P-value"] < 0.05 and r["x_overlap"] > 0:
                    w.writerow({"gene_set": name, **r})

    summary = {
        "sdp_causal_relations": len(causal),
        "sdp_unique_genes": len(gene_stats),
        "background_genes": len(background),
        "terms_tested": len(term_to_genes),
        "top_genes": [
            {
                "RAP_ID": g,
                "relations": s["relations"],
                "unique_pairs": len(s["pairs"]),
                "unique_pmids": len(s["pmids"]),
                "high_conf": s["high_conf"],
                "positive": s["positive"],
                "negative": s["negative"],
                "neutral": s["neutral"],
                "top_traits": s["traits"].most_common(5),
            }
            for g, s in ranked[:30]
        ],
        "sets": {},
        "top_significant": {},
    }
    for name, records in all_results.items():
        sig = [r for r in records if r["Adjusted P-value"] < 0.05 and r["x_overlap"] > 0]
        summary["sets"][name] = {
            "query_size": len(sets[name]),
            "mapped_to_background": len(set(sets[name]) & background),
            "hit_terms": len([r for r in records if r["x_overlap"] > 0]),
            "significant_terms_q_lt_0.05": len(sig),
        }
        summary["top_significant"][name] = [
            {k: r[k] for k in ["Term", "x_overlap", "k_query_mapped", "m_term_in_background",
                               "N_background", "fold_enrichment", "P-value", "Adjusted P-value", "Genes"]}
            for r in sig[:20]
        ]
    (RESULTS_DIR / "part4_enrichment_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"[4/4] 完成: {all_csv}")
    print(f"          {sig_csv}")
    print(f"          {RESULTS_DIR / 'part4_enrichment_summary.json'}")


if __name__ == "__main__":
    main()