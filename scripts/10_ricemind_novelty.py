# -*- coding: utf-8 -*-
"""
Step 10 (阶段四): RiceMind 交叉验证、候选新靶点排序与基因-性状网络
==========================================================================
输入:
  - results/sdp_relations.tsv    阶段三 SDP 因果关系
  - RiceMind API  http://lit-evi.hzau.edu.cn/ricemind-api/gene-profile/?gene=<RAP_ID>

输出:
  - results/part4_ricemind_crosscheck.csv
  - results/part4_novel_targets.csv
  - results/part4_network_nodes.csv
  - results/part4_network_edges.csv
  - results/part4_summary.json
"""

import csv
import json
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import networkx as nx

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"
CACHE_DIR = BASE_DIR / "data" / "ricemind_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
SDP_FILE = RESULTS_DIR / "sdp_relations.tsv"
CAUSAL = {"positive_regulation", "negative_regulation", "regulation"}
API = "http://lit-evi.hzau.edu.cn/ricemind-api/gene-profile/?gene="


def query_gene(rap_id):
    cache = CACHE_DIR / f"{rap_id}.json"
    if cache.exists():
        try:
            return json.loads(cache.read_text(encoding="utf-8"))
        except Exception:
            pass
    url = API + urllib.parse.quote(rap_id)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "HZAU-BioNLP-Course/1.0"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode("utf-8", errors="replace"))
        cache.write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")
        return data
    except Exception as e:
        return {"error": str(e)}


def main():
    rows = list(csv.DictReader(open(SDP_FILE, encoding="utf-8"), delimiter="\t"))
    causal = [r for r in rows if r["relation"] in CAUSAL]

    # 基因统计
    gene_stats = defaultdict(lambda: {
        "relations": 0, "pairs": set(), "traits": Counter(), "pmids": set(),
        "positive": 0, "negative": 0, "neutral": 0, "high": 0,
    })
    pair_counter = Counter()
    rap_to_symbol = {}
    for r in causal:
        g = r["RAP_ID"]
        rap_to_symbol.setdefault(g, r["gene_canonical"])
        s = gene_stats[g]
        s["relations"] += 1
        s["pairs"].add((g, r["RTO_ID"]))
        s["traits"][r["trait_canonical"]] += 1
        s["pmids"].add(r["PMID"])
        if r["relation"] == "positive_regulation":
            s["positive"] += 1
        elif r["relation"] == "negative_regulation":
            s["negative"] += 1
        else:
            s["neutral"] += 1
        if r["confidence"] == "high":
            s["high"] += 1
        pair_counter[(g, r["trait_canonical"], r["RTO_ID"])] += 1

    ranked = sorted(gene_stats.items(), key=lambda kv: (len(kv[1]["pairs"]), kv[1]["high"], kv[1]["relations"]), reverse=True)
    top = ranked[:100]

    print(f"[1/4] SDP 因果基因 {len(gene_stats):,} 个；查询 Top {len(top)} ...")
    cross_rows = []
    for i, (rap, s) in enumerate(top, 1):
        data = query_gene(rap)
        if "error" in data:
            found = False
            sym = ""
            rm_rap = ""
            trait_class = ""
            explanation = ""
            cgsnl = ""
            error = data["error"]
        else:
            info = data.get("official_info", {}) or {}
            ann = data.get("annotations", {}) or {}
            found = bool(info.get("standard_rap_id") or info.get("primary_id"))
            sym = info.get("primary_symbol", "") or ""
            rm_rap = info.get("standard_rap_id", "") or info.get("primary_id", "") or ""
            cgsnl = info.get("cgsnl_name", "") or ""
            trait_class = ann.get("trait_class", "") or ""
            explanation = ann.get("explanation", "") or ""
            error = ""
        text = f"{trait_class} {explanation}".lower()
        drought_mentioned = bool(any(k in text for k in ["drought", "dehydration", "water deficit", "water deprivation"]))
        aba_mentioned = bool(any(k in text for k in ["abscisic", "aba", "abscisic acid"]))
        stress_mentioned = bool(any(k in text for k in ["stress", "tolerance", "resistance"]))
        cross_rows.append({
            "RAP_ID": rap,
            "canonical_symbol": rap_to_symbol.get(rap, ""),
            "sdp_relations": s["relations"],
            "sdp_unique_pairs": len(s["pairs"]),
            "sdp_unique_pmids": len(s["pmids"]),
            "sdp_high_conf": s["high"],
            "sdp_positive": s["positive"],
            "sdp_negative": s["negative"],
            "sdp_neutral": s["neutral"],
            "top_traits": "; ".join(f"{t}({c})" for t, c in s["traits"].most_common(5)),
            "ricemind_found": int(found),
            "ricemind_symbol": sym,
            "ricemind_cgsnl_name": cgsnl,
            "ricemind_rap": rm_rap,
            "ricemind_trait_class": trait_class[:300],
            "ricemind_drought_mentioned": int(drought_mentioned),
            "ricemind_aba_mentioned": int(aba_mentioned),
            "ricemind_stress_mentioned": int(stress_mentioned),
            "ricemind_annotation_chars": len(explanation),
            "error": error,
        })
        if i % 20 == 0:
            print(f"    已查询 {i}/{len(top)} ...")
        time.sleep(0.05)

    # 写 crosscheck
    cross_csv = RESULTS_DIR / "part4_ricemind_crosscheck.csv"
    with open(cross_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(cross_rows[0].keys()))
        w.writeheader()
        w.writerows(cross_rows)

    print("[2/4] 计算覆盖率与候选新靶点 ...")
    found = [r for r in cross_rows if r["ricemind_found"] == 1]
    drought_ann = [r for r in found if r["ricemind_drought_mentioned"] == 1]
    stress_ann = [r for r in found if r["ricemind_stress_mentioned"] == 1]

    # 新靶点：SDP 证据 >=3 个唯一 pair，且 RiceMind 未收录或完全没有 stress/tolerance 注释
    core14 = {"Os09g0522200", "Os07g0154100", "Os02g0766700", "Os03g0815100",
              "Os01g0884300", "Os06g0211200", "Os09g0522000", "Os07g0448800",
              "Os05g0542500", "Os05g0455500", "Os01g0165000", "Os01g0859300",
              "Os03g0786400", "Os09g0439800"}
    candidates = []
    for r in cross_rows:
        if r["sdp_unique_pairs"] < 3 or r["RAP_ID"] in core14:
            continue
        if r["ricemind_found"] == 0:
            status = "not_in_ricemind"
        elif r["ricemind_stress_mentioned"] == 0:
            status = "in_ricemind_no_stress_annotation"
        else:
            continue
        candidates.append((status, r))
    candidates.sort(key=lambda x: (x[1]["sdp_unique_pairs"], x[1]["sdp_high_conf"], x[1]["sdp_relations"]), reverse=True)

    novel_rows = []
    for rank, (status, r) in enumerate(candidates[:20], 1):
        if status == "not_in_ricemind":
            rationale = "RiceMind 未检索到该 RAP ID；SDP 文献证据支持其与抗旱/农艺性状关联"
        else:
            rationale = "RiceMind 有该基因条目但 trait_class/annotation 无 stress/tolerance 注释；SDP 文献证据补足抗逆维度"
        novel_rows.append({
            "rank": rank,
            "RAP_ID": r["RAP_ID"],
            "ricemind_symbol": r["ricemind_symbol"],
            "status": status,
            "sdp_unique_pairs": r["sdp_unique_pairs"],
            "sdp_high_conf": r["sdp_high_conf"],
            "sdp_positive": r["sdp_positive"],
            "sdp_negative": r["sdp_negative"],
            "top_traits": r["top_traits"],
            "rationale": rationale,
        })
    novel_csv = RESULTS_DIR / "part4_novel_targets.csv"
    with open(novel_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(novel_rows[0].keys()) if novel_rows else ["rank"])
        w.writeheader()
        w.writerows(novel_rows)

    print("[3/4] 构建基因-性状网络并计算中心性 ...")
    G = nx.Graph()
    for (g, tname, tid), cnt in pair_counter.items():
        G.add_node(g, bipartite="gene")
        G.add_node(tid, bipartite="trait", name=tname)
        G.add_edge(g, tid, weight=cnt)
    deg = dict(G.degree(weight="weight"))
    gene_nodes = [n for n, d in G.nodes(data=True) if d.get("bipartite") == "gene"]
    trait_nodes = [n for n, d in G.nodes(data=True) if d.get("bipartite") == "trait"]
    top_genes = sorted(gene_nodes, key=lambda n: deg.get(n, 0), reverse=True)[:30]
    top_traits = sorted(trait_nodes, key=lambda n: deg.get(n, 0), reverse=True)[:30]

    nodes_csv = RESULTS_DIR / "part4_network_nodes.csv"
    with open(nodes_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["node", "type", "label", "weighted_degree", "degree"])
        w.writeheader()
        for n in top_genes:
            w.writerow({"node": n, "type": "gene", "label": n, "weighted_degree": deg[n], "degree": G.degree(n)})
        for n in top_traits:
            w.writerow({"node": n, "type": "trait", "label": G.nodes[n].get("name", n), "weighted_degree": deg[n], "degree": G.degree(n)})

    sub_nodes = set(top_genes) | set(top_traits)
    edges_csv = RESULTS_DIR / "part4_network_edges.csv"
    with open(edges_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["source", "target", "weight", "target_name"])
        w.writeheader()
        for u, v, d in G.edges(data=True):
            if u in sub_nodes and v in sub_nodes:
                w.writerow({"source": u, "target": v, "weight": d.get("weight", 1), "target_name": G.nodes[v].get("name", v)})

    print("[4/4] 写出汇总 ...")
    summary = {
        "sdp_causal_relations": len(causal),
        "sdp_unique_genes": len(gene_stats),
        "top100_queried": len(cross_rows),
        "ricemind_found": len(found),
        "ricemind_found_rate": len(found) / len(cross_rows) if cross_rows else 0.0,
        "ricemind_drought_annotated": len(drought_ann),
        "ricemind_aba_annotated": len([r for r in found if r["ricemind_aba_mentioned"] == 1]),
        "ricemind_stress_annotated": len(stress_ann),
        "novel_candidates": len(novel_rows),
        "network_nodes": G.number_of_nodes(),
        "network_edges": G.number_of_edges(),
        "top_genes_by_degree": [
            {"RAP_ID": n, "weighted_degree": deg[n], "degree": G.degree(n)} for n in top_genes[:15]
        ],
        "top_traits_by_degree": [
            {"RTO_ID": n, "name": G.nodes[n].get("name", n), "weighted_degree": deg[n], "degree": G.degree(n)}
            for n in top_traits[:15]
        ],
        "novel_targets_top10": novel_rows[:10],
        "crosscheck_top20": cross_rows[:20],
    }
    (RESULTS_DIR / "part4_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"    完成: {cross_csv}")
    print(f"          {novel_csv}")
    print(f"          {nodes_csv}")
    print(f"          {edges_csv}")
    print(f"          {RESULTS_DIR / 'part4_summary.json'}")


if __name__ == "__main__":
    main()