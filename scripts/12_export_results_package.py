# -*- coding: utf-8 -*-
"""
Step 12: 导出“可直接写文章”的结果表、关键数字、数据字典与图表索引
输出目录: results/tables/
"""

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"
TABLES_DIR = RESULTS_DIR / "tables"
TABLES_DIR.mkdir(parents=True, exist_ok=True)

CAUSAL = {"positive_regulation", "negative_regulation", "regulation"}


def load_json(name):
    return json.loads((RESULTS_DIR / name).read_text(encoding="utf-8"))


def load_tsv(name):
    return list(csv.DictReader(open(RESULTS_DIR / name, encoding="utf-8"), delimiter="\t"))


def write_csv(name, rows, fields=None):
    path = TABLES_DIR / name
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = fields or list(rows[0].keys())
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  {path.name}: {len(rows)} rows")


def main():
    print("[1/6] 读取汇总 ...")
    part1 = load_json("part1_lexical_and_tfidf_summary.json")
    part2 = load_json("part2_normalization_and_evidence_report.json")
    sdp = load_json("sdp_summary.json")
    audit = load_json("sdp_audit_summary.json")
    enr = load_json("part4_enrichment_summary.json")
    p4 = load_json("part4_summary.json")
    sdp_rows = load_tsv("sdp_relations.tsv")
    labels = load_json("sdp_audit_labels.json")

    print("[2/6] 发表级表格 ...")
    # Table: corpus stats
    corpus = [
        {"Metric": "PubMed_hits", "Value": 3311, "Note": "esearch 命中"},
        {"Metric": "PubTator_documents_parsed", "Value": part1["num_documents"], "Note": "5 篇 2026 新文献未被 PubTator 索引"},
        {"Metric": "Sentences", "Value": part1["num_sentences"], "Note": "NLTK sent_tokenize"},
        {"Metric": "Tokens", "Value": part1["ttr_metrics"]["total_tokens_N"], "Note": "小写词元"},
        {"Metric": "Types", "Value": part1["ttr_metrics"]["unique_types_V"], "Note": "独立词型"},
        {"Metric": "Raw_TTR", "Value": round(part1["ttr_metrics"]["raw_ttr"], 6), "Note": "V/N"},
        {"Metric": "PubTator_entities", "Value": part1["pubtator_entities_total"], "Note": ""},
        {"Metric": "PubTator_offset_verified", "Value": part1["pubtator_entities_offset_verified"], "Note": "text[b:e]==mention"},
        {"Metric": "Total_entity_mentions", "Value": part2["total_entity_mentions"], "Note": "基因+性状"},
        {"Metric": "RAP_mapped_gene_mentions", "Value": part2["rap_mapped_gene_mentions"], "Note": "水稻基因"},
        {"Metric": "Trait_mentions", "Value": part2["total_trait_mentions"], "Note": "RTO 归一化"},
        {"Metric": "PubTator_only_gene_mentions", "Value": part2["total_gene_mentions"] - part2["rap_mapped_gene_mentions"], "Note": "非水稻/未映射"},
        {"Metric": "Rice_candidate_pairs", "Value": 7868, "Note": "RAP ID 以 Os 开头"},
        {"Metric": "Rice_unique_pairs", "Value": 4678, "Note": ""},
        {"Metric": "Rice_unique_genes", "Value": 1344, "Note": ""},
        {"Metric": "Rice_unique_traits", "Value": 190, "Note": ""},
        {"Metric": "SDP_causal_relations", "Value": sdp["supported_relations_rice"], "Note": "positive/negative/regulation"},
        {"Metric": "SDP_unique_causal_pairs", "Value": sdp["supported_unique_pairs_rice"], "Note": ""},
        {"Metric": "SDP_association_relations", "Value": sdp["relation_counts"].get("association", 0), "Note": ""},
        {"Metric": "SDP_audit_precision", "Value": audit["sdp_precision"], "Note": "100 条分层审计"},
        {"Metric": "SDP_audit_recall", "Value": audit["sdp_recall_estimate"], "Note": "分层外推"},
        {"Metric": "SDP_audit_F1", "Value": audit["sdp_f1_estimate"], "Note": ""},
        {"Metric": "SDP_polarity_accuracy", "Value": audit["polarity_accuracy"], "Note": "34 条可判极性"},
        {"Metric": "Enrichment_significant_terms_top100", "Value": enr["sets"]["top100"]["significant_terms_q_lt_0.05"], "Note": "BH q<0.05"},
        {"Metric": "RiceMind_coverage_top100", "Value": p4["ricemind_found_rate"], "Note": ""},
        {"Metric": "RiceMind_stress_annotated_top100", "Value": p4["ricemind_stress_annotated"], "Note": ""},
        {"Metric": "RiceMind_drought_annotated_top100", "Value": p4["ricemind_drought_annotated"], "Note": ""},
        {"Metric": "Novel_annotation_gap_candidates", "Value": p4["novel_candidates"], "Note": ""},
        {"Metric": "Network_nodes", "Value": p4["network_nodes"], "Note": ""},
        {"Metric": "Network_edges", "Value": p4["network_edges"], "Note": ""},
    ]
    write_csv("table1_corpus_and_key_stats.csv", corpus)

    # Table: TTR
    ttr_rows = []
    for k, v in part1["ttr_metrics"]["standardized_samples"].items():
        ttr_rows.append({
            "Window": k,
            "Head_types": v["head_types"],
            "Head_TTR": round(v["head_ttr"], 6),
            "Shuf_TTR": round(v["shuf_mean_ttr"], 6),
            "Delta_pp": round((v["shuf_mean_ttr"] - v["head_ttr"]) * 100, 2),
        })
    ttr_rows.append({
        "Window": "Full",
        "Head_types": part1["ttr_metrics"]["unique_types_V"],
        "Head_TTR": round(part1["ttr_metrics"]["raw_ttr"], 6),
        "Shuf_TTR": round(part1["ttr_metrics"]["raw_ttr"], 6),
        "Delta_pp": 0.0,
    })
    write_csv("table2_ttr.csv", ttr_rows)

    # Table: TF-IDF top terms
    tfidf = list(csv.DictReader(open(RESULTS_DIR / "part1_tfidf_top_terms.tsv", encoding="utf-8"), delimiter="\t"))
    write_csv("table3_tfidf_top_terms.csv", tfidf)

    # Table: top genes / traits / pairs
    top_genes = []
    for i, r in enumerate(part2["top_genes"][:20], 1):
        top_genes.append({"Rank": i, "Gene": r["canonical_symbol"], "RAP_ID": r["rap_id"], "Mentions": r["count"]})
    write_csv("table4_top_genes.csv", top_genes)
    top_traits = []
    for i, r in enumerate(part2["top_traits"][:20], 1):
        top_traits.append({"Rank": i, "Trait": r["canonical_trait"], "RTO_ID": r["rto_id"], "Mentions": r["count"]})
    write_csv("table5_top_traits.csv", top_traits)
    top_pairs = []
    for i, r in enumerate(part2["top_gene_trait_pairs"][:20], 1):
        top_pairs.append({
            "Rank": i, "Gene": r["gene_symbol"], "RAP_ID": r["rap_id"],
            "Trait": r["trait_name"], "RTO_ID": r["rto_id"],
            "Evidence_sentences": r["evidence_sentences_count"],
        })
    write_csv("table6_top_gene_trait_pairs.csv", top_pairs)

    # Table: SDP relation counts
    rel = {k: v for k, v in sdp["relation_counts"].items() if not k.startswith("rice_")}
    total = sum(rel.values())
    rel_rows = [{"Relation": k, "Count": v, "Share_pct": round(v / total * 100, 2)} for k, v in sorted(rel.items(), key=lambda x: -x[1])]
    write_csv("table7_sdp_relation_counts.csv", rel_rows)

    # Table: confidence
    conf = sdp["rice_confidence_counts"]
    conf_total = sum(conf.values())
    write_csv("table8_sdp_confidence.csv", [
        {"Confidence": k, "Count": v, "Share_pct": round(v / conf_total * 100, 2)} for k, v in conf.items()
    ])

    # Table: top predicates
    write_csv("table9_sdp_top_predicates.csv", [
        {"Rank": i, "Predicate": p, "Count": c} for i, (p, c) in enumerate(sdp["predicate_top"][:25], 1)
    ])

    # Table: audit 100 with labels
    audit_rows = list(csv.DictReader(open(RESULTS_DIR / "sdp_audit_100.csv", encoding="utf-8")))
    out_audit = []
    for r in audit_rows:
        lab = labels.get(r["audit_id"], {})
        gold = str(lab.get("label", "uncertain"))
        sdp_causal = r["sdp_relation"] in CAUSAL
        gold_causal = gold in {"1", "-1", "reg"}
        polarity_ok = ""
        if gold in {"1", "-1"} and r["sdp_polarity"] in {"positive", "negative"}:
            polarity_ok = int((gold == "1") == (r["sdp_polarity"] == "positive"))
        out_audit.append({
            "audit_id": r["audit_id"], "stratum": r["stratum"],
            "PMID": r["PMID"], "sent_id": r["sent_id"],
            "gene_mention": r["gene_mention"], "gene_canonical": r["gene_canonical"], "RAP_ID": r["RAP_ID"],
            "trait_mention": r["trait_mention"], "trait_canonical": r["trait_canonical"], "RTO_ID": r["RTO_ID"],
            "sdp_relation": r["sdp_relation"], "sdp_polarity": r["sdp_polarity"],
            "sdp_confidence": r["sdp_confidence"], "predicate": r["predicate"], "sdp_path": r["sdp_path"],
            "gold_label": gold, "sdp_causal": int(sdp_causal), "gold_causal": int(gold_causal),
            "polarity_correct": polarity_ok, "note": lab.get("note", ""), "sentence": r["sentence"],
        })
    write_csv("table10_sdp_audit_100.csv", out_audit)

    # Error categories
    err = Counter()
    for r in out_audit:
        if r["sdp_causal"] == 1 and r["gold_causal"] == 0:
            err["SDP_false_positive"] += 1
        if r["sdp_causal"] == 0 and r["gold_causal"] == 1:
            err["SDP_false_negative"] += 1
        if r["polarity_correct"] == 0:
            err["polarity_error"] += 1
    write_csv("table11_sdp_error_categories.csv", [{"Category": k, "Count": v} for k, v in err.most_common()])

    # Enrichment top 25
    enr_rows = []
    for i, r in enumerate(enr["top_significant"]["top100"][:25], 1):
        enr_rows.append({
            "Rank": i, "Term": r["Term"], "x_overlap": r["x_overlap"], "k_query": r["k_query_mapped"],
            "m_background": r["m_term_in_background"], "N_background": r["N_background"],
            "Fold_enrichment": round(r["fold_enrichment"], 4), "P_value": r["P-value"],
            "Adjusted_P": r["Adjusted P-value"], "Genes": r["Genes"],
        })
    write_csv("table12_enrichment_top25.csv", enr_rows)

    # RiceMind coverage
    write_csv("table13_ricemind_coverage.csv", [
        {"Metric": "Top100_queried", "Value": p4["top100_queried"]},
        {"Metric": "Found_in_RiceMind", "Value": p4["ricemind_found"]},
        {"Metric": "Coverage_rate", "Value": p4["ricemind_found_rate"]},
        {"Metric": "Stress_annotated", "Value": p4["ricemind_stress_annotated"]},
        {"Metric": "Drought_annotated", "Value": p4["ricemind_drought_annotated"]},
        {"Metric": "ABA_annotated", "Value": p4["ricemind_aba_annotated"]},
        {"Metric": "Novel_candidates", "Value": p4["novel_candidates"]},
    ])

    # Network top
    net = list(csv.DictReader(open(RESULTS_DIR / "part4_network_nodes.csv", encoding="utf-8")))
    net_rows = [{"Rank": i, "Node": n["node"], "Type": n["type"], "Label": n["label"],
                 "Weighted_degree": n["weighted_degree"], "Degree": n["degree"]}
                for i, n in enumerate(net[:30], 1)]
    write_csv("table14_network_top30.csv", net_rows)

    # Novel targets copy
    novel = list(csv.DictReader(open(RESULTS_DIR / "part4_novel_targets.csv", encoding="utf-8")))
    write_csv("table15_novel_targets.csv", novel)

    # SDP top causal pairs
    pair_counter = Counter()
    for r in sdp_rows:
        if r["relation"] in CAUSAL:
            pair_counter[(r["gene_canonical"], r["RAP_ID"], r["trait_canonical"], r["RTO_ID"], r["relation"])] += 1
    pair_rows = []
    for i, ((g, rid, t, tid, rel), c) in enumerate(pair_counter.most_common(30), 1):
        pair_rows.append({"Rank": i, "Gene": g, "RAP_ID": rid, "Trait": t, "RTO_ID": tid,
                          "Relation": rel, "Evidence_sentences": c})
    write_csv("table16_sdp_top_causal_pairs.csv", pair_rows)

    print("[3/6] key_numbers.json 与 rice-only / causal 细分统计 ...")
    # Rice-only co-occurrence stats
    ev_rows = load_tsv("rice_drought_sentence_evidence.tsv")
    rice_ev = [r for r in ev_rows if r["RAP_ID"].startswith("Os")]
    rice_pairs = Counter((r["RAP_ID"], r["RTO_ID"]) for r in rice_ev)
    rice_genes = Counter((r["gene_canonical"], r["RAP_ID"]) for r in rice_ev)
    rice_traits = Counter((r["trait_canonical"], r["RTO_ID"]) for r in rice_ev)
    rice_stats = {
        "all_candidate_pairs": len(ev_rows),
        "non_rice_candidate_pairs": len(ev_rows) - len(rice_ev),
        "rice_candidate_pairs": len(rice_ev),
        "rice_unique_pairs": len(rice_pairs),
        "rice_unique_genes": len(rice_genes),
        "rice_unique_traits": len(rice_traits),
        "top_rice_genes": [
            {"gene": g, "RAP_ID": rid, "mentions": c}
            for (g, rid), c in rice_genes.most_common(20)
        ],
        "top_rice_traits": [
            {"trait": t, "RTO_ID": tid, "mentions": c}
            for (t, tid), c in rice_traits.most_common(20)
        ],
        "top_rice_pairs": [
            {"RAP_ID": rid, "RTO_ID": tid, "evidence_sentences": c}
            for (rid, tid), c in rice_pairs.most_common(20)
        ],
    }
    (RESULTS_DIR / "part2_rice_only_stats.json").write_text(
        json.dumps(rice_stats, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # Causal SDP stats
    causal_rows = [r for r in sdp_rows if r["relation"] in CAUSAL]
    causal_stats = {
        "causal_relations": len(causal_rows),
        "unique_pairs": len({(r["RAP_ID"], r["RTO_ID"]) for r in causal_rows}),
        "unique_genes": len({r["RAP_ID"] for r in causal_rows}),
        "unique_traits": len({r["RTO_ID"] for r in causal_rows}),
        "unique_pmids": len({r["PMID"] for r in causal_rows}),
        "unique_sentences": len({(r["PMID"], r["sent_id"]) for r in causal_rows}),
        "confidence": dict(Counter(r["confidence"] for r in causal_rows)),
        "patterns": dict(Counter(r["relation_pattern"] for r in causal_rows)),
        "constructs": dict(Counter(r["construct"] for r in causal_rows)),
        "relation_counts": dict(Counter(r["relation"] for r in causal_rows)),
        "top_causal_pairs": [
            {"gene": g, "RAP_ID": rid, "trait": t, "RTO_ID": tid,
             "relation": rel, "evidence_sentences": c}
            for (g, rid, t, tid, rel), c in pair_counter.most_common(30)
        ],
    }
    (RESULTS_DIR / "sdp_causal_stats.json").write_text(
        json.dumps(causal_stats, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("  part2_rice_only_stats.json")
    print("  sdp_causal_stats.json")

    key = {
        "generated_at": "2026-11-25",
        "corpus": {
            "pubmed_hits": 3311, "pubtator_docs": part1["num_documents"],
            "sentences": part1["num_sentences"], "tokens": part1["ttr_metrics"]["total_tokens_N"],
            "types": part1["ttr_metrics"]["unique_types_V"], "raw_ttr": part1["ttr_metrics"]["raw_ttr"],
        },
        "entities": {
            "total_mentions": part2["total_entity_mentions"],
            "rap_mapped_gene_mentions": part2["rap_mapped_gene_mentions"],
            "trait_mentions": part2["total_trait_mentions"],
            "pubtator_gene_mentions": part2["pubtator_gene_mentions"],
            "pubtator_mapped_to_rap": part2["pubtator_mapped_to_rap"],
            "rice_candidate_pairs": 7868,
            "rice_unique_pairs": 4678,
            "rice_unique_genes": 1344,
            "rice_unique_traits": 190,
        },
        "sdp": {
            "causal": sdp["supported_relations_rice"],
            "unique_causal_pairs": sdp["supported_unique_pairs_rice"],
            "positive": sdp["positive_relations_rice"],
            "negative": sdp["negative_relations_rice"],
            "neutral": sdp["relation_counts"].get("regulation", 0),
            "association": sdp["relation_counts"].get("association", 0),
            "excluded_direction": sdp["relation_counts"].get("excluded_direction", 0),
            "excluded_respectively": sdp["relation_counts"].get("excluded_respectively", 0),
            "no_effect": sdp["relation_counts"].get("no_effect", 0),
            "unsupported": sdp["relation_counts"].get("unsupported", 0),
            "confidence": sdp["rice_confidence_counts"],
        },
        "audit": {
            "sample_size": 100,
            "precision": audit["sdp_precision"],
            "precision_95ci": audit["sdp_precision_95ci"],
            "recall": audit["sdp_recall_estimate"],
            "f1": audit["sdp_f1_estimate"],
            "baseline_precision": audit["cooccurrence_baseline_precision_weighted"],
            "polarity_accuracy": audit["polarity_accuracy"],
            "polarity_n": audit["polarity_n"],
            "confusion": audit["polarity_confusion"],
        },
        "enrichment": {
            "top100_significant": enr["sets"]["top100"]["significant_terms_q_lt_0.05"],
            "top50_significant": enr["sets"]["top50"]["significant_terms_q_lt_0.05"],
            "top20_significant": enr["sets"]["top20"]["significant_terms_q_lt_0.05"],
            "core14_significant": enr["sets"]["core14"]["significant_terms_q_lt_0.05"],
        },
        "ricemind": {
            "coverage": p4["ricemind_found_rate"],
            "stress_annotated": p4["ricemind_stress_annotated"],
            "drought_annotated": p4["ricemind_drought_annotated"],
            "aba_annotated": p4["ricemind_aba_annotated"],
            "novel_candidates": p4["novel_candidates"],
            "network_nodes": p4["network_nodes"],
            "network_edges": p4["network_edges"],
        },
    }
    (RESULTS_DIR / "key_numbers.json").write_text(json.dumps(key, ensure_ascii=False, indent=2), encoding="utf-8")
    print("  key_numbers.json")

    print("[4/6] 图表索引 ...")
    fig_rows = [
        {"Figure": "Fig. 1", "File": "fig1_ttr_heaps.pdf", "Panels": "a,b",
         "Content": "Heaps law / TTR head vs shuf", "Source": "part1_lexical_and_tfidf_summary.json"},
        {"Figure": "Fig. 2", "File": "fig3_pipeline_ablation.pdf", "Panels": "a,b",
         "Content": "PubTator vs dual-dictionary ablation; KG scale", "Source": "part2_normalization_and_evidence_report.json"},
        {"Figure": "Fig. 3", "File": "fig4_top_genes_traits.pdf", "Panels": "a,b",
         "Content": "Top14 genes and traits", "Source": "part2_normalization_and_evidence_report.json"},
        {"Figure": "Fig. 4", "File": "fig5_sdp_framework.pdf", "Panels": "a,b,c",
         "Content": "SDP path example; filtering funnel; polarity composition", "Source": "sdp_relations.tsv, sdp_summary.json"},
        {"Figure": "Fig. 5", "File": "fig6_evaluation_enrichment.pdf", "Panels": "a,b",
         "Content": "Audited P/R/F1; TO/RTO enrichment", "Source": "sdp_audit_summary.json, part4_enrichment_summary.json"},
        {"Figure": "Fig. 6", "File": "fig7_gene_trait_network.pdf", "Panels": "single",
         "Content": "Top gene-trait bipartite network", "Source": "part4_network_nodes.csv, part4_network_edges.csv"},
    ]
    write_csv("table17_figure_index.csv", fig_rows)

    print("[5/6] 数据清单 ...")
    manifest = []
    for p in sorted(RESULTS_DIR.glob("*")):
        if p.is_file() and p.suffix in {".tsv", ".csv", ".json"}:
            try:
                if p.suffix == ".json":
                    obj = json.loads(p.read_text(encoding="utf-8"))
                    rows = len(obj) if isinstance(obj, list) else "JSON object"
                else:
                    with open(p, encoding="utf-8") as f:
                        rows = sum(1 for _ in f) - 1
                manifest.append({"File": p.name, "Rows": rows, "Size_KB": round(p.stat().st_size / 1024, 1), "Type": p.suffix[1:]})
            except Exception:
                manifest.append({"File": p.name, "Rows": "?", "Size_KB": round(p.stat().st_size / 1024, 1), "Type": p.suffix[1:]})
    write_csv("table18_data_manifest.csv", manifest)

    print("[6/6] 完成。")


if __name__ == "__main__":
    main()