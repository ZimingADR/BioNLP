# -*- coding: utf-8 -*-
"""
Step 8 (阶段三评测): SDP 关系抽取的分层人工审计与 Precision/Recall/F1 估计
==========================================================================
抽样设计:
  A 层: 从 SDP 判为因果调控 (positive/negative/regulation) 的 639 条中随机抽 50 条;
  B 层: 从其余 7,229 条候选 (association / excluded / no_effect / unsupported) 中随机抽 50 条。
  A 层用于估计 SDP precision; B 层用于估计被 SDP 漏掉的 gold positives, 进而估计 recall。

标注定义 (写入 sdp_audit_labels.json):
  label =  1 : 句子明确表明该基因正向调控该性状 (含 loss-of-function 表型反证)
  label = -1 : 句子明确表明该基因负向调控该性状
  label =  0 : 仅共现 / 关联 / 表达相关 / 调控的是其他性状 / 否定或无效应 / 实体错配
  label = "uncertain" : 语义确实无法判断 (从 metrics 中剔除)

输出:
  - results/sdp_audit_100.csv       分层审计表 (含空 gold_label 列)
  - results/sdp_audit_labels.json   人工标签 (由审计者填写)
  - results/sdp_audit_summary.json  指标与错误分析
"""

import argparse
import csv
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"
SDP_FILE = RESULTS_DIR / "sdp_relations.tsv"
AUDIT_CSV = RESULTS_DIR / "sdp_audit_100.csv"
LABELS_FILE = RESULTS_DIR / "sdp_audit_labels.json"
SUMMARY_FILE = RESULTS_DIR / "sdp_audit_summary.json"

CAUSAL = {"positive_regulation", "negative_regulation", "regulation"}
SEED = 42
N_A = 50
N_B = 50


def load_sdp():
    return list(csv.DictReader(open(SDP_FILE, encoding="utf-8"), delimiter="\t"))


def stratified_sample(rows):
    causal = [r for r in rows if r["relation"] in CAUSAL]
    noncausal = [r for r in rows if r["relation"] not in CAUSAL]
    rng = random.Random(SEED)
    a = rng.sample(causal, min(N_A, len(causal)))
    # B 层按 relation 分层，保证 association/excluded/no_effect 都有代表
    by_rel = defaultdict(list)
    for r in noncausal:
        by_rel[r["relation"]].append(r)
    quota = {
        "association": 20,
        "excluded_direction": 15,
        "excluded_respectively": 5,
        "no_effect": 5,
        "unsupported": 5,
    }
    b = []
    for rel, k in quota.items():
        pool = by_rel.get(rel, [])
        b.extend(rng.sample(pool, min(k, len(pool))))
    # 若某些层不足，从剩余池补齐到 50
    remaining = [r for r in noncausal if r not in b]
    if len(b) < N_B:
        b.extend(rng.sample(remaining, min(N_B - len(b), len(remaining))))
    return a[:N_A], b[:N_B]


def write_audit_csv(a, b):
    rows = []
    for stratum, items in [("A_sdp_positive", a), ("B_sdp_negative", b)]:
        for i, r in enumerate(items, 1):
            rows.append({
                "audit_id": f"{stratum[0]}{i:02d}",
                "stratum": stratum,
                "PMID": r["PMID"],
                "sent_id": r["sent_id"],
                "gene_mention": r["gene_mention"],
                "gene_canonical": r["gene_canonical"],
                "RAP_ID": r["RAP_ID"],
                "trait_mention": r["trait_mention"],
                "trait_canonical": r["trait_canonical"],
                "RTO_ID": r["RTO_ID"],
                "sdp_relation": r["relation"],
                "sdp_polarity": r["polarity"],
                "sdp_pattern": r["relation_pattern"],
                "sdp_confidence": r["confidence"],
                "predicate": r["predicate"],
                "sdp_path": r["sdp_path"],
                "sentence": r["sentence"],
                "gold_label": "",
                "gold_polarity": "",
                "note": "",
            })
    with open(AUDIT_CSV, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"写出审计表: {AUDIT_CSV} ({len(rows)} 条)")


def wilson_ci(k, n, z=1.96):
    if n == 0:
        return [0.0, 0.0]
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * n)) / n) / denom
    return [max(0.0, center - half), min(1.0, center + half)]


def evaluate():
    rows = list(csv.DictReader(open(AUDIT_CSV, encoding="utf-8")))
    labels = json.loads(LABELS_FILE.read_text(encoding="utf-8"))
    by_id = {r["audit_id"]: r for r in rows}
    for r in rows:
        lab = labels.get(r["audit_id"], {})
        r["label"] = lab.get("label", "uncertain")
        r["note"] = lab.get("note", "")

    A = [r for r in rows if r["stratum"] == "A_sdp_positive" and r["label"] != "uncertain"]
    B = [r for r in rows if r["stratum"] == "B_sdp_negative" and r["label"] != "uncertain"]
    # 总体规模
    total = len(load_sdp())
    n_causal = sum(1 for r in load_sdp() if r["relation"] in CAUSAL)
    n_noncausal = total - n_causal

    def is_pos(x):
        return x["label"] in {1, -1, "1", "-1", "reg"}

    posA = sum(1 for r in A if is_pos(r))
    posB = sum(1 for r in B if is_pos(r))

    # 估计量
    precision = posA / len(A) if A else 0.0
    missed_rate = posB / len(B) if B else 0.0
    tp_total_est = precision * n_causal
    missed_total_est = missed_rate * n_noncausal
    recall = tp_total_est / (tp_total_est + missed_total_est) if (tp_total_est + missed_total_est) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    baseline_precision = (posA + posB) / (len(A) + len(B)) if (A or B) else 0.0
    # 分层加权估计的全池共现 precision
    baseline_precision_weighted = (
        (posA / len(A) * n_causal + posB / len(B) * n_noncausal) / total
        if A and B else baseline_precision
    )

    # 极性准确率: 只在 A 层 gold positive 上比较
    pol_ok, pol_n = 0, 0
    polar_conf = Counter()
    for r in A:
        if r["label"] not in {1, -1, "1", "-1"}:
            continue
        gold_sign = 1 if r["label"] in {1, "1"} else -1
        pred_sign = 1 if r["sdp_polarity"] == "positive" else (-1 if r["sdp_polarity"] == "negative" else 0)
        if pred_sign == 0:
            continue
        pol_n += 1
        if pred_sign == gold_sign:
            pol_ok += 1
        polar_conf[(gold_sign, pred_sign)] += 1

    # 分类错误分析
    errors = []
    for r in A:
        if not is_pos(r):
            errors.append({
                "audit_id": r["audit_id"], "type": "SDP_false_positive",
                "gene": r["gene_canonical"], "trait": r["trait_canonical"],
                "relation": r["sdp_relation"], "sentence": r["sentence"],
                "note": r["note"],
            })
    for r in B:
        if is_pos(r):
            errors.append({
                "audit_id": r["audit_id"], "type": "SDP_false_negative",
                "gene": r["gene_canonical"], "trait": r["trait_canonical"],
                "relation": r["sdp_relation"], "sentence": r["sentence"],
                "note": r["note"],
            })

    summary = {
        "sample_design": {
            "seed": SEED, "A_population": n_causal, "B_population": n_noncausal,
            "A_sample": len(A), "B_sample": len(B),
        },
        "gold_positive_in_A": posA,
        "gold_positive_in_B": posB,
        "sdp_precision": precision,
        "sdp_precision_95ci": wilson_ci(posA, len(A)),
        "sdp_recall_estimate": recall,
        "sdp_f1_estimate": f1,
        "estimated_tp_total": tp_total_est,
        "estimated_missed_total": missed_total_est,
        "cooccurrence_baseline_precision_unweighted": baseline_precision,
        "cooccurrence_baseline_precision_weighted": baseline_precision_weighted,
        "polarity_accuracy": (pol_ok / pol_n) if pol_n else None,
        "polarity_n": pol_n,
        "polarity_confusion": {f"gold{g}_pred{p}": c for (g, p), c in polar_conf.items()},
        "error_cases": errors,
    }
    SUMMARY_FILE.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", action="store_true")
    parser.add_argument("--evaluate", action="store_true")
    args = parser.parse_args()
    if args.sample:
        rows = load_sdp()
        a, b = stratified_sample(rows)
        write_audit_csv(a, b)
    elif args.evaluate:
        evaluate()
    else:
        print("请使用 --sample 或 --evaluate")


if __name__ == "__main__":
    main()