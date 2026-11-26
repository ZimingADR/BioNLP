# -*- coding: utf-8 -*-
"""
Step 11: 生成阶段三/四论文图
  fig5_sdp_framework.pdf   (a) SDP 路径与 loss-of-function 极性翻转示例
                           (b) 7,868 候选过滤漏斗
                           (c) 884 条因果关系的极性构成
  fig6_evaluation_enrichment.pdf
                           (a) 共现基线 vs SDP 的 Precision/Recall/F1
                           (b) Top100 SDP 基因集的 TO/RTO 富集气泡图
  fig7_gene_trait_network.pdf
                           Top 基因-性状二分网络
"""

import csv
import json
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"
FIG_DIR = RESULTS_DIR / "figures"
BIBM_FIG_DIR = BASE_DIR / "paper_bibm" / "BIBM-Latex" / "figure"
FIG_DIR.mkdir(parents=True, exist_ok=True)
BIBM_FIG_DIR.mkdir(parents=True, exist_ok=True)

PALETTE = {
    "navy": "#0F4D92",
    "blue": "#1F65A6",
    "sky": "#3775BA",
    "steel": "#6A96C8",
    "grey": "#99A9BF",
    "dark": "#272727",
    "mid": "#6E7785",
    "red": "#B64342",
    "green": "#2E7D5B",
    "amber": "#C98A1B",
    "purple": "#7B5AA6",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "Liberation Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.9,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "legend.frameon": False,
    "figure.dpi": 300,
})


def save(fig, stem):
    for fmt in ["pdf", "svg", "png"]:
        fig.savefig(BIBM_FIG_DIR / f"{stem}.{fmt}", format=fmt, dpi=300, bbox_inches="tight")
        shutil.copyfile(BIBM_FIG_DIR / f"{stem}.{fmt}", FIG_DIR / f"{stem}.{fmt}")
    plt.close(fig)


def panel_label(ax, label, x=-0.14, y=1.03, fs=10):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=fs, fontweight="bold",
            ha="left", va="bottom", color=PALETTE["dark"])


def box(ax, xy, text, color, width=0.25, height=0.14, fs=8.4):
    x, y = xy
    p = FancyBboxPatch((x - width / 2, y - height / 2), width, height,
                       boxstyle="round,pad=0.012,rounding_size=0.02",
                       linewidth=1.0, edgecolor=color, facecolor=color + "1A")
    ax.add_patch(p)
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=PALETTE["dark"])


def arrow(ax, start, end, label, color, rad=0.0):
    a = FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=9,
                        linewidth=1.1, color=color,
                        connectionstyle=f"arc3,rad={rad}")
    ax.add_patch(a)
    mx, my = (start[0] + end[0]) / 2, (start[1] + end[1]) / 2
    ax.text(mx, my + 0.035, label, ha="center", va="bottom", fontsize=7.0,
            color=color, fontstyle="italic")


def fig5():
    sdp = json.loads((RESULTS_DIR / "sdp_summary.json").read_text(encoding="utf-8"))
    rel = sdp["relation_counts"]
    counts = {k: v for k, v in rel.items() if not k.startswith("rice_")}

    fig = plt.figure(figsize=(7.2, 2.45))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.05, 0.75], wspace=0.45)

    # (a) SDP 路径示意（纵向链条，避免重叠）
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    nodes = [
        (0.36, 0.88, "osbZIP23", "gene / pobj", PALETTE["sky"]),
        (0.36, 0.66, "mutations", "nsubj", PALETTE["steel"]),
        (0.36, 0.44, "weakened", "ccomp / negative", PALETTE["red"]),
        (0.36, 0.22, "drought tolerance", "dobj / trait", PALETTE["green"]),
    ]
    for x, y, label, sub, color in nodes:
        box(ax, (x, y), label, color, width=0.58, height=0.115, fs=7.4)
        ax.text(0.68, y, sub, fontsize=6.3, va="center", color=color, style="italic")
    for y0, y1, lab in [(0.805, 0.725, "pobj"), (0.585, 0.505, "nsubj"), (0.365, 0.285, "dobj")]:
        arrow(ax, (0.36, y0), (0.36, y1), "", PALETTE["mid"], rad=0)
        ax.text(0.40, (y0 + y1) / 2, lab, fontsize=6.2, va="center", ha="left",
                color=PALETTE["mid"], style="italic")
    ax.text(0.5, 0.045, "LOSS × NEG → gene effect POSITIVE",
            ha="center", fontsize=6.6, color=PALETTE["navy"], fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.22", facecolor="#EAF1F8", edgecolor="none"))
    panel_label(ax, "a", x=-0.16, y=1.02)

    # (b) 过滤漏斗
    ax = fig.add_subplot(gs[0, 1])
    stages = ["Candidates", "Causal", "Association", "Filtered"]
    vals = [7868, counts.get("positive_regulation", 0) + counts.get("negative_regulation", 0) + counts.get("regulation", 0),
            counts.get("association", 0),
            counts.get("excluded_direction", 0) + counts.get("excluded_respectively", 0) +
            counts.get("no_effect", 0) + counts.get("unsupported", 0)]
    colors = [PALETTE["grey"], PALETTE["navy"], PALETTE["sky"], PALETTE["steel"]]
    bars = ax.barh(range(len(stages))[::-1], vals, color=colors, height=0.60)
    ax.set_yticks(range(len(stages))[::-1])
    ax.set_yticklabels(stages, fontsize=7.0)
    ax.set_xscale("log")
    ax.set_xlim(20, 15000)
    ax.set_xlabel("Candidate pairs (log)", fontsize=7.0)
    for b, v in zip(bars, vals):
        ax.text(v * 1.15, b.get_y() + b.get_height() / 2, f"{v:,}", va="center", fontsize=6.4,
                color=PALETTE["dark"])
    panel_label(ax, "b", x=-0.30, y=1.02)

    # (c) 因果极性构成
    ax = fig.add_subplot(gs[0, 2])
    pol = [counts.get("positive_regulation", 0), counts.get("negative_regulation", 0), counts.get("regulation", 0)]
    labels = ["Pos.", "Neg.", "Neutral"]
    colors = [PALETTE["navy"], PALETTE["red"], PALETTE["steel"]]
    bars = ax.bar(labels, pol, color=colors, width=0.62)
    for b, v in zip(bars, pol):
        ax.text(b.get_x() + b.get_width() / 2, v + 10, f"{v}", ha="center", fontsize=6.5)
    ax.set_ylim(0, max(pol) * 1.20)
    ax.set_ylabel("Causal relations", fontsize=6.8)
    ax.tick_params(axis="x", labelsize=6.6)
    panel_label(ax, "c", x=-0.34, y=1.02)

    fig.suptitle("SDP relation extraction and polarity-aware filtering",
                 fontsize=8.2, fontweight="bold", y=1.04)
    save(fig, "fig5_sdp_framework")


def fig6():
    audit = json.loads((RESULTS_DIR / "sdp_audit_summary.json").read_text(encoding="utf-8"))
    enr = json.loads((RESULTS_DIR / "part4_enrichment_summary.json").read_text(encoding="utf-8"))

    fig = plt.figure(figsize=(7.2, 2.55))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.82, 1.18], wspace=0.48)

    # (a) P/R/F1
    ax = fig.add_subplot(gs[0, 0])
    metrics = ["Precision", "Recall", "F1"]
    baseline = [audit["cooccurrence_baseline_precision_weighted"], 1.0, 0.0]
    sdp = [audit["sdp_precision"], audit["sdp_recall_estimate"], audit["sdp_f1_estimate"]]
    x = np.arange(len(metrics))
    w = 0.34
    b1 = ax.bar(x - w / 2, baseline, w, color=PALETTE["grey"], label="Co-occurrence")
    b2 = ax.bar(x + w / 2, sdp, w, color=PALETTE["navy"], label="SDP")
    for i, b in enumerate(list(b1) + list(b2)):
        h = b.get_height()
        ytxt = h + 0.03
        if abs(h - sdp[0]) < 1e-6 and i >= len(baseline):
            ytxt = h + 0.10
        ax.text(b.get_x() + b.get_width() / 2, ytxt, f"{h:.2f}", ha="center", fontsize=6.2)
    lo, hi = audit["sdp_precision_95ci"]
    ax.errorbar(x[0] + w / 2, sdp[0], yerr=[[sdp[0] - lo], [hi - sdp[0]]], fmt="none",
                ecolor=PALETTE["dark"], capsize=2.5, linewidth=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics, fontsize=6.8)
    ax.set_ylim(0, 1.24)
    ax.set_ylabel("Score", fontsize=6.8)
    ax.legend(fontsize=6.0, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.02))
    panel_label(ax, "a", x=-0.18, y=1.03)

    # (b) enrichment bubble
    ax = fig.add_subplot(gs[0, 1])
    terms = enr["top_significant"]["top100"][:12][::-1]
    y = np.arange(len(terms))
    folds = [t["fold_enrichment"] for t in terms]
    overlaps = [t["x_overlap"] for t in terms]
    qvals = [max(t["Adjusted P-value"], 1e-300) for t in terms]
    colors = -np.log10(qvals)
    sizes = np.array(overlaps) * 4.8
    sc = ax.scatter(folds, y, s=sizes, c=colors, cmap="Blues", edgecolors=PALETTE["navy"],
                    linewidths=0.45, alpha=0.9, zorder=3)
    ax.set_yticks(y)
    labels = [t["Term"].split(" (")[0] for t in terms]
    ax.set_yticklabels(labels, fontsize=5.9)
    ax.set_ylim(-0.9, len(terms) + 0.7)
    ax.set_xlabel("Fold enrichment", fontsize=6.6)
    ax.set_xlim(0, max(folds) * 1.9)
    for i, (f, o) in enumerate(zip(folds, overlaps)):
        ax.text(f + 0.9, i, f"x={o}", va="center", fontsize=5.0, color=PALETTE["mid"])
    cbar = fig.colorbar(sc, ax=ax, fraction=0.032, pad=0.015)
    cbar.set_label(r"$-\log_{10}q$", fontsize=5.8, labelpad=2)
    cbar.ax.tick_params(labelsize=5.0)
    panel_label(ax, "b", x=-0.28, y=1.03)

    fig.suptitle("Audited SDP performance and TO/RTO enrichment",
                 fontsize=8.2, fontweight="bold", y=1.04)
    save(fig, "fig6_evaluation_enrichment")


def fig7():
    nodes = list(csv.DictReader(open(RESULTS_DIR / "part4_network_nodes.csv", encoding="utf-8")))
    edges = list(csv.DictReader(open(RESULTS_DIR / "part4_network_edges.csv", encoding="utf-8")))
    genes = [n for n in nodes if n["type"] == "gene"][:10]
    traits = [n for n in nodes if n["type"] == "trait"][:8]
    # RAP -> symbol mapping for readable labels
    rap_to_symbol = {}
    for r in csv.DictReader(open(RESULTS_DIR / "part4_ricemind_crosscheck.csv", encoding="utf-8")):
        if r["canonical_symbol"]:
            rap_to_symbol[r["RAP_ID"]] = r["canonical_symbol"]

    fig, ax = plt.subplots(figsize=(3.4, 2.85))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    gene_y = {n["node"]: 0.90 - i * (0.80 / max(len(genes) - 1, 1)) for i, n in enumerate(genes)}
    trait_y = {n["node"]: 0.90 - i * (0.80 / max(len(traits) - 1, 1)) for i, n in enumerate(traits)}
    # edges behind labels
    for e in edges:
        u, v = e["source"], e["target"]
        if u in gene_y and v in trait_y:
            ax.plot([0.08, 0.92], [gene_y[u], trait_y[v]], color=PALETTE["steel"],
                    linewidth=0.35 + 0.20 * float(e["weight"]), alpha=0.40, zorder=1)
        elif v in gene_y and u in trait_y:
            ax.plot([0.08, 0.92], [gene_y[v], trait_y[u]], color=PALETTE["steel"],
                    linewidth=0.35 + 0.20 * float(e["weight"]), alpha=0.40, zorder=1)
    for n in genes:
        y = gene_y[n["node"]]
        ax.scatter([0.06], [y], s=18, color=PALETTE["navy"], zorder=3)
        label = rap_to_symbol.get(n["node"], n["node"])
        ax.text(0.09, y, label, fontsize=5.4, ha="left", va="center", color=PALETTE["navy"],
                fontweight="bold", zorder=4,
                bbox=dict(boxstyle="round,pad=0.08", facecolor="white", edgecolor="none", alpha=0.85))
    for n in traits:
        y = trait_y[n["node"]]
        label = n["label"]
        if len(label) > 22:
            label = label[:20] + "…"
        ax.scatter([0.94], [y], s=18, color=PALETTE["green"], zorder=3)
        ax.text(0.90, y, label, fontsize=5.2, ha="right", va="center", color=PALETTE["green"],
                zorder=4,
                bbox=dict(boxstyle="round,pad=0.08", facecolor="white", edgecolor="none", alpha=0.85))
    ax.text(0.02, 1.0, "Gene", fontsize=5.8, color=PALETTE["navy"], fontweight="bold")
    ax.text(0.98, 1.0, "Trait", fontsize=5.8, color=PALETTE["green"], fontweight="bold", ha="right")
    save(fig, "fig7_gene_trait_network")


if __name__ == "__main__":
    print("[1/3] fig5 ...")
    fig5()
    print("[2/3] fig6 ...")
    fig6()
    print("[3/3] fig7 ...")
    fig7()
    print("done")