# -*- coding: utf-8 -*-
"""
Regenerate the 6 data-driven publication figures for paper_final/BIBM-Latex/figure/
directly from workspace/results/ using a cohesive high-contrast viridis palette.
"""

import csv
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from nltk.tokenize import sent_tokenize, word_tokenize

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"
OUT_FIG_DIR = BASE_DIR / "paper_final" / "BIBM-Latex" / "figure"
PREVIEW_DIR = Path(r"C:\Users\ADR\.gemini\antigravity\brain\c53bbf5c-0823-4df1-a425-bffaf6ab54e1\scratch\fig_previews_viridis")
OUT_FIG_DIR.mkdir(parents=True, exist_ok=True)
PREVIEW_DIR.mkdir(parents=True, exist_ok=True)

# High-contrast Viridis stops legible on white background
PALETTE = {
    "v_dark": "#440154",    # viridis 0.00 - primary hero / strong highlights
    "v_indigo": "#3B528B",  # viridis 0.25 - secondary primary
    "v_teal": "#21918C",    # viridis 0.50 - mid-tier / positive / traits
    "v_green": "#5EC962",   # viridis 0.75 - accent / soft highlight
    "v_lime": "#9FDA3A",    # viridis 0.85 - bright fill (use dark text)
    "muted": "#A8B0BC",     # neutral grey for baseline comparisons
    "dark": "#222222",
    "mid": "#5C6673",
    "pale": "#F2EEF6",
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


def save_fig(fig, stem: str):
    pdf_path = OUT_FIG_DIR / f"{stem}.pdf"
    png_path = PREVIEW_DIR / f"{stem}.png"
    fig.savefig(pdf_path, format="pdf", dpi=300, bbox_inches="tight", pad_inches=0.05)
    fig.savefig(png_path, format="png", dpi=200, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"  saved {pdf_path.name}")


def panel_tag(ax, letter: str, x=-0.12, y=1.04, fs=11):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=fs, fontweight="bold",
            color=PALETTE["dark"], ha="left", va="bottom")


def load_tokens():
    step2 = SourceFileLoader("step2", str(BASE_DIR / "scripts" / "02_lexical_stats_ttr_tfidf.py")).load_module()
    docs = step2.parse_pubtator_documents()
    tokens = []
    for d in docs:
        for s in sent_tokenize(d["text"]):
            tokens.extend([t.lower() for t in word_tokenize(s) if any(c.isalnum() for c in t)])
    return tokens


def make_fig1_lexical_tfidf(tokens, part1):
    """Fig. 1: 3-panel lexical scaling + multi-scale TTR + Top 10 TF-IDF bigrams (viridis)."""
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(7.2, 2.45),
                                        gridspec_kw={"width_ratios": [1.0, 1.0, 1.0], "wspace": 0.85})
    fig.patch.set_facecolor("#FFFFFF")

    # (a) Heaps' Law
    n_total = len(tokens)
    step = max(500, n_total // 55)
    ns, vs, ttrs = [], [], []
    seen = set()
    next_cp = step
    for idx, tok in enumerate(tokens, start=1):
        seen.add(tok)
        if idx == next_cp or idx == n_total:
            ns.append(idx / 1000.0)
            v_cnt = len(seen)
            vs.append(v_cnt / 1000.0)
            ttrs.append(v_cnt / idx * 100.0)
            next_cp += step

    ax1_t = ax1.twinx()
    ax1_t.spines["top"].set_visible(False)
    ax1_t.spines["right"].set_visible(True)
    ax1_t.spines["right"].set_color(PALETTE["v_teal"])
    l1 = ax1.plot(ns, vs, color=PALETTE["v_dark"], lw=1.6, marker="o", ms=2.2, markevery=6, label="Types V (k)")
    l2 = ax1_t.plot(ns, ttrs, color=PALETTE["v_teal"], lw=1.5, ls="--", marker="s", ms=2.0, markevery=6, label="TTR (%)")
    ax1.set_xlabel("Tokens N (x10^3)", fontsize=7.0)
    ax1.set_ylabel("Types V (x10^3)", fontsize=7.0, color=PALETTE["v_dark"], labelpad=2)
    ax1_t.set_ylabel("Raw TTR (%)", fontsize=7.0, color=PALETTE["v_teal"], labelpad=2)
    ax1.tick_params(labelsize=6.2)
    ax1_t.tick_params(labelsize=6.2, colors=PALETTE["v_teal"])
    ax1.set_title("Heaps' Law & TTR", fontsize=7.6, fontweight="bold", pad=5)
    lines = l1 + l2
    ax1.legend(lines, [l.get_label() for l in lines], loc="center right", fontsize=6.0)
    panel_tag(ax1, "a", x=-0.26, y=1.05, fs=10)

    # (b) Multi-scale TTR
    samples = part1["ttr_metrics"]["standardized_samples"]
    keys = ["N=1000", "N=5000", "N=10000", "N=50000", "N=100000", "N=500000"]
    xlabels = ["1k", "5k", "10k", "50k", "100k", "500k", "All"]
    head_v = [samples[k]["head_ttr"] * 100 for k in keys] + [part1["ttr_metrics"]["raw_ttr"] * 100]
    shuf_v = [samples[k]["shuf_mean_ttr"] * 100 for k in keys] + [part1["ttr_metrics"]["raw_ttr"] * 100]
    x = np.arange(len(xlabels))
    w = 0.37
    ax2.bar(x - w / 2, head_v, width=w, color=PALETTE["muted"], label="head")
    b2 = ax2.bar(x + w / 2, shuf_v, width=w, color=PALETTE["v_indigo"], label="shuf")
    for b in b2:
        h = b.get_height()
        ax2.text(b.get_x() + b.get_width() / 2, h + 1.0, f"{h:.1f}", ha="center", va="bottom",
                 fontsize=5.3, color=PALETTE["v_dark"], fontweight="bold")
    ax2.set_xticks(x)
    ax2.set_xticklabels(xlabels, fontsize=6.0)
    ax2.set_xlabel("Window N", fontsize=7.0)
    ax2.set_ylabel("TTR (%)", fontsize=7.0, labelpad=2)
    ax2.tick_params(axis="y", labelsize=6.2)
    ax2.set_ylim(0, 66)
    ax2.set_title("Multi-scale TTR", fontsize=7.6, fontweight="bold", pad=5)
    ax2.legend(loc="upper right", fontsize=6.0)
    panel_tag(ax2, "b", x=-0.24, y=1.05, fs=10)

    # (c) Top 10 TF-IDF bigrams
    bi = part1["top_bigrams_tfidf"][:10][::-1]
    blabels = [f"{t['term']} ({t['doc_freq']:,})" for t in bi]
    bscores = [t["tfidf_score"] for t in bi]
    hi_terms = {"drought stress", "drought tolerance", "drought resistance"}
    bcols = [PALETTE["v_dark"] if t["term"] in hi_terms else PALETTE["v_teal"] for t in bi]
    bars = ax3.barh(blabels, bscores, color=bcols, height=0.68)
    max_b = max(bscores)
    for b in bars:
        val = b.get_width()
        ax3.text(val + max_b * 0.02, b.get_y() + b.get_height() / 2, f"{val:.1f}",
                 va="center", fontsize=5.8, color=PALETTE["dark"])
    ax3.set_xlim(0, max_b * 1.24)
    ax3.set_xlabel("TF-IDF weight", fontsize=7.0)
    ax3.tick_params(labelsize=5.9)
    ax3.set_title("Top 10 bigrams (DF)", fontsize=7.6, fontweight="bold", pad=5)
    panel_tag(ax3, "c", x=-0.52, y=1.05, fs=10)

    save_fig(fig, "fig1_lexical_tfidf")


def make_fig3_ablation(part2, key_nums):
    """Fig. 2: PubTator vs Dual-Dict uplift (a) and Rice-only curated scale (b) (viridis)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.55), gridspec_kw={"wspace": 0.28})
    fig.patch.set_facecolor("#FFFFFF")

    cats = ["Gene\nmentions", "RAP-DB\nmapped", "RTO trait\nmentions", "Rice-only\nco-occur."]
    pub_vals = [
        part2["pubtator_gene_mentions"],
        part2["pubtator_mapped_to_rap"],
        0,
        0,
    ]
    our_vals = [
        part2["total_gene_mentions"],
        part2["rap_mapped_gene_mentions"],
        part2["total_trait_mentions"],
        key_nums["entities"]["rice_candidate_pairs"],
    ]
    x = np.arange(len(cats))
    w = 0.36
    b1 = ax1.bar(x - w / 2, pub_vals, width=w, color=PALETTE["muted"], label="PubTator only")
    b2 = ax1.bar(x + w / 2, our_vals, width=w, color=PALETTE["v_indigo"], label="Ours (Oryzabase+RTO)")
    max_v = max(our_vals)
    for b in b1:
        h = int(b.get_height())
        ax1.text(b.get_x() + b.get_width() / 2, h + max_v * 0.02, f"{h:,}",
                 ha="center", va="bottom", fontsize=7.0, color=PALETTE["mid"])
    for b in b2:
        h = int(b.get_height())
        ax1.text(b.get_x() + b.get_width() / 2, h + max_v * 0.02, f"{h:,}",
                 ha="center", va="bottom", fontsize=7.2, color=PALETTE["v_dark"], fontweight="bold")
    ax1.annotate("22.6x", xy=(1 + w / 2, 10708), xytext=(1.35, 17500),
                 arrowprops=dict(arrowstyle="->", color=PALETTE["v_teal"], lw=1.2),
                 fontsize=7.6, fontweight="bold", color=PALETTE["v_teal"])
    ax1.set_xticks(x)
    ax1.set_xticklabels(cats, fontsize=7.6)
    ax1.set_ylabel("Count (3,306 docs)", fontsize=8.0)
    ax1.tick_params(axis="y", labelsize=7.4)
    ax1.set_ylim(0, max_v * 1.18)
    ax1.set_title("Recall & RAP normalization gain", fontsize=8.5, fontweight="bold", pad=6)
    ax1.legend(loc="upper left", fontsize=7.2)
    panel_tag(ax1, "a", x=-0.16, y=1.04, fs=10.5)

    kg_cats = ["Unique\nrice genes", "Unique RTO\ntraits", "Unique\nG-T pairs", "Rice-only\ncandidates"]
    kg_vals = [
        key_nums["entities"]["rice_unique_genes"],       # 1,344
        key_nums["entities"]["rice_unique_traits"],      # 190
        key_nums["entities"]["rice_unique_pairs"],       # 4,678
        key_nums["entities"]["rice_candidate_pairs"],    # 7,868
    ]
    kg_cols = [PALETTE["v_dark"], PALETTE["v_indigo"], PALETTE["v_teal"], PALETTE["v_green"]]
    bars_kg = ax2.bar(np.arange(len(kg_cats)), kg_vals, width=0.52, color=kg_cols)
    max_kg = max(kg_vals)
    for b in bars_kg:
        h = int(b.get_height())
        ax2.text(b.get_x() + b.get_width() / 2, h + max_kg * 0.02, f"{h:,}",
                 ha="center", va="bottom", fontsize=7.5, fontweight="bold", color=PALETTE["dark"])
    ax2.text(0.03, 0.96, "Pre-filter total: 9,030\n(1,162 non-rice removed)",
             transform=ax2.transAxes, ha="left", va="top", fontsize=6.6, color=PALETTE["v_dark"],
             bbox=dict(boxstyle="round,pad=0.25", facecolor=PALETTE["pale"], edgecolor=PALETTE["v_indigo"], lw=0.6))
    ax2.set_xticks(np.arange(len(kg_cats)))
    ax2.set_xticklabels(kg_cats, fontsize=7.6)
    ax2.set_ylabel("Rice-only entity / pair count", fontsize=8.0)
    ax2.tick_params(axis="y", labelsize=7.4)
    ax2.set_ylim(0, max_kg * 1.20)
    ax2.set_title("Rice-only co-occurrence scale (Os*)", fontsize=8.5, fontweight="bold", pad=6)
    panel_tag(ax2, "b", x=-0.16, y=1.04, fs=10.5)

    save_fig(fig, "fig3_pipeline_ablation")


def make_fig4_top_entities(part2):
    """Fig. 3: Top 14 normalized rice genes & Top 14 RTO traits (viridis)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.05), gridspec_kw={"wspace": 0.90})
    fig.patch.set_facecolor("#FFFFFF")

    top_g = part2["top_genes"][:14][::-1]
    g_labels = [f"{g['canonical_symbol']} ({g['rap_id']})" for g in top_g]
    g_cnts = [g["count"] for g in top_g]
    g_cols = [PALETTE["v_dark"] if i >= len(top_g) - 5 else PALETTE["v_indigo"] for i in range(len(top_g))]
    bg = ax1.barh(g_labels, g_cnts, color=g_cols, height=0.68)
    mg = max(g_cnts)
    for b in bg:
        w = int(b.get_width())
        ax1.text(w + mg * 0.02, b.get_y() + b.get_height() / 2, f"{w}",
                 va="center", fontsize=6.6, fontweight="bold", color=PALETTE["dark"])
    ax1.set_xlim(0, mg * 1.22)
    ax1.set_xlabel("Sentence-level mentions", fontsize=7.6)
    ax1.tick_params(labelsize=6.6)
    ax1.set_title("Top 14 normalized rice genes", fontsize=8.2, fontweight="bold", pad=6)
    panel_tag(ax1, "a", x=-0.52, y=1.03, fs=10.5)

    top_t = part2["top_traits"][:14][::-1]
    t_labels = [f"{t['canonical_trait']} ({t['rto_id']})" for t in top_t]
    t_cnts = [t["count"] for t in top_t]
    t_cols = [PALETTE["v_dark"] if i >= len(top_t) - 3 else PALETTE["v_teal"] if i >= len(top_t) - 8 else PALETTE["v_green"]
              for i in range(len(top_t))]
    bt = ax2.barh(t_labels, t_cnts, color=t_cols, height=0.68)
    mt = max(t_cnts)
    for b in bt:
        w = int(b.get_width())
        ax2.text(w + mt * 0.02, b.get_y() + b.get_height() / 2, f"{w:,}",
                 va="center", fontsize=6.6, color=PALETTE["dark"])
    ax2.set_xlim(0, mt * 1.24)
    ax2.set_xlabel("Sentence-level mentions", fontsize=7.6)
    ax2.tick_params(labelsize=6.6)
    ax2.set_title("Top 14 normalized RTO traits", fontsize=8.2, fontweight="bold", pad=6)
    panel_tag(ax2, "b", x=-0.56, y=1.03, fs=10.5)

    save_fig(fig, "fig4_top_genes_traits")


def make_fig5_sdp(sdp):
    """Fig. 4: SDP example, filtering funnel (7,868), and polarity breakdown (884) (viridis)."""
    rel = sdp["relation_counts"]
    fig = plt.figure(figsize=(7.0, 2.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.08, 1.02, 0.72], wspace=0.45)

    # (a) LOF SDP chain
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    nodes = [
        (0.32, 0.88, "osbZIP23", "gene / pobj", PALETTE["v_indigo"]),
        (0.32, 0.66, "mutations", "nsubj (LOF)", PALETTE["v_teal"]),
        (0.32, 0.44, "weakened", "ccomp (NEG)", PALETTE["v_dark"]),
        (0.32, 0.22, "drought tolerance", "dobj / trait", PALETTE["v_green"]),
    ]
    for x, y, lab, sub, col in nodes:
        p = FancyBboxPatch((x - 0.28, y - 0.065), 0.56, 0.13,
                           boxstyle="round,pad=0.012,rounding_size=0.02",
                           linewidth=1.1, edgecolor=col, facecolor=col + "1A")
        ax.add_patch(p)
        ax.text(x, y, lab, ha="center", va="center", fontsize=7.4, color=PALETTE["dark"])
        ax.text(0.64, y, sub, fontsize=6.8, va="center", color=col, style="italic")
    for y0, y1, dep in [(0.805, 0.735, "pobj"), (0.585, 0.515, "nsubj"), (0.365, 0.295, "dobj")]:
        ax.add_patch(FancyArrowPatch((0.32, y0), (0.32, y1), arrowstyle="-|>",
                                     mutation_scale=9, linewidth=1.1, color=PALETTE["mid"]))
        ax.text(0.35, (y0 + y1) / 2, dep, fontsize=6.6, va="center", color=PALETTE["mid"], style="italic")
    ax.text(0.50, 0.045, "LOSS x NEG -> gene effect POSITIVE",
            ha="center", fontsize=7.0, color=PALETTE["v_dark"], fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.22", facecolor=PALETTE["pale"], edgecolor="none"))
    panel_tag(ax, "a", x=-0.10, y=1.02, fs=10.5)

    # (b) Filtering funnel
    ax = fig.add_subplot(gs[0, 1])
    stages = ["Candidates", "Causal", "Association", "Filtered"]
    causal_tot = rel["positive_regulation"] + rel["negative_regulation"] + rel["regulation"]
    filt_tot = rel["excluded_direction"] + rel["excluded_respectively"] + rel["no_effect"] + rel["unsupported"]
    vals = [7868, causal_tot, rel["association"], filt_tot]
    cols = [PALETTE["muted"], PALETTE["v_dark"], PALETTE["v_teal"], PALETTE["v_indigo"]]
    bars = ax.barh(range(4)[::-1], vals, color=cols, height=0.58)
    ax.set_yticks(range(4)[::-1])
    ax.set_yticklabels(stages, fontsize=7.4)
    ax.set_xscale("log")
    ax.set_xlim(150, 22000)
    ax.set_xlabel("Candidate pairs (log)", fontsize=7.5)
    ax.tick_params(axis="x", labelsize=7.0)
    for b, v in zip(bars, vals):
        ax.text(v * 1.15, b.get_y() + b.get_height() / 2, f"{v:,}", va="center", fontsize=7.0, color=PALETTE["dark"])
    panel_tag(ax, "b", x=-0.30, y=1.02, fs=10.5)

    # (c) Polarity
    ax = fig.add_subplot(gs[0, 2])
    pol = [rel["positive_regulation"], rel["negative_regulation"], rel["regulation"]]
    labs = ["Pos.", "Neg.", "Neutral"]
    bars = ax.bar(labs, pol, color=[PALETTE["v_teal"], PALETTE["v_dark"], PALETTE["v_green"]], width=0.60)
    for b, v in zip(bars, pol):
        ax.text(b.get_x() + b.get_width() / 2, v + 12, f"{v}", ha="center", fontsize=7.2, fontweight="bold")
    ax.set_ylim(0, max(pol) * 1.22)
    ax.set_ylabel("Causal relations", fontsize=7.5)
    ax.tick_params(labelsize=7.2)
    panel_tag(ax, "c", x=-0.34, y=1.02, fs=10.5)

    fig.suptitle("SDP relation extraction and polarity-aware filtering", fontsize=8.8, fontweight="bold", y=1.03)
    save_fig(fig, "fig5_sdp_framework")


def make_fig6_eval_enrichment(audit, enr):
    """Fig. 5: Audit metrics + TO/RTO enrichment bubble plot (viridis)."""
    fig = plt.figure(figsize=(7.0, 2.7))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.95, 1.15], wspace=0.50)

    # (a) Audit metrics
    ax = fig.add_subplot(gs[0, 0])
    metrics = ["Precision", "Recall", "F1", "Polarity\nAcc."]
    p_base = audit["cooccurrence_baseline_precision_weighted"]
    r_base = 1.0
    f1_base = 2 * p_base * r_base / (p_base + r_base)  # 0.5132 -> 0.51
    baseline = [p_base, r_base, f1_base, 0.0]
    sdp_vals = [
        audit["sdp_precision"],
        audit["sdp_recall_estimate"],
        audit["sdp_f1_estimate"],
        audit["polarity_accuracy"],
    ]
    x = np.arange(len(metrics))
    w = 0.35
    b1 = ax.bar(x - w / 2, baseline, w, color=PALETTE["muted"], label="Co-occurrence")
    b2 = ax.bar(x + w / 2, sdp_vals, w, color=PALETTE["v_dark"], label="SDP")
    for i, b in enumerate(b1):
        h = b.get_height()
        if i == 3:
            ax.text(b.get_x() + b.get_width() / 2, 0.03, "N/A", ha="center", fontsize=6.6, color=PALETTE["mid"])
        else:
            ax.text(b.get_x() + b.get_width() / 2, h + 0.025, f"{h:.2f}", ha="center", fontsize=6.8, color=PALETTE["mid"])
    for i, b in enumerate(b2):
        h = b.get_height()
        ytxt = h + 0.10 if i == 0 else h + 0.025
        ax.text(b.get_x() + b.get_width() / 2, ytxt, f"{h:.2f}", ha="center", fontsize=7.0,
                fontweight="bold", color=PALETTE["v_dark"])
    lo, hi = audit["sdp_precision_95ci"]
    ax.errorbar(x[0] + w / 2, sdp_vals[0], yerr=[[sdp_vals[0] - lo], [hi - sdp_vals[0]]],
                fmt="none", ecolor=PALETTE["v_teal"], capsize=2.8, linewidth=1.1)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics, fontsize=7.3)
    ax.set_ylim(0, 1.24)
    ax.set_ylabel("Score", fontsize=7.8)
    ax.tick_params(axis="y", labelsize=7.2)
    ax.legend(fontsize=7.0, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.03))
    panel_tag(ax, "a", x=-0.18, y=1.04, fs=10.5)

    # (b) Enrichment bubble with viridis colormap
    ax = fig.add_subplot(gs[0, 1])
    terms = enr["top_significant"]["top100"][:12][::-1]
    y = np.arange(len(terms))
    folds = [t["fold_enrichment"] for t in terms]
    overlaps = [t["x_overlap"] for t in terms]
    qvals = [max(t["Adjusted P-value"], 1e-300) for t in terms]
    colors = -np.log10(qvals)
    sizes = np.array(overlaps) * 3.8 + 18
    sc = ax.scatter(folds, y, s=sizes, c=colors, cmap="viridis", edgecolors=PALETTE["dark"],
                    linewidths=0.55, alpha=0.92, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([t["Term"].split(" (")[0] for t in terms], fontsize=6.9)
    ax.set_ylim(-0.9, len(terms) + 0.6)
    ax.set_xlabel("Fold enrichment", fontsize=7.8)
    ax.set_xlim(0, max(folds) * 1.48)
    ax.tick_params(axis="x", labelsize=7.2)
    for i, (f, o) in enumerate(zip(folds, overlaps)):
        offset = 1.65 if o >= 60 else 1.05
        ax.text(f + offset, i, f"x={o}", va="center", fontsize=6.4, color=PALETTE["dark"])
    cbar = fig.colorbar(sc, ax=ax, fraction=0.036, pad=0.02)
    cbar.set_label(r"$-\log_{10}q$", fontsize=7.0, labelpad=2)
    cbar.ax.tick_params(labelsize=6.5)
    panel_tag(ax, "b", x=-0.34, y=1.04, fs=10.5)

    fig.suptitle("Audited SDP performance and TO/RTO functional enrichment", fontsize=8.8, fontweight="bold", y=1.04)
    save_fig(fig, "fig6_evaluation_enrichment")


def make_fig7_network():
    """Fig. 6: Legible Top 10 genes x Top 8 traits bipartite network (viridis)."""
    nodes = list(csv.DictReader(open(RESULTS_DIR / "part4_network_nodes.csv", encoding="utf-8")))
    edges = list(csv.DictReader(open(RESULTS_DIR / "part4_network_edges.csv", encoding="utf-8")))
    genes = [n for n in nodes if n["type"] == "gene"][:10]
    traits = [n for n in nodes if n["type"] == "trait"][:8]
    rap_to_symbol = {}
    for r in csv.DictReader(open(RESULTS_DIR / "part4_ricemind_crosscheck.csv", encoding="utf-8")):
        if r["canonical_symbol"]:
            rap_to_symbol[r["RAP_ID"]] = r["canonical_symbol"]

    fig, ax = plt.subplots(figsize=(3.45, 2.85))
    fig.patch.set_facecolor("#FFFFFF")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.06, 1.03)
    ax.axis("off")

    gene_y = {n["node"]: 0.89 - i * (0.78 / max(len(genes) - 1, 1)) for i, n in enumerate(genes)}
    trait_y = {n["node"]: 0.89 - i * (0.78 / max(len(traits) - 1, 1)) for i, n in enumerate(traits)}

    for e in edges:
        u, v, w = e["source"], e["target"], float(e["weight"])
        if u in gene_y and v in trait_y:
            ax.plot([0.08, 0.92], [gene_y[u], trait_y[v]], color=PALETTE["v_indigo"],
                    linewidth=0.45 + 0.30 * w, alpha=0.38, zorder=1)
        elif v in gene_y and u in trait_y:
            ax.plot([0.08, 0.92], [gene_y[v], trait_y[u]], color=PALETTE["v_indigo"],
                    linewidth=0.45 + 0.30 * w, alpha=0.38, zorder=1)

    for n in genes:
        y = gene_y[n["node"]]
        ax.scatter([0.06], [y], s=26, color=PALETTE["v_dark"], zorder=3)
        lab = rap_to_symbol.get(n["node"], n["node"])
        ax.text(0.10, y, lab, fontsize=7.0, ha="left", va="center", color=PALETTE["v_dark"],
                fontweight="bold", zorder=4,
                bbox=dict(boxstyle="round,pad=0.10", facecolor="white", edgecolor="none", alpha=0.88))

    for n in traits:
        y = trait_y[n["node"]]
        lab = n["label"]
        if len(lab) > 21:
            lab = lab[:19] + "..."
        ax.scatter([0.94], [y], s=26, color=PALETTE["v_teal"], zorder=3)
        ax.text(0.90, y, lab, fontsize=6.8, ha="right", va="center", color=PALETTE["v_teal"],
                fontweight="bold", zorder=4,
                bbox=dict(boxstyle="round,pad=0.10", facecolor="white", edgecolor="none", alpha=0.88))

    ax.text(0.02, 0.98, "Top 10 Genes", fontsize=7.6, color=PALETTE["v_dark"], fontweight="bold")
    ax.text(0.98, 0.98, "Top 8 Traits", fontsize=7.6, color=PALETTE["v_teal"], fontweight="bold", ha="right")

    leg_handles = [
        Line2D([0], [0], color=PALETTE["v_indigo"], lw=0.75, alpha=0.6, label="1 sent."),
        Line2D([0], [0], color=PALETTE["v_indigo"], lw=1.35, alpha=0.6, label="3 sent."),
        Line2D([0], [0], color=PALETTE["v_indigo"], lw=1.95, alpha=0.6, label="5+ sent."),
    ]
    ax.legend(handles=leg_handles, loc="lower center", ncol=3, fontsize=6.4,
              handlelength=1.6, columnspacing=1.0, bbox_to_anchor=(0.5, -0.04))

    save_fig(fig, "fig7_gene_trait_network")


def main():
    part1 = json.loads((RESULTS_DIR / "part1_lexical_and_tfidf_summary.json").read_text(encoding="utf-8"))
    part2 = json.loads((RESULTS_DIR / "part2_normalization_and_evidence_report.json").read_text(encoding="utf-8"))
    sdp = json.loads((RESULTS_DIR / "sdp_summary.json").read_text(encoding="utf-8"))
    audit = json.loads((RESULTS_DIR / "sdp_audit_summary.json").read_text(encoding="utf-8"))
    enr = json.loads((RESULTS_DIR / "part4_enrichment_summary.json").read_text(encoding="utf-8"))
    key_nums = json.loads((RESULTS_DIR / "key_numbers.json").read_text(encoding="utf-8"))
    tokens = load_tokens()

    # Clean up removed pipeline figure if present
    old_fig0 = OUT_FIG_DIR / "fig0_pipeline_overview.pdf"
    if old_fig0.exists():
        old_fig0.unlink()

    print("Generating 6 viridis figures...")
    make_fig1_lexical_tfidf(tokens, part1)
    make_fig3_ablation(part2, key_nums)
    make_fig4_top_entities(part2)
    make_fig5_sdp(sdp)
    make_fig6_eval_enrichment(audit, enr)
    make_fig7_network()
    print("All 6 viridis figures regenerated.")


if __name__ == "__main__":
    main()
