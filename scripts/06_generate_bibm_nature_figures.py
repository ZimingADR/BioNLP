# -*- coding: utf-8 -*-
"""
Step 6: 按照 Nature Figure 规范 + BIBM 经典深海蓝配色 (#1F65A6 / #0F4D92 / #3775BA / #99A9BF)
生成阶段性课程大论文的 4 张矢量 PDF 配图 (同时输出 .pdf, .svg, .png 供 LaTeX 编译与视觉质检)
"""

import json
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from nltk.tokenize import sent_tokenize, word_tokenize

# ── MANDATORY Nature-Figure font & vector export rules ────────────────────────
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"   # 保持 SVG 文本为可编辑 <text> 节点
plt.rcParams["pdf.fonttype"] = 42       # 保持 PDF 文本为 TrueType 可编辑字体 (Type 42)
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.linewidth"] = 0.9
plt.rcParams["xtick.major.width"] = 0.8
plt.rcParams["ytick.major.width"] = 0.8
plt.rcParams["legend.frameon"] = False
plt.rcParams["figure.dpi"] = 300

# ── BIBM Classic Deep Ocean Blue + Nature Neutral Palette ────────────────────
PALETTE_BIBM = {
    "navy_hero":    "#0F4D92",  # Nature / BIBM 主强调深海蓝
    "bibm_blue":    "#1F65A6",  # BIBM 模板 fig_example.pdf 标准柱状图蓝
    "sky_mid":      "#3775BA",  # 次级信号蓝
    "steel_soft":   "#6A96C8",  # 柔和过渡蓝
    "slate_grey":   "#99A9BF",  # BIBM 模板 fig_example.pdf 对照组蓝灰 ("Others" 同款色)
    "neutral_dark": "#272727",  # 文本与坐标轴深黑
    "neutral_mid":  "#6E7785",  # 辅助虚线/次要文字灰
    "accent_red":   "#B64342",  # Nature 强调红 (仅用于稀释下降趋势线或关键阈值)
}

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
FIG_DIR = RESULTS_DIR / "figures"
BIBM_FIG_DIR = BASE_DIR / "paper_bibm" / "BIBM-Latex" / "figure"
ARTIFACT_DIR = Path(r"C:\Users\ADR\.gemini\antigravity\brain\4cc9992a-5a00-4b89-a726-bd6d3a2c6b76")

FIG_DIR.mkdir(parents=True, exist_ok=True)
BIBM_FIG_DIR.mkdir(parents=True, exist_ok=True)


def add_panel_label(ax, label, x=-0.12, y=1.04, fontsize=11):
    """Nature-style bold lowercase panel label (a, b)."""
    ax.text(
        x, y, label,
        transform=ax.transAxes,
        fontsize=fontsize,
        fontweight="bold",
        color=PALETTE_BIBM["neutral_dark"],
        ha="left",
        va="bottom",
    )


def save_publication_figure(fig, stem_name: str):
    """Save figure as editable vector PDF + SVG + 300 DPI PNG to both BIBM and results folders."""
    fig.tight_layout(pad=1.2)
    for fmt in ["pdf", "svg", "png"]:
        out_bibm = BIBM_FIG_DIR / f"{stem_name}.{fmt}"
        out_res = FIG_DIR / f"{stem_name}.{fmt}"
        fig.savefig(out_bibm, format=fmt, dpi=300, bbox_inches="tight")
        shutil.copyfile(out_bibm, out_res)
        if fmt == "png" and ARTIFACT_DIR.exists():
            shutil.copyfile(out_bibm, ARTIFACT_DIR / f"{stem_name}.png")
    plt.close(fig)


def load_tokens():
    from importlib.machinery import SourceFileLoader
    step2 = SourceFileLoader("step2", str(BASE_DIR / "scripts" / "02_lexical_stats_ttr_tfidf.py")).load_module()
    docs = step2.parse_pubtator_documents()
    tokens = []
    for d in docs:
        for s in sent_tokenize(d["text"]):
            tokens.extend([t.lower() for t in word_tokenize(s) if any(c.isalnum() for c in t)])
    return tokens


def make_fig1_ttr_heaps(tokens, part1_summary):
    """
    Figure 1: Lexical scaling under Heaps' Law and multi-scale TTR standardization (3,306 articles).
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 2.85))
    fig.patch.set_facecolor("#FFFFFF")

    # Panel a: Heaps' Law V(N) growth and Raw TTR decay
    n_total = len(tokens)
    step = max(500, n_total // 55)
    ns, vs, ttrs = [], [], []
    seen = set()
    next_cp = step
    for idx, tok in enumerate(tokens, start=1):
        seen.add(tok)
        if idx == next_cp or idx == n_total:
            ns.append(idx / 1000.0)  # in thousands (k tokens)
            v_cnt = len(seen)
            vs.append(v_cnt / 1000.0)  # in thousands (k types)
            ttrs.append(v_cnt / idx * 100.0)
            next_cp += step

    ax1_twin = ax1.twinx()
    ax1_twin.spines["top"].set_visible(False)
    ax1_twin.spines["right"].set_visible(True)
    ax1_twin.spines["right"].set_color(PALETTE_BIBM["accent_red"])
    ax1_twin.spines["right"].set_linewidth(0.8)

    l1 = ax1.plot(
        ns, vs, color=PALETTE_BIBM["navy_hero"], lw=1.8,
        marker="o", markersize=2.6, markevery=4, label="Vocabulary size V (×10³)"
    )
    l2 = ax1_twin.plot(
        ns, ttrs, color=PALETTE_BIBM["accent_red"], lw=1.6, linestyle="--",
        marker="s", markersize=2.5, markevery=4, label="Raw TTR V/N (%)"
    )

    ax1.set_xlabel("Cumulative corpus size N (×10³ tokens)", fontsize=8.5)
    ax1.set_ylabel("Unique word types V (×10³)", fontsize=8.5, color=PALETTE_BIBM["navy_hero"])
    ax1_twin.set_ylabel("Raw TTR (%)", fontsize=8.5, color=PALETTE_BIBM["accent_red"])
    ax1.tick_params(axis="both", labelsize=7.5)
    ax1_twin.tick_params(axis="y", labelsize=7.5, colors=PALETTE_BIBM["accent_red"])
    ax1.set_title("Heaps' Law & TTR dilution (3,306 docs)", fontsize=8.8, fontweight="bold", pad=6)

    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc="center right", fontsize=7.0)
    add_panel_label(ax1, "a", x=-0.16, y=1.04)

    # Panel b: Fixed-window sequential truncation (head) vs global random sampling (shuf)
    samples = part1_summary["ttr_metrics"]["standardized_samples"]
    scale_keys = ["N=1000", "N=5000", "N=10000", "N=50000", "N=100000", "N=500000"]
    scale_labels = ["1k", "5k", "10k", "50k", "100k", "500k", "All\n(795k)"]
    head_vals = [samples[k]["head_ttr"] * 100 for k in scale_keys] + [part1_summary["ttr_metrics"]["raw_ttr"] * 100]
    shuf_vals = [samples[k]["shuf_mean_ttr"] * 100 for k in scale_keys] + [part1_summary["ttr_metrics"]["raw_ttr"] * 100]

    x = np.arange(len(scale_labels))
    w = 0.36
    b1 = ax2.bar(
        x - w / 2, head_vals, width=w,
        color=PALETTE_BIBM["slate_grey"], edgecolor="none",
        label="Sequential cut (head)"
    )
    b2 = ax2.bar(
        x + w / 2, shuf_vals, width=w,
        color=PALETTE_BIBM["bibm_blue"], edgecolor="none",
        label="Random sample (shuf)"
    )

    for b in b2:
        h = b.get_height()
        ax2.text(
            b.get_x() + b.get_width() / 2, h + 1.0,
            f"{h:.1f}", ha="center", va="bottom", fontsize=6.5,
            color=PALETTE_BIBM["navy_hero"], fontweight="bold"
        )

    ax2.set_xticks(x)
    ax2.set_xticklabels(scale_labels, fontsize=7.5)
    ax2.set_xlabel("Standardized token window size (N)", fontsize=8.5)
    ax2.set_ylabel("Type-Token Ratio (%)", fontsize=8.5)
    ax2.tick_params(axis="y", labelsize=7.5)
    ax2.set_ylim(0, 66)
    ax2.set_title("Multi-scale TTR standardization", fontsize=8.8, fontweight="bold", pad=6)
    ax2.legend(loc="upper right", fontsize=7.0)
    add_panel_label(ax2, "b", x=-0.15, y=1.04)

    save_publication_figure(fig, "fig1_ttr_heaps")


def make_fig2_tfidf_terms(part1_summary):
    """
    Figure 2: Top 12 domain unigrams and bigrams by cumulative TF-IDF across 3,306 articles.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.0))
    fig.patch.set_facecolor("#FFFFFF")

    uni = part1_summary["top_unigrams_tfidf"][:12][::-1]
    bi = part1_summary["top_bigrams_tfidf"][:12][::-1]

    # Panel a: Unigrams
    u_labels = [f"{x['term']} ({x['doc_freq']:,})" for x in uni]
    u_scores = [x["tfidf_score"] for x in uni]
    u_colors = [
        PALETTE_BIBM["navy_hero"] if x["term"] in {"drought", "stress", "genes", "tolerance"}
        else PALETTE_BIBM["bibm_blue"] if i >= 6 else PALETTE_BIBM["steel_soft"]
        for i, x in enumerate(uni)
    ]

    bars1 = ax1.barh(u_labels, u_scores, color=u_colors, height=0.68)
    max_u = max(u_scores)
    for b in bars1:
        w = b.get_width()
        ax1.text(w + max_u * 0.02, b.get_y() + b.get_height() / 2, f"{w:.1f}", va="center", fontsize=6.8)

    ax1.set_xlim(0, max_u * 1.20)
    ax1.set_xlabel("Cumulative TF-IDF weight", fontsize=8.5)
    ax1.tick_params(axis="both", labelsize=7.5)
    ax1.set_title("Top 12 unigrams (with doc freq)", fontsize=8.8, fontweight="bold", pad=6)
    add_panel_label(ax1, "a", x=-0.32, y=1.03)

    # Panel b: Bigrams (highlighting drought stress vs drought tolerance vs drought resistance)
    b_labels = [f"{x['term']} ({x['doc_freq']:,})" for x in bi]
    b_scores = [x["tfidf_score"] for x in bi]
    drought_variants = {"drought stress", "drought tolerance", "drought resistance"}
    b_colors = [
        PALETTE_BIBM["navy_hero"] if x["term"] in drought_variants
        else PALETTE_BIBM["bibm_blue"] if "stress" in x["term"] or "gene" in x["term"]
        else PALETTE_BIBM["slate_grey"]
        for x in bi
    ]

    bars2 = ax2.barh(b_labels, b_scores, color=b_colors, height=0.68)
    max_b = max(b_scores)
    for b in bars2:
        w = b.get_width()
        ax2.text(w + max_b * 0.02, b.get_y() + b.get_height() / 2, f"{w:.1f}", va="center", fontsize=6.8)

    ax2.set_xlim(0, max_b * 1.22)
    ax2.set_xlabel("Cumulative TF-IDF weight", fontsize=8.5)
    ax2.tick_params(axis="both", labelsize=7.5)
    ax2.set_title("Top 12 bigrams (query validation)", fontsize=8.8, fontweight="bold", pad=6)
    add_panel_label(ax2, "b", x=-0.42, y=1.03)

    save_publication_figure(fig, "fig2_tfidf_terms")


def make_fig3_pipeline_ablation(part2_report):
    """
    Figure 3: Entity recognition & normalization uplift over PubTator baseline and knowledge graph scale.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 2.85))
    fig.patch.set_facecolor("#FFFFFF")

    # Panel a: PubTator baseline vs Our Pipeline (Oryzabase + RTO)
    cats = ["Gene\nmentions", "RAP-DB\nnormalized", "RTO trait\nmentions", "Sentence\nevidence"]
    pub_vals = [
        part2_report.get("pubtator_gene_mentions", 1767),
        part2_report.get("pubtator_mapped_to_rap", 473),
        0,
        0,
    ]
    our_vals = [
        part2_report.get("total_gene_mentions", 12651),
        part2_report.get("rap_mapped_gene_mentions", 10708),
        part2_report.get("total_trait_mentions", 29966),
        part2_report.get("sentence_cooccurrence_evidence_count", 9030),
    ]

    x = np.arange(len(cats))
    w = 0.36
    b1 = ax1.bar(x - w / 2, pub_vals, width=w, color=PALETTE_BIBM["slate_grey"], label="PubTator only")
    b2 = ax1.bar(x + w / 2, our_vals, width=w, color=PALETTE_BIBM["bibm_blue"], label="Ours (Oryzabase+RTO)")

    max_v = max(our_vals)
    for b in b1:
        h = int(b.get_height())
        ax1.text(b.get_x() + b.get_width() / 2, h + max_v * 0.02, f"{h:,}", ha="center", va="bottom", fontsize=6.5, color=PALETTE_BIBM["neutral_mid"])
    for b in b2:
        h = int(b.get_height())
        ax1.text(b.get_x() + b.get_width() / 2, h + max_v * 0.02, f"{h:,}", ha="center", va="bottom", fontsize=6.8, color=PALETTE_BIBM["navy_hero"], fontweight="bold")

    ax1.set_xticks(x)
    ax1.set_xticklabels(cats, fontsize=7.5)
    ax1.set_ylabel("Extracted count (3,306 articles)", fontsize=8.5)
    ax1.tick_params(axis="y", labelsize=7.5)
    ax1.set_ylim(0, max_v * 1.18)
    ax1.set_title("Recall & normalization gain", fontsize=8.8, fontweight="bold", pad=6)
    ax1.legend(loc="upper left", fontsize=6.8)
    add_panel_label(ax1, "a", x=-0.18, y=1.04)

    # Panel b: Normalized Rice Drought Knowledge Graph Scale
    kg_cats = ["Unique\ngenes", "Unique RTO\ntraits", "Unique\nG-T pairs", "Evidence\nsentences"]
    kg_vals = [
        part2_report.get("unique_genes_count", 2557),
        part2_report.get("unique_traits_count", 388),
        part2_report.get("unique_gene_trait_pairs_count", 5575),
        part2_report.get("sentence_cooccurrence_evidence_count", 9030),
    ]
    kg_colors = [
        PALETTE_BIBM["navy_hero"],
        PALETTE_BIBM["sky_mid"],
        PALETTE_BIBM["bibm_blue"],
        PALETTE_BIBM["steel_soft"],
    ]

    bars_kg = ax2.bar(np.arange(len(kg_cats)), kg_vals, width=0.52, color=kg_colors)
    max_kg = max(kg_vals)
    for b in bars_kg:
        h = int(b.get_height())
        ax2.text(b.get_x() + b.get_width() / 2, h + max_kg * 0.02, f"{h:,}", ha="center", va="bottom", fontsize=7.0, fontweight="bold", color=PALETTE_BIBM["neutral_dark"])

    ax2.set_xticks(np.arange(len(kg_cats)))
    ax2.set_xticklabels(kg_cats, fontsize=7.5)
    ax2.set_ylabel("Normalized entity / edge count", fontsize=8.5)
    ax2.tick_params(axis="y", labelsize=7.5)
    ax2.set_ylim(0, max_kg * 1.18)
    ax2.set_title("Curated knowledge base scale", fontsize=8.8, fontweight="bold", pad=6)
    add_panel_label(ax2, "b", x=-0.18, y=1.04)

    save_publication_figure(fig, "fig3_pipeline_ablation")


def make_fig4_top_genes_traits(part2_report):
    """
    Figure 4: Top 14 normalized rice drought genes (with RAP-DB IDs) and Top 14 RTO traits.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.15))
    fig.patch.set_facecolor("#FFFFFF")

    top_genes = part2_report["top_genes"][:14][::-1]
    g_labels = [f"{g['canonical_symbol']} ({g['rap_id']})" for g in top_genes]
    g_counts = [g["count"] for g in top_genes]
    g_colors = [
        PALETTE_BIBM["navy_hero"] if i >= len(top_genes) - 5 else PALETTE_BIBM["bibm_blue"]
        for i in range(len(top_genes))
    ]

    bars_g = ax1.barh(g_labels, g_counts, color=g_colors, height=0.68)
    max_g = max(g_counts)
    for b in bars_g:
        w = int(b.get_width())
        ax1.text(w + max_g * 0.02, b.get_y() + b.get_height() / 2, f"{w}", va="center", fontsize=6.8, fontweight="bold")

    ax1.set_xlim(0, max_g * 1.18)
    ax1.set_xlabel("Sentence-level mentions", fontsize=8.5)
    ax1.tick_params(axis="both", labelsize=7.0)
    ax1.set_title("Top 14 normalized rice genes", fontsize=8.8, fontweight="bold", pad=6)
    add_panel_label(ax1, "a", x=-0.45, y=1.03)

    top_traits = part2_report["top_traits"][:14][::-1]
    t_labels = [f"{t['canonical_trait']} ({t['rto_id']})" for t in top_traits]
    t_counts = [t["count"] for t in top_traits]
    t_colors = [
        PALETTE_BIBM["navy_hero"] if i >= len(top_traits) - 3
        else PALETTE_BIBM["sky_mid"] if i >= len(top_traits) - 8
        else PALETTE_BIBM["slate_grey"]
        for i in range(len(top_traits))
    ]

    bars_t = ax2.barh(t_labels, t_counts, color=t_colors, height=0.68)
    max_t = max(t_counts)
    for b in bars_t:
        w = int(b.get_width())
        ax2.text(w + max_t * 0.02, b.get_y() + b.get_height() / 2, f"{w:,}", va="center", fontsize=6.8)

    ax2.set_xlim(0, max_t * 1.22)
    ax2.set_xlabel("Sentence-level mentions", fontsize=8.5)
    ax2.tick_params(axis="both", labelsize=7.0)
    ax2.set_title("Top 14 normalized RTO traits", fontsize=8.8, fontweight="bold", pad=6)
    add_panel_label(ax2, "b", x=-0.48, y=1.03)

    save_publication_figure(fig, "fig4_top_genes_traits")


def main():
    part1_summary = json.loads((RESULTS_DIR / "part1_lexical_and_tfidf_summary.json").read_text(encoding="utf-8"))
    part2_report = json.loads((RESULTS_DIR / "part2_normalization_and_evidence_report.json").read_text(encoding="utf-8"))
    tokens = load_tokens()

    make_fig1_ttr_heaps(tokens, part1_summary)
    make_fig2_tfidf_terms(part1_summary)
    make_fig3_pipeline_ablation(part2_report)
    make_fig4_top_genes_traits(part2_report)
    print("Done! All 4 BIBM + Nature-style vector PDF/SVG/PNG figures have been exported.")


if __name__ == "__main__":
    main()
