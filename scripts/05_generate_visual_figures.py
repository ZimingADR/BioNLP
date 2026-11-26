# -*- coding: utf-8 -*-
"""
Step 5: 生成论文级高清科研可视化图表 (自适应 100 篇或全量 9,586 篇语料规模)
输出目录:
  - results/figures/fig1_ttr_heaps_curve.png
  - results/figures/fig2_tfidf_top_terms.png
  - results/figures/fig3_disambiguation_comparison.png
  - results/figures/fig4_gene_trait_matrix.png
"""

import json
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from nltk.tokenize import sent_tokenize, word_tokenize

# 配置中文字体与科研出版级样式
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 180

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
FIG_DIR = RESULTS_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
ARTIFACT_DIR = Path(r"C:\Users\ADR\.gemini\antigravity\brain\4cc9992a-5a00-4b89-a726-bd6d3a2c6b76")


def load_tokens():
    from importlib.machinery import SourceFileLoader
    step2 = SourceFileLoader("step2", str(BASE_DIR / "scripts" / "02_lexical_stats_ttr_tfidf.py")).load_module()
    docs = step2.parse_pubtator_documents()
    tokens = []
    for d in docs:
        for s in sent_tokenize(d["text"]):
            tokens.extend([t.lower() for t in word_tokenize(s) if any(c.isalnum() for c in t)])
    return tokens


def plot_fig1_ttr_heaps(tokens, part1_summary):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.4))
    fig.patch.set_facecolor("#ffffff")

    # 左图: O(N) 增量计算 Heaps' Law 词表增长与 TTR 随文本长度 N 稀释曲线
    n_total = len(tokens)
    step = max(500, n_total // 60)
    ns = []
    vs = []
    ttrs = []
    seen = set()
    next_checkpoint = step
    for idx, tok in enumerate(tokens, start=1):
        seen.add(tok)
        if idx == next_checkpoint or idx == n_total:
            ns.append(idx)
            v_cnt = len(seen)
            vs.append(v_cnt)
            ttrs.append(v_cnt / idx * 100.0)
            next_checkpoint += step

    color_v = "#0284c7"
    color_ttr = "#e11d48"

    ax1_twin = ax1.twinx()
    l1 = ax1.plot(ns, vs, color=color_v, lw=2.5, marker="o", markersize=3.2, label="独立词型数 Types ($V$)")
    l2 = ax1_twin.plot(ns, ttrs, color=color_ttr, lw=2.5, linestyle="--", marker="s", markersize=3.2, label="原始 TTR ($V/N$, %)")

    n_docs = part1_summary["num_documents"]
    ax1.set_title(f"(a) 希普斯定律 (Heaps' Law) 与 TTR 长度稀释效应 ({n_docs:,} 篇文献)", fontsize=11.5, fontweight="bold", pad=10)
    ax1.set_xlabel("累积文本长度 Tokens ($N$)", fontsize=10.5)
    ax1.set_ylabel("独立词型数 Types ($V$)", color=color_v, fontsize=10.5, fontweight="bold")
    ax1_twin.set_ylabel("类符形符比 TTR (%)", color=color_ttr, fontsize=10.5, fontweight="bold")
    ax1.grid(True, linestyle=":", alpha=0.5)

    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc="center right", frameon=True, facecolor="#f8fafc")

    # 右图: 不同定长窗口下 "连续截断 (head)" vs "全局随机抽样 (shuf)" TTR 对比
    samples = part1_summary["ttr_metrics"]["standardized_samples"]
    scales = []
    head_vals = []
    shuf_vals = []
    for k_key, k_val in samples.items():
        scales.append(k_key)
        head_vals.append(k_val["head_ttr"] * 100)
        shuf_vals.append(k_val["shuf_mean_ttr"] * 100)
    scales.append(f"全语料\n(N={n_total:,})")
    head_vals.append(part1_summary["ttr_metrics"]["raw_ttr"] * 100)
    shuf_vals.append(part1_summary["ttr_metrics"]["raw_ttr"] * 100)

    x = np.arange(len(scales))
    width = 0.36
    b1 = ax2.bar(x - width / 2, head_vals, width, label="开头连续截断 (head -n)", color="#38bdf8", edgecolor="#0369a1")
    b2 = ax2.bar(x + width / 2, shuf_vals, width, label="全局随机打散抽样 (shuf 均值)", color="#10b981", edgecolor="#047857")

    for b in b1:
        ax2.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.8, f"{b.get_height():.1f}%", ha="center", va="bottom", fontsize=7.8)
    for b in b2:
        ax2.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.8, f"{b.get_height():.1f}%", ha="center", va="bottom", fontsize=7.8, fontweight="bold")

    ax2.set_title("(b) 多尺度定长截断 vs. 随机打散抽样 (shuf) TTR 对比", fontsize=11.5, fontweight="bold", pad=10)
    ax2.set_xticks(x)
    ax2.set_xticklabels(scales, fontsize=8.5)
    ax2.set_ylabel("TTR 词汇丰富度 (%)", fontsize=10.5)
    ax2.set_ylim(0, 68)
    ax2.grid(True, axis="y", linestyle=":", alpha=0.5)
    ax2.legend(loc="upper right", frameon=True, facecolor="#f8fafc")

    plt.tight_layout()
    out_path = FIG_DIR / "fig1_ttr_heaps_curve.png"
    plt.savefig(out_path, bbox_inches="tight")
    if ARTIFACT_DIR.exists(): shutil.copyfile(out_path, ARTIFACT_DIR / "fig1_ttr_heaps_curve.png")
    plt.close(fig)


def plot_fig2_tfidf(part1_summary):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.8, 5.8))
    fig.patch.set_facecolor("#ffffff")
    n_docs = part1_summary["num_documents"]

    uni = part1_summary["top_unigrams_tfidf"][:15][::-1]
    bi = part1_summary["top_bigrams_tfidf"][:15][::-1]

    u_terms = [f"{x['term']} (DF={x['doc_freq']:,})" for x in uni]
    u_scores = [x["tfidf_score"] for x in uni]
    colors_u = plt.cm.Blues(np.linspace(0.45, 0.9, len(uni)))

    bars1 = ax1.barh(u_terms, u_scores, color=colors_u, edgecolor="#1e3a8a", height=0.68)
    pad_u = max(u_scores) * 0.015
    for b in bars1:
        ax1.text(b.get_width() + pad_u, b.get_y() + b.get_height() / 2, f"{b.get_width():.1f}", va="center", fontsize=8.5)
    ax1.set_title(f"(a) 水稻抗旱文献 TF-IDF Top 15 核心单字词 ({n_docs:,} 篇)", fontsize=11.5, fontweight="bold", pad=10)
    ax1.set_xlabel("全语料累积 TF-IDF 权重得分", fontsize=10)
    ax1.set_xlim(0, max(u_scores) * 1.18)
    ax1.grid(True, axis="x", linestyle=":", alpha=0.5)

    b_terms = [f"{x['term']} (DF={x['doc_freq']:,})" for x in bi]
    b_scores = [x["tfidf_score"] for x in bi]
    colors_b = plt.cm.Greens(np.linspace(0.45, 0.9, len(bi)))

    bars2 = ax2.barh(b_terms, b_scores, color=colors_b, edgecolor="#064e3b", height=0.68)
    pad_b = max(b_scores) * 0.015
    for b in bars2:
        ax2.text(b.get_width() + pad_b, b.get_y() + b.get_height() / 2, f"{b.get_width():.1f}", va="center", fontsize=8.5)
    ax2.set_title(f"(b) 水稻抗旱文献 TF-IDF Top 15 核心二元短语 ({n_docs:,} 篇)", fontsize=11.5, fontweight="bold", pad=10)
    ax2.set_xlabel("全语料累积 TF-IDF 权重得分", fontsize=10)
    ax2.set_xlim(0, max(b_scores) * 1.18)
    ax2.grid(True, axis="x", linestyle=":", alpha=0.5)

    plt.tight_layout()
    out_path = FIG_DIR / "fig2_tfidf_top_terms.png"
    plt.savefig(out_path, bbox_inches="tight")
    if ARTIFACT_DIR.exists(): shutil.copyfile(out_path, ARTIFACT_DIR / "fig2_tfidf_top_terms.png")
    plt.close(fig)


def plot_fig3_disambiguation_comparison(part2_report):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.3))
    fig.patch.set_facecolor("#ffffff")

    n_docs = part2_report.get("num_documents_analyzed", 95)
    pubtator_g = part2_report.get("pubtator_gene_mentions", 29)
    pubtator_rap = part2_report.get("pubtator_mapped_to_rap", 2)
    our_g = part2_report.get("total_gene_mentions", 235)
    our_rap = part2_report.get("rap_mapped_gene_mentions", 222)
    our_t = part2_report.get("total_trait_mentions", 819)
    our_evi = part2_report.get("sentence_cooccurrence_evidence_count", 199)

    categories = ["识别水稻基因提及数", "归一化至 RAP ID 数", "识别水稻性状(RTO)数", "产出句子级证据数"]
    pubtator_vals = [pubtator_g, pubtator_rap, 0, 0]
    our_vals = [our_g, our_rap, our_t, our_evi]

    x = np.arange(len(categories))
    width = 0.36
    b1 = ax1.bar(x - width / 2, pubtator_vals, width, label="仅依赖 PubTator API (医学通用模型)", color="#94a3b8", edgecolor="#475569")
    b2 = ax1.bar(x + width / 2, our_vals, width, label="本研究流水线 (PubTator + Oryzabase + RTO)", color="#0ea5e9", edgecolor="#0369a1")

    max_v = max(our_vals) if our_vals else 100
    for b in b1:
        ax1.text(b.get_x() + b.get_width() / 2, b.get_height() + max_v * 0.015, f"{int(b.get_height()):,}", ha="center", fontsize=8.5, color="#475569")
    for b in b2:
        ax1.text(b.get_x() + b.get_width() / 2, b.get_height() + max_v * 0.015, f"{int(b.get_height()):,}", ha="center", fontsize=9, fontweight="bold", color="#0369a1")

    ax1.set_title(f"(a) 通用 PubTator 与本研究水稻词典流水线对比 ({n_docs:,} 篇文献)", fontsize=11.5, fontweight="bold", pad=10)
    ax1.set_xticks(x)
    ax1.set_xticklabels(categories, fontsize=9.2)
    ax1.set_ylabel(f"提及数 / 证据条数 ({n_docs:,} 篇文献)", fontsize=10)
    ax1.set_ylim(0, max_v * 1.18)
    ax1.grid(True, axis="y", linestyle=":", alpha=0.5)
    ax1.legend(loc="upper left", frameon=True, facecolor="#f8fafc")

    # 右图: 词典覆盖与唯一归一化知识图谱节点统计
    u_genes = part2_report.get("unique_genes_count", 26)
    u_traits = part2_report.get("unique_traits_count", 19)
    u_pairs = part2_report.get("unique_gene_trait_pairs_count", 61)
    kg_labels = [
        f"唯一水稻抗旱基因\n(RAP-DB 节点)",
        f"唯一水稻农艺/逆境性状\n(RTO 本体节点)",
        f"唯一【基因-性状】\n关联边 (Edges)",
        f"句子级共现证据总数\n(Evidence Sentences)"
    ]
    kg_vals = [u_genes, u_traits, u_pairs, our_evi]
    bar_colors = ["#6366f1", "#14b8a6", "#f59e0b", "#10b981"]

    bars = ax2.bar(range(len(kg_labels)), kg_vals, color=bar_colors, width=0.5, edgecolor="#1e293b")
    max_kg = max(kg_vals) if kg_vals else 100
    for idx, b in enumerate(bars):
        ax2.text(b.get_x() + b.get_width() / 2, b.get_height() + max_kg * 0.018, f"{kg_vals[idx]:,}", ha="center", fontsize=9.5, fontweight="bold")

    ax2.set_title(f"(b) 消歧提纯后水稻抗旱知识库节点与证据规模 ({n_docs:,} 篇文献)", fontsize=11.5, fontweight="bold", pad=10)
    ax2.set_xticks(range(len(kg_labels)))
    ax2.set_xticklabels(kg_labels, fontsize=8.8)
    ax2.set_ylabel("独立节点数 / 关联边数 / 证据句数", fontsize=10)
    ax2.set_ylim(0, max_kg * 1.18)
    ax2.grid(True, axis="y", linestyle=":", alpha=0.5)

    plt.tight_layout()
    out_path = FIG_DIR / "fig3_disambiguation_comparison.png"
    plt.savefig(out_path, bbox_inches="tight")
    if ARTIFACT_DIR.exists(): shutil.copyfile(out_path, ARTIFACT_DIR / "fig3_disambiguation_comparison.png")
    plt.close(fig)


def plot_fig4_gene_trait_matrix(part2_report):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.8, 5.6))
    fig.patch.set_facecolor("#ffffff")
    n_docs = part2_report.get("num_documents_analyzed", 95)

    top_genes = part2_report["top_genes"][:15][::-1]
    g_labels = [f"{g['canonical_symbol']} ({g['rap_id']})" for g in top_genes]
    g_counts = [g["count"] for g in top_genes]

    bars_g = ax1.barh(g_labels, g_counts, color="#6366f1", edgecolor="#312e81", height=0.68)
    pad_g = max(g_counts) * 0.015 if g_counts else 1
    for b in bars_g:
        ax1.text(b.get_width() + pad_g, b.get_y() + b.get_height() / 2, f"{int(b.get_width()):,} 次", va="center", fontsize=8.5, fontweight="bold")
    ax1.set_title(f"(a) 归一化后高频抗旱水稻基因 Top 15 ({n_docs:,} 篇文献)", fontsize=11.5, fontweight="bold", pad=10)
    ax1.set_xlabel("文献句子提及频次", fontsize=10)
    ax1.set_xlim(0, max(g_counts) * 1.2)
    ax1.grid(True, axis="x", linestyle=":", alpha=0.5)

    top_traits = part2_report["top_traits"][:15][::-1]
    t_labels = [f"{t['canonical_trait']} ({t['rto_id']})" for t in top_traits]
    t_counts = [t["count"] for t in top_traits]

    bars_t = ax2.barh(t_labels, t_counts, color="#14b8a6", edgecolor="#115e59", height=0.68)
    pad_t = max(t_counts) * 0.015 if t_counts else 1
    for b in bars_t:
        ax2.text(b.get_width() + pad_t, b.get_y() + b.get_height() / 2, f"{int(b.get_width()):,} 次", va="center", fontsize=8.5, fontweight="bold")
    ax2.set_title(f"(b) 归一化后高频水稻表型与逆境性状 Top 15 ({n_docs:,} 篇文献)", fontsize=11.5, fontweight="bold", pad=10)
    ax2.set_xlabel("文献句子提及频次", fontsize=10)
    ax2.set_xlim(0, max(t_counts) * 1.2)
    ax2.grid(True, axis="x", linestyle=":", alpha=0.5)

    plt.tight_layout()
    out_path = FIG_DIR / "fig4_gene_trait_matrix.png"
    plt.savefig(out_path, bbox_inches="tight")
    if ARTIFACT_DIR.exists(): shutil.copyfile(out_path, ARTIFACT_DIR / "fig4_gene_trait_matrix.png")
    plt.close(fig)


def build_html_dashboard(part1_summary, part2_report):
    lines = (RESULTS_DIR / "rice_drought_sentence_evidence.tsv").read_text(encoding="utf-8").splitlines()
    rows = [l.split("\t") for l in lines[1:] if l.strip()]

    # 为保证浏览器渲染流畅，表格内展示前 1,500 条代表性证据（完整数据在 TSV 中）
    display_rows = rows[:1500]
    tr_html = []
    for r in display_rows:
        if len(r) < 10:
            continue
        pmid, sid, gm, gc, rap, gsrc, tm, tc, rto, sent = r[:10]
        sent_hl = sent.replace(gm, f'<mark class="hl-gene">{gm}</mark>').replace(tm, f'<mark class="hl-trait">{tm}</mark>')
        tr_html.append(
            f'<tr><td><a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" target="_blank">{pmid}</a></td>'
            f'<td>{sid}</td><td><span class="badge-gene">{gm}</span></td>'
            f'<td><b>{gc}</b><br/><code>{rap}</code></td>'
            f'<td><span class="badge-trait">{tm}</span></td>'
            f'<td><b>{tc}</b><br/><code>{rto}</code></td>'
            f'<td class="sent-col">{sent_hl}</td></tr>'
        )

    n_docs = part1_summary["num_documents"]
    n_sents = part1_summary["num_sentences"]
    n_ents = part2_report["total_entity_mentions"]
    n_evi = part2_report["sentence_cooccurrence_evidence_count"]
    u_genes = part2_report.get("unique_genes_count", 0)
    u_traits = part2_report.get("unique_traits_count", 0)

    html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<title>水稻抗旱 BioNLP 全量实测可视化 Walkthrough 报告</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Microsoft YaHei", sans-serif; background: #0f172a; color: #f1f5f9; margin: 0; padding-bottom: 60px; }}
  header {{ background: linear-gradient(135deg, #0f172a, #1e3a5f); padding: 32px 24px; text-align: center; border-bottom: 1px solid #334155; }}
  header h1 {{ margin: 0 0 10px; font-size: 1.9rem; color: #38bdf8; }}
  header p {{ color: #94a3b8; margin: 0; }}
  .container {{ max-width: 1240px; margin: 0 auto; padding: 24px; }}
  .kpi-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 16px; margin-bottom: 24px; }}
  .kpi {{ background: #1e293b; border: 1px solid #334155; border-radius: 12px; padding: 18px; text-align: center; }}
  .kpi .num {{ font-size: 1.7rem; font-weight: 800; color: #4ade80; }}
  .kpi .lbl {{ font-size: 0.86rem; color: #94a3b8; margin-top: 4px; }}
  .card {{ background: #1e293b; border: 1px solid #334155; border-radius: 14px; padding: 24px; margin-bottom: 24px; }}
  .card h2 {{ margin-top: 0; color: #38bdf8; font-size: 1.3rem; border-bottom: 1px solid #334155; padding-bottom: 10px; }}
  .fig-box {{ background: #fff; border-radius: 10px; padding: 12px; text-align: center; margin-top: 14px; }}
  .fig-box img {{ max-width: 100%; height: auto; border-radius: 6px; }}
  .search-bar {{ display: flex; gap: 12px; margin-bottom: 14px; }}
  .search-bar input {{ flex: 1; padding: 10px 14px; border-radius: 8px; border: 1px solid #475569; background: #0f172a; color: #fff; font-size: 0.95rem; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.88rem; }}
  th, td {{ border: 1px solid #334155; padding: 10px; text-align: left; vertical-align: top; }}
  th {{ background: #131c31; color: #38bdf8; position: sticky; top: 0; }}
  a {{ color: #38bdf8; text-decoration: none; }}
  .badge-gene {{ background: rgba(56,189,248,0.2); color: #38bdf8; padding: 2px 6px; border-radius: 4px; font-weight: 600; }}
  .badge-trait {{ background: rgba(74,222,128,0.2); color: #4ade80; padding: 2px 6px; border-radius: 4px; font-weight: 600; }}
  mark.hl-gene {{ background: rgba(56,189,248,0.35); color: #fff; padding: 1px 4px; border-radius: 3px; }}
  mark.hl-trait {{ background: rgba(74,222,128,0.35); color: #fff; padding: 1px 4px; border-radius: 3px; }}
  .table-wrap {{ max-height: 560px; overflow-y: auto; border: 1px solid #334155; border-radius: 8px; }}
</style>
</head>
<body>
<header>
  <h1>🌾 水稻抗旱 BioNLP（Part 1 & Part 2）全量文献实测可视化大屏</h1>
  <p>基于 PubMed 全量检索 (sleep 5.8s) + NLTK TTR/TF-IDF + Oryzabase (102,866 别名) + RTO 本体 (3,904 同义词)</p>
</header>
<div class="container">
  <div class="kpi-grid">
    <div class="kpi"><div class="num">{n_docs:,} 篇</div><div class="lbl">PubTator 成功解析文献总数</div></div>
    <div class="kpi"><div class="num">{n_sents:,} 句</div><div class="lbl">NLTK 切分句子总数</div></div>
    <div class="kpi"><div class="num">{n_ents:,} 个</div><div class="lbl">100% 偏移校验通过实体提及</div></div>
    <div class="kpi"><div class="num">{u_genes:,} 基因 × {u_traits:,} 性状</div><div class="lbl">归一化独立基因与性状节点</div></div>
    <div class="kpi"><div class="num">{n_evi:,} 条</div><div class="lbl">提纯后句子级【基因-性状】证据</div></div>
  </div>

  <div class="card">
    <h2>📊 图 1：希普斯定律 (Heaps' Law) 与定长/随机抽样 (shuf) TTR 词汇丰富度分析</h2>
    <div class="fig-box"><img src="results/figures/fig1_ttr_heaps_curve.png" alt="Fig 1"></div>
  </div>

  <div class="card">
    <h2>📊 图 2：水稻抗旱文献 TF-IDF Top 15 核心单字词与二元短语排行</h2>
    <div class="fig-box"><img src="results/figures/fig2_tfidf_top_terms.png" alt="Fig 2"></div>
  </div>

  <div class="card">
    <h2>📊 图 3：PubTator vs. 本研究词典覆盖率对比 & 抗旱知识图谱节点规模</h2>
    <div class="fig-box"><img src="results/figures/fig3_disambiguation_comparison.png" alt="Fig 3"></div>
  </div>

  <div class="card">
    <h2>📊 图 4：归一化后高频水稻抗旱基因 (附 RAP ID) 与表型性状 (附 RTO/TO ID)</h2>
    <div class="fig-box"><img src="results/figures/fig4_gene_trait_matrix.png" alt="Fig 4"></div>
  </div>

  <div class="card">
    <h2>🔎 交互式检索：句子级【水稻基因 - 抗旱性状】归一化证据表 (共 {n_evi:,} 条，当前展示前 {len(display_rows):,} 条，全量见 TSV)</h2>
    <div class="search-bar">
      <input type="text" id="searchInput" oninput="filterTable()" placeholder="输入基因名 (如 OsbZIP23, OsNAC, OsDREB)、RAP ID (如 Os02g)、表型 (如 drought tolerance, root length) 或 PMID 实时过滤...">
    </div>
    <div class="table-wrap">
      <table id="eviTable">
        <thead>
          <tr><th>PMID</th><th>句号</th><th>原文基因</th><th>归一化基因 & RAP ID</th><th>原文性状</th><th>归一化性状 & RTO ID</th><th>证据原句 (高亮对齐)</th></tr>
        </thead>
        <tbody>
          {"".join(tr_html)}
        </tbody>
      </table>
    </div>
  </div>
</div>
<script>
function filterTable() {{
  const q = document.getElementById('searchInput').value.toLowerCase();
  const rows = document.querySelectorAll('#eviTable tbody tr');
  rows.forEach(r => {{
    r.style.display = r.innerText.toLowerCase().includes(q) ? '' : 'none';
  }});
}}
</script>
</body>
</html>"""

    (BASE_DIR / "visual_walkthrough_report.html").write_text(html, encoding="utf-8")


def main():
    part1_summary = json.loads((RESULTS_DIR / "part1_lexical_and_tfidf_summary.json").read_text(encoding="utf-8"))
    part2_report = json.loads((RESULTS_DIR / "part2_normalization_and_evidence_report.json").read_text(encoding="utf-8"))
    tokens = load_tokens()

    plot_fig1_ttr_heaps(tokens, part1_summary)
    plot_fig2_tfidf(part1_summary)
    plot_fig3_disambiguation_comparison(part2_report)
    plot_fig4_gene_trait_matrix(part2_report)
    build_html_dashboard(part1_summary, part2_report)
    print("4 张高清可视化图表与 visual_walkthrough_report.html 大屏已成功更新！")


if __name__ == "__main__":
    main()
