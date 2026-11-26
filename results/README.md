# Results 目录索引

本目录是水稻抗旱 BioNLP 实验的完整结果包。**先读 `../RESULTS_REPORT.md`**，它包含全部数字、表格、图、局限与写作映射。

## 快速入口

| 需求 | 文件 |
|---|---|
| 全部关键数字 | `key_numbers.json` |
| 发表级 CSV 表 | `tables/`（18 张表 + 索引） |
| 论文图（PDF/SVG/PNG） | `figures/` |
| 水稻-only 统计（避免 9,030/7,868 混淆） | `part2_rice_only_stats.json` |
| SDP 因果细分统计 | `sdp_causal_stats.json` |
| 100 条人工审计 | `sdp_audit_100.csv` + `sdp_audit_labels.json` + `tables/table10_sdp_audit_100.csv` |
| SDP 全部判定 | `sdp_relations.tsv`（7,868 行） |
| 阶段二全部候选 | `rice_drought_sentence_evidence.tsv`（9,030 行） |
| 全部实体 span | `all_sentence_entities.tsv`（42,617 行） |
| 富集全部结果 | `part4_enrichment_all_terms.csv` / `part4_enrichment_significant.csv` |
| RiceMind 交叉验证 | `part4_ricemind_crosscheck.csv` |
| 注释缺口候选 | `part4_novel_targets.csv` |
| 网络 | `part4_network_nodes.csv` / `part4_network_edges.csv` |
| 端到端日志 | `run_all.log` |

## 关键口径

- 9,030 = 阶段二全部候选（含 1,162 条非水稻 PubTator-only）。
- 7,868 = 水稻-only 候选（`RAP_ID` 以 `Os` 开头），SDP/富集/网络都用这个口径。
- 884 = SDP 保留的因果候选，不是 884 条全部正确；审计 precision 86.0%。
- 13 = RiceMind 注释缺口候选，不是“新基因”。

## 表索引（`tables/`）

| 表 | 内容 |
|---|---|
| `table1_corpus_and_key_stats.csv` | 30 个核心指标 |
| `table2_ttr.csv` | 多尺度 TTR |
| `table3_tfidf_top_terms.csv` | 60 个 TF-IDF 词/短语 |
| `table4_top_genes.csv` | Top15 基因 |
| `table5_top_traits.csv` | Top15 性状 |
| `table6_top_gene_trait_pairs.csv` | Top15 共现对 |
| `table7_sdp_relation_counts.csv` | 8 类关系 |
| `table8_sdp_confidence.csv` | 置信度 |
| `table9_sdp_top_predicates.csv` | Top25 谓词 |
| `table10_sdp_audit_100.csv` | 100 条审计 |
| `table11_sdp_error_categories.csv` | 错误类别 |
| `table12_enrichment_top25.csv` | Top25 富集 |
| `table13_ricemind_coverage.csv` | RiceMind 覆盖 |
| `table14_network_top30.csv` | Top30 网络节点 |
| `table15_novel_targets.csv` | 13 个注释缺口 |
| `table16_sdp_top_causal_pairs.csv` | Top30 因果对 |
| `table17_figure_index.csv` | 图表索引 |
| `table18_data_manifest.csv` | 结果文件清单 |

## 图索引（`figures/`）

- `fig1_ttr_heaps.*`：Heaps 定律 + TTR head/shuf
- `fig3_pipeline_ablation.*`：PubTator vs 双词典消融 + 知识库规模
- `fig4_top_genes_traits.*`：Top14 基因/性状
- `fig5_sdp_framework.*`：SDP 路径 + 过滤漏斗 + 极性
- `fig6_evaluation_enrichment.*`：审计 P/R/F1 + TO/RTO 富集
- `fig7_gene_trait_network.*`：Top 基因-性状网络
- `fig2_tfidf_terms.*`：TF-IDF 图（保留，未进入最终 8 页草稿）