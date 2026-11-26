# 数据字段字典

仅列出写作者最可能直接使用的文件。所有文本文件均为 UTF-8；TSV 为制表符分隔；CSV 为逗号分隔（`results/tables/` 下带 UTF-8 BOM，便于 Excel 打开）。

## 1. `results/rice_drought_sentence_evidence.tsv`（9,030 行）

阶段二全部候选共现。**包含 1,162 条非水稻 PubTator-only 实体**；做水稻统计时先过滤 `RAP_ID` 以 `Os` 开头。

| 列 | 类型 | 说明 |
|---|---|---|
| `PMID` | str | PubMed ID |
| `sent_id` | int | 句号（1-based，NLTK `sent_tokenize` 顺序） |
| `gene_mention` | str | 基因原文提及 |
| `gene_canonical` | str | 归一化标准符号 |
| `RAP_ID` | str | `Os...` 或 `NCBIGene:...` / `Oryzabase:...` |
| `gene_source` | str | `Oryzabase_Dict` / `PubTator_Only` / `PubTator+Oryzabase` |
| `trait_mention` | str | 性状原文提及 |
| `trait_canonical` | str | RTO 标准名 |
| `RTO_ID` | str | `TO:` / `WTO:` |
| `sentence` | str | 证据句原文 |

## 2. `results/all_sentence_entities.tsv`（42,617 行）

所有基因/性状实体的精确字符偏移。

| 列 | 类型 | 说明 |
|---|---|---|
| `PMID`, `sent_id` | str/int | 文献与句号 |
| `begin`, `end` | int | 全文绝对字符偏移；已断言 `text[begin:end] == mention` |
| `mention` | str | 原文提及 |
| `type` | str | `Gene` / `Trait` |
| `canonical` | str | 标准名 |
| `norm_id` | str | RAP / RTO / NCBIGene |
| `source` | str | 词典/掩码来源 |
| `sentence` | str | 所属句子 |

## 3. `results/sdp_relations.tsv`（7,868 行，32 列）

水稻-only 候选的 SDP 判定。**一行 = 一个基因-性状候选对**（同一句多对会有多行；同句同对只保留一次）。

| 列 | 说明 |
|---|---|
| `PMID`, `sent_id` | 文献/句号 |
| `gene_mention`, `gene_canonical`, `RAP_ID`, `gene_source` | 基因信息 |
| `is_rice_gene` | 1（本文件全为 1） |
| `gene_begin`, `gene_end` | 基因绝对偏移 |
| `gene_head` | 基因中心词/依存标签 |
| `trait_mention`, `trait_canonical`, `RTO_ID` | 性状信息 |
| `trait_begin`, `trait_end` | 性状绝对偏移 |
| `trait_span_expanded` | 1 = 性状 span 从短 mention 扩展到完整标准短语（如 drought → drought tolerance） |
| `trait_head` | 性状中心词/依存标签 |
| `predicate` | SDP 选中的谓词原文 |
| `predicate_lemma` | 谓词 lemma |
| `predicate_reason` | `ok` / `no_argument_pattern` / `no_verb_on_path` / `no_path` |
| `relation_pattern` | `active` / `passive` / `nominalized` / `chain` / `light_verb_adjective` / `required_for` / `subtree_only` / `no_path` |
| `direction` | `gene_to_trait` / 空 / 其他 |
| `construct` | `direct_gene_subject` / `gain_of_function_construct` / `loss_of_function_construct` |
| `verb_polarity` | positive / negative / neutral / no_effect / unknown |
| `polarity` | 最终基因效应极性（已做 loss-of-function 反转） |
| `polarity_reason` | 极性判定依据 |
| `relation` | positive_regulation / negative_regulation / regulation / association / no_effect / excluded_respectively / excluded_direction / unsupported |
| `confidence` | high / medium / low |
| `path_len` | SDP 边数；-1 = 无路径 |
| `sdp_path` | 路径字符串 `token/dep -> ...` |
| `respectively_reason` | `no_respectively` / `respectively_aligned` / `respectively_mismatch_gX_tY` / `respectively_insufficient_lists` |
| `sentence` | 证据句 |

## 4. `results/sdp_audit_100.csv`（100 行）与 `sdp_audit_labels.json`

`sdp_audit_100.csv` 是抽样审计原表，包含 SDP 输出、句子、空标签列；`sdp_audit_labels.json` 是实际标签：

```json
{
  "A01": {"label": "0", "note": "..."},
  "A02": {"label": "-1", "note": "..."},
  "B05": {"label": "reg", "note": "..."}
}
```

label 取值：`1` 正向、`-1` 负向、`reg` 因果但方向未知、`0` 非因果、`uncertain`（未使用）。

## 5. `results/sdp_audit_summary.json`

| 字段 | 说明 |
|---|---|
| `sample_design` | A/B 总体、样本、seed |
| `gold_positive_in_A/B` | 各层 gold 因果数 |
| `sdp_precision` | A 层精度 |
| `sdp_precision_95ci` | Wilson CI |
| `sdp_recall_estimate` | 分层外推 recall |
| `sdp_f1_estimate` | F1 |
| `cooccurrence_baseline_precision_unweighted` | 100 条简单平均（不建议写） |
| `cooccurrence_baseline_precision_weighted` | 按总体规模加权（建议写 34.5%） |
| `polarity_accuracy` / `polarity_n` | 极性准确率与可判样本数 |
| `polarity_confusion` | gold/pred 混淆矩阵 |
| `estimated_tp_total` / `estimated_missed_total` | 外推 TP/漏检 |
| `error_cases` | 逐条错误案例 |

## 6. `results/part4_enrichment_all_terms.csv` / `part4_enrichment_significant.csv`

长表；`gene_set` 列区分 `top100` / `top50` / `top20` / `core14`。

| 列 | 说明 |
|---|---|
| `gene_set` | 查询集 |
| `Term` | TO/RTO 条目 |
| `N_background` | 背景基因数（2,218） |
| `m_term_in_background` | 背景中命中该条目的基因数 |
| `k_query_mapped` | 查询集映射到背景的基因数 |
| `x_overlap` | 查询集命中数 |
| `expected` | m*k/N |
| `fold_enrichment` | x / expected |
| `P-value` | 超几何上尾概率 |
| `Adjusted P-value` | BH q |
| `Genes` | 命中的 RAP ID，分号分隔 |

## 7. `results/part4_ricemind_crosscheck.csv`（100 行）

| 列 | 说明 |
|---|---|
| `RAP_ID`, `canonical_symbol` | 基因 |
| `sdp_relations`, `sdp_unique_pairs`, `sdp_unique_pmids`, `sdp_high_conf` | SDP 证据 |
| `sdp_positive`, `sdp_negative`, `sdp_neutral` | 极性构成 |
| `top_traits` | Top5 性状（频次） |
| `ricemind_found` | 0/1 |
| `ricemind_symbol`, `ricemind_cgsnl_name`, `ricemind_rap` | RiceMind 标准信息 |
| `ricemind_trait_class` | RiceMind `trait_class`（截断 300 字符） |
| `ricemind_drought_mentioned` | 0/1；annotation/trait_class 是否出现 drought/dehydration/water deficit |
| `ricemind_aba_mentioned` | 0/1 |
| `ricemind_stress_mentioned` | 0/1；是否出现 stress/tolerance/resistance |
| `ricemind_annotation_chars` | annotation 长度 |
| `error` | API 错误（正常为空） |

## 8. `results/part4_novel_targets.csv`（13 行）

| 列 | 说明 |
|---|---|
| `rank` | 按唯一对/高置信排序 |
| `RAP_ID`, `ricemind_symbol` | 基因 |
| `status` | `not_in_ricemind` / `in_ricemind_no_stress_annotation` |
| `sdp_unique_pairs`, `sdp_high_conf`, `sdp_positive`, `sdp_negative` | SDP 证据 |
| `top_traits` | 主要性状（频次） |
| `rationale` | 为什么是注释缺口候选（不是“新基因”） |

## 9. `results/part4_network_nodes.csv` / `part4_network_edges.csv`

- nodes：`node`, `type`（gene/trait）, `label`, `weighted_degree`, `degree`。
- edges：`source`, `target`, `weight`（证据句数）, `target_name`。
- 网络是无向二分图；边权不是因果强度，只是证据数。

## 10. `results/tables/table18_data_manifest.csv`

自动生成的结果文件清单（文件、行数、大小、类型）。写作者可据此核对数据完整性。  
`results/run_all.log` 记录 02→12 每一步的 exit code 与耗时。  
`results/part2_rice_only_stats.json` 与 `results/sdp_causal_stats.json` 是避免口径混淆的机器可读统计。