# 水稻抗旱 BioNLP 实验：完整结果报告与数据交付说明

- **生成日期**：2026-11-25
- **工作目录**：`C:\Users\ADR\Desktop\2026-2027-1\BioNLP\workspace`
- **代码入口**：`run_all.py`
- **结果目录**：`results/`；发表级表格：`results/tables/`；图：`results/figures/`
- **说明**：本报告只提供实验事实、数字、表格、图和文件说明，不是论文正文；后续写作者可直接取用。

---

## 0. 实验完成度

| 阶段 | 内容 | 状态 | 主要产物 |
|---|---|---|---|
| 阶段一 | PubMed 检索、PubTator 下载、NLTK 分句分词、TTR/Heaps、TF-IDF | ✅ 完成 | `result_Pubtator_all.txt`、`part1_lexical_and_tfidf_summary.json`、`part1_tfidf_top_terms.tsv`、图 1 |
| 阶段二 | Oryzabase/RTO 词典、五重消歧、句子级实体与共现证据 | ✅ 完成 | `rice_gene_synonym_to_rap.json`、`rto_trait_synonym_to_id.json`、`rice_drought_sentence_evidence.tsv`、图 2–3 |
| 阶段三 | spaCy 依存树、最短依存路径（SDP）、极性/否定/respectively 消歧、100 条人工审计 | ✅ 完成 | `sdp_relations.tsv`、`sdp_summary.json`、`sdp_audit_100.csv`、`sdp_audit_summary.json`、图 4–5 |
| 阶段四 | TO/RTO 超几何富集、RiceMind 交叉验证、注释缺口候选、基因-性状网络 | ✅ 完成 | `part4_enrichment_*.csv/json`、`part4_ricemind_crosscheck.csv`、`part4_novel_targets.csv`、图 5–6 |
| 复现 | 一键脚本 02→12，端到端约 1.5 分钟（不含首次下载） | ✅ 通过 | `run_all.py`、`results/run_all.log` |

**重要口径说明**

1. 阶段二共输出 **9,030 条**候选共现；其中 **1,162 条**是 PubTator-only 的非水稻/未映射基因（如大豆 `GmGolS`）。阶段三按“RAP ID 以 `Os` 开头”过滤，实际输入 **7,868 条水稻候选**。写文章时请勿把 9,030 直接当作水稻关系数。
2. 阶段二实体总数 **42,617** = RAP 归一化水稻基因 **10,708** + RTO 性状 **29,966** + PubTator-only 非水稻基因 **1,943**。
3. SDP 的 **884 条因果**是“保留的候选因果”，不是 884 条全部正确；100 条审计估计 precision **86.0%**，约 **760 条**为真实调控（见 §6）。

---

## 1. 语料与检索

### 1.1 检索式

```text
("Oryza sativa"[Title/Abstract] OR "Oryza sativa"[MeSH Terms] OR "rice"[Title/Abstract])
AND "drought"[Title/Abstract]
```

- PubMed esearch 命中：**3,311** 篇（`data/pmid_raw.xml`、`data/pmid_list_all.txt`）。
- PubTator 成功解析：**3,306** 篇；**5 篇 2026 年新文献尚未被 PubTator 索引**。
- 下载策略：单篇循环，每篇 `time.sleep(5.8)`；批量模式断点续传。脚本：`scripts/01_fetch_pubmed_and_pubtator.py`、`.sh`。
- 5 篇未索引文献已用 PubMed efetch 补齐题录/摘要，存于 `data/unindexed_pmids.tsv`（不参与主统计，避免与阶段二已提交数字冲突）。

| 未索引 PMID | 标题（截断） |
|---|---|
| 42880016 | Mild alternate wetting and drying irrigation improves stomatal and mesophyll conductance... |
| 42881188 | WLG coordinates multi-stress tolerance in rice by regulating DST-dependent... |
| 42884139 | Proteomic and Molecular Biological Analyses of Response Mechanism in Rice... |
| 42884770 | Specialized Metabolites in Cereal Crops Under Abiotic Stress... |
| 42884918 | Crop-Level Screening and Physiological Evaluation Reveal Drought-Protective Effects... |

### 1.2 语料规模

| 指标 | 数值 |
|---|---:|
| 文献数 | 3,306 |
| 句子数（NLTK `sent_tokenize`） | 35,345 |
| 词元数（小写、含标点过滤后的词） | 795,401 |
| 独立词型数 | 34,247 |
| 原始 TTR（V/N） | 4.31% |
| PubTator 原始实体 | 47,157 |
| 偏移校验通过 | 47,156（99.998%） |

唯一未通过校验的实体：PMID `40588031` 的 `all-trans-β-carotene`（Unicode β 与 PubTator 原始编码差异，属 Chemical，不影响基因/性状）。

---

## 2. 词典与归一化

### 2.1 词典规模

| 资源 | 原始规模 | 清洗后 | 说明 |
|---|---:|---:|---|
| Oryzabase 基因记录 | 22,023 | 22,023 | 15,862 有 RAP ID；6,386 回退 Oryzabase ID |
| 基因别名映射 | — | **102,860** | 任意别名 → `rap_id` + 标准符号 |
| RTO 标准性状节点 | — | **2,022** | 过滤 `is_obsolete` 与顶层抽象词 |
| RTO 可匹配同义词 | — | **3,904** | 同义词/精确名 → `TO:`/`WTO:` ID |

### 2.2 五重农学消歧规则（写作时可直接作为 Methods 条目）

1. **缺省占位符回退**：CGSNL 首列为 `_` 的 6,000 余条基因，剔除 `_`，回退首个有效同义词作为展示名。
2. **暂定名 `(t)` 剥离**：正则 `\(t\)\*?` 剥离 `Hd6(t)*` 的 `t`，避免单字母别名污染。
3. **激素/生理缩写黑名单**：`ABA`、`AsA`、`JA`、`ROS`、`MDA`、`CAT`、`SOD` 等与突变体同形缩写。
4. **育种世代/泛家族名**：`M1`–`M5`、`BC1`、`B2`、`B3`；单独出现的 `MYB`、`WRKY`、`NAC`。
5. **产量单位/品种名**：`t ha-1` 中的 `ha-1`（与 `AHA8` 别名冲突，50 次）、对照品种 `Apo`（与穗型基因 `APO` 冲突，47 次）。

### 2.3 归一化自检样例

| 输入别名 | 标准符号 | RAP ID |
|---|---|---|
| `OsDREB1A` | DREB1A | `Os09g0522200` |
| `OsbZIP23` | BZIP23 | `Os02g0766700` |
| `SNAC1` | ONAC2 | `Os03g0815100` |
| `OsCKA2-1` / `HD6` | HD6 | `Os03g0762000` |
| `drought resistance` | drought tolerance | `TO:0000276` |
| `leaf rolling` | leaf rolling | `TO:0000085` |
| `stomatal closure` | stomatal resistance | `TO:0000522` |

---

## 3. 词汇统计（阶段一）

### 3.1 多尺度 TTR

| 窗口 N | head 词型数 | head TTR | shuf 平均 TTR | 差值（百分点） |
|---:|---:|---:|---:|---:|
| 1,000 | 542 | 54.20% | 56.19% | +1.99 |
| 5,000 | 1,736 | 34.72% | 36.39% | +1.67 |
| 10,000 | 2,648 | 26.48% | 28.90% | +2.42 |
| 50,000 | 7,015 | 14.03% | 15.65% | +1.62 |
| 100,000 | 10,484 | 10.48% | 11.74% | +1.26 |
| 500,000 | 26,221 | 5.24% | 5.47% | +0.23 |
| 全量 795,401 | 34,247 | 4.31% | 4.31% | 0 |

**可用结论**：原始 TTR 随 N 幂律衰减，符合 Heaps 定律；固定窗口下全局随机抽样（shuf）高于顺序截断（head）：1k–50k 窗口约高 1.6–2.4 个百分点，100k 高 1.26 个百分点，500k 仅高 0.23 个百分点。跨文献随机抽样在中小尺度上更能反映领域词汇多样性。

### 3.2 TF-IDF Top 12 单字词

| 排名 | 词 | 语料频次 | 文档频次 DF | TF-IDF |
|---:|---|---:|---:|---:|
| 1 | drought | 10,819 | 3,079 | 64.8793 |
| 2 | stress | 8,951 | 2,535 | 61.1784 |
| 3 | genes | 6,272 | 1,829 | 51.6604 |
| 4 | tolerance | 4,580 | 1,731 | 43.0203 |
| 5 | expression | 4,099 | 1,684 | 37.1190 |
| 6 | gene | 3,974 | 1,678 | 36.0633 |
| 7 | abiotic | 2,643 | 1,312 | 27.2110 |
| 8 | water | 2,222 | 899 | 26.5615 |
| 9 | yield | 2,133 | 896 | 25.3363 |
| 10 | root | 1,851 | 605 | 25.2574 |
| 11 | stresses | 2,241 | 1,173 | 24.9027 |
| 12 | aba | 1,551 | 589 | 24.1021 |

### 3.3 TF-IDF Top 12 二元短语

| 排名 | 短语 | 语料频次 | 文档频次 DF | TF-IDF |
|---:|---|---:|---:|---:|
| 1 | drought stress | 2,752 | 1,253 | 30.4537 |
| 2 | drought tolerance | 1,878 | 883 | 25.4562 |
| 3 | abiotic stress | 1,224 | 755 | 16.2018 |
| 4 | stress tolerance | 1,049 | 669 | 15.0275 |
| 5 | abiotic stresses | 1,039 | 719 | 14.0910 |
| 6 | gene family | 733 | 365 | 12.0124 |
| 7 | drought resistance | 582 | 281 | 11.4595 |
| 8 | grain yield | 589 | 297 | 10.2297 |
| 9 | transcription factor | 533 | 391 | 10.0633 |
| 10 | gene expression | 593 | 421 | 9.5702 |
| 11 | abscisic acid | 586 | 505 | 9.5539 |
| 12 | transcription factors | 554 | 403 | 9.3701 |

**检索反证**：`drought stress` 覆盖 1,253 篇，`drought tolerance` 覆盖 883 篇，前者多覆盖 **370 篇（41.9%）**；`drought resistance` 额外覆盖 281 篇。窄口径 `rice AND drought tolerance` 会漏检约三分之二核心文献。

---

## 4. 实体识别与共现基线（阶段二）

### 4.1 总量

| 指标 | 数值 |
|---|---:|
| 实体提及总数 | 42,617 |
| RAP 归一化水稻基因 | 10,708 |
| RTO 性状 | 29,966 |
| PubTator-only 非水稻/未映射基因 | 1,943 |
| PubTator 原生 Gene 提及 | 1,767 |
| PubTator 映射到 RAP | 473 |
| 水稻候选共现对（RAP=Os） | 7,868 |
| 水稻唯一基因-性状对 | 4,678 |
| 水稻唯一基因 | 1,344 |
| 水稻唯一性状 | 190 |
| 全部候选（含非水稻） | 9,030 |

**对比结论**：PubTator 仅 473 处映射到 RAP；双词典流水线 RAP 归一化 **10,708**，提升 **22.6 倍**；性状从 0 到 29,966，补齐了 PubTator 完全缺失的维度。

### 4.2 Top 15 基因（按提及数）

| 排名 | 基因 | RAP ID | 提及数 |
|---:|---|---|---:|
| 1 | DREB1A | Os09g0522200 | 103 |
| 2 | NCED3 | Os07g0154100 | 83 |
| 3 | BZIP23 | Os02g0766700 | 80 |
| 4 | ONAC2 | Os03g0815100 | 70 |
| 5 | NAC6 | Os01g0884300 | 70 |
| 6 | BZIP46 | Os06g0211200 | 68 |
| 7 | DREB1B | Os09g0522000 | 67 |
| 8 | PIP2A | Os07g0448800 | 65 |
| 9 | LEA3-1 | Os05g0542500 | 62 |
| 10 | P5CS | Os05g0455500 | 59 |
| 11 | DREB2A | Os01g0165000 | 57 |
| 12 | ABI5 | Os01g0859300 | 53 |
| 13 | DST | Os03g0786400 | 53 |
| 14 | DRO1 | Os09g0439800 | 51 |
| 15 | D1 | Os05g0333200 | 51 |

**注意**：D1（`Os05g0333200`）是赤霉素信号基因，也常被写作 `RGA`；它与“root growth angle (RGA)”缩写冲突，是阶段二/三的已知误差源（见 §6 的 A01）。

### 4.3 Top 15 性状（按提及数）

| 排名 | 性状 | RTO/TO ID | 提及数 |
|---:|---|---|---:|
| 1 | drought tolerance | TO:0000276 | 5,824 |
| 2 | drought | WTO:0000218 | 5,422 |
| 3 | salt in the soil | WTO:0000581 | 2,452 |
| 4 | growth | WTO:0000021 | 2,081 |
| 5 | development | WTO:0000011 | 1,635 |
| 6 | water quantity | WTO:0000117 | 1,408 |
| 7 | abiotic stress | WTO:0000062 | 988 |
| 8 | cold | WTO:0000129 | 773 |
| 9 | crop yield | WTO:0000071 | 712 |
| 10 | stress trait | TO:0000164 | 669 |
| 11 | heat | WTO:0000152 | 642 |
| 12 | soil property | WTO:0000109 | 641 |
| 13 | quality | WTO:0000048 | 366 |
| 14 | water stress | WTO:0000203 | 279 |
| 15 | salt tolerance | TO:0006001 | 252 |

### 4.4 Top 10 基因-性状共现对

| 排名 | 基因（RAP） | 性状（RTO/TO） | 证据句数 |
|---:|---|---|---:|
| 1 | BZIP23 (Os02g0766700) | drought tolerance (TO:0000276) | 22 |
| 2 | DST (Os03g0786400) | drought (WTO:0000218) | 20 |
| 3 | ONAC2 (Os03g0815100) | drought tolerance (TO:0000276) | 20 |
| 4 | DREB2A (Os01g0165000) | drought (WTO:0000218) | 16 |
| 5 | DREB1A (Os09g0522200) | drought tolerance (TO:0000276) | 16 |
| 6 | ONAC2 (Os03g0815100) | drought (WTO:0000218) | 16 |
| 7 | P5CS (Os05g0455500) | drought (WTO:0000218) | 15 |
| 8 | DST (Os03g0786400) | salt tolerance (TO:0006001) | 15 |
| 9 | DRO1 (Os09g0439800) | drought (WTO:0000218) | 14 |
| 10 | BZIP46 (Os06g0211200) | drought tolerance (TO:0000276) | 14 |

### 4.5 别名汇聚实证（Top 10 基因）

| 基因 | RAP | 语料中的异构写法（频次） |
|---|---|---|
| DREB1A | Os09g0522200 | DREB1A (37), OsDREB1A (31), CBF3 (24), DREB1a (9), Cbf3 (2) |
| NCED3 | Os07g0154100 | OsNCED3 (36), OsNCED4 (27), NCED4 (10), NCED3 (8) |
| BZIP23 | Os02g0766700 | OsbZIP23 (70), osbZIP23 (4), OsZIP23 (2), OsAREB2 (2) |
| ONAC2 | Os03g0815100 | SNAC1 (51), OsNAC9 (12), OsSNAC1 (6), ONAC2 (1) |
| NAC6 | Os01g0884300 | OsNAC6 (40), OsNAC2 (14), SNAC2 (10), NAC6 (3) |
| BZIP46 | Os06g0211200 | OsbZIP46 (41), OsABF2 (8), OsbZIP46CA1 (7), ABF3 (5) |
| DREB1B | Os09g0522000 | OsDREB1B (38), DREB1B (12), CBF1 (9), OsDREB1C (7) |
| PIP2A | Os07g0448800 | OsPIP2 (49), PIP2 (12), OsPIP2a (3), OsPIP2-1 (1) |
| LEA3-1 | Os05g0542500 | OsLEA3 (34), OsLEA3-1 (11), LEA3 (7), Lea3 (4) |
| P5CS | Os05g0455500 | P5CS (33), OsP5CS (12), OsP5CS1 (9), P5CS1 (5) |

---

## 5. SDP 因果抽取（阶段三）

### 5.1 算法要点（可直接写入 Methods）

1. 对每个候选句用 `spaCy en_core_web_sm` 构建依存树。
2. 实体字符 span → token span → 名词短语中心词。
3. 在依存图上求基因中心词与性状中心词的最短路径 `P`，长度为 `L`。
4. 沿 `P` 及其支配动词选择谓词 `v`，要求满足以下论元模式之一：
   - 主动：基因在 `v` 的 `nsubj/nsubjpass` 子树，性状在 `dobj/attr` 或 `to/of/for` 介词宾语；
   - 被动：性状为 `nsubjpass`，基因在 `by`-agent；
   - 名词化：`overexpression/knockout of GENE` 作主语，性状作宾语；
   - 链式：`v1` 主语含基因，`v2` 宾语为性状，且 `v2` 在 `v1` 的补足语/关系从句中。
5. 极性：positive/negative/neutral 谓词词典 + 程度副词 + 否定词；若主语含 loss-of-function 构建（knockout、silencing、mutation、inactivation、repression 等），极性取反。
6. `respectively`：按基因/性状出现顺序建立配对，索引不一致则标记 `excluded_respectively`。
7. 置信度：路径长度 + 论元模式 + 极性明确度 → high/medium/low。

### 5.2 总体关系分布

| 关系 | 数量 | 占 7,868 |
|---|---:|---:|
| positive_regulation | 575 | 7.31% |
| negative_regulation | 150 | 1.91% |
| regulation（中性） | 159 | 2.02% |
| **因果小计** | **884** | **11.24%** |
| association | 1,379 | 17.53% |
| excluded_direction | 5,489 | 69.76% |
| excluded_respectively | 48 | 0.61% |
| no_effect | 35 | 0.44% |
| unsupported（无路径） | 33 | 0.42% |

- 884 条因果覆盖 **385 个基因、72 个性状、693 个唯一基因-性状对、414 篇 PMID、647 个句子**。
- 正/负极性比 **3.83:1**（575/150）。
- association 覆盖 1,189 个唯一对（与因果层有重叠）。

### 5.3 置信度与句法模式

| 置信度 | 全部 7,868 条 | 884 条因果中 |
|---|---:|---:|
| high | 640 | **511** |
| medium | 957 | **200** |
| low | 6,271 | **173** |

884 条因果的句法模式：active 635、chain 182、light_verb_adjective 64、passive 3。构建类型：direct_gene_subject 557、gain_of_function 254、loss_of_function 73。可见 loss-of-function 反转规则实际覆盖 73 条，是极性判别的重要组成。

### 5.4 Top 15 谓词

| 排名 | 谓词 lemma | 次数 |
|---:|---|---:|
| 1 | regulate | 611 |
| 2 | induce | 516 |
| 3 | be | 419 |
| 4 | enhance | 376 |
| 5 | increase | 290 |
| 6 | improve | 280 |
| 7 | show | 251 |
| 8 | play | 191 |
| 9 | exhibit | 190 |
| 10 | confer | 154 |
| 11 | involve | 142 |
| 12 | upregulate | 107 |
| 13 | reduce | 105 |
| 14 | result | 101 |
| 15 | decrease | 95 |

### 5.5 `respectively` 处理统计

| 情况 | 数量 |
|---|---:|
| 无 `respectively` | 7,742 |
| 有 `respectively` 但并列列表不足 | 53 |
| 对齐正确 | 24 |
| 对齐错误（被排除） | 48（g0_t1=11, g1_t0=7, g2_t0=3 …） |

**可直接写的结论**：SDP 显式排除 48 条 `respectively` 交叉错配；阶段二的两个经典错误（`OsBIHD1` 对应 blast defense 而非 drought；`SAMs` 与 `SAMS1` 同形）中，前者已被消除，后者仍属实体词典问题。

---

## 6. 100 条分层人工审计（阶段三评测）

### 6.1 抽样设计

| 层 | 总体 | 抽样 | 说明 |
|---|---:|---:|---|
| A：SDP 判为因果 | 884 | 50 | 估计 precision |
| B：SDP 判为非因果 | 6,984 | 50 | 估计漏检率 → recall |
| 合计 | 7,868 | 100 | 随机种子 42 |

B 层按 relation 分层配额：association 20、excluded_direction 15、excluded_respectively 5、no_effect 5、unsupported 5。

### 6.2 标注定义

- `1`：句子明确表明该基因正向调控该性状（含 loss-of-function 表型反证）。
- `-1`：句子明确表明该基因负向调控该性状。
- `reg`：明确因果调控，但方向未知（如 “regulates drought tolerance”）。
- `0`：仅共现/关联/表达相关/调控其他性状/否定或无效应/实体错配。

标签文件：`results/sdp_audit_labels.json`；审计表：`results/sdp_audit_100.csv`；合并表：`results/tables/table10_sdp_audit_100.csv`。

### 6.3 指标

| 指标 | 数值 | 说明 |
|---|---:|---|
| SDP precision | **86.0%** | 43/50 A 层 gold 因果 |
| Precision 95% CI | **73.8%–93.0%** | Wilson |
| 共现基线 precision（加权） | **34.5%** | 按 A/B 总体规模加权 |
| 共现基线 precision（未加权） | 57.0% | 仅作对照，不建议写 |
| SDP recall | **28.0%** | 分层外推 |
| SDP F1 | **42.2%** | 由 P/R 计算 |
| 估计 TP | 760.2 | 86.0% × 884 |
| 估计漏检 | 1,955.5 | 28.0% × 6,984 |
| 极性准确率 | **88.2%** | 34 条可判极性样本 |

### 6.4 极性混淆矩阵

| Gold \ Pred | Pred negative | Pred positive |
|---|---:|---:|
| Gold negative | 7 | 2 |
| Gold positive | 2 | 23 |

- 正确 30/34；错误 4 条：A04、A07、A39、A42。
- 另有 A09（SDP 中性、gold 正向）未被计入“极性错误”，属于方向漏判。

### 6.5 错误类别

| 类别 | 数量 | 审计 ID |
|---|---:|---|
| SDP false positive（A 层 gold=0） | 7 | A01, A03, A17, A18, A22, A45, A48 |
| SDP false negative（B 层 gold 因果） | 14 | B05, B06, B09, B12, B16, B19, B22, B26, B28, B29, B31, B32, B33, B38 |
| 极性冲突（gold ±1 与 SDP 相反） | 4 | A04, A07, A39, A42 |

**错误原因归类**

1. **实体/缩写错配**：A01（`RGA` 根生长角 vs 基因 `D1`）、A48（转录组列表无表型）。
2. **间接调控/下游机制**：A03（OsCaM1.1 调控 CNGC 蛋白）、A17（调控分蘖以稳定产量）、A18（抑制耐旱负调控因子）、A22（影响 salt 下渗透势）、A45（诱导 ABA/气孔关闭但未测耐旱）。
3. **谓词/论元选择错误**：A04（increased sensitivity to drought → 实际负向）、A42（decreased ABA sensitivity 与 increased yield 并列，SDP 取错分句）。
4. **loss-of-function 反转未命中**：A07（knocking out 未触发反转）、A39（loss of function 未触发）。
5. **名词语义/轻动词未覆盖**：B 层 14 条，典型为 `showed enhanced tolerance`、`is required for development`、`mutant with increased tolerance`、`positively involved in`、`were markedly tolerant to`。

### 6.6 代表性审计案例（可直接做论文表）

| ID | 基因/性状 | SDP 输出 | Gold | 诊断 | 原文（节选） |
|---|---|---|---|---|---|
| A02 | SST / salt | negative | -1 | 正确：直接负调控 | OsSPL10 ... negatively regulate salt and drought stress tolerance in rice. |
| A05 | SST / drought tolerance | positive | +1 | 正确：直接正调控 | OsSPL10 confers drought tolerance by regulating OsNAC2 expression. |
| A46 | PRX130 / drought tolerance | positive | +1 | 正确：knockout 降低 tolerance → 基因正向 | Knockout of OsPRX130 ... compromising the drought tolerance in rice. |
| A36 | BHLH148 / drought tolerance | regulation | reg | 正确：respectively 对齐 | OsbHLH148 and OsBIHD1, which regulate drought tolerance and blast defense, respectively. |
| A13 | GRX22 / drought tolerance | negative | -1 | 正确：silencing 提高 tolerance → 基因负向 | Silencing of OsGRXS17 in rice improves drought stress tolerance. |
| A01 | D1 / water quantity | regulation | 0 | 错误：实体错配 | The root growth angle (RGA) ... affects the area in which roots capture water. |
| A04 | PFA-DSP1 / drought tolerance | positive | -1 | 错误：overexpression 增加 drought sensitivity | Ectopic overexpression of OsPFA-DSP1 ... increased sensitivity to drought stress. |
| A07 | GSTU17 / drought tolerance | negative | +1 | 错误：knockout 降低 tolerance 未反转 | knocking out OsGSTU17 significantly reduced the drought tolerance. |
| A39 | LG1 / salt | positive | -1 | 错误：loss of function 未反转 | the loss of OsSPL8 function enhanced tolerance to drought and salt stresses. |
| A42 | PP2C49 / crop yield | negative | +1 | 错误：取错并列分句 | OsPP2C49 overexpression led to decreased ABA sensitivity and increased rice grain yield. |

---

## 7. TO/RTO 富集（阶段四）

### 7.1 设计

- 查询集：SDP 因果证据排序的 **Top100 基因**（全部 100 个均能映射到背景库）。
- 背景：`lab03/task3_to_enrichment/results/RiceMind_TO_RTO_background.gmt`，**N = 2,218** 个基因，**272** 个 TO/RTO 条目。
- 检验：`scipy.stats.hypergeom.sf(x-1, N, m, k)`，BH 校正。
- 结果：151 个条目命中，**53 个 q<0.05**。
- 对照：Top50 46 个显著；Top20 40 个；core14 68 个（core14 有选择偏倚，建议只作补充）。

### 7.2 Top 15 显著条目（Top100 查询集）

| 排名 | 条目 | x | k | m | N | Fold | q |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | drought tolerance (TO:0000276) | 84 | 100 | 803 | 2,218 | 2.32 | 2.58e-21 |
| 2 | drought (WTO:0000218) | 84 | 100 | 867 | 2,218 | 2.15 | 5.51e-19 |
| 3 | salt in the soil (WTO:0000581) | 66 | 100 | 531 | 2,218 | 2.76 | 5.99e-18 |
| 4 | growth (WTO:0000021) | 43 | 100 | 315 | 2,218 | 3.03 | 3.76e-11 |
| 5 | crop yield (WTO:0000071) | 24 | 100 | 97 | 2,218 | 5.49 | 4.31e-11 |
| 6 | abscisic acid sensitivity (TO:0002667) | 13 | 100 | 27 | 2,218 | 10.68 | 7.86e-10 |
| 7 | abiotic stress trait (TO:0000168) | 21 | 100 | 86 | 2,218 | 5.42 | 1.40e-09 |
| 8 | salt tolerance (TO:0006001) | 21 | 100 | 91 | 2,218 | 5.12 | 3.93e-09 |
| 9 | development (WTO:0000011) | 35 | 100 | 262 | 2,218 | 2.96 | 1.24e-08 |
| 10 | stress trait (TO:0000164) | 29 | 100 | 187 | 2,218 | 3.44 | 1.60e-08 |
| 11 | cold (WTO:0000129) | 32 | 100 | 231 | 2,218 | 3.07 | 2.96e-08 |
| 12 | stomatal resistance (TO:0000522) | 17 | 100 | 68 | 2,218 | 5.54 | 5.22e-08 |
| 13 | water quantity (WTO:0000117) | 23 | 100 | 137 | 2,218 | 3.72 | 2.36e-07 |
| 14 | grain size (TO:0000397) | 8 | 100 | 20 | 2,218 | 8.87 | 2.02e-05 |
| 15 | development (RTO:0001542) | 8 | 100 | 22 | 2,218 | 8.07 | 4.44e-05 |

**可用结论**：ABA 敏感性与气孔抗性富集倍数最高（10.68×、5.54×），说明 SDP 因果基因集在功能层面高度内聚；`drought tolerance`、`drought`、`salt` 的高命中是领域预期。不要写成“富集证明因果关系”，只能写“功能一致性/生物学保真度”。

---

## 8. RiceMind 交叉验证与注释缺口（阶段四）

### 8.1 覆盖率

| 指标 | 数值 |
|---|---:|
| 查询 Top100 基因 | 100 |
| 在 RiceMind 找到标准 RAP 记录 | 99（99%） |
| 有 stress/tolerance 注释 | 85 |
| 明确提到 drought/water deficit | **13** |
| 提到 ABA/abscisic acid | 5 |
| 注释缺口候选（≥3 唯一对且无 stress 注释） | **13** |

**可用结论**：RiceMind 对 Top100 基因的覆盖很高（99%），但“一般 stress 注释”与“drought 特异注释”之间存在明显落差（85 vs 13）；这说明文献挖掘的价值主要在补足干旱特异维度，而不是重新发现基因。

### 8.2 注释缺口候选 Top 5

| 排名 | 基因 | RAP ID | 唯一对 | 主要 SDP 性状（频次） |
|---:|---|---|---:|---|
| 1 | DRAP1 | Os11g0544700 | 5 | drought tolerance (2), salt tolerance (1) |
| 2 | ANN1 | Os01g0497400 | 4 | stress trait (1), abiotic stress (1), heat (1) |
| 3 | CALA | Os05g0491100 | 4 | abiotic stress (1), heat (1), growth (1) |
| 4 | BHC | Os04g0624500 | 4 | stomatal opening (2), salt (2) |
| 5 | — | Os01g0178600 | 4 | drought (1), drought tolerance (1), cold (1) |

**写作注意**：这些不是“新基因”，而是 **RiceMind 当前注释缺口候选**。CALA、BHC 等已有文献报道，不能宣称首次发现；应写成“需人工审编/实验验证的候选”。

---

## 9. 基因-性状网络

- 节点：**457**；边：**693**（与 693 个唯一因果对一致）。
- 网络由 SDP 因果层构建，边权为证据句数。

### 9.1 Top 10 基因（按加权度）

| 排名 | 基因 | RAP ID | 加权度 | 度 |
|---:|---|---|---:|---:|
| 1 | DSR2 | Os01g0839200 | 15 | 5 |
| 2 | ONAC2 | Os03g0815100 | 10 | 6 |
| 3 | WOX11 | Os07g0684900 | 9 | 3 |
| 4 | HDA704 | Os07g0164100 | 9 | 6 |
| 5 | OsMYB4 | Os04g0517100 | 9 | 5 |
| 6 | HSP18.0 | Os03g0267000 | 7 | 5 |
| 7 | CDPK9 | Os03g0688300 | 7 | 5 |
| 8 | RRS1 | Os04g0594100 | 7 | 5 |
| 9 | DST | Os03g0786400 | 7 | 5 |
| 10 | OFP4 | Os01g0823500 | 7 | 6 |

### 9.2 Top 10 性状（按加权度）

| 排名 | 性状 | RTO/TO ID | 加权度 | 度 |
|---:|---|---|---:|---:|
| 1 | drought tolerance | TO:0000276 | 286 | 196 |
| 2 | drought | WTO:0000218 | 119 | 93 |
| 3 | salt in the soil | WTO:0000581 | 102 | 71 |
| 4 | growth | WTO:0000021 | 47 | 42 |
| 5 | stress trait | TO:0000164 | 29 | 27 |
| 6 | development | WTO:0000011 | 28 | 23 |
| 7 | crop yield | WTO:0000071 | 24 | 19 |
| 8 | stomatal resistance | TO:0000522 | 22 | 17 |
| 9 | cold | WTO:0000129 | 17 | 14 |
| 10 | abscisic acid sensitivity | TO:0002667 | 15 | 15 |

**可用结论**：`drought tolerance` 是绝对核心 hub；`salt in the soil`、`growth`、`crop yield`、`stomatal resistance`、`ABA sensitivity` 构成第二层。不要把网络中心性写成因果方向，边是无向共现/调控关系。

---

## 10. 图表清单

| 图 | 文件 | 面板 | 内容 | 数据来源 |
|---|---|---|---|---|
| Fig. 1 | `results/figures/fig1_ttr_heaps.pdf` | a,b | Heaps 定律 / TTR head vs shuf | `part1_lexical_and_tfidf_summary.json` |
| Fig. 2 | `results/figures/fig3_pipeline_ablation.pdf` | a,b | PubTator vs 双词典消融；知识库规模 | `part2_normalization_and_evidence_report.json` |
| Fig. 3 | `results/figures/fig4_top_genes_traits.pdf` | a,b | Top14 基因/性状 | `part2_normalization_and_evidence_report.json` |
| Fig. 4 | `results/figures/fig5_sdp_framework.pdf` | a,b,c | SDP 路径示例；过滤漏斗；极性构成 | `sdp_relations.tsv`、`sdp_summary.json` |
| Fig. 5 | `results/figures/fig6_evaluation_enrichment.pdf` | a,b | 审计 P/R/F1；TO/RTO 富集 | `sdp_audit_summary.json`、`part4_enrichment_summary.json` |
| Fig. 6 | `results/figures/fig7_gene_trait_network.pdf` | 单 | Top 基因-性状网络 | `part4_network_nodes.csv`、`part4_network_edges.csv` |

每张图同时有 `.pdf`/`.svg`/`.png`；LaTeX 用 PDF，网页/预览用 PNG。`fig2_tfidf_terms.pdf`（TF-IDF 图）仍保留在 `results/figures/`，但最终 8 页论文未使用；若写作者需要，可加入。

---

## 11. 数据文件清单与字段

### 11.1 核心数据

| 文件 | 行数 | 说明 |
|---|---:|---|
| `data/result_Pubtator_all.txt` | 3,306 篇 | PubTator 全文+实体 |
| `data/unindexed_pmids.tsv` | 5 | 未索引文献题录/摘要 |
| `results/rice_drought_sentence_evidence.tsv` | 9,030 | 全部候选共现（含 1,162 非水稻） |
| `results/all_sentence_entities.tsv` | 42,617 | 所有实体 span 与归一化 |
| `results/sdp_relations.tsv` | 7,868 | 水稻候选 SDP 判定（32 列） |
| `results/sdp_error_cases.tsv` | 6,680 | 低置信/排除/否定/并列错配案例 |
| `results/sdp_audit_100.csv` | 100 | 分层审计原表（含空标签列） |
| `results/sdp_audit_labels.json` | 100 | 人工标签与备注 |
| `results/part4_enrichment_all_terms.csv` | 272×4 组 | 四个查询集的全部富集结果 |
| `results/part4_enrichment_significant.csv` | — | q<0.05 子集 |
| `results/part4_ricemind_crosscheck.csv` | 100 | Top100 基因 RiceMind 画像 |
| `results/part4_novel_targets.csv` | 13 | 注释缺口候选 |
| `results/part4_network_nodes.csv` | 30 | Top 基因/性状节点 |
| `results/part4_network_edges.csv` | — | Top 子网络边 |

### 11.2 关键字段（`sdp_relations.tsv`）

| 字段 | 含义 |
|---|---|
| `PMID`, `sent_id` | 文献与句子编号 |
| `gene_mention`, `gene_canonical`, `RAP_ID` | 基因原文提及/标准符号/RAP |
| `gene_begin`, `gene_end` | 基因绝对字符偏移 |
| `trait_mention`, `trait_canonical`, `RTO_ID` | 性状原文/标准名/本体 ID |
| `trait_begin`, `trait_end`, `trait_span_expanded` | 性状偏移；是否扩展到完整短语 |
| `predicate`, `predicate_lemma` | SDP 选中的谓词 |
| `relation_pattern` | active/passive/nominalized/chain/light_verb_adjective/required_for |
| `direction` | gene_to_trait / 其他 |
| `construct` | gain/loss/direct |
| `verb_polarity`, `polarity` | 谓词极性 vs 最终基因效应极性 |
| `relation` | positive/negative/regulation/association/no_effect/excluded_* |
| `confidence` | high/medium/low |
| `path_len`, `sdp_path` | 路径长度与路径字符串 |
| `respectively_reason` | `respectively` 对齐诊断 |
| `sentence` | 原始证据句 |

### 11.3 发表级表格（`results/tables/`）

| 表 | 内容 |
|---|---|
| `table1_corpus_and_key_stats.csv` | 30 个核心指标 |
| `table2_ttr.csv` | 多尺度 TTR |
| `table3_tfidf_top_terms.csv` | 60 个 TF-IDF 词/短语 |
| `table4_top_genes.csv` | Top15 基因 |
| `table5_top_traits.csv` | Top15 性状 |
| `table6_top_gene_trait_pairs.csv` | Top15 共现对 |
| `table7_sdp_relation_counts.csv` | 8 类关系计数/占比 |
| `table8_sdp_confidence.csv` | 置信度分布 |
| `table9_sdp_top_predicates.csv` | Top25 谓词 |
| `table10_sdp_audit_100.csv` | 100 条审计 + gold + 极性正确性 |
| `table11_sdp_error_categories.csv` | 错误类别计数 |
| `table12_enrichment_top25.csv` | Top25 富集条目 |
| `table13_ricemind_coverage.csv` | RiceMind 覆盖指标 |
| `table14_network_top30.csv` | Top30 网络节点 |
| `table15_novel_targets.csv` | 13 个注释缺口候选 |
| `table16_sdp_top_causal_pairs.csv` | Top30 因果对 |
| `table17_figure_index.csv` | 图表索引 |
| `table18_data_manifest.csv` | 结果文件清单 |

---

## 12. 复现方式

### 12.1 环境

- Python 3.14.7
- spaCy 3.8.18 + `en_core_web_sm`
- NLTK 3.10.3
- networkx 3.7
- numpy 2.5.2 / scipy 1.18.1
- matplotlib 3.11.1
- XeLaTeX（MiKTeX 25.12；原 TeX Live 2026 的 E: 盘已不可用）

安装：
```powershell
pip install -r requirements.txt
python -m spacy download en_core_web_sm
```

### 12.2 一键运行

```powershell
cd C:\Users\ADR\Desktop\2026-2027-1\BioNLP\workspace
python run_all.py
```

- 默认复用 `data/result_Pubtator_all.txt`，不重新下载。
- 10 步实测耗时（本机）：02 18.9s；03 0.9s；04 11.6s；06 12.7s；07 33.9s；08 0.3s；09 1.3s；10 6.6s；11 5.9s；12 0.3s。
- 从零下载 3,306 篇约 5.3 小时（5.8s/篇），无需重复。
- 审计样本已冻结：`run_all.py` 默认只跑 `08 --evaluate`，不会覆盖 `sdp_audit_100.csv`；如需重新抽样，使用 `python run_all.py --with-audit-sample`，之后必须重新标注。

### 12.3 关键随机性

- TTR shuf：`random.seed(42)`，M=10。
- 审计抽样：`random.Random(42)`，A/B 各 50，B 层按 relation 配额。
- SDP 无随机性（spaCy 确定性解析）。
- BH 校正为确定性。

---

## 13. 局限与写作注意（重要）

1. **9,030 vs 7,868**：阶段二候选含 1,162 条非水稻；SDP/富集/网络用 7,868 条水稻-only。论文中必须写清楚，否则“水稻知识图谱”口径不成立。
2. **884 条因果 ≠ 884 条正确**：审计 precision 86.0%，估计 TP 约 760；不要把 884 写成“高置信 884 条全部正确”。
3. **Recall 仅 28.0%**：SDP 是 precision-first 过滤器，漏掉大量名词语义/轻动词表达；应把 1,379 条 association 作为高覆盖补充层，而不是宣称 SDP 全面优于共现。
4. **100 条审计是分层样本，不是全量随机**：A 层过采样因果，B 层过采样非因果；precision/recall 是分层外推，需在 Methods 写明。
5. **人工标签由作者/Agent 预标注**：`sdp_audit_labels.json` 可直接复核；写作者应至少抽查表 10 中的 10 条，避免把标签当“金标准”不加说明。
6. **RiceMind 注释缺口候选不是新基因**：13 个候选是“数据库未覆盖/未注释抗逆关联”，不是首次发现；需人工审编与实验验证。
7. **5 篇未索引文献**：未进入主统计；`data/unindexed_pmids.tsv` 提供题录/摘要，写限制时可用。
8. **实体同形**：`RGA`/`D1`、`SAMs`/`SAMS1`、`Apo`/`APO`、`ha-1`/`AHA8` 是已知误差源；`A01` 是典型错误案例。
9. **富集选择偏倚**：core14 的 68 个显著条目来自高频抗旱基因，存在选择偏倚；建议主文用 Top100 的 53 个条目，core14 作补充。
10. **网络边是无向关系**：度中心性不等于因果方向；SDP 方向已给出，但网络图未画箭头。

---

## 14. 建议的论文结果小节映射（只给事实，不写正文）

| 论文章节 | 可直接用的表/图 | 核心事实 |
|---|---|---|
| 数据与语料 | `table1`、`table2` | 3,306 篇/35,345 句/795,401 词元；TTR 4.31%；shuf 高于 head 1.6–2.4 pp |
| 检索策略 | `table3`、`table6` | `drought stress` 覆盖 1,253 篇 vs `drought tolerance` 883 篇（+41.9%） |
| 实体归一化 | `table4`、`table5`、Fig. 2–3 | 42,617 提及；10,708 RAP vs PubTator 473（22.6×）；1,344 基因/190 性状/4,678 对 |
| SDP 方法 | Fig. 4 | 7,868 → 884 因果 + 1,379 association；48 条 respectively 排除 |
| SDP 评测 | Fig. 5a、`table10`、`table11` | precision 86.0%（95% CI 73.8–93.0）；基线 34.5%；recall 28.0%；F1 42.2%；极性 88.2% |
| 富集 | Fig. 5b、`table12` | Top100：53 个 q<0.05；ABA sensitivity 10.68×；stomatal resistance 5.54× |
| RiceMind | `table13`、`table15` | 覆盖 99%；85 有 stress 注释；13 提 drought；13 个注释缺口候选 |
| 网络 | Fig. 6、`table14` | 457 节点/693 边；drought tolerance 加权度 286；DSR2/ONAC2/WOX11 高连接 |

---

## 15. 版本与文件位置

- 本报告：`RESULTS_REPORT.md`
- 关键数字：`results/key_numbers.json`
- 一键复现：`run_all.py`、`results/run_all.log`
- 结果表：`results/tables/`
- 图：`results/figures/`（PDF/SVG/PNG）
- 论文草稿（非本次交付重点）：`paper_bibm/BIBM-Latex/main.tex`、`main.pdf`（7 页，已启用模板左下角 AI 页脚；后续写作者可参考排版与数据引用，不必沿用文字）

---

## 附：最小可用事实清单（写作时可直接复制）

- 3,311 PubMed hits；3,306 PubTator docs；35,345 sentences；795,401 tokens；34,247 types；raw TTR 4.31%。
- 47,157 PubTator entities；47,156 offset-verified（99.998%）。
- 42,617 entity mentions = 10,708 RAP genes + 29,966 traits + 1,943 PubTator-only non-rice genes.
- Oryzabase 22,023 records → 102,860 aliases；RTO 2,022 terms → 3,904 synonyms.
- Rice-only candidates 7,868 → 4,678 unique pairs → 1,344 genes / 190 traits.
- SDP 7,868 → 884 causal (575 positive, 150 negative, 159 neutral) + 1,379 association；693 unique causal pairs；385 genes；414 PMIDs；647 sentences.
- Audit 100 (50+50, seed 42): precision 86.0% (95% CI 73.8–93.0), baseline 34.5%, recall 28.0%, F1 42.2%, polarity 88.2% (30/34).
- Error counts: 7 false positives, 14 false negatives, 4 polarity conflicts.
- Enrichment Top100: 151 hit, 53 q<0.05; ABA sensitivity 10.68×; stomatal resistance 5.54×; grain size 8.87×.
- RiceMind Top100: 99 found (99%), 85 stress-annotated, 13 drought-annotated, 13 annotation-gap candidates.
- Network: 457 nodes, 693 edges; drought tolerance weighted degree 286; top genes DSR2 (15), ONAC2 (10), WOX11/HDA704/OsMYB4 (9).