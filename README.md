# BioNLP — 基于 PubMed 文献挖掘与水稻性状本体的水稻抗旱基因–性状调控关系发现

> 课程：《生物文本挖掘与知识发现》（华中农业大学 信息学院）
> 作者：易子茗（2024317240212，生信 2402）
> 论文终稿：[`paper_final/BIBM-Latex/main.pdf`](paper_final/BIBM-Latex/main.pdf)

本仓库是课程论文的**全部代码、结果与论文源文件**，对应论文中的
Data and Code Availability 一节。

---

## 1. 这个项目在做什么

水稻抗旱的遗传证据分散在几万篇 PubMed 摘要里，通用生物医学文本挖掘管线
（BioGRID / RelEx 一类）直接迁移到作物上会遇到三个具体的坑：

1. **基因符号体系不同**：医学 NER 词典里没有 RAP-DB / Oryzabase 的水稻基因位点，
   PubTator 原生映射只覆盖 473 处提及；
2. **表型词表缺失**：农艺性状（分蘖数、结实率、气孔导度…）不在通用医学本体里，
   必须换用水稻性状本体 RTO/TO；
3. **同形词冲突严重**：`ha-1`（产量单位）与 `Ha-1`（基因）、`Apo`（基因）与普通英文词、
   `red`/`dull`/`lazy` 这类既是性状词又是基因别名的词，会造成假阳性爆炸。

本项目搭建了一条从 PubMed 检索到知识网络验证的完整管线，并在每一步给出可审计的数字。

```
PubMed 检索 (3,311) ─► PubTator 标注 (3,306) ─► NLTK 分句分词 (35,345 句 / 795,401 词元)
        │
        ├─► 词汇统计：TTR / Heaps 定律 / TF-IDF 二元短语
        │
        ├─► 农学五重消歧 + 双词典归一化（Oryzabase 基因 + RTO 性状）
        │        └─► 42,617 处实体提及 = 10,708 RAP 基因 + 29,966 性状 + 1,943 非水稻基因
        │
        ├─► 句子级共现：9,030 条候选 ─► 水稻-only 过滤 ─► 7,868 条（4,678 唯一对）
        │
        ├─► spaCy 依存树 + 最短依存路径（SDP）：884 条因果候选 + 1,379 条 association
        │        └─► 100 条分层人工审计：precision 86.0%，polarity 88.2%
        │
        └─► 超几何富集 / RiceMind 交叉验证 / 基因–性状网络（457 节点 693 边）
```

## 2. 目录结构

```
.
├─ run_all.py                     # 一键复现（阶段一~四 + 结果包导出）
├─ requirements.txt
├─ RESULTS_REPORT.md              # 主报告：全部实验数字、方法、表格、图、局限
├─ scripts/
│  ├─ 01_fetch_pubmed_and_pubtator.py   # PubMed esearch + PubTator 批量下载（断点续传）
│  ├─ 02_lexical_stats_ttr_tfidf.py     # TTR / Heaps / TF-IDF
│  ├─ 03_build_oryza_and_rto_dicts.py   # Oryzabase + RTO 词典构建与五重消歧
│  ├─ 04_extract_sentence_evidence.py   # 句子级实体与共现证据
│  ├─ 05_generate_visual_figures.py     # 阶段一/二图
│  ├─ 06_generate_bibm_nature_figures.py
│  ├─ 07_sdp_relation_extraction.py     # spaCy 依存树 + SDP 因果/极性抽取
│  ├─ 08_sdp_audit_eval.py              # 100 条分层审计评测
│  ├─ 09_enrichment_validation.py       # TO/RTO 超几何富集
│  ├─ 10_ricemind_novelty.py            # RiceMind 交叉验证 + 网络 + 注释缺口
│  ├─ 11_generate_part34_figures.py
│  └─ 12_export_results_package.py      # 导出表格 / key_numbers / 数据清单
├─ data/
│  ├─ pmid_raw.xml / pmid_list_all.txt  # 检索结果（3,311 PMID）
│  ├─ result_Pubtator_all.txt           # PubTator 标注全文（3,306 篇，约 7.5 MB）
│  ├─ result_Pubtator_100.txt           # 100 篇子集（调试用）
│  └─ unindexed_pmids.tsv               # 5 篇尚未被 PubTator 索引的文献题录
├─ dicts/
│  ├─ OryzabaseGeneListEn.txt           # Oryzabase 基因命名词典原始文件
│  └─ RTO-1.0.obo                       # 水稻性状本体 Rice Trait Ontology
├─ results/
│  ├─ key_numbers.json                  # 机器可读关键数字
│  ├─ README.md / DATA_DICTIONARY.md     # 结果索引与字段字典
│  ├─ tables/                           # 18 张发表级 CSV
│  ├─ figures/                          # 图（PDF / SVG / PNG）
│  ├─ rice_drought_sentence_evidence.tsv  # 阶段二 9,030 条候选
│  ├─ sdp_relations.tsv                 # 阶段三 7,868 条 SDP 判定
│  ├─ sdp_audit_100.csv + sdp_audit_labels.json  # 100 条人工审计
│  └─ part4_*.csv/json                  # 富集 / RiceMind / 网络 / 缺口候选
├─ paper_final/
│  ├─ regenerate_final_figures.py       # 论文终稿 6 张图（viridis 配色）
│  └─ BIBM-Latex/                       # 论文终稿（IEEEtran 模板）
│     ├─ main.tex / main.pdf            # 10 页终稿
│     ├─ reference.bib / IEEEtran.cls
│     └─ figure/                        # 论文内嵌矢量图
└─ paper_bibm/BIBM-Latex/               # 阶段二参考草稿（历史版本，仅供追溯）
```

## 3. 快速开始

### 3.1 环境

```powershell
pip install -r requirements.txt
python -m spacy download en_core_web_sm
```

`requirements.txt` 覆盖 `nltk`、`spacy`、`networkx`、`numpy`、`scipy`、`matplotlib`、`pymupdf`。

### 3.2 复现全部结果

仓库已包含语料与词典，**不需要重新下载**：

```powershell
python run_all.py
```

本机实测端到端约 **1.5 分钟**（复用 `data/result_Pubtator_all.txt`）。
产出写入 `results/`，日志写入 `results/run_all.log`。

常用开关：

```powershell
python run_all.py --skip-figures        # 跳过绘图，只算数字
python run_all.py --with-audit-sample   # 重新抽 100 条审计样本（会覆盖已冻结的标签表，慎用）
```

### 3.3 从零重新抓取文献（可选）

完整下载 3,306 篇 PubTator 标注按课件要求限速 `sleep 5.8s`，约需 **5.3 小时**。
单篇缓存写入 `data/pubtator_cache/`，支持断点续传。

```powershell
python scripts/01_fetch_pubmed_and_pubtator.py
```

### 3.4 编译论文

需要 TeX Live 2026（含 `xelatex`、`bibtex`）：

```powershell
cd paper_final/BIBM-Latex
xelatex main.tex
bibtex main
xelatex main.tex
xelatex main.tex
```

## 4. 关键结果

| 环节 | 结果 |
|---|---|
| 语料 | 3,311 篇命中 → 3,306 篇 PubTator 解析；35,345 句 / 795,401 词元 / 34,247 词型；raw TTR 4.31% |
| 实体 | 47,157 PubTator 原始实体，47,156 偏移校验通过（99.998%）；42,617 处提及 |
| 归一化 | 基因提及较 PubTator 原生 RAP 映射（473）提升 **22.6 倍**；Oryzabase 22,023 条记录 → 102,860 条别名 |
| 共现 | 9,030 条原始候选 → 剔除 1,162 条外源基因 → **7,868 条水稻候选**（4,678 唯一对 / 1,344 基因 / 190 性状） |
| SDP | 884 条因果候选（575 正向 / 150 负向 / 159 中性）+ 1,379 条 association；693 个唯一因果对 / 285 基因 / 414 PMID / 947 句 |
| 审计 | 100 条分层抽样（50+50，seed 42）：precision **86.0%**（95% CI 73.8–93.0），baseline 34.5%，recall 28.0%，F1 42.2%，polarity 88.2% |
| 富集 | Top100 中 151 个条目命中、53 个 q<0.05；ABA sensitivity 10.68×、stomatal resistance 5.54×、grain size 8.87× |
| 验证 | RiceMind Top100 覆盖 99%，85 个有 stress 注释、13 个有 drought 注释、13 个注释缺口候选 |
| 网络 | 457 节点 / 693 边；drought tolerance 加权度 286；高连接基因 DSR2(15)、ONAC2(10)、WOX11/HDA704/OsMYB4(9) |

**口径提醒**（写作与引用时务必区分）：

- `9,030` = 阶段二全部候选（含 1,162 条非水稻 PubTator-only 基因）；
- `7,868` = 水稻-only 候选（`RAP_ID` 以 `Os` 开头），SDP / 富集 / 网络均用这个口径；
- `884` = SDP **保留的候选**因果，不是 884 条全部正确；按 86.0% precision 估计约 760 条为真实调控；
- `13` = RiceMind **注释缺口候选**，不是“新基因”。

完整数字与图表见 [`RESULTS_REPORT.md`](RESULTS_REPORT.md)，机器可读版本见
[`results/key_numbers.json`](results/key_numbers.json)。

## 5. 数据与许可

| 资源 | 来源 | 说明 |
|---|---|---|
| 文献标注 | [PubTator](https://www.ncbi.nlm.nih.gov/research/pubtator/) | 公开 API，遵守 `sleep 5.8s` 限速 |
| 基因词典 | [Oryzabase](https://shigen.nig.ac.jp/rice/oryzabase/) | `dicts/OryzabaseGeneListEn.txt` |
| 性状本体 | Rice Trait Ontology (RTO-1.0) | `dicts/RTO-1.0.obo` |
| 外部验证 | [RiceMind](http://ricemind.gao-lab.org/) | 交叉验证与注释缺口分析 |

代码以 MIT 许可发布，供课程学习与后续研究使用；第三方词典与文献标注的版权归各自来源所有。

## 6. 已知限制

- 静态词典 + 规则匹配对未见别名与跨句指代无能为力；
- SDP 精确率优先，代价是 recall 仅 28.0%（1,379 条 association 留作补充检索池，未计入主干）；
- 5 篇 2026 年新文献尚未进入 PubTator 索引，已单列于 `data/unindexed_pmids.tsv`，不参与主统计；
- 富集与网络结论建立在 3,306 篇干旱相关摘要之上，属于文献层面的证据聚合，不等同于实验验证。
