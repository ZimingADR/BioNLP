import re

with open('paper_bibm/BIBM-Latex/main.tex', 'r', encoding='utf-8') as f:
    tex = f.read()

# Ensure dbltop parameters in preamble
if 'dbltopnumber' not in tex:
    tex = tex.replace(
        '\\begin{document}',
        '\\setcounter{dbltopnumber}{2}\n\\renewcommand{\\dbltopfraction}{0.9}\n\\renewcommand{\\dblfloatpagefraction}{0.8}\n\\begin{document}'
    )

# 1. Clean up any existing Fig 1 and Fig 2
fig1_match = re.search(r'(\\begin\{figure\*\}\[t\].*?\\includegraphics.*?fig1_ttr_heaps\.pdf.*?\\end\{figure\*\})', tex, re.DOTALL)
fig2_match = re.search(r'(\\begin\{figure\*\}\[t\].*?\\includegraphics.*?fig2_tfidf_terms\.pdf.*?\\end\{figure\*\})', tex, re.DOTALL)
fig1_code = fig1_match.group(1)
fig2_code = fig2_match.group(1)

tex = tex.replace(fig1_code, '').replace(fig2_code, '')

# Place Fig 1 and Fig 2 at the beginning of Section IV: Results
results_hdr = r'\section{结果 Results}'
tex = tex.replace(results_hdr, results_hdr + '\n\n' + fig1_code + '\n\n' + fig2_code)

# 2. Extract Table II and Table III
tab2_match = re.search(r'(\\begin\{table\*\}\[t\].*?\\caption\{全量文献挖掘出的.*?\\end\{table\*\})', tex, re.DOTALL)
tab3_match = re.search(r'(\\begin\{table\*\}\[t\].*?\\caption\{符合讲义验证规范.*?\\end\{table\*\})', tex, re.DOTALL)
tab2_code = tab2_match.group(1)
tab3_code = tab3_match.group(1)

tex = tex.replace(tab2_code, '').replace(tab3_code, '')

# Create optimized, compact Table II (p{6.6cm} for col 5 so it wraps in 2 lines, fontsize 7.8pt/9.2pt)
optimized_tab2 = r'''\begin{table*}[t]
  \centering
  \fontsize{7.8pt}{9.2pt}\selectfont
  \renewcommand{\arraystretch}{0.85}
  \setlength{\tabcolsep}{3.5pt}
  \caption{全量文献挖掘出的 Top 10 水稻抗旱核心基因及其多源异构别名向唯一 RAP-DB ID 的自动汇聚实证}
  \label{tab:top_genes}
  \begin{tabularx}{\textwidth}{c l l c p{6.8cm} X}
    \toprule
    排名 & 标准基因符号 & RAP-DB 位点编号 & 提及数 & 语料中自动归并的异构文献写法 (提及频次) & 基因全称与抗旱分子调控机制 \\
    \midrule
    1 & \textbf{DREB1A} & \texttt{Os09g0522200} & 103 & \texttt{DREB1A} (37), \texttt{OsDREB1A} (31), \texttt{CBF3} (24), \texttt{DREB1a} (9), \texttt{Cbf3} (2) & AP2/ERF 转录因子，结合 DRE 顺式元件激活抗逆基因表达 \\
    2 & \textbf{NCED3} & \texttt{Os07g0154100} & 83 & \texttt{OsNCED3} (36), \texttt{OsNCED4} (27), \texttt{NCED4} (10), \texttt{NCED3} (8) & 9-顺式-环氧类胡萝卜素双加氧酶，ABA 生物合成关键限速酶 \\
    3 & \textbf{BZIP23} & \texttt{Os02g0766700} & 80 & \texttt{OsbZIP23} (70), \texttt{osbZIP23} (4), \texttt{OsZIP23} (2), \texttt{OsAREB2} (2), \texttt{bZIP23} (1) & 华农克隆代表性抗旱 bZIP 转录因子，正向调控 ABA 敏感性与抗旱性 \\
    4 & \textbf{ONAC2} & \texttt{Os03g0815100} & 70 & \texttt{SNAC1} (51), \texttt{OsNAC9} (12), \texttt{OsSNAC1} (6), \texttt{ONAC2} (1) & 华农明星抗旱基因 \texttt{SNAC1}，诱导气孔关闭并改良根系构型提高抗旱产量 \\
    5 & \textbf{NAC6} & \texttt{Os01g0884300} & 70 & \texttt{OsNAC6} (40), \texttt{OsNAC2} (14), \texttt{SNAC2} (10), \texttt{OsSNAC2} (3), \texttt{NAC6} (3) & NAC 家族转录因子，激活过氧化物酶与根系发育通路增强耐旱性 \\
    6 & \textbf{BZIP46} & \texttt{Os06g0211200} & 68 & \texttt{OsbZIP46} (41), \texttt{OsABF2} (8), \texttt{OsbZIP46CA1} (7), \texttt{ABF3} (5), \texttt{OsAREB8} (2) & 华农克隆 AREB/ABF 亚家族因子，受干旱与 ABA 诱导正调控抗旱响应 \\
    7 & \textbf{DREB1B} & \texttt{Os09g0522000} & 67 & \texttt{OsDREB1B} (38), \texttt{DREB1B} (12), \texttt{CBF1} (9), \texttt{OsDREB1C} (7) & 非生物逆境诱导型 DREB 转录因子，提升植株脱水耐受性 \\
    8 & \textbf{PIP2A} & \texttt{Os07g0448800} & 65 & \texttt{OsPIP2} (49), \texttt{PIP2} (12), \texttt{OsPIP2a} (3), \texttt{OsPIP2-1} (1) & 质膜水通道蛋白，调控根系水力导度与干旱胁迫下细胞水分平衡 \\
    9 & \textbf{LEA3-1} & \texttt{Os05g0542500} & 62 & \texttt{OsLEA3} (34), \texttt{OsLEA3-1} (11), \texttt{LEA3} (7), \texttt{Lea3} (4), \texttt{OsLEA3-A} (2) & 华农克隆胚胎晚期丰富蛋白基因，过表达显著保障大田干旱条件下结实率 \\
    10 & \textbf{P5CS} & \texttt{Os05g0455500} & 59 & \texttt{P5CS} (33), \texttt{OsP5CS} (12), \texttt{OsP5CS1} (9), \texttt{P5CS1} (5) & $\Delta^1$-吡咯琳-5-羧酸合成酶，催化游离脯氨酸积累以维持细胞渗透势 \\
    \bottomrule
  \end{tabularx}
\end{table*}'''

# Create optimized, compact Table III (fontsize 7.8pt/9.2pt)
optimized_tab3 = r'''\begin{table*}[t]
  \centering
  \fontsize{7.8pt}{9.2pt}\selectfont
  \renewcommand{\arraystretch}{0.85}
  \setlength{\tabcolsep}{3.5pt}
  \caption{符合讲义验证规范的真实句子级证据输出样例与典型边界错误案例对比表}
  \label{tab:examples}
  \begin{tabularx}{\textwidth}{p{2.1cm} p{1.6cm} p{3.2cm} p{3.4cm} X}
    \toprule
    案例类型 & PMID (句号) & 原文提及 (偏移区间) & 归一化结果 (标准名 \& 唯一 ID) & 原始英文证据句 (高亮目标实体) \\
    \midrule
    \textbf{正面样例 1}\newline (无 Os 前缀蛋白归并) & \texttt{42824614}\newline (\texttt{sent\_id=8}) & Gene: \texttt{HSP18.0} $[1190,1197)$\newline Trait: \texttt{drought resistance} $[1238,1256)$ & Gene: \textbf{HSP18.0} (\texttt{Os03g0267000})\newline Trait: \textbf{drought tolerance} (\texttt{TO:0000276}) & \textit{``The results establish an OsMPK16-\textbf{HSP18.0} module that positively regulates \textbf{drought resistance} in rice.''} \\
    \midrule
    \textbf{正面样例 2}\newline (下划线缺省符号恢复) & \texttt{42731987}\newline (\texttt{sent\_id=9}) & Gene: \texttt{OsSET37} $[943,950)$\newline Trait: \texttt{drought tolerance} $[999,1016)$ & Gene: \textbf{OsSET37} (\texttt{Os02g0811600})\newline Trait: \textbf{drought tolerance} (\texttt{TO:0000276}) & \textit{``Overall, our findings demonstrate that \textbf{OsSET37} negatively regulates \textbf{drought tolerance} in rice.''} \\
    \midrule
    \textbf{局限案例 1}\newline (组织器官缩写同形异义) & \texttt{42739424}\newline (\texttt{sent\_id=2}) & Gene?: \texttt{SAMs} $[188,192)$\newline Trait: \texttt{drought} $[206,213)$ & 误判为 \textbf{SAMS1} (\texttt{Os05g0135700})\newline 实际语义: \textit{shoot apical meristems} & \textit{``Here, we Performed single-nucleus RNA sequencing on rice shoot apical meristems (\textbf{SAMs}) under \textbf{drought} conditions.''} \\
    \midrule
    \textbf{局限案例 2}\newline (并列分别调控共现误判) & \texttt{42768028}\newline (\texttt{sent\_id=6}) & Gene: \texttt{OsBIHD1} $[954,961)$\newline Trait: \texttt{drought tolerance} $[978,995)$ & Gene: \textbf{BIHD1} (\texttt{Os03g0680800})\newline 实际调控: \textit{blast defense} (稻瘟病) & \textit{``...by regulating \textbf{OsbHLH148} and \textbf{OsBIHD1}, which regulate \textbf{drought tolerance} and blast defense, respectively.''} \\
    \bottomrule
  \end{tabularx}
\end{table*}'''

# Place Table II and Table III right before Section IV-B text (after fig 4)
fig4_end = tex.find(r'\label{fig:top_entities}')
fig4_close = tex.find(r'\end{figure*}', fig4_end) + len(r'\end{figure*}')
tex = tex[:fig4_close] + '\n\n' + optimized_tab2 + '\n\n' + optimized_tab3 + tex[fig4_close:]

with open('paper_bibm/BIBM-Latex/main.tex', 'w', encoding='utf-8') as f:
    f.write(tex)

print('Updated main.tex with compact tables and float positions.')
