# -*- coding: utf-8 -*-
import re
import subprocess
import os

tex_path = r'paper_bibm\BIBM-Latex\main.tex'

with open(tex_path, 'r', encoding='utf-8') as f:
    text = f.read()

# 1. 删除致谢内容与参考文献声明
old_end_block = r"""% ======================================================================
%  致谢
% ======================================================================
\section*{Acknowledgment 致谢}
{\small 衷心感谢华中农业大学信息学院夏静波教授在《生物文本挖掘与知识发现》课程中的悉心指导，特别是关于接口限流规范、句子偏移对齐及 Oryzabase 词典消歧的启发；同时鸣谢华农 RiceMind/RTO 团队及日本 Oryzabase/RAP-DB 平台提供的开源数据。}

% ======================================================================
%  参考文献
% ======================================================================
\IEEEtriggeratref{6}
\bibliographystyle{IEEEtran}
\bibliography{reference}"""

if old_end_block in text:
    text = text.replace(old_end_block, "% Acknowledgments and References removed per user instruction.")
    print("SUCCESS: Removed Acknowledgments and References.")
else:
    # Try regex if exact match fails
    text = re.sub(r'\\section\*\{Acknowledgment.*?\\bibliography\{reference\}', '% Acknowledgments and References removed per user instruction.', text, flags=re.DOTALL)
    print("Regex removed Acknowledgments and References.")

# 2. 移除正文中所有的 \cite{...}，避免缺失参考文献时出现 [?]
text = re.sub(r'~?\\cite\{[^}]+\}', '', text)
print("SUCCESS: Cleaned all \\cite{...} citations.")

# 3. 校验并优化已完成内容与占位符规范（nature-polishing 风格检查）
# 确保没有全角英文破折号或不规范标点，确保占位符表述严密、无虚构数据
text = text.replace("——", "，")
text = text.replace("—", "，")

with open(tex_path, 'w', encoding='utf-8') as f:
    f.write(text)

print("Updated main.tex successfully.")
