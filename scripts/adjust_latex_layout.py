import re

with open('paper_bibm/BIBM-Latex/main.tex', 'r', encoding='utf-8') as f:
    text = f.read()

# 1. Compact quote in III-A
quote_pattern = r'\\begin\{quote\}\s*\\small\\ttfamily\s*\("Oryza sativa".*?\)\s*AND\s*"drought"\[Title/Abstract\]\s*\\end\{quote\}'
new_quote = r'{\\small\\ttfamily ("Oryza sativa"[Title/Abstract] OR "Oryza sativa"[MeSH Terms] OR "rice"[Title/Abstract]) AND "drought"[Title/Abstract]}'
text = re.sub(quote_pattern, new_quote, text, flags=re.DOTALL)

# 2. Table fonts to 8pt for professional IEEE table look
text = text.replace(
    '\\begin{table*}[t]\n  \\centering\n  \\footnotesize\n  \\renewcommand{\\arraystretch}{0.92}\n  \\caption{全量文献挖掘出的',
    '\\begin{table*}[t]\n  \\centering\n  \\fontsize{8pt}{9.6pt}\\selectfont\n  \\renewcommand{\\arraystretch}{0.88}\n  \\caption{全量文献挖掘出的'
)
text = text.replace(
    '\\begin{table*}[t]\n  \\centering\n  \\footnotesize\n  \\renewcommand{\\arraystretch}{0.92}\n  \\caption{符合讲义验证规范',
    '\\begin{table*}[t]\n  \\centering\n  \\fontsize{8pt}{9.6pt}\\selectfont\n  \\renewcommand{\\arraystretch}{0.88}\n  \\caption{符合讲义验证规范'
)

# 3. Compact lists
text = text.replace(
    r'\begin{enumerate}\setlength{\itemsep}{2pt}\setlength{\parsep}{0pt}',
    r'\begin{enumerate}\setlength{\itemsep}{0pt}\setlength{\parsep}{0pt}\setlength{\topsep}{1pt}\setlength{\partopsep}{0pt}'
)
text = text.replace(
    r'\begin{itemize}\setlength{\itemsep}{1pt}\setlength{\parsep}{0pt}',
    r'\begin{itemize}\setlength{\itemsep}{0pt}\setlength{\parsep}{0pt}\setlength{\topsep}{1pt}\setlength{\partopsep}{0pt}'
)

with open('paper_bibm/BIBM-Latex/main.tex', 'w', encoding='utf-8') as f:
    f.write(text)

print('Adjustments completed successfully.')
