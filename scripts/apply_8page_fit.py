import re

with open('paper_bibm/BIBM-Latex/main.tex', 'r', encoding='utf-8') as f:
    tex = f.read()

# 1. Figure widths to 0.75\textwidth
tex = tex.replace('width=0.82\\textwidth]{figure/fig1_ttr_heaps.pdf}', 'width=0.74\\textwidth]{figure/fig1_ttr_heaps.pdf}')
tex = tex.replace('width=0.82\\textwidth]{figure/fig2_tfidf_terms.pdf}', 'width=0.74\\textwidth]{figure/fig2_tfidf_terms.pdf}')
tex = tex.replace('width=0.82\\textwidth]{figure/fig3_pipeline_ablation.pdf}', 'width=0.74\\textwidth]{figure/fig3_pipeline_ablation.pdf}')
tex = tex.replace('width=0.82\\textwidth]{figure/fig4_top_genes_traits.pdf}', 'width=0.74\\textwidth]{figure/fig4_top_genes_traits.pdf}')

# 2. Quote to compact format
q_pattern = r'\\begin\{quote\}\s*\\small\\ttfamily\s*\("Oryza sativa".*?\)\s*AND\s*"drought"\[Title/Abstract\]\s*\\end\{quote\}'
new_q = r'{\\centering\\vspace{-2pt}\\small\\ttfamily ("Oryza sativa"[Title/Abstract] OR "Oryza sativa"[MeSH Terms] OR "rice"[Title/Abstract]) AND "drought"[Title/Abstract]\\par\\vspace{-2pt}}'
tex = re.sub(q_pattern, new_q, tex, flags=re.DOTALL)

# 3. Tighten list spacing
tex = tex.replace(
    r'\begin{enumerate}\setlength{\itemsep}{1pt}\setlength{\parsep}{0pt}\setlength{\topsep}{2pt}',
    r'\begin{enumerate}\setlength{\itemsep}{0pt}\setlength{\parsep}{0pt}\setlength{\topsep}{1pt}'
)
tex = tex.replace(
    r'\begin{itemize}\setlength{\itemsep}{1pt}\setlength{\parsep}{0pt}\setlength{\topsep}{1pt}',
    r'\begin{itemize}\setlength{\itemsep}{0pt}\setlength{\parsep}{0pt}\setlength{\topsep}{0.5pt}'
)

# 4. In references, trigger balance at ref 5
tex = tex.replace(r'\IEEEtriggeratref{6}', r'\IEEEtriggeratref{5}')

with open('paper_bibm/BIBM-Latex/main.tex', 'w', encoding='utf-8') as f:
    f.write(tex)

print('Applied 8-page adjustments.')
