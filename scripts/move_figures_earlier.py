import re

with open('paper_bibm/BIBM-Latex/main.tex', 'r', encoding='utf-8') as f:
    tex = f.read()

# 1. Extract Figure 1 and Figure 2
fig1_match = re.search(r'(\\begin\{figure\*\}\[t\].*?\\includegraphics.*?fig1_ttr_heaps\.pdf.*?\\end\{figure\*\})', tex, re.DOTALL)
fig2_match = re.search(r'(\\begin\{figure\*\}\[t\].*?\\includegraphics.*?fig2_tfidf_terms\.pdf.*?\\end\{figure\*\})', tex, re.DOTALL)

fig1_code = fig1_match.group(1)
fig2_code = fig2_match.group(1)

# Remove fig1 and fig2 from their current location
tex_clean = tex.replace(fig1_code, '').replace(fig2_code, '')

# Move Figure 1 and Figure 2 to the beginning of Section III (around line 140)
# Right before \subsection{数据获取与防封禁限流策略 Data Acquisition}
target_str = r'\subsection{数据获取与防封禁限流策略 Data Acquisition}'
tex_new = tex_clean.replace(target_str, fig1_code + '\n\n' + fig2_code + '\n\n' + target_str)

# Also let's check Table III position: keep Table II and Table III where they make sense
with open('paper_bibm/BIBM-Latex/main.tex', 'w', encoding='utf-8') as f:
    f.write(tex_new)

print('Moved Fig 1 and Fig 2 to Section III.')
