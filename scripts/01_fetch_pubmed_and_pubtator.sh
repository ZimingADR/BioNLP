#!/bin/bash
# ==============================================================================
# 课程：生物自然语言处理与知识发现 (BioNLP)
# 对应课件：JingboSlides_8.PubMedTermsNERandShellProgramming.pdf (Page 11-19 & Page 31)
# 功能：
#   1. 使用 E-utilities (esearch) 检索水稻抗旱文献 PMID 列表
#   2. 严格按照课件 Page 31 使用 while read 循环 + sleep 5.8s 从 PubTator 下载标注
# ==============================================================================

mkdir -p data

QUERY='("Oryza sativa"[Title/Abstract] OR "Oryza sativa"[MeSH Terms] OR "rice"[Title/Abstract]) AND "drought"[Title/Abstract]'
F_XML="data/pmid_raw.xml"
F_ALL_LIST="data/pmid_list_all.txt"
F_LIST="data/pmid_list_100.txt"
F_OUT="data/result_Pubtator.txt"

echo "Step 1: Fetching PMIDs from PubMed E-utilities..."
curl -s -G "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi" \
    --data-urlencode "db=pubmed" \
    --data-urlencode "term=$QUERY" \
    --data-urlencode "retmax=9999" > "$F_XML"

# 课件 Page 18: 使用 grep 和 sed 提取纯 PMID 列表
cat "$F_XML" | grep "<Id>" | sed 's/<Id>//g' | sed 's/<\/Id>//g' | sed 's/[[:blank:]]//g' > "$F_ALL_LIST"

# 截取前 100 篇用于首轮测试
head -n 100 "$F_ALL_LIST" > "$F_LIST"

# 课件 Page 31: 循环请求 PubTator API 并严格执行 sleep 5.8s 防止校园网 IP 被封
echo -e "\n Let us curate the results!\n"
echo -e "\n" > "$F_OUT"

total=$(wc -l < "$F_LIST")
i=1
while IFS= read -r line
do
    url="https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/pubtator?pmids=$line"
    curl -s "$url" >> "$F_OUT"
    printf "%d-th result out of %d is processing (PMID: %s)...\n" "$i" "$total" "$line"
    i=$((i+1))
    sleep 5.8s
done < "$F_LIST"

echo "Done! Saved to $F_OUT"
