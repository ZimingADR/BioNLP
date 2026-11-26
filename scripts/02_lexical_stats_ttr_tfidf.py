# -*- coding: utf-8 -*-
"""
Step 2 (Part 1 核心): 水稻抗旱文献的词汇计算、TTR 词汇丰富度分析与 TF-IDF 核心术语挖掘
对应课件:
  - Ch2_词汇计算基础_Slide.pdf & 实验1_TTR / 实验1_TFIDF
  - JingboSlides_8.PubMedTermsNERandShellProgramming.pdf (Page 38-49: NLTK 分句与分词)
输出:
  - results/part1_lexical_and_tfidf_summary.json
  - results/part1_tfidf_top_terms.tsv
"""

import json
import math
import random
import re
from collections import Counter
from pathlib import Path

import nltk
from nltk.tokenize import sent_tokenize, word_tokenize

# 确保 NLTK 所需资源可用
for res in ["punkt_tab", "stopwords"]:
    try:
        nltk.data.find(f"tokenizers/{res}" if "punkt" in res else f"corpora/{res}")
    except LookupError:
        nltk.download(res, quiet=True)

from nltk.corpus import stopwords

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
PUBTATOR_ALL = DATA_DIR / "result_Pubtator_all.txt"
PUBTATOR_100 = DATA_DIR / "result_Pubtator_100.txt"
PUBTATOR_FILE = PUBTATOR_ALL if (PUBTATOR_ALL.exists() and PUBTATOR_ALL.stat().st_size > 1000) else PUBTATOR_100


def parse_pubtator_documents(pubtator_path: Path = None):
    """
    按照课件 Page 42 的 parse_pubtator 规范解析多篇 PubTator 文献记录。
    同时执行自动校验: text[begin:end] == mention
    """
    if pubtator_path is None:
        pubtator_path = PUBTATOR_ALL if (PUBTATOR_ALL.exists() and PUBTATOR_ALL.stat().st_size > 1000) else PUBTATOR_100
    raw = pubtator_path.read_text(encoding="utf-8", errors="replace")
    docs = []
    current_pmid = None
    title = ""
    abstract = ""
    entities = []

    def flush_doc():
        nonlocal current_pmid, title, abstract, entities
        if current_pmid and (title or abstract):
            full_text = title + (" " if title and abstract else "") + abstract
            verified_cnt = 0
            for ent in entities:
                b, e = ent["begin"], ent["end"]
                slice_str = full_text[b:e]
                ent["offset_verified"] = (slice_str == ent["mention"])
                if ent["offset_verified"]:
                    verified_cnt += 1
            docs.append({
                "pmid": current_pmid,
                "title": title,
                "abstract": abstract,
                "text": full_text,
                "entities": entities,
                "verified_entities": verified_cnt,
            })
        current_pmid = None
        title = ""
        abstract = ""
        entities = []

    for line in raw.splitlines():
        line_str = line.strip()
        if not line_str:
            flush_doc()
            continue
        if "|t|" in line_str:
            pmid, t = line_str.split("|t|", 1)
            if current_pmid and current_pmid != pmid:
                flush_doc()
            current_pmid = pmid
            title = t
        elif "|a|" in line_str:
            pmid, a = line_str.split("|a|", 1)
            current_pmid = pmid
            abstract = a
        else:
            parts = line_str.split("\t")
            if len(parts) >= 5 and parts[1].isdigit() and parts[2].isdigit():
                entities.append({
                    "pmid": parts[0],
                    "begin": int(parts[1]),
                    "end": int(parts[2]),
                    "mention": parts[3],
                    "type": parts[4],
                    "identifier": parts[5] if len(parts) > 5 else "",
                })
    flush_doc()
    return docs


def compute_ttr_metrics(all_tokens: list):
    """
    计算原始 TTR、固定长度截断 TTR 以及多次随机抽样平均 TTR (消除 Heaps' Law 长度稀释效应)
    """
    n_tokens = len(all_tokens)
    n_types = len(set(all_tokens))
    raw_ttr = n_types / n_tokens if n_tokens else 0.0

    # 设定等长抽样窗口 (覆盖从 1,000 到 1,000,000 的多尺度窗口)
    random.seed(42)
    sample_metrics = {}
    for k in [1000, 5000, 10000, 50000, 100000, 500000, 1000000]:
        if n_tokens >= k:
            head_types = len(set(all_tokens[:k]))
            head_ttr = head_types / k
            # 10 次随机打散抽样 (对应课件 shuf 方案)
            shuf_ttrs = []
            for _ in range(10):
                sampled = random.sample(all_tokens, k)
                shuf_ttrs.append(len(set(sampled)) / k)
            sample_metrics[f"N={k}"] = {
                "head_types": head_types,
                "head_ttr": round(head_ttr, 4),
                "shuf_mean_ttr": round(sum(shuf_ttrs) / len(shuf_ttrs), 4),
            }
    return {
        "total_tokens_N": n_tokens,
        "unique_types_V": n_types,
        "raw_ttr": round(raw_ttr, 4),
        "standardized_samples": sample_metrics,
    }


def compute_tfidf_keywords(docs: list, top_k: int = 40):
    """
    实现 TF-IDF 核心术语提取（同时支持 1-gram 专业词与 2-gram 领域短语）
    采用稀疏累加 O(NNZ) 算法，毫秒级支持万篇文献规模
    TF(t, d) = count(t in d) / total_terms(d)
    IDF(t) = log((1 + N_docs) / (1 + df(t))) + 1
    """
    eng_stops = set(stopwords.words("english")) | {
        "rice", "oryza", "sativa", "plant", "plants", "study", "studies", "results",
        "showed", "using", "used", "also", "may", "however", "two", "one", "three",
        "under", "significantly", "high", "low", "different", "analysis", "important",
        "role", "roles", "effects", "effect", "via", "well", "new", "data", "found",
        "increased", "decreased", "enhanced", "improved", "reduced", "compared", "among",
        "title", "abstract", "article", "author", "authors", "rights", "reserved"
    }

    n_docs = len(docs)
    df_counter = Counter()
    corpus_tf_counter = Counter()
    sum_normalized_tf = Counter()

    for doc in docs:
        raw_tokens = [w.lower() for w in word_tokenize(doc["text"])]
        clean_unigrams = [
            w for w in raw_tokens
            if len(w) >= 3 and re.match(r"^[a-z][a-z0-9\-]*$", w) and w not in eng_stops
        ]
        clean_bigrams = []
        for i in range(len(raw_tokens) - 1):
            w1, w2 = raw_tokens[i], raw_tokens[i + 1]
            if (
                len(w1) >= 3 and len(w2) >= 3
                and w1.isalpha() and w2.isalpha()
                and w1 not in eng_stops and w2 not in eng_stops
            ):
                clean_bigrams.append(f"{w1} {w2}")

        combined_terms = clean_unigrams + clean_bigrams
        doc_len = len(combined_terms)
        if doc_len < 10:
            continue
        tf_cnt = Counter(combined_terms)
        corpus_tf_counter.update(tf_cnt)
        for term, cnt in tf_cnt.items():
            df_counter[term] += 1
            sum_normalized_tf[term] += cnt / doc_len

    min_df = max(2, int(n_docs * 0.005))
    tfidf_scores = {}
    for term, df in df_counter.items():
        if df < min_df:
            continue
        idf = math.log((1 + n_docs) / (1 + df)) + 1.0
        total_tfidf = sum_normalized_tf[term] * idf
        tfidf_scores[term] = {
            "term": term,
            "corpus_freq": corpus_tf_counter[term],
            "doc_freq": df,
            "idf": round(idf, 4),
            "tfidf_score": round(total_tfidf, 4),
        }

    ranked = sorted(tfidf_scores.values(), key=lambda x: x["tfidf_score"], reverse=True)
    unigrams_top = [r for r in ranked if " " not in r["term"]][:top_k]
    bigrams_top = [r for r in ranked if " " in r["term"]][:top_k]
    return unigrams_top, bigrams_top


def main():
    if not PUBTATOR_FILE.exists():
        print(f"未找到 {PUBTATOR_FILE}，请先运行 01_fetch_pubmed_and_pubtator.py")
        return

    docs = parse_pubtator_documents(PUBTATOR_FILE)
    print(f"=== Part 1 文献解析与偏移量校验 ===")
    total_ents = sum(len(d["entities"]) for d in docs)
    verified_ents = sum(d["verified_entities"] for d in docs)
    print(f"成功解析文献数: {len(docs)} 篇")
    print(f"PubTator 原始实体标注总数: {total_ents:,} 个, 字符偏移量 text[begin:end]==mention 验证通过: {verified_ents:,} 个 ({verified_ents/max(1,total_ents)*100:.1f}%)")

    # 句子切分与分词
    all_sentences = []
    all_alpha_tokens = []
    for d in docs:
        sents = sent_tokenize(d["text"])
        all_sentences.extend(sents)
        for s in sents:
            toks = [t.lower() for t in word_tokenize(s) if any(c.isalnum() for c in t)]
            all_alpha_tokens.extend(toks)

    ttr_stats = compute_ttr_metrics(all_alpha_tokens)
    print(f"\n=== 词汇丰富度 (TTR) 统计结果 ===")
    print(f"总句子数: {len(all_sentences):,}")
    print(f"总形符数 (Tokens, N): {ttr_stats['total_tokens_N']:,}")
    print(f"独立类符数 (Types, V): {ttr_stats['unique_types_V']:,}")
    print(f"原始全语料 TTR (V/N): {ttr_stats['raw_ttr']*100:.2f}%")
    for k, v in ttr_stats["standardized_samples"].items():
        print(f"  定长尺度 [{k}]: 开头截断 TTR = {v['head_ttr']*100:.2f}%, 随机抽样(shuf)平均 TTR = {v['shuf_mean_ttr']*100:.2f}%")

    unigrams_top, bigrams_top = compute_tfidf_keywords(docs, top_k=30)
    print(f"\n=== TF-IDF 核心术语 Top 15 (单字词 Unigrams) ===")
    for i, item in enumerate(unigrams_top[:15], 1):
        print(f"  {i:2d}. {item['term']:<18} | TF={item['corpus_freq']:3d} | DF={item['doc_freq']:2d} | TF-IDF={item['tfidf_score']:.4f}")

    print(f"\n=== TF-IDF 核心短语 Top 15 (二元词组 Bigrams) ===")
    for i, item in enumerate(bigrams_top[:15], 1):
        print(f"  {i:2d}. {item['term']:<26} | TF={item['corpus_freq']:3d} | DF={item['doc_freq']:2d} | TF-IDF={item['tfidf_score']:.4f}")

    # 保存汇总结果
    summary = {
        "num_documents": len(docs),
        "num_sentences": len(all_sentences),
        "pubtator_entities_total": total_ents,
        "pubtator_entities_offset_verified": verified_ents,
        "ttr_metrics": ttr_stats,
        "top_unigrams_tfidf": unigrams_top,
        "top_bigrams_tfidf": bigrams_top,
    }
    (RESULTS_DIR / "part1_lexical_and_tfidf_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    with open(RESULTS_DIR / "part1_tfidf_top_terms.tsv", "w", encoding="utf-8") as f:
        f.write("Type\tRank\tTerm\tCorpus_Freq\tDoc_Freq\tIDF\tTFIDF_Score\n")
        for idx, item in enumerate(unigrams_top, 1):
            f.write(f"Unigram\t{idx}\t{item['term']}\t{item['corpus_freq']}\t{item['doc_freq']}\t{item['idf']}\t{item['tfidf_score']}\n")
        for idx, item in enumerate(bigrams_top, 1):
            f.write(f"Bigram\t{idx}\t{item['term']}\t{item['corpus_freq']}\t{item['doc_freq']}\t{item['idf']}\t{item['tfidf_score']}\n")


if __name__ == "__main__":
    main()
