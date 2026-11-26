# -*- coding: utf-8 -*-
"""
Step 7 (阶段三): 基于 spaCy 依存句法树的最短依存路径 (SDP) 基因-性状关系抽取
==========================================================================
核心思想:
  先求基因与性状中心词的最短依存路径 (SDP)，再在路径/其支配动词中寻找
  "调控谓词"，最后用句法论元角色 (主语/宾语/被动施事/名词化构建) 判定
  关系方向与调控极性。只有满足明确句法模式的对才保留为 relation，
  其余降级为 association / cooccurrence_only，从源头抑制共现假阳性。

输出:
  - results/sdp_relations.tsv
  - results/sdp_summary.json
  - results/sdp_error_cases.tsv
"""

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from importlib.machinery import SourceFileLoader
from pathlib import Path

import networkx as nx
import spacy

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
PUBTATOR_FILE = DATA_DIR / "result_Pubtator_all.txt"
EVIDENCE_FILE = RESULTS_DIR / "rice_drought_sentence_evidence.tsv"
ENTITIES_FILE = RESULTS_DIR / "all_sentence_entities.tsv"

step2 = SourceFileLoader("step2", str(BASE_DIR / "scripts" / "02_lexical_stats_ttr_tfidf.py")).load_module()
step4 = SourceFileLoader("step4", str(BASE_DIR / "scripts" / "04_extract_sentence_evidence.py")).load_module()

# ---------------------------------------------------------------------------
# 调控谓词与构建类型词典
# ---------------------------------------------------------------------------
POSITIVE_VERBS = {
    "improve", "improves", "improved", "improving", "enhance", "enhances", "enhanced",
    "enhancing", "increase", "increases", "increased", "increasing", "promote", "promotes",
    "promoted", "promoting", "confer", "confers", "conferred", "conferring", "elevate",
    "elevates", "elevated", "raise", "raises", "raised", "boost", "boosts", "boosted",
    "strengthen", "strengthens", "strengthened", "activate", "activates", "activated",
    "induce", "induces", "induced", "inducing", "upregulate", "upregulates", "upregulated",
    "up-regulate", "up-regulates", "up-regulated", "overexpress", "overexpresses",
    "overexpressed", "facilitate", "facilitates", "facilitated", "maintain", "maintains",
    "maintained", "rescue", "rescues", "rescued", "restore", "restores", "restored",
    "accelerate", "accelerates", "accelerated", "stimulate", "stimulates", "stimulated",
    "improve tolerance", "enhance tolerance", "positively",
}

NEGATIVE_VERBS = {
    "reduce", "reduces", "reduced", "reducing", "decrease", "decreases", "decreased",
    "decreasing", "inhibit", "inhibits", "inhibited", "inhibiting", "suppress", "suppresses",
    "suppressed", "suppressing", "impair", "impairs", "impaired", "impairing", "negatively",
    "diminish", "diminishes", "diminished", "lower", "lowers", "lowered", "restrict",
    "restricts", "restricted", "limit", "limits", "limited", "compromise", "compromises",
    "compromised", "abolish", "abolishes", "abolished", "downregulate", "downregulates",
    "downregulated", "down-regulate", "down-regulates", "down-regulated", "attenuate",
    "attenuates", "attenuated", "weaken", "weakens", "weakened", "delay", "delays",
    "delayed", "repress", "represses", "repressed", "prevent", "prevents", "prevented",
    "block", "blocks", "blocked", "disrupt", "disrupts", "disrupted",
}

NEUTRAL_VERBS = {
    "regulate", "regulates", "regulated", "regulating", "modulate", "modulates",
    "modulated", "affect", "affects", "affected", "affecting", "influence", "influences",
    "influenced", "alter", "alters", "altered", "altering", "control", "controls",
    "controlled", "mediate", "mediates", "mediated", "determine", "determines",
    "determined", "change", "changes", "changed",
}

ASSOCIATION_VERBS = {
    "associate", "associated", "associates", "correlate", "correlated", "correlates",
    "involve", "involved", "involves", "contribute", "contributes", "contributed",
    "participate", "participates", "participated", "play", "plays", "played",
    "function", "functions", "functioned", "act", "acts", "acted",
}

# 明确不是因果调控的“报告/方法/轻动词”
LIGHT_VERBS = {
    "be", "is", "are", "was", "were", "been", "being", "have", "has", "had",
    "show", "shows", "showed", "shown", "reveal", "reveals", "revealed",
    "suggest", "suggests", "suggested", "indicate", "indicates", "indicated",
    "demonstrate", "demonstrates", "demonstrated", "report", "reports", "reported",
    "observe", "observes", "observed", "find", "finds", "found", "identify",
    "identifies", "identified", "detect", "detects", "detected", "include",
    "includes", "included", "including", "establish", "establishes", "established",
    "generate", "generates", "generated", "pinpoint", "pinpoints", "pinpointed",
    "analyze", "analyzes", "analyzed", "analyse", "analyses", "analysed",
    "characterize", "characterizes", "characterized", "study", "studies", "studied",
    "investigate", "investigates", "investigated", "compare", "compares", "compared",
    "evaluate", "evaluates", "evaluated", "assess", "assesses", "assessed",
    "test", "tests", "tested", "examine", "examines", "examined", "sequence",
    "sequences", "sequenced", "map", "maps", "mapped", "clone", "clones", "cloned",
    "isolate", "isolates", "isolated", "express", "expresses", "expressed",
}

LOSS_CONSTRUCT = {
    "knockout", "knock-out", "knockdown", "knock-down", "knock down", "silencing",
    "silenced", "silence", "mutation", "mutant", "mutants", "loss-of-function",
    "loss of function", "rnai", "antisense", "inactivation", "inactivated",
    "repression", "repressed", "deficiency", "deficient", "disruption", "disrupted", "null",
}

GAIN_CONSTRUCT = {
    "overexpression", "over-expression", "overexpress", "overexpressed",
    "over-expressed", "ectopic", "activation", "gain-of-function", "gain of function",
    "transgenic", "constitutive",
}

REGULATION_NOUNS = {
    "tolerance", "resistance", "response", "sensitivity", "accumulation", "content",
    "growth", "yield", "survival", "adaptation", "acclimation", "performance",
    "development", "architecture", "length", "weight", "rate", "index", "phenotype",
}

# 与性状/普通名词同形的基因符号：仅当 mention 是裸词时过滤（如 OsSALT1 保留）
GENE_MENTION_STOPLIST = {
    "salt", "drought", "heat", "cold", "stress", "tolerance", "resistance",
    "growth", "yield", "water", "rice", "plant", "gene", "genes", "protein",
    "mutant", "mutants", "wild", "type", "line", "lines", "root", "roots",
    "leaf", "leaves", "seed", "seeds", "grain", "grains", "height", "weight",
    "length", "content", "level", "levels", "rate", "index", "phenotype",
    "development", "architecture", "sensitivity", "accumulation", "response",
}

NEGATION_CUES = {"not", "no", "none", "neither", "nor", "without", "failed", "fail", "lacks", "lack"}

POSITIVE_ADJ = {
    "enhanced", "increased", "improved", "elevated", "higher", "greater", "promoted",
    "positive", "strong", "robust", "better", "upregulated", "up-regulated",
}
NEGATIVE_ADJ = {
    "reduced", "decreased", "lower", "less", "impaired", "diminished", "negative",
    "weaker", "suppressed", "downregulated", "down-regulated",
}
REQUIRED_ADJ = {"required", "essential", "necessary", "important", "critical"}
NO_EFFECT_PATTERNS = [
    r"\bno\s+(significant\s+)?(effect|change|difference|impact|influence)\b",
    r"\bdid\s+not\s+(alter|change|affect|influence|improve|enhance|reduce|increase)\b",
    r"\bdoes\s+not\s+(alter|change|affect|influence|improve|enhance|reduce|increase)\b",
    r"\bno\s+obvious\b",
    r"\bunaffected\b",
    r"\bwithout\s+(compromising|affecting|reducing|altering|changing)\b",
    r"\bdid\s+not\s+differ\b",
    r"\bunchanged\b",
]


def token_span_for_entity(doc, sent_begin: int, begin: int, end: int):
    rel_b, rel_e = begin - sent_begin, end - sent_begin
    return doc.char_span(rel_b, rel_e, alignment_mode="expand")


def entity_head_token(span):
    """取实体名词短语的句法中心词。"""
    if span is None or len(span) == 0:
        return None
    priority = {"nsubj", "nsubjpass", "dobj", "obj", "attr", "pobj", "conj", "ROOT", "appos"}
    for t in span:
        if t.dep_ in priority and (t.head.i < span.start or t.head.i >= span.end):
            return t
    candidates = [t for t in span if t.dep_ in priority]
    if candidates:
        return candidates[-1]
    nouns = [t for t in span if t.pos_ in {"NOUN", "PROPN", "PRON"}]
    return nouns[-1] if nouns else span[-1]


def build_dep_graph(doc):
    g = nx.Graph()
    for t in doc:
        g.add_node(t.i)
        if t.head.i != t.i:
            g.add_edge(t.i, t.head.i)
    return g


def path_tokens(doc, graph, a: int, b: int):
    try:
        nodes = nx.shortest_path(graph, a, b)
    except (nx.NetworkXNoPath, nx.NodeNotFound):
        return None
    return [doc[i] for i in nodes]


def subtree_text(token):
    return " ".join(t.text.lower() for t in token.subtree)


def is_verb(token):
    return token.pos_ in {"VERB", "AUX"} or token.dep_ == "ROOT"


def subject_phrase(verb):
    """返回动词的主语短语 token（含被动主语），可能为 None。"""
    for c in verb.children:
        if c.dep_ in {"nsubj", "nsubjpass", "csubj"}:
            return c
    return None


def agent_phrase(verb):
    for c in verb.children:
        if c.dep_ in {"agent", "obl"}:
            return c
    return None


def gene_in_subject(verb, gene_head):
    """基因（或其名词化构建）是否位于动词主语短语中。"""
    s = subject_phrase(verb)
    if s is not None and gene_head in s.subtree:
        return True, s
    # 关系从句: "the gene that regulates ..."
    if verb.dep_ == "relcl" and verb.head == gene_head:
        return True, gene_head
    # 并列主语: "A and B regulate ..."
    if s is not None:
        for c in s.children:
            if c.dep_ == "conj" and gene_head in c.subtree:
                return True, c
    return False, None


def trait_in_object(verb, trait_head):
    """性状是否位于动词宾语/表语中；返回 (bool, object_token, pattern)。"""
    for c in verb.children:
        if c.dep_ in {"dobj", "obj", "attr", "oprd"}:
            if trait_head in c.subtree:
                return True, c, "direct_object"
    # 被动: trait 是 nsubjpass
    for c in verb.children:
        if c.dep_ == "nsubjpass" and trait_head in c.subtree:
            return True, c, "passive_subject"
    # 介词宾语: "regulates tolerance to drought" / "increases sensitivity to salt"
    for p in verb.children:
        if p.dep_ == "prep" and p.lemma_ in {"to", "of", "for", "against", "on"}:
            for c in p.children:
                if c.dep_ == "pobj" and trait_head in c.subtree:
                    return True, c, "prepositional_object"
    # 宽松宾语: spaCy 把 "to salt stress" 误挂为 advcl/obl/nmod；
    # 仅当前置介词是 to/of/for/against/on 时接受，under/in/during/with 等条件状语拒绝。
    for c in verb.children:
        if c.dep_ in {"advcl", "obl", "nmod"} and c.head == verb:
            if trait_head in c.subtree:
                prev = trait_head.doc[trait_head.i - 1] if trait_head.i > 0 else None
                if prev is not None and prev.lemma_.lower() in {"to", "of", "for", "against", "on"}:
                    return True, c, "loose_prepositional_object"
    # 关联动词: "associated with drought tolerance"
    if (verb.lemma_ or "").lower() in ASSOCIATION_VERBS:
        for p in verb.children:
            if p.dep_ == "prep" and p.lemma_ in {"with", "to"}:
                for c in p.children:
                    if c.dep_ == "pobj" and trait_head in c.subtree:
                        return True, c, "association_object"
    return False, None, ""


def gene_in_agent(verb, gene_head):
    """被动句: 基因是否位于 by-agent 短语中。"""
    for c in verb.children:
        if c.dep_ == "agent":
            if gene_head in c.subtree:
                return True, c
        if c.dep_ == "prep" and c.lemma_ == "by":
            if gene_head in c.subtree:
                return True, c
    return False, None


def is_nominalized_subject(verb, subject_token):
    """主语是否为 overexpression/knockout of gene 这类名词化构建。"""
    if subject_token is None:
        return False
    return subject_token.lemma_.lower() in {"overexpression", "over-expression", "knockout",
                                             "knock-out", "knockdown", "knock-down",
                                             "silencing", "mutation", "activation",
                                             "induction", "expression"}


def construct_sign(subject_token, agent_token, gene_head):
    """判断基因构建是功能获得 (+1) 还是功能缺失 (-1) 还是未知 (0)。"""
    text = ""
    if subject_token is not None:
        text += " " + subtree_text(subject_token)
    if agent_token is not None:
        text += " " + subtree_text(agent_token)
    text = text.lower()
    if any(k in text for k in LOSS_CONSTRUCT):
        return -1, "loss_of_function_construct"
    if any(k in text for k in GAIN_CONSTRUCT):
        return 1, "gain_of_function_construct"
    # 基因本身作主语（无构建名词）视为野生型/直接调控
    return 1, "direct_gene_subject"


def find_light_adjective(obj_token):
    """在宾语/表语子树中找 enhanced/increased/reduced 等评价形容词。"""
    if obj_token is None:
        return None
    for t in obj_token.subtree:
        lem = t.lemma_.lower()
        if lem in POSITIVE_ADJ:
            return "positive", t
        if lem in NEGATIVE_ADJ:
            return "negative", t
    return None


def expand_trait_span(stext, sbegin, gene_ent, trait_ent, canonical_name):
    """
    若句子里实际写的是更完整的标准性状短语（如 mention='drought'，原文='drought tolerance'），
    则把性状 span 扩展到完整短语，避免把 'under drought stress' 误当成宾语。
    """
    if not canonical_name:
        return trait_ent["begin"], trait_ent["end"], False
    low = stext.lower()
    if canonical_name.lower() not in low or len(canonical_name) <= len(trait_ent["mention"]):
        return trait_ent["begin"], trait_ent["end"], False
    for m in re.finditer(re.escape(canonical_name), stext, re.I):
        b, e = sbegin + m.start(), sbegin + m.end()
        if not (b < gene_ent["end"] and gene_ent["begin"] < e):
            return b, e, True
    return trait_ent["begin"], trait_ent["end"], False


def choose_predicate(doc, path, gene_head, trait_head):
    """
    在路径上的动词及其支配动词中，选择同时连接基因与性状的最佳调控谓词。
    返回 (verb, score, reason, direction, subject_token, object_token, pattern)
    """
    verbs = [t for t in path if is_verb(t)]
    # 路径节点上方的动词 (如 ccomp/relcl 的支配动词) 也纳入候选
    extra = []
    for t in path:
        h = t.head
        if is_verb(h) and h not in path:
            extra.append(h)
    candidates = list(dict.fromkeys(verbs + extra))
    if not candidates:
        return None, 0.0, "no_verb_on_path", None, None, None, ""

    best = None
    for v in candidates:
        g_sub, subj = gene_in_subject(v, gene_head)
        t_obj, obj, obj_pattern = trait_in_object(v, trait_head)
        g_agent, agent = gene_in_agent(v, gene_head)
        # 被动: 性状作主语 + 基因在 by-agent
        passive_ok = obj_pattern == "passive_subject" and g_agent
        # 主动: 基因作主语 + 性状作宾语
        active_ok = g_sub and t_obj and obj_pattern != "passive_subject"
        # 名词化: overexpression of gene 作主语 + 性状作宾语
        nominal_ok = g_sub and t_obj and is_nominalized_subject(v, subj)

        score = 0.0
        direction = None
        if active_ok:
            score += 6.0
            direction = "gene_to_trait"
        if passive_ok:
            score += 5.5
            direction = "gene_to_trait"
        if nominal_ok:
            score += 1.0
        if g_sub:
            score += 1.0
        if t_obj:
            score += 1.0
        if g_agent:
            score += 0.5

        lemma = (v.lemma_ or v.text).lower()
        if lemma in POSITIVE_VERBS or lemma in NEGATIVE_VERBS or lemma in NEUTRAL_VERBS:
            score += 1.0
        elif lemma in ASSOCIATION_VERBS:
            score -= 0.5
        elif lemma in LIGHT_VERBS:
            score -= 2.0
        if v in path:
            score += 0.3
        score -= 0.1 * (abs(path.index(v)) if v in path else 2)

        if best is None or score > best[0]:
            pattern = "active" if active_ok else ("passive" if passive_ok else ("nominalized" if nominal_ok else "subtree_only"))
            best = (score, v, "ok", direction, subj, obj, pattern)

    # 链式结构: "gene has the potential to improve trait" / "gene is a MKP which decreased trait"
    # V1 的主语是基因，V2 的宾语是性状，且 V2 位于 V1 的补足语/关系从句子树中。
    if best is None or best[3] is None:
        for v1 in candidates:
            gs, subj1 = gene_in_subject(v1, gene_head)
            if not gs:
                continue
            for v2 in candidates:
                if v2 is v1:
                    continue
                to, obj2, pat2 = trait_in_object(v2, trait_head)
                if not to:
                    continue
                if v2 in v1.subtree or v1 in v2.subtree:
                    score = 4.5
                    if (v2.lemma_ or "").lower() in POSITIVE_VERBS or (v2.lemma_ or "").lower() in NEGATIVE_VERBS:
                        score += 0.5
                    if best is None or score > best[0]:
                        best = (score, v2, "chain", "gene_to_trait", subj1, obj2, "chain")

    assert best is not None
    if best[0] <= 0:
        return best[1], best[0], "no_argument_pattern", None, best[4], best[5], best[6]
    return best[1], best[0], best[2], best[3], best[4], best[5], best[6]


def verb_polarity(sent_text: str, verb, direction):
    """判定动词本身的方向（未考虑基因构建类型）。"""
    if verb is None:
        return "unknown", "no_verb"
    low = sent_text.lower()
    lemma = (verb.lemma_ or verb.text).lower() if verb is not None else ""
    text = verb.text.lower() if verb is not None else ""

    for pat in NO_EFFECT_PATTERNS:
        window = sent_text[max(0, verb.idx - 110): verb.idx + 180]
        if re.search(pat, window):
            return "no_effect", "no_effect_pattern"

    # 动词子节点中的否定词
    if any(c.dep_ == "neg" for c in verb.children):
        return "no_effect", "negated_predicate"

    modifiers = {c.text.lower() for c in verb.children if c.dep_ == "advmod"}
    if "negatively" in modifiers or "negative" in modifiers:
        return "negative", "negative_adverb"
    if "positively" in modifiers or "positive" in modifiers:
        return "positive", "positive_adverb"

    if lemma in POSITIVE_VERBS or text in POSITIVE_VERBS:
        return "positive", "positive_verb"
    if lemma in NEGATIVE_VERBS or text in NEGATIVE_VERBS:
        return "negative", "negative_verb"
    if lemma in NEUTRAL_VERBS or text in NEUTRAL_VERBS:
        return "neutral", "neutral_verb"
    return "unknown", "unknown_predicate"


def final_polarity(verb_pol: str, construct: int):
    if verb_pol == "no_effect":
        return "no_effect"
    if construct == -1:
        if verb_pol == "positive":
            return "negative"
        if verb_pol == "negative":
            return "positive"
        return verb_pol
    return verb_pol


def detect_respectively(sent_text, gene_entities, trait_entities, cur_gene, cur_trait):
    low = sent_text.lower()
    if "respectively" not in low:
        return True, "no_respectively", None, None
    if len(gene_entities) < 2 or len(trait_entities) < 2:
        return True, "respectively_insufficient_lists", None, None
    genes = sorted(gene_entities, key=lambda e: e["begin"])
    traits = sorted(trait_entities, key=lambda e: e["begin"])
    gi = next((i for i, e in enumerate(genes) if e["begin"] == cur_gene["begin"] and e["end"] == cur_gene["end"]), None)
    ti = next((i for i, e in enumerate(traits) if e["begin"] == cur_trait["begin"] and e["end"] == cur_trait["end"]), None)
    if gi is None or ti is None:
        return True, "respectively_entity_not_in_lists", None, None
    if gi == ti:
        return True, "respectively_aligned", gi, ti
    return False, f"respectively_mismatch_g{gi}_t{ti}", gi, ti


def entity_rows_for(evidence_row, ent_index):
    pmid, sid = evidence_row["PMID"], evidence_row["sent_id"]
    g_cands = [e for e in ent_index.get((pmid, sid, evidence_row["gene_mention"], "Gene"), [])
               if e["norm_id"] == evidence_row["RAP_ID"]]
    t_cands = [e for e in ent_index.get((pmid, sid, evidence_row["trait_mention"], "Trait"), [])
               if e["norm_id"] == evidence_row["RTO_ID"]]
    return g_cands, t_cands


def choose_pair(g_cands, t_cands):
    best = None
    for g in g_cands:
        for t in t_cands:
            if g["begin"] < t["end"] and t["begin"] < g["end"]:
                continue
            dist = abs(g["begin"] - t["begin"])
            if best is None or dist < best[0]:
                best = (dist, g, t)
    if best is None:
        return (g_cands[0] if g_cands else None, t_cands[0] if t_cands else None)
    return best[1], best[2]


def confidence_label(path_len, pred, pred_score, polarity, resp_reason, pattern, construct):
    score = 0.0
    if pattern == "active":
        score += 3.0
    elif pattern == "passive":
        score += 2.5
    elif pattern == "nominalized":
        score += 2.0
    elif pattern == "subtree_only":
        score += 0.5
    if path_len <= 4:
        score += 2.0
    elif path_len <= 6:
        score += 1.0
    if polarity in {"positive", "negative"}:
        score += 1.0
    elif polarity == "neutral":
        score += 0.3
    if pred is not None and (pred.lemma_ or "").lower() in LIGHT_VERBS:
        score -= 2.0
    if resp_reason.startswith("respectively_mismatch"):
        score -= 2.0
    if not resp_reason.startswith("no_respectively"):
        score -= 0.3
    if score >= 5.0:
        return "high"
    if score >= 3.0:
        return "medium"
    return "low"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--rice-only", action="store_true")
    args = parser.parse_args()

    print("[1/5] 读取语料、候选证据与实体偏移 ...")
    docs = step2.parse_pubtator_documents(PUBTATOR_FILE)
    doc_by_pmid = {d["pmid"]: d for d in docs}

    sent_begin = defaultdict(dict)
    for d in docs:
        for sp in step4.sentence_spans(d["text"]):
            sent_begin[d["pmid"]][str(sp["sent_id"])] = sp["begin"]

    ent_index = defaultdict(list)
    sent_entities = defaultdict(lambda: {"Gene": [], "Trait": []})
    for r in csv.DictReader(open(ENTITIES_FILE, encoding="utf-8"), delimiter="\t"):
        r["begin"], r["end"] = int(r["begin"]), int(r["end"])
        ent_index[(r["PMID"], r["sent_id"], r["mention"], r["type"])].append(r)
        sent_entities[(r["PMID"], r["sent_id"])][r["type"]].append(r)

    evidence_rows = list(csv.DictReader(open(EVIDENCE_FILE, encoding="utf-8"), delimiter="\t"))
    if args.rice_only:
        evidence_rows = [r for r in evidence_rows if r["RAP_ID"].startswith("Os")]
    if args.limit:
        evidence_rows = evidence_rows[: args.limit]
    print(f"    候选证据 {len(evidence_rows):,} 条；语料 {len(docs):,} 篇")

    print("[2/5] 加载 spaCy en_core_web_sm 并对候选句去重解析 ...")
    nlp = spacy.load("en_core_web_sm", disable=["ner"])
    sent_keys = list(dict.fromkeys((r["PMID"], r["sent_id"]) for r in evidence_rows))
    sent_text, sent_doc = {}, {}
    for pmid, sid in sent_keys:
        d = doc_by_pmid.get(pmid)
        if d is None:
            continue
        sb = sent_begin.get(pmid, {}).get(sid)
        if sb is None:
            continue
        spans = step4.sentence_spans(d["text"])
        sp = next((x for x in spans if str(x["sent_id"]) == sid), None)
        if sp is None:
            continue
        sent_text[(pmid, sid)] = (sp["text"], sb)
    parsed = list(nlp.pipe([v[0] for v in sent_text.values()], batch_size=64, n_process=1))
    for key, doc in zip(sent_text.keys(), parsed):
        sent_doc[key] = doc
    print(f"    解析句子 {len(sent_doc):,} 个")

    print("[3/5] 逐候选对执行 SDP、谓词、方向与极性判定 ...")
    out_rows, error_rows = [], []
    counter = Counter()
    polarity_counter = Counter()
    relation_counter = Counter()
    pattern_counter = Counter()

    for idx, er in enumerate(evidence_rows, 1):
        pmid, sid = er["PMID"], er["sent_id"]
        key = (pmid, sid)
        if key not in sent_doc:
            continue
        sdoc = sent_doc[key]
        stext, sbegin = sent_text[key]
        g_cands, t_cands = entity_rows_for(er, ent_index)
        if not g_cands or not t_cands:
            counter["missing_entity_span"] += 1
            continue
        g_ent, t_ent = choose_pair(g_cands, t_cands)
        if g_ent is None or t_ent is None:
            counter["missing_entity_span"] += 1
            continue

        # 校验实体行确实属于当前句子，且字符切片与 mention 完全一致
        g_rel_b, g_rel_e = g_ent["begin"] - sbegin, g_ent["end"] - sbegin
        t_rel_b, t_rel_e = t_ent["begin"] - sbegin, t_ent["end"] - sbegin
        if (
            g_ent.get("sentence") != stext
            or t_ent.get("sentence") != stext
            or g_rel_b < 0 or t_rel_e > len(stext)
            or stext[g_rel_b:g_rel_e] != er["gene_mention"]
            or stext[t_rel_b:t_rel_e] != er["trait_mention"]
        ):
            counter["entity_sentence_mismatch"] += 1
            continue

        # 过滤与性状/普通名词同形的裸基因符号（如 "salt" → SALT 基因）
        if er["gene_mention"].lower() in GENE_MENTION_STOPLIST:
            counter["entity_ambiguity_filtered"] += 1
            continue

        # 性状短语扩展: "drought" → 原文中的 "drought tolerance"
        t_ent_orig = t_ent
        t_exp_b, t_exp_e, t_expanded = expand_trait_span(
            stext, sbegin, g_ent, t_ent, er["trait_canonical"]
        )
        if t_expanded:
            t_ent = dict(t_ent)
            t_ent["begin"], t_ent["end"] = t_exp_b, t_exp_e

        g_span = token_span_for_entity(sdoc, sbegin, g_ent["begin"], g_ent["end"])
        t_span = token_span_for_entity(sdoc, sbegin, t_ent["begin"], t_ent["end"])
        if g_span is None or t_span is None:
            counter["char_span_failed"] += 1
            continue
        g_head, t_head = entity_head_token(g_span), entity_head_token(t_span)
        if g_head is None or t_head is None:
            counter["head_failed"] += 1
            continue

        keep_resp, resp_reason, gi, ti = detect_respectively(
            stext, sent_entities[key]["Gene"], sent_entities[key]["Trait"], g_ent, t_ent_orig
        )

        graph = build_dep_graph(sdoc)
        path = path_tokens(sdoc, graph, g_head.i, t_head.i)
        if path is None:
            pred = None
            path_len, path_str = -1, ""
            pred_reason, direction, pattern = "no_path", None, "no_path"
            verb_pol, verb_reason = "unknown", "no_path"
            construct, construct_reason = 0, "no_path"
            relation = "unsupported"
        else:
            pred, pred_score, pred_reason, direction, subj, obj, pattern = choose_predicate(
                sdoc, path, g_head, t_head
            )
            path_len = len(path) - 1
            path_str = " -> ".join(f"{t.text}/{t.dep_}" for t in path)
            verb_pol, verb_reason = verb_polarity(stext, pred, direction)
            # 轻动词 + 评价形容词: "showed enhanced drought tolerance"
            if pattern == "active" and pred is not None and (pred.lemma_ or "").lower() in LIGHT_VERBS:
                adj = find_light_adjective(obj)
                if adj is not None:
                    pattern = "light_verb_adjective"
                    verb_pol, verb_reason = adj[0], f"light_verb_adjective:{adj[1].text}"
            # required/essential for trait
            if pred is not None and (pred.lemma_ or "").lower() in REQUIRED_ADJ:
                pattern = "required_for"
                verb_pol, verb_reason = "positive", "required_for"
            construct, construct_reason = construct_sign(subj, obj if pattern == "passive" else None, g_head)
            polarity = final_polarity(verb_pol, construct)
            # "reduces growth retardation" = improves growth; 对负向名词宾语翻转
            if polarity in {"positive", "negative"} and t_head.lemma_.lower() in {
                "retardation", "inhibition", "impairment", "deficiency", "reduction"
            }:
                polarity = "negative" if polarity == "positive" else "positive"

            if not keep_resp:
                relation = "excluded_respectively"
            elif direction != "gene_to_trait":
                relation = "excluded_direction"
            elif verb_pol == "no_effect":
                relation = "no_effect"
            elif pred is not None and (pred.lemma_ or "").lower() in ASSOCIATION_VERBS:
                relation = "association"
            elif pattern == "chain" and path_len > 5:
                relation = "association"
            elif pattern in {"active", "passive", "nominalized", "chain", "light_verb_adjective", "required_for"} and pred is not None:
                if polarity == "positive":
                    relation = "positive_regulation"
                elif polarity == "negative":
                    relation = "negative_regulation"
                elif polarity == "neutral":
                    relation = "regulation"
                else:
                    relation = "association"
            elif pred is not None and pattern == "subtree_only":
                relation = "association"
            else:
                relation = "cooccurrence_only"

        has_negation = bool({t.text.lower() for t in path} & NEGATION_CUES) if path else False
        conf = confidence_label(max(path_len, 0), pred, 0.0, polarity, resp_reason, pattern, construct)

        is_rice = er["RAP_ID"].startswith("Os")
        counter["total"] += 1
        relation_counter[relation] += 1
        polarity_counter[polarity] += 1
        pattern_counter[pattern] += 1
        if is_rice:
            relation_counter["rice_" + relation] += 1
            polarity_counter["rice_" + polarity] += 1
            pattern_counter["rice_" + pattern] += 1

        row = {
            "PMID": pmid,
            "sent_id": sid,
            "gene_mention": er["gene_mention"],
            "gene_canonical": er["gene_canonical"],
            "RAP_ID": er["RAP_ID"],
            "gene_source": er["gene_source"],
            "is_rice_gene": int(is_rice),
            "gene_begin": g_ent["begin"],
            "gene_end": g_ent["end"],
            "gene_head": f"{g_head.text}/{g_head.dep_}",
            "trait_mention": er["trait_mention"],
            "trait_canonical": er["trait_canonical"],
            "RTO_ID": er["RTO_ID"],
            "trait_begin": t_ent["begin"],
            "trait_end": t_ent["end"],
            "trait_span_expanded": int(t_expanded),
            "trait_head": f"{t_head.text}/{t_head.dep_}",
            "predicate": pred.text if pred is not None else "",
            "predicate_lemma": (pred.lemma_ or pred.text).lower() if pred is not None else "",
            "predicate_reason": pred_reason,
            "relation_pattern": pattern,
            "direction": direction or "",
            "construct": construct_reason,
            "verb_polarity": verb_pol,
            "polarity": polarity,
            "polarity_reason": verb_reason,
            "relation": relation,
            "confidence": conf,
            "path_len": path_len,
            "sdp_path": path_str,
            "respectively_reason": resp_reason,
            "sentence": stext,
        }
        out_rows.append(row)

        if conf == "low" or relation in {"excluded_respectively", "excluded_direction", "no_effect", "unsupported"} or has_negation:
            error_rows.append(row)

        if idx % 1000 == 0:
            print(f"    已处理 {idx:,}/{len(evidence_rows):,} ...")

    print("[4/5] 写出 TSV 与汇总 JSON ...")
    fields = list(out_rows[0].keys()) if out_rows else []
    with open(RESULTS_DIR / "sdp_relations.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    with open(RESULTS_DIR / "sdp_error_cases.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(error_rows)

    examples = defaultdict(list)
    for r in out_rows:
        if len(examples[r["relation"]]) < 3:
            examples[r["relation"]].append({
                "PMID": r["PMID"], "sent_id": r["sent_id"],
                "gene": r["gene_canonical"], "trait": r["trait_canonical"],
                "predicate": r["predicate"], "polarity": r["polarity"],
                "relation": r["relation"], "confidence": r["confidence"],
                "pattern": r["relation_pattern"], "sdp_path": r["sdp_path"],
                "sentence": r["sentence"],
            })

    rice_rows = [r for r in out_rows if r["is_rice_gene"] == 1]
    supported = {"positive_regulation", "negative_regulation", "regulation"}
    rice_supported = [r for r in rice_rows if r["relation"] in supported]
    summary = {
        "input_candidate_pairs": len(evidence_rows),
        "processed_pairs": len(out_rows),
        "rice_only_pairs": len(rice_rows),
        "rice_only_unique_gene_trait_pairs": len({(r["RAP_ID"], r["RTO_ID"]) for r in rice_rows}),
        "supported_relations_total": sum(1 for r in out_rows if r["relation"] in supported),
        "supported_relations_rice": len(rice_supported),
        "supported_unique_pairs_rice": len({(r["RAP_ID"], r["RTO_ID"]) for r in rice_supported}),
        "positive_relations_rice": sum(1 for r in rice_supported if r["relation"] == "positive_regulation"),
        "negative_relations_rice": sum(1 for r in rice_supported if r["relation"] == "negative_regulation"),
        "relation_counts": dict(relation_counter),
        "polarity_counts": dict(polarity_counter),
        "pattern_counts": dict(pattern_counter),
        "confidence_counts": dict(Counter(r["confidence"] for r in out_rows)),
        "rice_confidence_counts": dict(Counter(r["confidence"] for r in rice_rows)),
        "predicate_top": Counter(r["predicate_lemma"] for r in out_rows if r["predicate_lemma"]).most_common(25),
        "respectively_stats": dict(Counter(r["respectively_reason"] for r in out_rows)),
        "error_cases": len(error_rows),
        "examples": dict(examples),
        "offset_alignment": {
            "missing_entity_span": counter["missing_entity_span"],
            "char_span_failed": counter["char_span_failed"],
            "head_failed": counter["head_failed"],
            "entity_sentence_mismatch": counter["entity_sentence_mismatch"],
            "entity_ambiguity_filtered": counter["entity_ambiguity_filtered"],
        },
    }
    (RESULTS_DIR / "sdp_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print("[5/5] 完成。")
    print(f"    候选对: {len(out_rows):,} (水稻 {len(rice_rows):,})")
    print(f"    支持关系: {summary['supported_relations_total']:,} (水稻 {summary['supported_relations_rice']:,})")
    print(f"    关系分布: {dict(relation_counter)}")
    print(f"    模式分布: {dict(pattern_counter)}")
    print(f"    置信度: {summary['confidence_counts']}")
    print(f"    输出: {RESULTS_DIR / 'sdp_relations.tsv'}")


if __name__ == "__main__":
    main()