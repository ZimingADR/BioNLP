# -*- coding: utf-8 -*-
"""
Step 3 (Part 2 核心): 构建水稻基因归一化词典 (Oryzabase) 与水稻性状本体词典 (RTO)
输入:
  - dicts/OryzabaseGeneListEn.txt (日本国立遗传学研究所 Oryzabase 基因命名词典)
  - dicts/RTO-1.0.obo (华中农业大学 水稻性状本体 Rice Trait Ontology)
输出:
  - dicts/rice_gene_synonym_to_rap.json (基因同义词 -> 归一化信息)
  - dicts/rice_gene_dict.tsv (可读 TSV 格式基因词典)
  - dicts/rto_trait_synonym_to_id.json (性状同义词 -> RTO/TO 标准 ID)
  - dicts/rto_trait_dict.tsv (可读 TSV 格式性状词典)
"""

import json
import re
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DICTS_DIR = BASE_DIR / "dicts"
ORYZA_FILE = DICTS_DIR / "OryzabaseGeneListEn.txt"
RTO_FILE = DICTS_DIR / "RTO-1.0.obo"

# 常见英语普通词/停用词黑名单（防止如基因名 "L", "S", "as", "in", "red", "anova" 等造成假阳性爆炸）
COMMON_ENGLISH_STOPWORDS = {
    "a", "an", "the", "in", "on", "at", "to", "for", "of", "with", "by", "from", "as", "is", "are",
    "was", "were", "be", "been", "being", "have", "has", "had", "do", "does", "did", "but", "and",
    "or", "if", "because", "until", "while", "about", "against", "between", "into", "through",
    "during", "before", "after", "above", "below", "up", "down", "out", "off", "over", "under",
    "again", "further", "then", "once", "here", "there", "when", "where", "why", "how", "all",
    "any", "both", "each", "few", "more", "most", "other", "some", "such", "no", "nor", "not",
    "only", "own", "same", "so", "than", "too", "very", "can", "will", "just", "should", "now",
    "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten",
    "high", "low", "long", "short", "large", "small", "big", "little", "old", "new", "good", "bad",
    "early", "late", "first", "last", "next", "young", "right", "left", "open", "closed",
    "red", "green", "white", "black", "brown", "yellow", "gold", "purple", "dark", "light",
    "leaf", "root", "stem", "seed", "grain", "panicle", "flower", "pollen", "shoot", "plant",
    "rice", "oryza", "sativa", "japonica", "indica", "wild", "mutant", "normal", "type", "line",
    "trait", "gene", "protein", "dna", "rna", "mrna", "cdna", "pcr", "qtl", "snp", "gwas",
    "drought", "salt", "cold", "heat", "stress", "water", "soil", "field", "yield", "height",
    "weight", "length", "width", "ratio", "rate", "level", "content", "index", "number", "area",
    "time", "day", "days", "stage", "phase", "growth", "development", "response", "tolerance",
    "resistance", "sensitivity", "sterility", "fertility", "lethal", "dwarf", "tall", "erect",
    "rolled", "lazy", "lax", "dense", "spreading", "narrow", "wide", "thick", "thin", "soft", "hard",
    "sweet", "waxy", "dull", "floury", "chalky", "brittle", "fragile", "shattering", "heading",
    "flowering", "ripening", "senescence", "germination", "tillering", "branching", "rooting",
    "t", "mt", "cp", "unknown", "none", "null", "na", "n/a", "anova", "fig", "table", "et", "al",
    # 农业与抗旱文献高频非基因词/植物激素/生化指标/实验用语/基因家族统称 (消歧关键黑名单)
    "aba", "asa", "ja", "sa", "ga", "ga3", "br", "ck", "et", "iaa", "ros", "mda", "h2o2", "pro",
    "peg", "rwc", "el", "wt", "oe", "ko", "rnai", "crispr", "cas9", "pot", "field", "soil", "dry",
    "wet", "air", "gas", "ion", "atp", "adp", "nad", "nadp", "nadph", "gsh", "gssg", "cat", "pod",
    "sod", "apx", "gpx", "gst", "pal", "cad", "myb", "myc", "wrky", "nac", "bzip", "bhlh", "erf",
    "dreb", "ap2", "mads", "arf", "aux", "iaa", "grf", "tcp", "spl", "hd", "zip", "hsp", "lea",
    "mapk", "cdpk", "cipk", "snrk", "pp2c", "pyr", "pyl", "rcar", "abscisic acid", "ascorbic acid",
    "salicylic acid", "jasmonic acid", "gibberellin", "brassinosteroid", "cytokinin", "ethylene",
    "proline", "malondialdehyde", "superoxide dismutase", "peroxidase", "catalase",
    "ha-1", "ha1", "apo", "gdp", "nam", "sap", "sams"
}


def clean_gene_synonym_candidates(raw_field: str):
    """
    拆解并清洗 Oryzabase 中的基因符号与同义词字段。
    专门处理水稻遗传学命名中的特殊噪音:
      1. '_' -> Oryzabase 中表示尚无正式 CGSNL 符号的占位符，直接丢弃
      2. 'Hd6(t)*' -> 水稻中 '(t)' 表示 tentative(暂定名)，需提取出 'Hd6'，绝不能把 't' 当成基因名！
      3. 'pms2*(Pms2,pgms,ms2P)' -> 拆出 'pms2', 'Pms2', 'pgms', 'ms2P'
      4. '[CMS-54257]*' -> 去掉末尾 '*'
    """
    if not raw_field or not raw_field.strip():
        return []

    # 先把表示 tentative 的 (t) 或 (t)* 替换为特殊占位或直接剥离
    cleaned_field = re.sub(r"\(\s*t\s*\)\*?", "", raw_field, flags=re.IGNORECASE)

    # 提取括号内的子同义词，例如 a(b, c) -> a, b, c
    parts = []
    for chunk in re.split(r"[,;/]", cleaned_field):
        chunk = chunk.strip()
        if not chunk or chunk == "_":
            continue
        # 检查是否有括号包裹的子名，如 pms2*(Pms2)
        paren_matches = re.findall(r"\(([^)]+)\)", chunk)
        base = re.sub(r"\([^)]*\)", "", chunk).strip()
        if base and base != "_":
            parts.append(base)
        for pm in paren_matches:
            for sub in re.split(r"[,;/]", pm):
                if sub.strip() and sub.strip() != "_":
                    parts.append(sub.strip())

    results = []
    for p in parts:
        # 去除首尾的星号、引号、方括号、下划线、空白
        p = p.strip(" _*'\"`[]()")
        p = re.sub(r"\*+$", "", p).strip()
        if not p or p == "_":
            continue
        results.append(p)
    return results


def is_valid_gene_alias(alias: str) -> bool:
    """
    判断一个基因别名是否适合放入自动文本匹配词典，防止英文常用词及育种世代编号误匹配。
    """
    if not alias or alias == "_":
        return False
    low = alias.lower().strip()
    # 过滤长度 <= 2 的所有短词（含 M1, M2, M3, F1, F2, B3, T1, D1 等在文献中极易与世代/处理组编号冲突的短码）
    if len(low) <= 2:
        return False
    if low.isdigit():
        return False
    if low in COMMON_ENGLISH_STOPWORDS:
        return False
    # 过滤类似 M12, F23, T0, CK1 等单字母+数字的常见实验处理组编号
    if re.match(r"^[mftbck]\d+$", low):
        return False
    # 过滤纯拉丁常见形容词/名词短语（若全为字母且长度<=3，仅保留全大写且不在黑名单中的基因符号如 DST）
    if low.isalpha() and len(low) <= 3 and not alias.isupper():
        return False
    # 过滤首字母大写但全小写的普通英文单词（如 Pot, Dry, Low 等）
    if alias[0].isupper() and alias[1:].islower() and alias.isalpha() and not low.startswith("os"):
        return False
    # 过滤包含过长描述性句子的条目
    if len(low.split()) > 5:
        return False
    return True


def parse_rto_obo(obo_path: Path):
    """
    解析 RTO-1.0.obo 本体文件，提取所有非废弃的性状词条及其同义词。
    """
    text = obo_path.read_text(encoding="utf-8", errors="replace")
    blocks = text.split("[Term]")
    trait_records = []
    synonym_to_rto = {}

    for block in blocks[1:]:
        lines = [line.strip() for line in block.strip().splitlines() if line.strip()]
        term_id = None
        name = None
        is_obsolete = False
        synonyms = []
        is_a_list = []

        for line in lines:
            if line.startswith("["):
                break
            if line.startswith("id:"):
                term_id = line.split("id:", 1)[1].strip()
            elif line.startswith("name:"):
                name = line.split("name:", 1)[1].strip()
            elif line.startswith("is_obsolete:") and "true" in line.lower():
                is_obsolete = True
            elif line.startswith("synonym:"):
                m = re.search(r'"([^"]+)"', line)
                if m:
                    syn = m.group(1).strip()
                    # 去掉 RTO 中常见的 " (related)" 后缀
                    syn = re.sub(r"\s*\((related|exact|broad|narrow)\)\s*$", "", syn, flags=re.IGNORECASE).strip()
                    if syn:
                        synonyms.append(syn)
            elif line.startswith("is_a:"):
                parent = line.split("is_a:", 1)[1].split("!")[0].strip()
                is_a_list.append(parent)

        if is_obsolete or not term_id or not name:
            continue

        all_names = [name] + synonyms
        valid_syns = []
        for nm in all_names:
            nm_clean = nm.strip()
            if not nm_clean or len(nm_clean) <= 2:
                continue
            # 跳过顶层过于宽泛的词如 "Rice Trait"
            if nm_clean.lower() in {"rice trait", "trait", "plant trait", "other"}:
                continue
            valid_syns.append(nm_clean)
            key = nm_clean.lower()
            if key not in synonym_to_rto:
                synonym_to_rto[key] = {
                    "rto_id": term_id,
                    "canonical_name": name,
                    "matched_synonym": nm_clean,
                }

        trait_records.append({
            "rto_id": term_id,
            "name": name,
            "synonyms": list(dict.fromkeys(valid_syns)),
            "is_a": is_a_list,
        })

    # 补充文献中极高频的抗旱核心表述（若本体中只有带后缀形式则映射到对应标准 TO/RTO ID）
    extra_drought_mappings = {
        "drought tolerance": ("TO:0000276", "drought tolerance"),
        "drought resistance": ("TO:0000276", "drought tolerance"),
        "drought stress": ("TO:0000276", "drought tolerance"),
        "drought tolerant": ("TO:0000276", "drought tolerance"),
        "drought resistant": ("TO:0000276", "drought tolerance"),
        "drought adaptability": ("TO:0000276", "drought tolerance"),
        "drought sensitivity": ("TO:0000277", "drought sensitivity"),
        "drought susceptibility": ("TO:0000277", "drought sensitivity"),
        "water deficit": ("TO:0000276", "drought tolerance"),
        "water stress": ("TO:0000276", "drought tolerance"),
        "dehydration tolerance": ("TO:0000276", "drought tolerance"),
        "dehydration stress": ("TO:0000276", "drought tolerance"),
        "stomatal closure": ("TO:0000522", "stomatal resistance"),
        "stomatal conductance": ("TO:0000522", "stomatal resistance"),
        "stomatal aperture": ("TO:0000522", "stomatal resistance"),
        "water loss rate": ("TO:0000276", "drought tolerance"),
        "relative water content": ("TO:0000136", "leaf relative water content"),
        "survival rate": ("TO:0000276", "drought tolerance"),
        "proline content": ("TO:0006002", "proline content"),
        "aba sensitivity": ("TO:0002667", "abscisic acid sensitivity"),
        "root length": ("TO:0000227", "root length"),
        "root depth": ("TO:0000519", "rooting depth"),
        "deep rooting": ("TO:0000519", "rooting depth"),
        "leaf rolling": ("TO:0000085", "leaf rolling"),
        "grain yield": ("TO:0000396", "grain yield"),
    }
    for k, (tid, cname) in extra_drought_mappings.items():
        if k not in synonym_to_rto:
            synonym_to_rto[k] = {
                "rto_id": tid,
                "canonical_name": cname,
                "matched_synonym": k,
            }

    return trait_records, synonym_to_rto


def parse_oryzabase_genes(oryza_path: Path, rto_trait_keys: set):
    """
    解析 OryzabaseGeneListEn.txt，构建 基因别名 -> 标准 RAP ID (或 Oryzabase 回退 ID) 的词典。
    """
    lines = oryza_path.read_text(encoding="utf-8", errors="replace").splitlines()
    header = lines[0].split("\t")
    print(f"Oryzabase 列字段 ({len(header)} 列): {header}")

    gene_records = []
    synonym_to_gene = {}
    rap_count = 0
    fallback_count = 0

    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < 11:
            continue

        trait_gene_id = cols[0].strip()
        cgsnl_symbol = cols[1].strip()
        symbol_syns = cols[2].strip()
        cgsnl_name = cols[3].strip()
        name_syns = cols[4].strip()
        allele_field = cols[6].strip() if len(cols) > 6 else ""
        explanation = cols[8].strip() if len(cols) > 8 else ""
        rap_id_raw = cols[10].strip() if len(cols) > 10 else ""
        msu_id_raw = cols[11].strip() if len(cols) > 11 else ""
        gramene_id = cols[12].strip() if len(cols) > 12 else ""
        go_annot = cols[15].strip() if len(cols) > 15 else ""
        to_annot = cols[16].strip() if len(cols) > 16 else ""
        po_annot = cols[17].strip() if len(cols) > 17 else ""

        # 提取标准 RAP ID (如 Os03g0762000)，若有多个取第一个作为主 ID
        rap_ids = re.findall(r"Os\d{2}g\d{7}", rap_id_raw)
        if rap_ids:
            norm_id = rap_ids[0]
            rap_count += 1
        else:
            norm_id = f"Oryzabase:{trait_gene_id}"
            fallback_count += 1

        # 收集所有符号级候选别名
        raw_candidates = []
        raw_candidates.extend(clean_gene_synonym_candidates(cgsnl_symbol))
        raw_candidates.extend(clean_gene_synonym_candidates(symbol_syns))
        raw_candidates.extend(clean_gene_synonym_candidates(allele_field))
        raw_candidates.extend(rap_ids)

        # 提取 MSU Locus ID (如 LOC_Os03g55389.1 -> LOC_Os03g55389)
        msu_loci = re.findall(r"LOC_Os\d{2}g\d+", msu_id_raw + " " + explanation)
        raw_candidates.extend(msu_loci)

        valid_aliases = []
        for cand in raw_candidates:
            if not is_valid_gene_alias(cand):
                continue
            # 防止某个基因别名恰好与 RTO 表型名完全一致（例如早期突变体直接用表型命名）
            if cand.lower() in rto_trait_keys:
                continue
            valid_aliases.append(cand)
            # 若基因以 Os 开头（如 OsDREB1A），同时登记无连字符版
            if "-" in cand:
                no_hyphen = cand.replace("-", "")
                if is_valid_gene_alias(no_hyphen):
                    valid_aliases.append(no_hyphen)

        valid_aliases = list(dict.fromkeys(valid_aliases))
        if not valid_aliases:
            continue

        canonical_symbol = clean_gene_synonym_candidates(cgsnl_symbol)
        canonical_symbol = canonical_symbol[0] if canonical_symbol else valid_aliases[0]

        rec = {
            "oryzabase_id": trait_gene_id,
            "canonical_symbol": canonical_symbol,
            "cgsnl_name": cgsnl_name,
            "rap_id": norm_id,
            "has_rap_id": bool(rap_ids),
            "msu_ids": list(dict.fromkeys(msu_loci)),
            "aliases": valid_aliases,
            "to_annotations": to_annot,
            "po_annotations": po_annot,
            "go_annotations": go_annot,
        }
        gene_records.append(rec)

        for alias in valid_aliases:
            key = alias.lower()
            # 优先保留拥有真实 RAP ID 的基因映射
            if key not in synonym_to_gene or (bool(rap_ids) and not synonym_to_gene[key]["has_rap_id"]):
                synonym_to_gene[key] = {
                    "rap_id": norm_id,
                    "canonical_symbol": canonical_symbol,
                    "cgsnl_name": cgsnl_name,
                    "oryzabase_id": trait_gene_id,
                    "has_rap_id": bool(rap_ids),
                }

    print(f"基因记录解析完成: 共 {len(gene_records):,} 个基因词条 (其中拥有 RAP ID: {rap_count:,}, 回退 Oryzabase ID: {fallback_count:,})")
    print(f"归一化基因同义词映射表规模: {len(synonym_to_gene):,} 个独立别名入口")
    return gene_records, synonym_to_gene


def main():
    print("=== Step 3.1: 解析水稻性状本体 (RTO-1.0.obo) ===")
    trait_records, synonym_to_rto = parse_rto_obo(RTO_FILE)
    print(f"RTO 性状词条数: {len(trait_records):,}, 可匹配性状同义词总数: {len(synonym_to_rto):,}")

    print("\n=== Step 3.2: 解析 Oryzabase 水稻基因命名词典 (OryzabaseGeneListEn.txt) ===")
    gene_records, synonym_to_gene = parse_oryzabase_genes(ORYZA_FILE, set(synonym_to_rto.keys()))

    # 保存 JSON 映射词典
    (DICTS_DIR / "rto_trait_synonym_to_id.json").write_text(
        json.dumps(synonym_to_rto, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (DICTS_DIR / "rice_gene_synonym_to_rap.json").write_text(
        json.dumps(synonym_to_gene, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # 保存可读 TSV 词典
    with open(DICTS_DIR / "rto_trait_dict.tsv", "w", encoding="utf-8") as f:
        f.write("RTO_ID\tCanonical_Name\tSynonyms\n")
        for r in trait_records:
            f.write(f"{r['rto_id']}\t{r['name']}\t{'; '.join(r['synonyms'])}\n")

    with open(DICTS_DIR / "rice_gene_dict.tsv", "w", encoding="utf-8") as f:
        f.write("RAP_ID_or_Fallback\tCanonical_Symbol\tCGSNL_Name\tAliases\tTrait_Ontology_in_Oryzabase\n")
        for g in gene_records:
            f.write(
                f"{g['rap_id']}\t{g['canonical_symbol']}\t{g['cgsnl_name']}\t"
                f"{', '.join(g['aliases'])}\t{g['to_annotations']}\n"
            )

    print("\n=== Step 3.3: 典型水稻抗旱明星基因与表型归一化自检 ===")
    test_genes = ["OsDREB1A", "OsbZIP23", "SNAC1", "DST", "OsNAC6", "OsCKA2-1", "HD6"]
    for tg in test_genes:
        info = synonym_to_gene.get(tg.lower())
        print(f"  [Gene]  '{tg}' -> {info}")

    test_traits = ["drought tolerance", "drought resistance", "leaf rolling", "root length", "stomatal closure"]
    for tt in test_traits:
        info = synonym_to_rto.get(tt.lower())
        print(f"  [Trait] '{tt}' -> {info}")


if __name__ == "__main__":
    main()
