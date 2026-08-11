#!/usr/bin/env python3
"""
元数据提取与整合脚本
从各数据集的 iSeq/ENA 元数据文件中提取技术字段，
并生成待手工填写临床字段的模板 TSV。

用法：
    python3 scripts/18_extract_metadata.py \
        --meta-dir data/OLP/raw/meta \
        --output results/OLP/metadata/combined_metadata.tsv \
        --col-report results/OLP/metadata/all_columns.txt

输出：
    combined_metadata.tsv  — 所有样本合并表（技术字段自动填，临床字段标 FILL_ME）
    all_columns.txt        — 各数据集完整列名（用于确认哪些字段有临床信息）
"""

import argparse
import csv
import os
import sys
from pathlib import Path
from collections import defaultdict

# ── 各数据集固定信息（来自文献）────────────────────────────────────
# amplicon_region: V4 / V3-V4 / V1-V2 / WGS
# sample_type: saliva / swab / tissue / plaque_supra / plaque_sub / mixed
DATASET_FIXED = {
    "CRA008410":   {"amplicon_region": "V3-V4", "sample_type": "swab",         "primers": "341F/806R"},
    "PRJDB12280":  {"amplicon_region": "V1-V2", "sample_type": "saliva",        "primers": "27Fmod/338R"},
    "PRJEB90477":  {"amplicon_region": "V1-V2", "sample_type": "plaque_supra",  "primers": "27F/338R"},
    "PRJNA1043432":{"amplicon_region": "V3-V4", "sample_type": "swab",          "primers": "341F/806R"},
    "PRJNA1049117":{"amplicon_region": "V3-V4", "sample_type": "mixed",         "primers": "343F/798R"},
    "PRJNA1201607":{"amplicon_region": "WGS",   "sample_type": "plaque_supra",  "primers": "NA"},
    "PRJNA306560": {"amplicon_region": "V4",    "sample_type": "saliva",        "primers": "515F/806R"},
    "PRJNA512581": {"amplicon_region": "V1-V3_or_V3-V4", "sample_type": "swab_or_tissue", "primers": "mixed"},
    "PRJNA542018": {"amplicon_region": "V4",    "sample_type": "saliva",        "primers": "515F/806R"},
    "PRJNA555458": {"amplicon_region": "V3-V4", "sample_type": "saliva",        "primers": "341F/806R"},
    "PRJNA556311": {"amplicon_region": "V3-V4", "sample_type": "tissue",        "primers": "341F/806R"},
    "PRJNA598825": {"amplicon_region": "V3-V4", "sample_type": "swab",          "primers": "V3-V4_Illumina"},
    "PRJNA690677": {"amplicon_region": "V3-V4", "sample_type": "plaque_sub",    "primers": "338F/806R"},
}

# PRJDB12280 已有处理好的 condition（来自14_filter_prjdb12280.py 输出）
PRJDB12280_META = "data/OLP/raw/meta/PRJDB12280/metadata.tsv"

OUTPUT_COLS = [
    "sample_id", "study_id", "biosample_id",
    "condition",         # OLP / Control / FILL_ME
    "sample_type",       # saliva / swab / tissue / plaque_supra / plaque_sub
    "amplicon_region",   # V4 / V3-V4 / V1-V2 / WGS
    "primers",
    "library_layout",    # PAIRED / SINGLE
    "platform",
    "instrument_model",
    "read_count",
    "sample_title",      # 原始标题，辅助判断 condition
    "library_name",      # 同上
    "age",               # FILL_ME if unknown
    "sex",               # FILL_ME if unknown
    "notes",
]

# 临床字段关键词映射（用于从 sample_title / library_name 推断 condition）
OLP_KEYWORDS   = ["olp", "oral lichen planus", "lichen", "patient", "case"]
CTRL_KEYWORDS  = ["control", "ctrl", "healthy", "normal", "hc"]


def infer_condition(text: str) -> str:
    if not text:
        return "FILL_ME"
    t = text.lower()
    if any(k in t for k in OLP_KEYWORDS):
        return "OLP"
    if any(k in t for k in CTRL_KEYWORDS):
        return "Control"
    return "FILL_ME"


def sniff_delim(path: Path) -> str:
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as f:
        s = f.read(2048)
    return "\t" if s.count("\t") > s.count(",") else ","


def read_tsv_csv(path: Path):
    delim = sniff_delim(path)
    rows = []
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as f:
        reader = csv.DictReader(f, delimiter=delim)
        for row in reader:
            rows.append(dict(row))
    return rows


def col(row: dict, *candidates) -> str:
    """按候选列名顺序取第一个存在且非空的值"""
    for c in candidates:
        for k, v in row.items():
            if k.lower().strip() == c.lower() and v.strip():
                return v.strip()
    return ""


def process_ena_rows(rows, study_id, col_report_lines):
    """处理 ENA 51列格式（iSeq 下载的 metadata.tsv）"""
    if rows:
        col_report_lines.append(f"  列名（共{len(rows[0])}列）：")
        col_report_lines.append("  " + ", ".join(rows[0].keys()))
        col_report_lines.append("")

    fixed = DATASET_FIXED.get(study_id, {})
    records = []
    for row in rows:
        sample_title = col(row, "sample_title", "experiment_title", "study_title")
        lib_name     = col(row, "library_name")
        condition    = infer_condition(sample_title) or infer_condition(lib_name)

        rec = {
            "sample_id":       col(row, "run_accession"),
            "study_id":        study_id,
            "biosample_id":    col(row, "sample_accession", "secondary_sample_accession"),
            "condition":       condition,
            "sample_type":     fixed.get("sample_type", "FILL_ME"),
            "amplicon_region": fixed.get("amplicon_region", "FILL_ME"),
            "primers":         fixed.get("primers", "FILL_ME"),
            "library_layout":  col(row, "library_layout"),
            "platform":        col(row, "instrument_platform"),
            "instrument_model":col(row, "instrument_model"),
            "read_count":      col(row, "read_count"),
            "sample_title":    sample_title,
            "library_name":    lib_name,
            "age":             "FILL_ME",
            "sex":             "FILL_ME",
            "notes":           "",
        }
        records.append(rec)
    return records


def process_cra008410(meta_dir: Path, col_report_lines):
    """CRA008410：GSA 格式 CSV"""
    csv_path = meta_dir / "CRA008410" / "CRA008410.metadata.csv"
    if not csv_path.exists():
        return []
    rows = read_tsv_csv(csv_path)
    if rows:
        col_report_lines.append(f"  列名（共{len(rows[0])}列）：")
        col_report_lines.append("  " + ", ".join(rows[0].keys()))
        col_report_lines.append("")

    fixed = DATASET_FIXED["CRA008410"]
    records = []
    for row in rows:
        sample_title = col(row, "Title", "SampleType", "LibraryName")
        condition    = infer_condition(sample_title)
        biosample    = col(row, "BioSample")
        rec = {
            "sample_id":       col(row, "Run"),
            "study_id":        "CRA008410",
            "biosample_id":    biosample,
            "condition":       condition,
            "sample_type":     fixed["sample_type"],
            "amplicon_region": fixed["amplicon_region"],
            "primers":         fixed["primers"],
            "library_layout":  col(row, "LibraryLayout"),
            "platform":        col(row, "Platform"),
            "instrument_model":"NovaSeq 6000",
            "read_count":      "",
            "sample_title":    sample_title,
            "library_name":    col(row, "LibraryName"),
            "age":             "FILL_ME",
            "sex":             "FILL_ME",
            "notes":           "",
        }
        records.append(rec)
    return records


def process_prjdb12280(meta_dir: Path, col_report_lines):
    """PRJDB12280：合并 ENA metadata + 已有 condition 标注"""
    ena_path  = meta_dir / "PRJDB12280" / "PRJDB12280.metadata.tsv"
    cond_path = meta_dir / "PRJDB12280" / "metadata.tsv"

    # 读已有 condition 表
    cond_map = {}  # sample_id -> {condition, biosample}
    if cond_path.exists():
        for row in read_tsv_csv(cond_path):
            sid = row.get("sample_id", "").strip()
            if sid:
                cond_map[sid] = row

    rows = read_tsv_csv(ena_path) if ena_path.exists() else []
    if rows:
        col_report_lines.append(f"  列名（共{len(rows[0])}列）：")
        col_report_lines.append("  " + ", ".join(rows[0].keys()))
        col_report_lines.append("")

    fixed = DATASET_FIXED["PRJDB12280"]
    records = []
    for row in rows:
        sid = col(row, "run_accession")
        known = cond_map.get(sid, {})
        rec = {
            "sample_id":       sid,
            "study_id":        "PRJDB12280",
            "biosample_id":    known.get("biosample", col(row, "sample_accession")),
            "condition":       known.get("condition", "FILL_ME"),
            "sample_type":     fixed["sample_type"],
            "amplicon_region": fixed["amplicon_region"],
            "primers":         fixed["primers"],
            "library_layout":  col(row, "library_layout"),
            "platform":        col(row, "instrument_platform"),
            "instrument_model":col(row, "instrument_model"),
            "read_count":      col(row, "read_count"),
            "sample_title":    col(row, "sample_title", "experiment_title"),
            "library_name":    col(row, "library_name"),
            "age":             "FILL_ME",
            "sex":             "FILL_ME",
            "notes":           "",
        }
        records.append(rec)
    return records


def process_prjna1049117(meta_dir: Path, col_report_lines):
    """PRJNA1049117：每个 run 单独一个 tsv，合并"""
    d = meta_dir / "PRJNA1049117"
    tsv_files = sorted(d.glob("SRR*.metadata.tsv"))
    if not tsv_files:
        return []

    sample_rows = []
    shown = False
    for f in tsv_files:
        rows = read_tsv_csv(f)
        if rows:
            if not shown:
                col_report_lines.append(f"  列名（共{len(rows[0])}列，以第一个文件为例）：")
                col_report_lines.append("  " + ", ".join(rows[0].keys()))
                col_report_lines.append("")
                shown = True
            sample_rows.extend(rows)

    return process_ena_rows(sample_rows, "PRJNA1049117", [])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--meta-dir", default="data/OLP/raw/meta")
    parser.add_argument("--output",   default="results/OLP/metadata/combined_metadata.tsv")
    parser.add_argument("--col-report", default="results/OLP/metadata/all_columns.txt")
    args = parser.parse_args()

    meta_root = Path(args.meta_dir)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    all_records = []
    col_lines = []

    # ── CRA008410 ─────────────────────────────────────────────────
    col_lines.append("=== CRA008410 ===")
    all_records.extend(process_cra008410(meta_root, col_lines))

    # ── PRJDB12280 ────────────────────────────────────────────────
    col_lines.append("=== PRJDB12280 ===")
    all_records.extend(process_prjdb12280(meta_root, col_lines))

    # ── PRJNA1049117（per-run tsv）────────────────────────────────
    col_lines.append("=== PRJNA1049117 ===")
    all_records.extend(process_prjna1049117(meta_root, col_lines))

    # ── 其余 ENA 格式数据集 ───────────────────────────────────────
    ena_datasets = [
        "PRJEB90477", "PRJNA1043432", "PRJNA1201607",
        "PRJNA306560", "PRJNA512581", "PRJNA542018",
        "PRJNA555458", "PRJNA556311", "PRJNA598825", "PRJNA690677",
    ]
    for ds in ena_datasets:
        # 找主 metadata 文件
        candidates = list((meta_root / ds).glob(f"{ds}.metadata.tsv")) + \
                     list((meta_root / ds).glob(f"{ds}.metadata.csv"))
        if not candidates:
            col_lines.append(f"=== {ds} ===\n  ⚠️  未找到主元数据文件\n")
            continue
        col_lines.append(f"=== {ds} ===")
        rows = read_tsv_csv(candidates[0])
        all_records.extend(process_ena_rows(rows, ds, col_lines))

    # ── 写输出 ────────────────────────────────────────────────────
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_records)

    with open(args.col_report, "w", encoding="utf-8") as f:
        f.write("\n".join(col_lines))

    # ── 统计 ──────────────────────────────────────────────────────
    from collections import Counter
    by_study = Counter(r["study_id"] for r in all_records)
    cond_fill = sum(1 for r in all_records if r["condition"] == "FILL_ME")
    cond_olp  = sum(1 for r in all_records if r["condition"] == "OLP")
    cond_ctrl = sum(1 for r in all_records if r["condition"] == "Control")

    print(f"[OK] 合并元数据 → {args.output}")
    print(f"     总样本数：{len(all_records)}")
    print(f"     各数据集：")
    for ds, n in sorted(by_study.items()):
        print(f"       {ds}: {n}")
    print(f"     condition 分布：OLP={cond_olp}, Control={cond_ctrl}, FILL_ME={cond_fill}")
    print(f"\n[OK] 完整列名报告 → {args.col_report}")
    if cond_fill > 0:
        print(f"\n[!]  有 {cond_fill} 个样本的 condition 需要手工填写")
        print(f"     请编辑 {args.output}，将 FILL_ME 替换为 OLP 或 Control")


if __name__ == "__main__":
    main()
