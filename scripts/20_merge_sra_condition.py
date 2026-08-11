#!/usr/bin/env python3
"""
SraRunTable 条件整合脚本
读取各数据集目录下的 SraRunTable.csv，提取 OLP/Control 分组信息，
合并到 combined_metadata.tsv 中。

用法：
    python3 scripts/20_merge_sra_condition.py \
        --meta-dir data/OLP/raw/meta \
        --input    results/OLP/metadata/combined_metadata.tsv \
        --output   results/OLP/metadata/combined_metadata_annotated.tsv

SraRunTable.csv 放置规则：
    data/OLP/raw/meta/{STUDY_ID}/SraRunTable.csv
    也可以命名为 SraRunTable.txt / RunTable.csv / RunSelector.csv

脚本会自动：
  1. 扫描所有可能包含分组信息的列（disease_state、condition、LibraryName 等）
  2. 用 Run accession（SRR/DRR/ERR）关联到 combined_metadata.tsv
  3. 推断 OLP/Control 并填入 condition 列
  4. 同时填入 age / sex（如果 SraRunTable 有这些字段）
  5. 输出更新后的 metadata 和一份整合报告
"""

import argparse
import csv
import sys
from pathlib import Path
from collections import defaultdict

# ── 关键词 ───────────────────────────────────────────────────────────────────
OLP_KEYWORDS  = [
    "olp", "oral lichen planus", "lichen planus", "patient", "case",
    "diseased", "affected", "disease",
]
CTRL_KEYWORDS = [
    "control", "ctrl", "healthy", "normal", "hc", "non-olp",
    "unaffected", "disease free", "health",
]
EXCLUDE_RULES = [
    (["genital", "vulvar", "penile", "vaginal", "esophageal"], "Exclude_nonoral"),
    (["treated", "treatment", "after therapy", "post-treatment"],   "Exclude_treated"),
    (["hypothyroid", "hypot", "hypo-t"],                            "Exclude_other_disease"),
    (["menstrual", "menstruation", "cycle", "follicular", "luteal"],"Exclude_other_disease"),
    (["periodontitis", "gingivitis", "caries", "cancer", "tumor"],  "Exclude_other_disease"),
]

# SraRunTable 中可能含分组信息的列名（优先级顺序）
CONDITION_COLS = [
    "disease_state", "disease", "health_state", "host_health_state",
    "condition", "status", "diagnosis", "phenotype",
    "patient_group", "group", "Sample Name", "sample_name",
    "LibraryName", "library_name",
    "source_name", "isolation_source", "tissue",
    "subject_is_affected", "affected",
]

AGE_COLS = ["age", "host_age", "Age"]
SEX_COLS = ["sex", "gender", "Sex", "Gender", "host_sex"]

# SraRunTable 中 Run accession 可能用的列名
RUN_COLS = ["Run", "run_accession", "run", "Run Accession", "SRR", "DRR", "ERR"]


def infer_from_text(text: str) -> str:
    if not text:
        return ""
    t = text.lower().strip()
    if t in ("missing", "not applicable", "na", "n/a", "unknown", ""):
        return ""
    for keywords, label in EXCLUDE_RULES:
        if any(k in t for k in keywords):
            return label
    if any(k in t for k in OLP_KEYWORDS):
        return "OLP"
    if any(k in t for k in CTRL_KEYWORDS):
        return "Control"
    return ""


def normalize_sex(val: str) -> str:
    v = val.lower().strip()
    if v in ("male", "m"):
        return "M"
    if v in ("female", "f"):
        return "F"
    if v in ("missing", "not applicable", "na", ""):
        return ""
    return val.strip()


def sniff_delim(path: Path) -> str:
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as f:
        sample = f.read(4096)
    return "\t" if sample.count("\t") > sample.count(",") else ","


def read_table(path: Path) -> list[dict]:
    delim = sniff_delim(path)
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as f:
        reader = csv.DictReader(f, delimiter=delim)
        return [dict(r) for r in reader]


def find_run_col(fieldnames: list[str]) -> str | None:
    for c in RUN_COLS:
        if c in fieldnames:
            return c
    # 模糊匹配
    for f in fieldnames:
        if f.lower() in ("run", "run_accession", "accession"):
            return f
    return None


def extract_condition_from_row(row: dict) -> tuple[str, str]:
    """返回 (condition, evidence)"""
    for col in CONDITION_COLS:
        val = row.get(col, "").strip()
        if not val:
            continue
        cond = infer_from_text(val)
        if cond:
            return cond, f"{col}={val!r}"

    # 全列扫描（跳过明显无关列）
    skip = {"run", "experiment", "submission", "biosample", "bioproject",
            "bytes", "md5", "ftp", "bytes", "read_count", "base_count",
            "first_public", "last_updated", "first_created",
            "collection_date", "geo_loc_name", "lat_lon"}
    for col, val in row.items():
        if col.lower() in skip or not val:
            continue
        cond = infer_from_text(val)
        if cond:
            return cond, f"{col}={val!r}"

    return "", "no_match"


def extract_age_sex(row: dict) -> tuple[str, str]:
    age = ""
    for c in AGE_COLS:
        v = row.get(c, "").strip()
        if v and v.lower() not in ("missing", "not applicable", "na", ""):
            age = v
            break
    sex = ""
    for c in SEX_COLS:
        v = row.get(c, "").strip()
        if v and v.lower() not in ("missing", "not applicable", "na", ""):
            sex = normalize_sex(v)
            break
    return age, sex


def find_sra_table(dataset_dir: Path) -> Path | None:
    """在数据集目录中找 SraRunTable 文件"""
    candidates = [
        "SraRunTable.csv", "SraRunTable.txt",
        "RunTable.csv", "RunTable.txt",
        "RunSelector.csv", "SRAtable.csv",
        "metadata.csv",   # 部分用户可能这样命名
    ]
    for name in candidates:
        p = dataset_dir / name
        if p.exists():
            return p
    # 模糊匹配
    for p in dataset_dir.glob("*.csv"):
        if "sra" in p.name.lower() or "run" in p.name.lower() or "table" in p.name.lower():
            return p
    for p in dataset_dir.glob("*.txt"):
        if "sra" in p.name.lower() or "run" in p.name.lower():
            return p
    return None


def load_sra_table(dataset_dir: Path, study_id: str) -> dict[str, dict]:
    """
    返回 {run_accession: {"cond": ..., "evidence": ..., "age": ..., "sex": ...}}
    """
    sra_path = find_sra_table(dataset_dir)
    if sra_path is None:
        return {}

    try:
        rows = read_table(sra_path)
    except Exception as e:
        print(f"  [WARN] 读取 {sra_path} 失败: {e}")
        return {}

    if not rows:
        return {}

    run_col = find_run_col(list(rows[0].keys()))
    if run_col is None:
        print(f"  [WARN] {sra_path}: 找不到 Run accession 列，跳过")
        return {}

    result = {}
    for row in rows:
        run_id = row.get(run_col, "").strip()
        if not run_id:
            continue
        cond, evidence = extract_condition_from_row(row)
        age, sex = extract_age_sex(row)
        result[run_id] = {"cond": cond, "evidence": evidence,
                          "age": age, "sex": sex}

    print(f"  [{study_id}] 读取 {sra_path.name}：{len(result)} 条记录")
    if result:
        # 打印前几条帮助用户确认
        examples = [(k, v) for k, v in list(result.items())[:3] if v["cond"]]
        if examples:
            for run, info in examples:
                print(f"    示例: {run} → {info['cond']}  ({info['evidence']})")
        else:
            # 没有推断出 condition，展示原始字段
            sample_row = rows[0]
            interesting = {k: v for k, v in sample_row.items()
                           if v and v.lower() not in ("missing", "not applicable")}
            print(f"    [!] 未能推断 condition，样例行字段：")
            for k, v in list(interesting.items())[:8]:
                print(f"        {k!r}: {v!r}")

    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--meta-dir", default="data/OLP/raw/meta",
                        help="各数据集目录根路径")
    parser.add_argument("--input",    default="results/OLP/metadata/combined_metadata.tsv")
    parser.add_argument("--output",   default="results/OLP/metadata/combined_metadata_annotated.tsv")
    parser.add_argument("--report",   default="results/OLP/metadata/sra_merge_report.txt")
    args = parser.parse_args()

    meta_root = Path(args.meta_dir)
    in_path   = Path(args.input)

    if not in_path.exists():
        sys.exit(f"[ERROR] 找不到输入文件：{in_path}")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    # 读取现有 metadata
    with open(in_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames
        rows = list(reader)

    # 找出所有数据集
    study_ids = sorted({r["study_id"] for r in rows if r.get("study_id")})
    print(f"[INFO] combined_metadata 包含 {len(rows)} 条记录，{len(study_ids)} 个数据集\n")

    # 加载各数据集 SraRunTable
    run_lookup: dict[str, dict] = {}   # run_accession -> info
    tables_found = 0
    for study_id in study_ids:
        dataset_dir = meta_root / study_id
        if not dataset_dir.exists():
            continue
        info = load_sra_table(dataset_dir, study_id)
        if info:
            tables_found += 1
            run_lookup.update(info)

    print(f"\n[INFO] 找到 SraRunTable：{tables_found} 个，共 {len(run_lookup)} 条 Run 记录\n")

    if not run_lookup:
        print("[!] 未找到任何 SraRunTable。")
        print("    请将各数据集的 SraRunTable.csv 放入对应目录：")
        print("    data/OLP/raw/meta/{STUDY_ID}/SraRunTable.csv")
        sys.exit(1)

    # 更新字段
    cond_updated = age_updated = sex_updated = 0
    still_fill = 0
    report_lines = ["sample_id\tstudy_id\tcond_before\tcond_after\tevidence\tage\tsex"]

    for row in rows:
        run_id = row.get("sample_id", "").strip()
        info   = run_lookup.get(run_id, {})

        cond_before = row.get("condition", "").strip()

        # condition
        if cond_before == "FILL_ME":
            new_cond = info.get("cond", "")
            if new_cond:
                row["condition"] = new_cond
                cond_updated += 1
            else:
                still_fill += 1
            report_lines.append(
                f"{run_id}\t{row['study_id']}\t{cond_before}\t"
                f"{row['condition']}\t{info.get('evidence','')}\t"
                f"{info.get('age','')}\t{info.get('sex','')}"
            )

        # age
        if row.get("age", "").strip() in ("FILL_ME", "") and info.get("age"):
            row["age"] = info["age"]
            age_updated += 1

        # sex
        if row.get("sex", "").strip() in ("FILL_ME", "") and info.get("sex"):
            row["sex"] = info["sex"]
            sex_updated += 1

    # 写输出
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    with open(args.report, "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))

    # 统计
    by_cond = {}
    for row in rows:
        c = row.get("condition", "?")
        by_cond[c] = by_cond.get(c, 0) + 1

    print(f"[OK] → {args.output}")
    print(f"     condition 更新：{cond_updated} 条")
    print(f"     age 更新：      {age_updated} 条")
    print(f"     sex 更新：      {sex_updated} 条")
    print(f"     仍为 FILL_ME：  {still_fill} 条")
    print(f"\n     分组统计：")
    for k, v in sorted(by_cond.items()):
        print(f"       {k}: {v}")
    print(f"\n[OK] 整合报告 → {args.report}")
    if still_fill:
        print(f"\n[!]  {still_fill} 条无法自动标注，请手工编辑 {args.output}")
        print(f"     可参考各数据集论文的 Supplementary Table")


if __name__ == "__main__":
    main()
