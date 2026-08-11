#!/usr/bin/env python3
"""
数据集专属分组规则应用脚本
根据各数据集 SraRunTable 的具体字段，精确标注 condition。

规则来源：逐数据集人工核查 SraRunTable，确认分组编码方式。
每条规则都有文字说明，便于论文方法部分描述。

用法：
    python3 scripts/21_apply_dataset_rules.py \
        --meta-dir data/OLP/raw/meta \
        --input    results/OLP/metadata/combined_metadata.tsv \
        --output   results/OLP/metadata/combined_metadata_annotated.tsv

注意：
  - 已完成标注的数据集（CRA008410、PRJDB12280）不会被覆盖
  - PRJNA1201607（WGS）暂不处理
  - 每个数据集需在 data/OLP/raw/meta/{STUDY_ID}/ 下放置 SraRunTable.csv
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


# ── 数据集专属规则 ─────────────────────────────────────────────────────────────
#
# 每个函数接收 SraRunTable 某行的字典，返回 (condition, olp_subtype, notes)
#   condition  : OLP / Control / Exclude_* / FILL_ME
#   olp_subtype: erosive / non-erosive / '' （写入 notes 列）
#   notes      : 额外备注
#
# 规则说明：
#   PRJNA306560 — Title 列首字母编码
#       H* → Control（健康对照）
#       R* → OLP 非糜烂型（Reticular）
#       E* → OLP 糜烂型（Erosive）
#
#   PRJNA512581 — Title 列关键词（纯 OLP 队列，无健康对照）
#       含 Tissue  → OLP，组织活检
#       含 Mucosa  → OLP，黏膜拭子
#       含 OLP     → OLP，来源待确认
#       其他       → FILL_ME
#
#   PRJNA542018 — Sample Alias 列前缀编码
#       NE* → OLP 非糜烂型
#       E*  → OLP 糜烂型
#       C*  → Control
#       U*  → Exclude_other（溃疡型，不纳入）
#
#   PRJNA555458 — Sample Alias 列首字母
#       H* → Control
#       O* → OLP
#
#   PRJNA556311 — Sample Alias 列前缀
#       HC*  → Control
#       OLP* → OLP
#
#   PRJNA598825 — Title 列
#       含 OLP（包括 H-OLP）→ OLP
#       含 H（不含 OLP）    → Control
#
#   PRJNA690677 — Sample Alias 列
#       含 OLP → OLP 糜烂型
#       含 CP  → Control
#
#   PRJNA1043432 — Sample Alias 列
#       含 ITS → Exclude_other（真菌 ITS 测序，非 16S）
#       含 NE  → OLP 非糜烂型
#       含 E   → OLP 糜烂型
#       含 C   → Control
#
#   PRJNA1049117 — Sample Alias 列
#       含 HC  → Control
#       含 OLP → OLP
# ──────────────────────────────────────────────────────────────────────────────

def _get(row: dict, *colnames: str) -> str:
    """按优先级取第一个非空列值，不区分大小写"""
    lower_row = {k.lower(): v for k, v in row.items()}
    for c in colnames:
        v = lower_row.get(c.lower(), "").strip()
        if v:
            return v
    return ""


def rule_PRJNA306560(row: dict) -> tuple[str, str, str]:
    title = _get(row, "Title", "title", "Library Name", "sample_title")
    prefix = title[0].upper() if title else ""
    if prefix == "H":
        return "Control", "", ""
    if prefix == "R":
        return "OLP", "non-erosive", "reticular OLP"
    if prefix == "E":
        return "OLP", "erosive", "erosive OLP"
    return "FILL_ME", "", f"unrecognized title={title!r}"


def rule_PRJNA512581(row: dict) -> tuple[str, str, str]:
    title = _get(row, "Title", "title", "Library Name", "sample_title").lower()
    if "tissue" in title:
        return "OLP", "", "OLP tissue biopsy"
    if "mucosa" in title:
        return "OLP", "", "OLP mucosal swab"
    if "olp" in title:
        return "OLP", "", ""
    return "FILL_ME", "", f"unrecognized title={title!r}"


def rule_PRJNA542018(row: dict) -> tuple[str, str, str]:
    alias = _get(row, "Sample Alias", "sample_alias", "Sample Name", "Library Name")
    # NE 必须先于 E 检查
    if re.search(r"NE", alias, re.IGNORECASE):
        return "OLP", "non-erosive", ""
    if re.search(r"\bE\b|^E[_\-\d]|^E$", alias) or \
       (alias and alias[0].upper() == "E" and not alias.upper().startswith("NE")):
        return "OLP", "erosive", ""
    if re.search(r"\bC\b|^C[_\-\d]|^C$", alias) or \
       (alias and alias[0].upper() == "C"):
        return "Control", "", ""
    if re.search(r"\bU\b|^U[_\-\d]|^U$", alias) or \
       (alias and alias[0].upper() == "U"):
        return "Exclude_other", "", "ulcerative type excluded"
    return "FILL_ME", "", f"unrecognized alias={alias!r}"


def rule_PRJNA555458(row: dict) -> tuple[str, str, str]:
    alias = _get(row, "Sample Alias", "sample_alias", "Sample Name", "Library Name")
    prefix = alias[0].upper() if alias else ""
    if prefix == "H":
        return "Control", "", ""
    if prefix == "O":
        return "OLP", "", ""
    return "FILL_ME", "", f"unrecognized alias={alias!r}"


def rule_PRJNA556311(row: dict) -> tuple[str, str, str]:
    alias = _get(row, "Sample Alias", "sample_alias", "Sample Name", "Library Name")
    a = alias.upper()
    if a.startswith("HC"):
        return "Control", "", ""
    if a.startswith("OLP"):
        return "OLP", "", ""
    return "FILL_ME", "", f"unrecognized alias={alias!r}"


def rule_PRJNA598825(row: dict) -> tuple[str, str, str]:
    title = _get(row, "Title", "title", "Library Name", "sample_title")
    t = title.upper()
    if "OLP" in t:          # 包括 H-OLP
        return "OLP", "", ""
    if "H" in t:
        return "Control", "", ""
    return "FILL_ME", "", f"unrecognized title={title!r}"


def rule_PRJNA690677(row: dict) -> tuple[str, str, str]:
    alias = _get(row, "Sample Alias", "sample_alias", "Sample Name", "Library Name")
    a = alias.upper()
    if "OLP" in a:
        return "OLP", "erosive", "erosive OLP"
    if "CP" in a:
        return "Control", "", ""
    return "FILL_ME", "", f"unrecognized alias={alias!r}"


def rule_PRJNA1043432(row: dict) -> tuple[str, str, str]:
    alias = _get(row, "Sample Alias", "sample_alias", "Sample Name", "Library Name")
    a = alias.upper()
    # ITS 先排除（真菌测序）
    if "ITS" in a:
        return "Exclude_other", "", "ITS fungal sequencing excluded"
    # NE 先于 E
    if "NE" in a:
        return "OLP", "non-erosive", ""
    if re.search(r"(?<![A-Z])E(?![A-Z])|^E[^A-Z]|^E$", a):
        return "OLP", "erosive", ""
    if re.search(r"(?<![A-Z])C(?![A-Z])|^C[^A-Z]|^C$", a):
        return "Control", "", ""
    return "FILL_ME", "", f"unrecognized alias={alias!r}"


def rule_PRJNA1049117(row: dict) -> tuple[str, str, str]:
    alias = _get(row, "Sample Alias", "sample_alias", "Sample Name", "Library Name")
    a = alias.upper()
    if "HC" in a:
        return "Control", "", ""
    if "OLP" in a:
        return "OLP", "", ""
    return "FILL_ME", "", f"unrecognized alias={alias!r}"


# 规则注册表
RULES: dict[str, callable] = {
    "PRJNA306560":  rule_PRJNA306560,
    "PRJNA512581":  rule_PRJNA512581,
    "PRJNA542018":  rule_PRJNA542018,
    "PRJNA555458":  rule_PRJNA555458,
    "PRJNA556311":  rule_PRJNA556311,
    "PRJNA598825":  rule_PRJNA598825,
    "PRJNA690677":  rule_PRJNA690677,
    "PRJNA1043432": rule_PRJNA1043432,
    "PRJNA1049117": rule_PRJNA1049117,
    # PRJNA1201607: WGS 鸟枪法，暂不处理
    # CRA008410:    已由 script 18 处理，不覆盖
    # PRJDB12280:   已由 script 14 处理，不覆盖
}

# 已完成标注、不覆盖的数据集
DO_NOT_OVERWRITE = {"CRA008410", "PRJDB12280"}


# ── 工具函数 ───────────────────────────────────────────────────────────────────

def sniff_delim(path: Path) -> str:
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as f:
        s = f.read(4096)
    return "\t" if s.count("\t") > s.count(",") else ","


def read_table(path: Path) -> list[dict]:
    delim = sniff_delim(path)
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as f:
        return [dict(r) for r in csv.DictReader(f, delimiter=delim)]


def find_sra_table(dataset_dir: Path) -> Path | None:
    for name in ["SraRunTable.csv", "SraRunTable.txt", "RunTable.csv",
                 "RunSelector.csv", "SRAtable.csv"]:
        p = dataset_dir / name
        if p.exists():
            return p
    for p in list(dataset_dir.glob("*.csv")) + list(dataset_dir.glob("*.txt")):
        n = p.name.lower()
        if "sra" in n or ("run" in n and "table" in n):
            return p
    return None


def find_run_col(fieldnames: list[str]) -> str | None:
    for c in ["Run", "run_accession", "run", "Accession"]:
        if c in fieldnames:
            return c
    for f in fieldnames:
        if f.lower() in ("run", "run_accession", "accession"):
            return f
    return None


# ── 主流程 ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--meta-dir", default="data/OLP/raw/meta")
    parser.add_argument("--input",    default="results/OLP/metadata/combined_metadata.tsv")
    parser.add_argument("--output",   default="results/OLP/metadata/combined_metadata_annotated.tsv")
    args = parser.parse_args()

    meta_root = Path(args.meta_dir)
    in_path   = Path(args.input)
    if not in_path.exists():
        sys.exit(f"[ERROR] 找不到输入文件：{in_path}")
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    # 读取现有 metadata
    with open(in_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = list(reader.fieldnames)
        rows = list(reader)

    # 确保 olp_subtype 列存在（写入 notes）
    if "notes" not in fieldnames:
        fieldnames.append("notes")

    # 为每个有规则的数据集构建 run_id → (condition, subtype, note) 映射
    run_map: dict[str, tuple[str, str, str]] = {}

    for study_id, rule_fn in RULES.items():
        sra_path = find_sra_table(meta_root / study_id)
        if sra_path is None:
            print(f"[SKIP] {study_id}: 未找到 SraRunTable，跳过")
            continue

        sra_rows = read_table(sra_path)
        if not sra_rows:
            print(f"[WARN] {study_id}: SraRunTable 为空")
            continue

        run_col = find_run_col(list(sra_rows[0].keys()))
        if run_col is None:
            print(f"[WARN] {study_id}: 找不到 Run 列，跳过")
            continue

        counts = {"OLP": 0, "Control": 0, "Exclude": 0, "FILL_ME": 0}
        for sra_row in sra_rows:
            run_id = sra_row.get(run_col, "").strip()
            if not run_id:
                continue
            cond, subtype, note = rule_fn(sra_row)
            run_map[run_id] = (cond, subtype, note)
            key = cond if not cond.startswith("Exclude") else "Exclude"
            counts[key] = counts.get(key, 0) + 1

        print(f"[OK]  {study_id}: "
              f"OLP={counts['OLP']}  Control={counts['Control']}  "
              f"Exclude={counts['Exclude']}  FILL_ME={counts['FILL_ME']}")

    print()

    # 更新 combined_metadata
    cond_updated = subtype_updated = skipped_protected = still_fill = 0

    for row in rows:
        study_id = row.get("study_id", "")
        run_id   = row.get("sample_id", "").strip()

        # 保护已完成标注的数据集
        if study_id in DO_NOT_OVERWRITE:
            skipped_protected += 1
            continue

        # 只更新 FILL_ME（不覆盖已有的非 FILL_ME 标注）
        if row.get("condition", "").strip() != "FILL_ME":
            continue

        if run_id not in run_map:
            still_fill += 1
            continue

        cond, subtype, note = run_map[run_id]
        row["condition"] = cond
        cond_updated += 1

        # 追加 OLP 亚型到 notes
        existing_notes = row.get("notes", "").strip()
        parts = []
        if subtype:
            parts.append(subtype)
            subtype_updated += 1
        if note:
            parts.append(note)
        if parts:
            row["notes"] = (existing_notes + "; " if existing_notes else "") + "; ".join(parts)

    # 写输出
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    # 统计
    by_cond: dict[str, int] = {}
    for row in rows:
        c = row.get("condition", "?")
        by_cond[c] = by_cond.get(c, 0) + 1

    print(f"[OK] → {args.output}")
    print(f"     condition 更新：{cond_updated} 条")
    print(f"     亚型写入 notes：{subtype_updated} 条")
    print(f"     受保护不覆盖：  {skipped_protected} 条 (CRA008410/PRJDB12280)")
    print(f"     仍为 FILL_ME：  {still_fill} 条")
    print(f"\n     最终分组统计：")
    for k, v in sorted(by_cond.items()):
        print(f"       {k:30s}: {v}")


if __name__ == "__main__":
    main()
