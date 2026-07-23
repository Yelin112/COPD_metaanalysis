#!/usr/bin/env python3
"""
PRJDB12280 样本筛选脚本
从 SraRunTable.csv 中提取 OLP 患者和健康对照的 Run accession，
生成 iSeq 下载列表。

PRJDB12280 样本分组说明：
  OLP 患者：
    - OLP-1 ~ OLP-18       （命名直接含 OLP）
    - OLP_11               （命名含 OLP_）
    - O_OLP-L2/L3/L4       （OLP 病损部位）
  健康对照：
    - Cont-1 ~ Cont-14     （Cont 前缀）
    - O_C_Control_1~30     （口腔健康对照）
  排除组：
    - O_M0_Control_*       （月经周期相关，非普通健康对照）
    - O_M1_Control_*       （月经周期相关，非普通健康对照）
    - O_C_HypoT_*          （甲状腺功能减退，不是 OLP 对照）
    - O_M0_HypoT_*         （甲状腺功能减退）
    - O_M1_HypoT_*         （甲状腺功能减退）
    - 76008 等纯数字/字母ID （来源不明）
    - A* 系列               （来源不明）

用法：
  python3 scripts/14_filter_prjdb12280.py \\
      --sra-table /path/to/SraRunTable.csv \\
      --output data/OLP/raw/meta/PRJDB12280/accession_list.txt \\
      [--include-menstrual]   # 加此参数则同时保留 M0/M1 对照
"""

import argparse
import csv
import re
import sys
from pathlib import Path


def classify_sample(name: str, include_menstrual: bool = False) -> str:
    """返回 'OLP'、'Control' 或 None（排除）"""
    n = name.strip()

    # ── OLP 患者 ──────────────────────────────────────────────────────
    if re.match(r"^OLP[-_]", n, re.I):
        return "OLP"
    if re.match(r"^O_OLP", n, re.I):
        return "OLP"

    # ── 健康对照 ──────────────────────────────────────────────────────
    if re.match(r"^Cont[-_]", n, re.I):
        return "Control"
    if re.match(r"^O_C_Control", n, re.I):
        return "Control"

    # 月经周期对照（可选保留）
    if include_menstrual and re.match(r"^O_M[01]_Control", n, re.I):
        return "Control"

    # ── 排除 ──────────────────────────────────────────────────────────
    return None


def main():
    parser = argparse.ArgumentParser(description="Filter PRJDB12280 SraRunTable.csv")
    parser.add_argument("--sra-table", required=True,
                        help="SraRunTable.csv 路径")
    parser.add_argument("--output", required=True,
                        help="输出 accession list（每行一个 Run ID）")
    parser.add_argument("--metadata-out", default=None,
                        help="可选：同时输出带 condition 标注的 metadata TSV")
    parser.add_argument("--include-menstrual", action="store_true",
                        help="保留 O_M0/M1_Control 月经周期对照样本")
    args = parser.parse_args()

    sra_path = Path(args.sra_table)
    if not sra_path.exists():
        sys.exit(f"[ERROR] 找不到文件：{sra_path}")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    rows_olp, rows_ctrl, rows_skip = [], [], []

    with open(sra_path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        # 兼容 Run / run_accession 两种列名
        run_col   = next((c for c in reader.fieldnames if c.lower() in ("run", "run_accession")), None)
        # 优先匹配 Sample_name（下划线，含分组名），避免误选 Sample Name（空格，biosample ID）
        name_col  = next((c for c in reader.fieldnames if c == "Sample_name"), None) or \
                    next((c for c in reader.fieldnames if c.lower() in ("sample_name", "sample name")), None)
        if not run_col or not name_col:
            sys.exit(f"[ERROR] 找不到 Run 或 Sample_name 列。列名：{reader.fieldnames}")

        for row in reader:
            run  = row[run_col].strip()
            name = row[name_col].strip()
            cond = classify_sample(name, args.include_menstrual)
            if cond == "OLP":
                rows_olp.append((run, name, cond))
            elif cond == "Control":
                rows_ctrl.append((run, name, cond))
            else:
                rows_skip.append((run, name))

    kept = rows_olp + rows_ctrl

    print(f"[筛选结果]")
    print(f"  OLP 患者  : {len(rows_olp)}")
    print(f"  健康对照  : {len(rows_ctrl)}")
    print(f"  排除      : {len(rows_skip)}")
    print(f"  合计保留  : {len(kept)}")

    if rows_skip:
        print("\n[排除样本列举（前20）]")
        for run, name in rows_skip[:20]:
            print(f"  {run}\t{name}")
        if len(rows_skip) > 20:
            print(f"  ... 共 {len(rows_skip)} 个")

    # 写 accession list
    with open(args.output, "w") as fh:
        for run, name, cond in kept:
            fh.write(run + "\n")
    print(f"\n[OK] accession list → {args.output}")

    # 可选：写 metadata TSV
    if args.metadata_out:
        Path(args.metadata_out).parent.mkdir(parents=True, exist_ok=True)
        with open(args.metadata_out, "w") as fh:
            fh.write("sample_id\tcondition\tsample_name\n")
            for run, name, cond in kept:
                fh.write(f"{run}\t{cond}\t{name}\n")
        print(f"[OK] metadata     → {args.metadata_out}")


if __name__ == "__main__":
    main()
