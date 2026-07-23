#!/usr/bin/env python3
"""
PRJDB12280 样本筛选脚本
基于论文 Supplementary Table 2，仅保留 Presort fraction（全唾液，分选前）
样本，排除 IgA-enriched / IgA-nonenriched 分选后样本。

数据集设计说明（Kikuchi et al.）：
  该数据集为 IgA 分选测序研究，每受试者有 3 份平行：
    ① Presort fraction      ← 分选前全唾液，代表真实口腔微生物组，保留
    ② IgA-enriched fraction ← IgA 包被菌群，排除
    ③ IgA-nonenriched fraction ← 非 IgA 包被菌群，排除

保留样本（来自论文 Supplementary Table 2）：
  OLP（21 例）：SAMD01578002 – SAMD01578022
  Control（56 例）：SAMD01578047 – SAMD01578102

筛选方法：
  以 SraRunTable.csv 中 BioSample 列查找 Run accession（DRR*）

用法：
  python3 scripts/14_filter_prjdb12280.py \\
      --sra-table /path/to/SraRunTable.csv \\
      --output data/OLP/raw/meta/PRJDB12280/accession_list.txt \\
      --metadata-out data/OLP/raw/meta/PRJDB12280/metadata.tsv
"""

import argparse
import csv
import sys
from pathlib import Path

# ── 白名单：Presort fraction BioSample IDs（来自 Supplementary Table 2）──
OLP_BIOSAMPLES = {f"SAMD015780{n:02d}" for n in range(2, 23)}      # 02–22, n=21
CTRL_BIOSAMPLES = {f"SAMD015780{n:02d}" for n in range(47, 103)}   # 47–102, n=56
# SAMD01578100–102 的编号超过两位，需要单独生成
OLP_BIOSAMPLES  = {f"SAMD01578{n:03d}" for n in range(2, 23)}
CTRL_BIOSAMPLES = {f"SAMD01578{n:03d}" for n in range(47, 103)}


def main():
    parser = argparse.ArgumentParser(description="Filter PRJDB12280 to Presort fraction only")
    parser.add_argument("--sra-table", required=True,
                        help="SraRunTable.csv 路径")
    parser.add_argument("--output", required=True,
                        help="输出 accession list（每行一个 Run ID）")
    parser.add_argument("--metadata-out", default=None,
                        help="可选：同时输出带 condition 标注的 metadata TSV")
    args = parser.parse_args()

    sra_path = Path(args.sra_table)
    if not sra_path.exists():
        sys.exit(f"[ERROR] 找不到文件：{sra_path}")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    rows_olp, rows_ctrl, rows_skip = [], [], []

    with open(sra_path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        run_col  = next((c for c in reader.fieldnames if c.lower() == "run"), None)
        bio_col  = next((c for c in reader.fieldnames if c.lower() == "biosample"), None)
        if not run_col or not bio_col:
            sys.exit(f"[ERROR] 找不到 Run 或 BioSample 列。列名：{reader.fieldnames}")

        for row in reader:
            run  = row[run_col].strip()
            bio  = row[bio_col].strip()
            if bio in OLP_BIOSAMPLES:
                rows_olp.append((run, bio, "OLP"))
            elif bio in CTRL_BIOSAMPLES:
                rows_ctrl.append((run, bio, "Control"))
            else:
                rows_skip.append((run, bio))

    kept = rows_olp + rows_ctrl

    print("[筛选结果] 仅保留 Presort fraction（全唾液）")
    print(f"  OLP 患者  : {len(rows_olp)}  （SAMD01578002–022）")
    print(f"  健康对照  : {len(rows_ctrl)}  （SAMD01578047–102）")
    print(f"  排除      : {len(rows_skip)}  （IgA-enriched / IgA-nonenriched）")
    print(f"  合计保留  : {len(kept)}")

    # 检查白名单是否全部找到
    found_bio = {bio for _, bio, _ in kept}
    missing = (OLP_BIOSAMPLES | CTRL_BIOSAMPLES) - found_bio
    if missing:
        print(f"\n[WARNING] 以下 BioSample 未在 CSV 中找到（共 {len(missing)} 个）：")
        for b in sorted(missing)[:10]:
            print(f"  {b}")

    with open(args.output, "w") as fh:
        for run, bio, cond in kept:
            fh.write(run + "\n")
    print(f"\n[OK] accession list → {args.output}")

    if args.metadata_out:
        Path(args.metadata_out).parent.mkdir(parents=True, exist_ok=True)
        with open(args.metadata_out, "w") as fh:
            fh.write("sample_id\tcondition\tbiosample\n")
            for run, bio, cond in kept:
                fh.write(f"{run}\t{cond}\t{bio}\n")
        print(f"[OK] metadata     → {args.metadata_out}")


if __name__ == "__main__":
    main()
