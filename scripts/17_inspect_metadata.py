#!/usr/bin/env python3
"""
元数据文件探查脚本
扫描各数据集目录，列出元数据文件及其列名和前几行样例，
为后续整合脚本提供字段映射依据。

用法：
    python3 scripts/17_inspect_metadata.py \
        --meta-dir data/OLP/raw/meta \
        --output results/OLP/metadata/inspect_report.txt
"""

import argparse
import csv
import os
import sys
from pathlib import Path


META_EXTENSIONS = {".csv", ".tsv", ".txt", ".xlsx"}
SKIP_DIRS = {"__pycache__"}

# iSeq 下载时通常生成的元数据文件名模式
METADATA_KEYWORDS = [
    "metadata", "runinfo", "srarunTable", "SraRunTable",
    "run_info", "sample_info", "info", "meta"
]


def sniff_delimiter(path: Path) -> str:
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as f:
        sample = f.read(4096)
    if "\t" in sample:
        return "\t"
    return ","


def read_table(path: Path):
    """返回 (列名列表, 前3行数据列表)，失败返回 (None, None)"""
    suffix = path.suffix.lower()
    if suffix == ".xlsx":
        try:
            import openpyxl
            wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
            ws = wb.active
            rows = list(ws.iter_rows(values_only=True))
            wb.close()
            if not rows:
                return None, None
            headers = [str(c) if c is not None else "" for c in rows[0]]
            data = [[str(c) if c is not None else "" for c in r] for r in rows[1:4]]
            return headers, data
        except Exception as e:
            return None, [str(e)]

    delim = sniff_delimiter(path)
    try:
        with open(path, newline="", encoding="utf-8-sig", errors="replace") as f:
            reader = csv.reader(f, delimiter=delim)
            rows = []
            for i, row in enumerate(reader):
                rows.append(row)
                if i >= 3:
                    break
        if not rows:
            return None, None
        return rows[0], rows[1:4]
    except Exception as e:
        return None, [str(e)]


def is_metadata_file(path: Path) -> bool:
    name_lower = path.stem.lower()
    return any(kw.lower() in name_lower for kw in METADATA_KEYWORDS)


def inspect_dir(dataset_dir: Path, out):
    out.write(f"\n{'='*60}\n")
    out.write(f"数据集：{dataset_dir.name}\n")
    out.write(f"路径  ：{dataset_dir}\n")
    out.write(f"{'='*60}\n")

    found = []
    for f in sorted(dataset_dir.rglob("*")):
        if f.suffix.lower() in META_EXTENSIONS:
            found.append(f)

    if not found:
        out.write("  ⚠️  未找到任何表格文件\n")
        return

    for fpath in found:
        rel = fpath.relative_to(dataset_dir)
        marker = "★ " if is_metadata_file(fpath) else "  "
        out.write(f"\n{marker}文件：{rel}  ({fpath.stat().st_size // 1024} KB)\n")

        headers, rows = read_table(fpath)
        if headers is None:
            out.write(f"    [读取失败: {rows}]\n")
            continue

        out.write(f"    列数：{len(headers)}\n")
        out.write(f"    列名：{', '.join(headers[:20])}")
        if len(headers) > 20:
            out.write(f" ... (+{len(headers)-20}列)")
        out.write("\n")

        if rows:
            out.write(f"    样例（前{len(rows)}行）：\n")
            for row in rows:
                # 只显示前8列防止过长
                preview = row[:8]
                if len(row) > 8:
                    preview.append(f"...(+{len(row)-8}列)")
                out.write(f"      {' | '.join(str(v)[:30] for v in preview)}\n")


def main():
    parser = argparse.ArgumentParser(description="扫描各数据集元数据文件")
    parser.add_argument("--meta-dir", default="data/OLP/raw/meta",
                        help="微生物组数据根目录")
    parser.add_argument("--output", default=None,
                        help="输出报告路径（默认打印到屏幕）")
    args = parser.parse_args()

    meta_root = Path(args.meta_dir)
    if not meta_root.exists():
        sys.exit(f"[ERROR] 目录不存在：{meta_root}")

    datasets = sorted([d for d in meta_root.iterdir()
                       if d.is_dir() and d.name not in SKIP_DIRS])

    if not datasets:
        sys.exit(f"[ERROR] 未找到子目录：{meta_root}")

    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        out = open(args.output, "w", encoding="utf-8")
    else:
        out = sys.stdout

    out.write(f"元数据探查报告\n")
    out.write(f"扫描目录：{meta_root.resolve()}\n")
    out.write(f"数据集数：{len(datasets)}\n")
    out.write(f"★ = 文件名含 metadata/runinfo 等关键词（可能是主元数据文件）\n")

    for d in datasets:
        inspect_dir(d, out)

    out.write(f"\n\n{'='*60}\n")
    out.write("扫描完成\n")

    if args.output:
        out.close()
        print(f"[OK] 报告已保存：{args.output}")


if __name__ == "__main__":
    main()
