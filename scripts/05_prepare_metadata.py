"""
步骤5：从 GEO series_matrix.txt 中自动提取样本元数据
生成流程所需的 {GSE}_metadata.txt（含 sample_id 和 condition 列）

用法：python3 scripts/05_prepare_metadata.py
      python3 scripts/05_prepare_metadata.py --datadir data/COPD/raw/host
              --disease-keyword COPD --control-keyword normal

⚠ 重要：GEO 数据集的 condition 标注格式各不相同，脚本会尽力自动识别，
  但建议运行后人工检查每个数据集的 metadata 文件，确认 condition 列正确。
"""

import os
import re
import logging
import argparse
from pathlib import Path


# 常见的 COPD 相关关键词（不区分大小写）
DISEASE_KEYWORDS = ["copd", "chronic obstructive", "emphysema", "smoker with airflow"]
CONTROL_KEYWORDS = ["control", "normal", "healthy", "non-copd", "nonsmoker", "non-smoker",
                    "smoker without", "at-risk"]


def parse_series_matrix_meta(filepath: Path):
    """
    解析 series_matrix.txt 中 ! 开头的元数据行
    两遍读取：第一遍获取样本ID，第二遍获取 title 和 characteristics
    （GEO 文件中 !Sample_geo_accession 不一定在 !Sample_title 之前）
    """
    # 第一遍：只收集所有 ! 行
    meta_lines = []
    with open(filepath, encoding="latin-1") as f:
        for line in f:
            line = line.rstrip("\n").replace('"', "")
            if line.startswith("!"):
                meta_lines.append(line)

    # 获取样本 ID 列表
    sample_ids = []
    for line in meta_lines:
        if line.startswith("!Sample_geo_accession"):
            sample_ids = line.split("\t")[1:]
            break

    # 第二遍：按列对齐解析 title 和 characteristics
    titles = {}
    characteristics = {}
    for line in meta_lines:
        if line.startswith("!Sample_title"):
            parts = line.split("\t")
            for i, sid in enumerate(sample_ids):
                titles[sid] = parts[i + 1] if i + 1 < len(parts) else ""
        elif line.startswith("!Sample_characteristics_ch1"):
            parts = line.split("\t")
            for i, sid in enumerate(sample_ids):
                val = parts[i + 1] if i + 1 < len(parts) else ""
                characteristics.setdefault(sid, []).append(val)

    return sample_ids, titles, characteristics


def infer_condition(title: str, chars: list[str],
                    disease_kw: list[str], control_kw: list[str],
                    disease_label: str, control_label: str) -> str:
    """
    从样本标题和 characteristics 中推断 condition
    返回 disease_label、control_label 或 'Unknown'
    """
    text = " ".join([title] + chars).lower()

    for kw in disease_kw:
        if kw.lower() in text:
            return disease_label
    for kw in control_kw:
        if kw.lower() in text:
            return control_label
    return "Unknown"


def process_dataset(matrix_file: Path, disease_kw: list[str],
                    control_kw: list[str], disease_label: str,
                    control_label: str, force: bool = False) -> Path:
    gse = matrix_file.stem.replace("_series_matrix", "")
    outfile = matrix_file.parent / f"{gse}_metadata.txt"

    if outfile.exists() and not force:
        logging.info(f"跳过（已存在）: {gse}_metadata.txt  （用 --force 强制重新生成）")
        return outfile

    sample_ids, titles, chars = parse_series_matrix_meta(matrix_file)

    rows = []
    unknown_count = 0
    for sid in sample_ids:
        condition = infer_condition(
            titles.get(sid, ""),
            chars.get(sid, []),
            disease_kw, control_kw,
            disease_label, control_label
        )
        if condition == "Unknown":
            unknown_count += 1
        rows.append({"sample_id": sid, "condition": condition})

    with open(outfile, "w") as f:
        f.write("sample_id\tcondition\n")
        for row in rows:
            f.write(f"{row['sample_id']}\t{row['condition']}\n")

    status = "⚠ 需检查" if unknown_count > 0 else "✓"
    logging.info(
        f"{status} {gse}: {len(rows)} 样本，"
        f"{disease_label}={sum(1 for r in rows if r['condition']==disease_label)}, "
        f"{control_label}={sum(1 for r in rows if r['condition']==control_label)}, "
        f"Unknown={unknown_count}"
    )

    if unknown_count > 0:
        logging.warning(
            f"  {gse} 有 {unknown_count} 个样本未能自动识别 condition，"
            f"请手动编辑 {outfile}"
        )

    return outfile


def main():
    parser = argparse.ArgumentParser(description="从 GEO series_matrix 提取样本元数据")
    parser.add_argument("--datadir", default="data/COPD/raw/host")
    parser.add_argument("--disease-keyword", nargs="+", default=DISEASE_KEYWORDS)
    parser.add_argument("--control-keyword", nargs="+", default=CONTROL_KEYWORDS)
    parser.add_argument("--disease-label", default="COPD",
                        help="疾病样本的 condition 标签（默认 COPD，OLP 研究传入 OLP）")
    parser.add_argument("--control-label", default="Control",
                        help="对照样本的 condition 标签（默认 Control）")
    parser.add_argument("--force", action="store_true",
                        help="强制重新生成已存在的 metadata 文件")
    args = parser.parse_args()

    datadir = Path(args.datadir)
    Path("logs").mkdir(exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        datefmt="%H:%M:%S",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler("logs/prepare_metadata.log"),
        ],
    )

    matrix_files = sorted(datadir.glob("*_series_matrix.txt"))
    if not matrix_files:
        logging.error(f"在 {datadir} 中未找到 series_matrix.txt 文件，请先运行步骤1")
        return

    logging.info(f"处理 {len(matrix_files)} 个数据集...")
    logging.info(f"疾病关键词: {args.disease_keyword}  → 标签: {args.disease_label}")
    logging.info(f"对照关键词: {args.control_keyword}  → 标签: {args.control_label}")

    for mf in matrix_files:
        process_dataset(mf, args.disease_keyword, args.control_keyword,
                        args.disease_label, args.control_label, args.force)

    logging.info("\n✅ 元数据提取完成")
    logging.info("⚠  请务必检查标注为 'Unknown' 的样本，手动修改 condition 列")
    logging.info("   condition 列的值必须与 study_config.yaml 中完全一致：")
    logging.info(f"   disease_label: {args.disease_label}    control_label: {args.control_label}")


if __name__ == "__main__":
    main()
