"""
步骤2：下载 GEO 平台注释文件并解析为 probe_id → gene_symbol 映射表
输出：databases/platform_mappings/GPL*_mapping.txt（供流程步骤1使用）

涉及平台（来自 config/COPD/study_config.yaml）：
  GPL570   Affymetrix Human Genome U133 Plus 2.0
  GPL96    Affymetrix Human Genome U133A
  GPL14550 Agilent-028004 SurePrint G3 Human GE 8x60K Microarray
  GPL6244  Affymetrix Human Gene 1.0 ST Array
  GPL10558 Illumina HumanHT-12 V4.0 expression beadchip

用法：python3 scripts/02_download_platform_annotations.py
      python3 scripts/02_download_platform_annotations.py --outdir /your/db/path
"""

import os
import re
import gzip
import shutil
import logging
import argparse
import urllib.request
from pathlib import Path

GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/platforms"

# 平台列表：平台号 → 描述（⚠️ 请根据实际用到的数据集核实平台号）
PLATFORMS = {
    "GPL570":   "Affymetrix HG-U133_Plus_2",
    "GPL96":    "Affymetrix HG-U133A",
    "GPL14550": "Agilent-028004 SurePrint G3 Human 8x60K",
    "GPL6244":  "Affymetrix HuGene-1_0-st-v1",
    "GPL10558": "Illumina HumanHT-12 V4.0",
}


def gpl_ftp_prefix(gpl: str) -> str:
    """GPL96 → GPL0nnn，GPL570 → GPL0nnn，GPL14550 → GPL14nnn"""
    num = gpl[3:]
    return f"GPL{num[:-3]}nnn" if len(num) > 3 else f"GPL{num.zfill(4)[:-3]}nnn"


def download_and_parse_platform(gpl: str, outdir: Path) -> bool:
    """
    下载 GPL soft 文件，解析 probe_id → gene_symbol 映射
    Affymetrix/Agilent/Illumina 的 soft 文件格式不同但列名规律可循
    """
    outfile = outdir / f"{gpl}_mapping.txt"
    if outfile.exists():
        logging.info(f"跳过（已存在）: {gpl}")
        return True

    prefix = gpl_ftp_prefix(gpl)
    url = f"{GEO_FTP}/{prefix}/{gpl}/soft/{gpl}_family.soft.gz"
    gz_path = outdir / f"{gpl}_family.soft.gz"
    soft_path = outdir / f"{gpl}_family.soft"

    try:
        logging.info(f"下载 {gpl} soft 文件...")
        urllib.request.urlretrieve(url, gz_path)
        with gzip.open(gz_path, "rb") as f_in, open(soft_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        gz_path.unlink(missing_ok=True)
    except Exception as e:
        logging.error(f"✗ {gpl} 下载失败: {e}")
        return False

    # 解析 soft 文件：找到数据表头后提取 ID 和 Gene Symbol 列
    probe_gene = {}
    in_table = False
    id_col = gene_col = -1

    with open(soft_path, encoding="latin-1") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table:
                continue

            cols = line.split("\t")

            # 第一行是表头，找 ID 列和基因符号列
            if id_col == -1:
                headers_lower = [c.lower() for c in cols]
                # ID 列
                for candidate in ["id", "probe_id", "probeset_id"]:
                    if candidate in headers_lower:
                        id_col = headers_lower.index(candidate)
                        break
                if id_col == -1:
                    id_col = 0  # 默认第一列

                # Gene Symbol 列（不同平台命名不同）
                gene_candidates = [
                    "gene symbol", "gene_symbol", "genesymbol",
                    "symbol", "gene_assignment", "mrna_assignment",
                ]
                for candidate in gene_candidates:
                    if candidate in headers_lower:
                        gene_col = headers_lower.index(candidate)
                        break
                if gene_col == -1:
                    logging.warning(f"  {gpl}: 未找到基因符号列，headers={cols[:10]}")
                continue

            if len(cols) <= max(id_col, gene_col):
                continue

            probe_id = cols[id_col].strip()
            gene_raw = cols[gene_col].strip()

            if not probe_id or not gene_raw or gene_raw in ("---", "", "NA"):
                continue

            # Affymetrix gene_assignment 格式：
            # "NM_001005484 // OR4F5 // ..." → 提取基因符号
            if "//" in gene_raw:
                parts = gene_raw.split("//")
                gene = parts[1].strip() if len(parts) > 1 else parts[0].strip()
            else:
                gene = gene_raw

            # 跳过多基因映射（含 /// 的行）
            if "///" in gene:
                continue

            probe_gene[probe_id] = gene

    soft_path.unlink(missing_ok=True)

    if not probe_gene:
        logging.error(f"✗ {gpl}: 解析结果为空，请手动检查 soft 文件格式")
        return False

    with open(outfile, "w") as f:
        f.write("probe_id\tgene_symbol\n")
        for probe, gene in probe_gene.items():
            f.write(f"{probe}\t{gene}\n")

    logging.info(f"✓ {gpl}: {len(probe_gene):,} 个探针映射 → {outfile.name}")
    return True


def main():
    parser = argparse.ArgumentParser(description="下载并解析 GEO 平台注释文件")
    parser.add_argument("--outdir", default="databases/platform_mappings",
                        help="输出目录（需与 config/databases.yaml 中路径一致）")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    Path("logs").mkdir(exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        datefmt="%H:%M:%S",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler("logs/download_platform.log"),
        ],
    )

    logging.info(f"下载 {len(PLATFORMS)} 个平台注释文件 → {outdir}")

    failed = []
    for gpl, desc in PLATFORMS.items():
        logging.info(f"\n── {gpl}: {desc}")
        if not download_and_parse_platform(gpl, outdir):
            failed.append(gpl)

    if failed:
        logging.error(f"\n✗ 失败: {failed}")
        logging.error("请手动从 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL* 下载")
    else:
        logging.info(f"\n✅ 全部 {len(PLATFORMS)} 个平台注释解析完成")
        logging.info(f"请将 databases.yaml 中 probe_mappings.base_dir 设为: {outdir.resolve()}")


if __name__ == "__main__":
    main()
