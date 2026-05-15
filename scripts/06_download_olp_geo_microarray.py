"""
步骤6（OLP专用）：下载 OLP 相关 GEO 宿主转录组数据集

微阵列（series_matrix.txt）：
  GSE52130  Illumina HumanHT-12 V4.0 (GPL10558)
            7 口腔 OLP 上皮 vs 7 健康对照（含生殖器 LP，需过滤）
            Vo et al. PLoS ONE 2021 (PMID 34506598)

  GSE38616  Affymetrix Human Gene 1.0 ST Array (GPL6244)
            7 OLP 黏膜 vs 7 健康对照
            Gong et al. Arch Oral Biol 2017 (PMID 28237528)

Bulk RNA-seq（GEO 补充文件 GSE213346_raw_counts.txt.gz）：
  GSE213346 Illumina NovaSeq 6000
            40 OLP 口腔黏膜 vs 10 健康对照
            Cillo et al. eLife 2023 (PMID 37874637)

用法：python3 scripts/06_download_olp_geo_microarray.py
      python3 scripts/06_download_olp_geo_microarray.py --outdir data/OLP/raw/host
"""

import re
import gzip
import shutil
import logging
import argparse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

MICROARRAY_DATASETS = [
    "GSE52130",   # Illumina HumanHT-12 V4.0 — 口腔 OLP 上皮
    "GSE38616",   # Affymetrix HuGene-1_0-st — OLP 黏膜
]

# GSE213346 GEO 补充文件（RNA-seq 原始计数矩阵）
RNASEQ_SUPPLEMENTARY = {
    "GSE213346": "GSE213346_count.txt.gz",
}

GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"


def geo_ftp_prefix(gse: str) -> str:
    num = gse[3:]
    return f"GSE{num[:-3]}nnn" if len(num) > 3 else f"GSE{num.zfill(4)[:-3]}nnn"


def download_series_matrix(gse: str, outdir: Path, retries: int = 3) -> bool:
    prefix = geo_ftp_prefix(gse)
    base_url = f"{GEO_FTP}/{prefix}/{gse}/matrix"
    outfile = outdir / f"{gse}_series_matrix.txt"

    if outfile.exists():
        logging.info(f"跳过（已存在）: {gse}")
        return True

    candidate_urls = [f"{base_url}/{gse}_series_matrix.txt.gz"]

    try:
        with urllib.request.urlopen(f"{base_url}/", timeout=20) as resp:
            html = resp.read().decode()
        found = re.findall(r'href="([^"]*series_matrix[^"]*\.gz)"', html)
        for f in found:
            url = f"{base_url}/{f}" if not f.startswith("http") else f
            if url not in candidate_urls:
                candidate_urls.append(url)
    except Exception:
        pass

    for url in candidate_urls:
        gz_path = outfile.with_suffix(".txt.gz")
        for attempt in range(1, retries + 1):
            try:
                logging.info(f"下载 {gse} (尝试 {attempt}/{retries}): {url}")
                urllib.request.urlretrieve(url, gz_path)
                with gzip.open(gz_path, "rb") as f_in, open(outfile, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                gz_path.unlink(missing_ok=True)
                logging.info(f"✓ {gse} → {outfile.name}")
                return True
            except Exception as e:
                logging.warning(f"  失败: {e}")
                gz_path.unlink(missing_ok=True)

    logging.error(f"✗ {gse} 所有尝试均失败")
    return False


def download_rnaseq_supplementary(gse: str, filename: str,
                                   outdir: Path, retries: int = 3) -> bool:
    """下载 GEO 补充文件（RNA-seq 原始计数矩阵）"""
    prefix = geo_ftp_prefix(gse)
    base_url = f"{GEO_FTP}/{prefix}/{gse}/suppl"
    url = f"{base_url}/{filename}"

    # 输出文件名去掉 .gz 后缀（自动解压）
    stem = filename[:-3] if filename.endswith(".gz") else filename
    outfile = outdir / stem

    if outfile.exists():
        logging.info(f"跳过（已存在）: {gse} 补充文件")
        return True

    gz_path = outdir / filename
    for attempt in range(1, retries + 1):
        try:
            logging.info(f"下载 {gse} 补充文件 (尝试 {attempt}/{retries}): {url}")
            urllib.request.urlretrieve(url, gz_path)
            with gzip.open(gz_path, "rb") as f_in, open(outfile, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            gz_path.unlink(missing_ok=True)
            logging.info(f"✓ {gse} 补充文件 → {outfile.name}")
            return True
        except Exception as e:
            logging.warning(f"  失败: {e}")
            gz_path.unlink(missing_ok=True)

    logging.error(f"✗ {gse} 补充文件所有尝试均失败")
    return False


def main():
    parser = argparse.ArgumentParser(description="下载 OLP GEO 宿主转录组数据集")
    parser.add_argument("--threads", type=int, default=3)
    parser.add_argument("--outdir", default="data/OLP/raw/host")
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
            logging.FileHandler("logs/download_olp_host.log"),
        ],
    )

    total = len(MICROARRAY_DATASETS) + len(RNASEQ_SUPPLEMENTARY)
    logging.info(f"下载 {total} 个 OLP 宿主转录组数据集 → {outdir}")

    failed = []

    # 微阵列 series_matrix 并发下载
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(download_series_matrix, gse, outdir): gse
                   for gse in MICROARRAY_DATASETS}
        futures.update({
            executor.submit(download_rnaseq_supplementary, gse, fn, outdir): gse
            for gse, fn in RNASEQ_SUPPLEMENTARY.items()
        })
        for future in as_completed(futures):
            gse = futures[future]
            if not future.result():
                failed.append(gse)

    if failed:
        fail_log = Path("logs/download_olp_failed.txt")
        fail_log.write_text("\n".join(failed) + "\n")
        logging.error(f"✗ {len(failed)} 个数据集失败，见 {fail_log}")
    else:
        logging.info(f"✅ 全部 {total} 个 OLP 宿主转录组数据集下载成功")

    logging.info("\n后续步骤（提取元数据）：")
    logging.info("  # GSE52130：含口腔+生殖器样本，必须精确匹配 'oral lichen planus'")
    logging.info("  python3 scripts/05_prepare_metadata.py \\")
    logging.info("    --matrix  data/OLP/raw/host/GSE52130_series_matrix.txt \\")
    logging.info("    --output  data/OLP/raw/host/GSE52130_metadata.txt \\")
    logging.info('    --disease-keyword "oral lichen planus" \\')
    logging.info('    --control-keyword "control oral"')
    logging.info("")
    logging.info("  # GSE38616：全部为口腔样本，关键词较宽松")
    logging.info("  python3 scripts/05_prepare_metadata.py \\")
    logging.info("    --matrix  data/OLP/raw/host/GSE38616_series_matrix.txt \\")
    logging.info("    --output  data/OLP/raw/host/GSE38616_metadata.txt \\")
    logging.info('    --disease-keyword "oral lichen planus" \\')
    logging.info('    --control-keyword "healthy"')
    logging.info("")
    logging.info("  # GSE213346 RNA-seq metadata（从 GEO SRA RunTable 获取）：")
    logging.info("  python3 scripts/09_prepare_rnaseq_metadata.py \\")
    logging.info("    --geo    GSE213346 \\")
    logging.info("    --output data/OLP/raw/host/GSE213346_metadata.txt \\")
    logging.info('    --disease-keyword "oral lichen planus" "OLP" \\')
    logging.info('    --control-keyword "healthy" "control"')


if __name__ == "__main__":
    main()
