"""
步骤6（OLP专用）：下载 OLP 相关 GEO 微阵列数据集
从 GEO FTP 直接下载 series_matrix.txt 文件

数据集说明：
  GSE52130  Illumina HumanHT-12 V4.0 (GPL10558)
            7 OLP 上皮 vs 7 健康对照
            Vo et al. PLoS ONE 2021 (PMID 34506598)

  GSE38616  Affymetrix Human Gene 1.0 ST Array (GPL6244)
            7 OLP 黏膜 vs 7 健康对照
            Gong et al. Arch Oral Biol 2017 (PMID 28237528)

  GSE23558  OLP 黏膜 vs 正常对照（platform 待核实）
            Soltaninezhad et al. Curr Pharm Des 2024 (PMID 38310566)

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

DATASETS = [
    "GSE52130",   # Illumina HumanHT-12 V4.0 — OLP 上皮
    "GSE38616",   # Affymetrix HuGene-1_0-st — OLP 黏膜
    "GSE23558",   # ⚠ platform 待核实
]

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


def main():
    parser = argparse.ArgumentParser(description="下载 OLP GEO 微阵列数据集")
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
            logging.FileHandler("logs/download_olp_microarray.log"),
        ],
    )

    logging.info(f"下载 {len(DATASETS)} 个 OLP GEO 数据集 → {outdir}")

    failed = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(download_series_matrix, gse, outdir): gse
                   for gse in DATASETS}
        for future in as_completed(futures):
            gse = futures[future]
            if not future.result():
                failed.append(gse)

    if failed:
        fail_log = Path("logs/download_olp_failed.txt")
        fail_log.write_text("\n".join(failed) + "\n")
        logging.error(f"✗ {len(failed)} 个数据集失败，见 {fail_log}")
    else:
        logging.info(f"✅ 全部 {len(DATASETS)} 个 OLP 数据集下载成功")

    logging.info("\n后续步骤：")
    logging.info("  1. 提取样本元数据（⚠ OLP 关键词与 COPD 不同）：")
    logging.info("     python3 scripts/05_prepare_metadata.py \\")
    logging.info("       --datadir data/OLP/raw/host \\")
    logging.info('       --disease-keyword "lichen planus" "OLP" "erosive" \\')
    logging.info('       --control-keyword "control" "normal" "healthy"')
    logging.info("  2. 人工检查生成的 *_metadata.txt，确认 condition 列正确")
    logging.info("  3. 核实 GSE23558 的平台号，更新 config/OLP/study_config.yaml")
    logging.info("  4. 下载平台注释（若尚未下载）：")
    logging.info("     python3 scripts/02_download_platform_annotations.py")


if __name__ == "__main__":
    main()
