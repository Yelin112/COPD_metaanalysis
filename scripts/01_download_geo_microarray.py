"""
步骤1：下载 GEO 微阵列 series_matrix.txt 文件
从 GEO FTP 直接下载（无需 iSeq，iSeq 下载的是原始 CEL 文件）

用法：python3 scripts/01_download_geo_microarray.py
      python3 scripts/01_download_geo_microarray.py --threads 4 --outdir data/COPD/raw/host
"""

import os
import re
import gzip
import shutil
import logging
import argparse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# ── 配置 ──────────────────────────────────────────────────────────
DATASETS = [
    "GSE103174", "GSE106986", "GSE112260", "GSE11784", "GSE119040",
    "GSE11906",  "GSE12472",  "GSE13896",  "GSE16972", "GSE37147",
    "GSE37768",  "GSE38974",  "GSE47460",  "GSE56341", "GSE57148",
    "GSE73395",  "GSE76925",  "GSE8581",   "GSE86064",
]

GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"


def geo_ftp_prefix(gse: str) -> str:
    """GSE11906 → GSE11nnn，GSE103174 → GSE103nnn"""
    num = gse[3:]
    return f"GSE{num[:-3]}nnn" if len(num) > 3 else f"GSE{num[:-len(num)]}nnn"


def download_series_matrix(gse: str, outdir: Path, retries: int = 3) -> bool:
    prefix = geo_ftp_prefix(gse)
    base_url = f"{GEO_FTP}/{prefix}/{gse}/matrix"
    outfile = outdir / f"{gse}_series_matrix.txt"

    if outfile.exists():
        logging.info(f"跳过（已存在）: {gse}")
        return True

    # 先尝试标准文件名，再列目录找备用文件名
    candidate_urls = [f"{base_url}/{gse}_series_matrix.txt.gz"]

    # 列目录获取实际文件名（应对多平台数据集）
    try:
        with urllib.request.urlopen(f"{base_url}/", timeout=20) as resp:
            html = resp.read().decode()
        found = re.findall(r'href="([^"]*series_matrix[^"]*\.gz)"', html)
        for f in found:
            url = f"{base_url}/{f}" if not f.startswith("http") else f
            if url not in candidate_urls:
                candidate_urls.append(url)
    except Exception:
        pass  # 使用默认候选 URL

    for url in candidate_urls:
        gz_path = outfile.with_suffix(".txt.gz")
        for attempt in range(1, retries + 1):
            try:
                logging.info(f"下载 {gse} (尝试 {attempt}/{retries}): {url}")
                urllib.request.urlretrieve(url, gz_path)
                # 解压
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
    parser = argparse.ArgumentParser(description="下载 GEO series_matrix 文件")
    parser.add_argument("--threads", type=int, default=4, help="并发下载线程数（默认 4）")
    parser.add_argument("--outdir", default="data/COPD/raw/host", help="输出目录")
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
            logging.FileHandler("logs/download_geo_microarray.log"),
        ],
    )

    logging.info(f"开始下载 {len(DATASETS)} 个 GEO 微阵列数据集，并发数 {args.threads}")

    failed = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(download_series_matrix, gse, outdir): gse
                   for gse in DATASETS}
        for future in as_completed(futures):
            gse = futures[future]
            if not future.result():
                failed.append(gse)

    if failed:
        fail_log = Path("logs/download_failed_geo.txt")
        fail_log.write_text("\n".join(failed) + "\n")
        logging.error(f"✗ {len(failed)} 个数据集下载失败，详见 {fail_log}")
    else:
        logging.info(f"✅ 全部 {len(DATASETS)} 个数据集下载成功")


if __name__ == "__main__":
    main()
