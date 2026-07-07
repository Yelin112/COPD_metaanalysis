#!/bin/bash
# PICRUSt2 功能预测脚本（QIIME2 输出 → EC 丰度表）
#
# 用法：
#   conda activate picrust2
#   cd /home/usr/yel/exp-OLP_meta/pipeline
#   bash scripts/13_picrust2.sh PRJNA542018
#   bash scripts/13_picrust2.sh CRA008410

set -euo pipefail

PROJ="${1:?请传入数据集名称，如：bash scripts/13_picrust2.sh PRJNA542018}"
QIIME_OUT="data/OLP/qiime2/${PROJ}/exported"
OUTDIR="data/OLP/picrust2/${PROJ}"
THREADS=8

SEQS="${QIIME_OUT}/rep-seqs/dna-sequences.fasta"
TABLE="${QIIME_OUT}/table/feature-table.biom"

mkdir -p "${OUTDIR}" logs

echo "[$(date '+%H:%M:%S')] ── ${PROJ} PICRUSt2 开始 ──"

# ── PICRUSt2 一步式流程 ───────────────────────────────────────────
picrust2_pipeline.py \
    --study_fasta  "${SEQS}" \
    --input_table  "${TABLE}" \
    --output       "${OUTDIR}" \
    --threads      ${THREADS} \
    --stratified \
    --verbose

# ── 提取 EC 丰度表 ────────────────────────────────────────────────
EC_FILE="${OUTDIR}/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz"
EC_OUT="data/OLP/raw/meta/${PROJ}.tsv"

echo "[$(date '+%H:%M:%S')] 提取 EC 丰度表..."
gunzip -c "${EC_FILE}" > "${EC_OUT}"

echo "[$(date '+%H:%M:%S')] ✅ ${PROJ} PICRUSt2 完成"
echo "  EC 丰度表 → ${EC_OUT}"
echo "  （可直接用于 Snakemake meta_metagenomics 流程）"
