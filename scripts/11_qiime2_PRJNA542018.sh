#!/bin/bash
# PRJNA542018 QIIME2 处理脚本
# 数据类型：merged paired-end（以单端格式提交），V4，515F/806R 引物需截除
#
# 用法：
#   conda activate qiime2-amplicon-2023
#   cd /home/usr/yel/exp-OLP_meta/pipeline
#   nohup bash scripts/11_qiime2_PRJNA542018.sh > logs/qiime2_PRJNA542018.log 2>&1 &

set -euo pipefail

PROJ="PRJNA542018"
RAWDIR="data/OLP/raw/meta/${PROJ}"
OUTDIR="data/OLP/qiime2/${PROJ}"
THREADS=8

PRIMER_F="GTGYCAGCMGCCGCGGTAA"   # 515F
PRIMER_R="GGACTACNVGGGTWTCTAAT"  # 806R

mkdir -p "${OUTDIR}" logs

echo "[$(date '+%H:%M:%S')] ── ${PROJ} QIIME2 流程开始 ──"

# ── 步骤1：生成 manifest 文件 ─────────────────────────────────────
echo "[$(date '+%H:%M:%S')] 生成 manifest..."
MANIFEST="${OUTDIR}/manifest.tsv"
echo -e "sample-id\tabsolute-filepath" > "${MANIFEST}"
for fq in "${RAWDIR}"/SRR*.fastq.gz; do
    sample=$(basename "${fq}" .fastq.gz)
    echo -e "${sample}\t$(realpath ${fq})"
done >> "${MANIFEST}"
echo "  样本数：$(tail -n +2 ${MANIFEST} | wc -l)"

# ── 步骤2：导入数据 ───────────────────────────────────────────────
echo "[$(date '+%H:%M:%S')] 导入 FASTQ..."
qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path "${MANIFEST}" \
    --input-format SingleEndFastqManifestPhred33V2 \
    --output-path "${OUTDIR}/demux.qza"

# ── 步骤3：cutadapt 截除引物 ──────────────────────────────────────
echo "[$(date '+%H:%M:%S')] cutadapt 截除引物..."
qiime cutadapt trim-single \
    --i-demultiplexed-sequences "${OUTDIR}/demux.qza" \
    --p-front "${PRIMER_F}" \
    --p-adapter "${PRIMER_R}" \
    --p-discard-untrimmed \
    --p-minimum-length 100 \
    --p-cores ${THREADS} \
    --o-trimmed-sequences "${OUTDIR}/demux_trimmed.qza" \
    --verbose 2>&1 | tail -20

# ── 步骤4：DADA2 去噪 ────────────────────────────────────────────
# trunc-len 0：不截断（merged reads 长度可变，截断会丢失信息）
echo "[$(date '+%H:%M:%S')] DADA2 去噪..."
qiime dada2 denoise-single \
    --i-demultiplexed-seqs "${OUTDIR}/demux_trimmed.qza" \
    --p-trunc-len 0 \
    --p-trim-left 0 \
    --p-max-ee 2.0 \
    --p-n-threads ${THREADS} \
    --o-table "${OUTDIR}/table.qza" \
    --o-representative-sequences "${OUTDIR}/rep-seqs.qza" \
    --o-denoising-stats "${OUTDIR}/dada2-stats.qza"

# ── 步骤5：导出 feature table 和代表序列 ──────────────────────────
echo "[$(date '+%H:%M:%S')] 导出结果..."
qiime tools export --input-path "${OUTDIR}/table.qza" \
    --output-path "${OUTDIR}/exported/table"
qiime tools export --input-path "${OUTDIR}/rep-seqs.qza" \
    --output-path "${OUTDIR}/exported/rep-seqs"

# biom 转 tsv
biom convert \
    -i "${OUTDIR}/exported/table/feature-table.biom" \
    -o "${OUTDIR}/exported/table/feature-table.tsv" \
    --to-tsv

echo "[$(date '+%H:%M:%S')] ✅ ${PROJ} QIIME2 完成"
echo "  feature table : ${OUTDIR}/exported/table/feature-table.tsv"
echo "  rep seqs      : ${OUTDIR}/exported/rep-seqs/dna-sequences.fasta"
echo "  下一步：bash scripts/13_picrust2.sh ${PROJ}"
