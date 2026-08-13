#!/bin/bash
# PRJNA306560 QIIME2 流程：双端 2×251bp，V4，515F/806R，唾液
# 数据类型：双端 paired-end，引物未截除
#
# 用法：
#   conda activate qiime2-amplicon-2023
#   cd /home/usr/yel/exp-OLP_meta/pipeline
#   nohup bash scripts/22_qiime2_PRJNA306560.sh > logs/qiime2_PRJNA306560.log 2>&1 &

set -euo pipefail

CONDA_R_LIB="$(conda info --base 2>/dev/null)/envs/qiime2-amplicon-2023/lib/R/library"
[[ -d "${CONDA_R_LIB}" ]] && export R_LIBS_SITE="${CONDA_R_LIB}"

PROJ="PRJNA306560"
RAWDIR="data/OLP/raw/meta/${PROJ}"
OUTDIR="data/OLP/qiime2/${PROJ}"
THREADS=8

PRIMER_F="GTGYCAGCMGCCGCGGTAA"
PRIMER_R="GGACTACNVGGGTWTCTAAT"

mkdir -p "${OUTDIR}" logs

echo "[$(date '+%H:%M:%S')] ── ${PROJ} QIIME2 流程开始 ──"

# Step 1: manifest
echo "[$(date '+%H:%M:%S')] 生成 manifest..."
MANIFEST="${OUTDIR}/manifest.tsv"
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "${MANIFEST}"
for fwd in "${RAWDIR}"/SRR*_1.fastq.gz; do
    sample=$(basename "${fwd}" _1.fastq.gz)
    rev="${RAWDIR}/${sample}_2.fastq.gz"
    if [[ -f "${rev}" ]]; then
        echo -e "${sample}\t$(realpath ${fwd})\t$(realpath ${rev})"
    fi
done >> "${MANIFEST}"
echo "  样本数：$(tail -n +2 ${MANIFEST} | wc -l)"

# Step 2: import
echo "[$(date '+%H:%M:%S')] 导入 FASTQ..."
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "${MANIFEST}" \
    --input-format PairedEndFastqManifestPhred33V2 \
    --output-path "${OUTDIR}/demux.qza"

# Step 3: cutadapt
echo "[$(date '+%H:%M:%S')] cutadapt 截除引物..."
qiime cutadapt trim-paired \
    --i-demultiplexed-sequences "${OUTDIR}/demux.qza" \
    --p-front-f "${PRIMER_F}" \
    --p-front-r "${PRIMER_R}" \
    --p-discard-untrimmed \
    --p-minimum-length 100 \
    --p-cores ${THREADS} \
    --o-trimmed-sequences "${OUTDIR}/demux_trimmed.qza" \
    --verbose 2>&1 | tail -20

# Step 4: DADA2
# V4 扩增子 ~253bp（去引物后），251bp reads overlap 充足（230+150-253=127bp）
# 原始 reads 开头质量偏低（Phred 16-24），max-ee 由 2.0 放宽至 5.0 提升保留率
echo "[$(date '+%H:%M:%S')] DADA2 去噪（双端）..."
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "${OUTDIR}/demux_trimmed.qza" \
    --p-trunc-len-f 230 \
    --p-trunc-len-r 150 \
    --p-trim-left-f 0 \
    --p-trim-left-r 0 \
    --p-max-ee-f 5.0 \
    --p-max-ee-r 5.0 \
    --p-n-threads ${THREADS} \
    --o-table "${OUTDIR}/table.qza" \
    --o-representative-sequences "${OUTDIR}/rep-seqs.qza" \
    --o-denoising-stats "${OUTDIR}/dada2-stats.qza"

# Step 5: export
echo "[$(date '+%H:%M:%S')] 导出结果..."
qiime tools export --input-path "${OUTDIR}/table.qza" \
    --output-path "${OUTDIR}/exported/table"
qiime tools export --input-path "${OUTDIR}/rep-seqs.qza" \
    --output-path "${OUTDIR}/exported/rep-seqs"
biom convert \
    -i "${OUTDIR}/exported/table/feature-table.biom" \
    -o "${OUTDIR}/exported/table/feature-table.tsv" \
    --to-tsv

echo "[$(date '+%H:%M:%S')] ✅ ${PROJ} QIIME2 完成"
echo "  feature table : ${OUTDIR}/exported/table/feature-table.tsv"
echo "  rep seqs      : ${OUTDIR}/exported/rep-seqs/dna-sequences.fasta"
echo "  下一步：bash scripts/13_picrust2.sh ${PROJ}"
