#!/bin/bash
# PRJNA1049117 QIIME2 流程：双端 2×250bp，V3-V4，343F/798R，混合样本类型
# 数据类型：双端 paired-end，引物未截除
# 注意：798R 序列请以原文引物信息核对
#
# 用法：
#   conda activate qiime2-amplicon-2023
#   cd /home/usr/yel/exp-OLP_meta/pipeline
#   nohup bash scripts/28_qiime2_PRJNA1049117.sh > logs/qiime2_PRJNA1049117.log 2>&1 &

set -euo pipefail

CONDA_R_LIB="$(conda info --base 2>/dev/null)/envs/qiime2-amplicon-2023/lib/R/library"
[[ -d "${CONDA_R_LIB}" ]] && export R_LIBS_SITE="${CONDA_R_LIB}"

PROJ="PRJNA1049117"
RAWDIR="data/OLP/raw/meta/${PROJ}"
OUTDIR="data/OLP/qiime2/${PROJ}"
THREADS=8

PRIMER_F="TACGGRAGGCAGCAG"     # 343F
PRIMER_R="AGGGTATCTAATCCT"     # 798R（请核对原文）

mkdir -p "${OUTDIR}" logs

echo "[$(date '+%H:%M:%S')] ── ${PROJ} QIIME2 流程开始 ──"

# ── 步骤1：生成 manifest ───────────────────────────────────────────
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

# ── 步骤2：导入数据 ───────────────────────────────────────────────
echo "[$(date '+%H:%M:%S')] 导入 FASTQ..."
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "${MANIFEST}" \
    --input-format PairedEndFastqManifestPhred33V2 \
    --output-path "${OUTDIR}/demux.qza"

# ── 步骤3：cutadapt 截除引物 ──────────────────────────────────────
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

# ── 步骤4：DADA2 去噪 ────────────────────────────────────────────
# V3-V4 amplicon ~428bp（去引物后），250bp 双端，overlap 较紧
# trunc-f=240 + trunc-r=200 → overlap = 440-428 = 12bp（最低可用）
echo "[$(date '+%H:%M:%S')] DADA2 去噪（双端）..."
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "${OUTDIR}/demux_trimmed.qza" \
    --p-trunc-len-f 240 \
    --p-trunc-len-r 200 \
    --p-trim-left-f 0 \
    --p-trim-left-r 0 \
    --p-max-ee-f 2.0 \
    --p-max-ee-r 2.0 \
    --p-n-threads ${THREADS} \
    --o-table "${OUTDIR}/table.qza" \
    --o-representative-sequences "${OUTDIR}/rep-seqs.qza" \
    --o-denoising-stats "${OUTDIR}/dada2-stats.qza"

# ── 步骤5：导出 ───────────────────────────────────────────────────
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
