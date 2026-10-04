#!/bin/bash
# CRA008410 QIIME2 处理脚本
# 数据类型：双端 2×250bp，V3-V4，引物已截除
#
# 用法：
#   conda activate qiime2-amplicon-2023
#   cd /home/usr/yel/exp-OLP_meta/pipeline
#   nohup bash scripts/12_qiime2_CRA008410.sh > logs/qiime2_CRA008410.log 2>&1 &

set -euo pipefail

# 让系统 R 也能找到 conda 环境中的 R 包（如 optparse），不改变 PATH
CONDA_R_LIB="$(conda info --base 2>/dev/null)/envs/qiime2-amplicon-2023/lib/R/library"
[[ -d "${CONDA_R_LIB}" ]] && export R_LIBS_SITE="${CONDA_R_LIB}"

PROJ="CRA008410"
RAWDIR="data/OLP/raw/meta/${PROJ}"
OUTDIR="data/OLP/qiime2/${PROJ}"
THREADS=8

mkdir -p "${OUTDIR}" logs

echo "[$(date '+%H:%M:%S')] ── ${PROJ} QIIME2 流程开始 ──"

# ── 步骤1：生成 paired-end manifest ──────────────────────────────
echo "[$(date '+%H:%M:%S')] 生成 manifest..."
MANIFEST="${OUTDIR}/manifest.tsv"
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "${MANIFEST}"

for fwd in "${RAWDIR}"/CRR*_f1.fq.gz; do
    sample=$(basename "${fwd}" _f1.fq.gz)
    rev="${RAWDIR}/${sample}_r2.fq.gz"
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

# ── 步骤3：DADA2 去噪（引物已截除，直接去噪）─────────────────────
# V3-V4 amplicon 去引物后 ~428bp；reads 固定 250bp，无需 cutadapt
# 原参数 230+200=430 - 428 = 2bp overlap < 12bp 最低要求 → 几乎全部 merge 失败
# 修正：240+240=480 - 428 = 52bp overlap ✓，均在 250bp read 长度以内
echo "[$(date '+%H:%M:%S')] DADA2 去噪（双端）..."
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "${OUTDIR}/demux.qza" \
    --p-trunc-len-f 240 \
    --p-trunc-len-r 240 \
    --p-trim-left-f 0 \
    --p-trim-left-r 0 \
    --p-max-ee-f 2.0 \
    --p-max-ee-r 2.0 \
    --p-n-threads ${THREADS} \
    --o-table "${OUTDIR}/table.qza" \
    --o-representative-sequences "${OUTDIR}/rep-seqs.qza" \
    --o-denoising-stats "${OUTDIR}/dada2-stats.qza"

# ── 步骤4：导出 ───────────────────────────────────────────────────
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
