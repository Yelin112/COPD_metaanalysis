#!/usr/bin/env bash
# =============================================================================
# 步骤4：QIIME2 质控 + PICRUSt2 功能预测
# 输入：步骤3 下载的 FASTQ 文件
# 输出：taxonomy.tsv + pred_metagenome_strat.tsv（流程直接使用）
#
# 依赖安装：
#   conda create -n qiime2 --file https://data.qiime2.org/distro/core/qiime2-2024.2-py39-linux-conda.yml
#   conda install -n qiime2 -c bioconda picrust2
#
# 用法：bash scripts/04_run_qiime2_picrust2.sh SRP066375
# =============================================================================

set -euo pipefail

STUDY=${1:?"用法: bash $0 <STUDY_ID>，如 SRP066375"}
FASTQ_DIR="data/COPD/raw/meta/fastq/${STUDY}"
META_OUTDIR="data/COPD/raw/meta"
THREADS=${2:-16}
LOG="logs/qiime2_picrust2_${STUDY}.log"

mkdir -p "$META_OUTDIR" logs

log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

log "═══ QIIME2 + PICRUSt2 处理：$STUDY ═══"

# ── QIIME2：导入数据 ──────────────────────────────────────────────
log "1/5 导入 FASTQ 数据..."
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$FASTQ_DIR" \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt \
    --output-path "tmp/${STUDY}_demux.qza"

# ── QIIME2：质控 + 去噪（DADA2）────────────────────────────────
log "2/5 DADA2 去噪（此步骤较慢，约 1-4 小时）..."
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "tmp/${STUDY}_demux.qza" \
    --p-trim-left-f 0 \
    --p-trim-left-r 0 \
    --p-trunc-len-f 230 \
    --p-trunc-len-r 230 \
    --p-n-threads "$THREADS" \
    --o-table        "tmp/${STUDY}_table.qza" \
    --o-representative-sequences "tmp/${STUDY}_rep_seqs.qza" \
    --o-denoising-stats "tmp/${STUDY}_stats.qza"

# ── QIIME2：分类注释 ──────────────────────────────────────────────
log "3/5 物种分类注释（需要预训练的 SILVA 分类器）..."
# 分类器下载：https://docs.qiime2.org/2024.2/data-resources/
CLASSIFIER="databases/silva-138-99-seqs-515-806-classifier.qza"
if [ ! -f "$CLASSIFIER" ]; then
    log "⚠ 未找到分类器：$CLASSIFIER"
    log "  请下载：https://data.qiime2.org/2024.2/common/silva-138-99-seqs-515-806.qza"
    exit 1
fi

qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER" \
    --i-reads "tmp/${STUDY}_rep_seqs.qza" \
    --p-n-jobs "$THREADS" \
    --o-classification "tmp/${STUDY}_taxonomy.qza"

# 导出 taxonomy.tsv
qiime tools export \
    --input-path "tmp/${STUDY}_taxonomy.qza" \
    --output-path "tmp/${STUDY}_taxonomy_export"

cp "tmp/${STUDY}_taxonomy_export/taxonomy.tsv" \
   "${META_OUTDIR}/${STUDY}.taxonomy.tsv"
log "  ✓ taxonomy.tsv → ${META_OUTDIR}/${STUDY}.taxonomy.tsv"

# ── PICRUSt2：功能预测 ────────────────────────────────────────────
log "4/5 PICRUSt2 功能预测..."

# 导出 ASV 特征表（biom 格式）
qiime tools export \
    --input-path "tmp/${STUDY}_table.qza" \
    --output-path "tmp/${STUDY}_table_export"

# 导出代表序列
qiime tools export \
    --input-path "tmp/${STUDY}_rep_seqs.qza" \
    --output-path "tmp/${STUDY}_seqs_export"

picrust2_pipeline.py \
    -s "tmp/${STUDY}_seqs_export/dna-sequences.fasta" \
    -i "tmp/${STUDY}_table_export/feature-table.biom" \
    -o "tmp/${STUDY}_picrust2_out" \
    -p "$THREADS" \
    --in_traits EC

log "5/5 整理输出文件..."

# 将 stratified EC 贡献表（流程需要的文件）复制到目标目录
cp "tmp/${STUDY}_picrust2_out/EC_metagenome_out/pred_metagenome_strat.tsv.gz" \
   "${META_OUTDIR}/${STUDY}.tsv.gz"
gunzip -f "${META_OUTDIR}/${STUDY}.tsv.gz"
log "  ✓ EC贡献表 → ${META_OUTDIR}/${STUDY}.tsv"

log "═══ $STUDY 处理完成 ═══"
log "输出文件："
log "  ${META_OUTDIR}/${STUDY}.tsv           ← 流程步骤4需要"
log "  ${META_OUTDIR}/${STUDY}.taxonomy.tsv  ← LOGO分析需要"
