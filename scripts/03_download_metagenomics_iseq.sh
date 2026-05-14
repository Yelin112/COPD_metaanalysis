#!/usr/bin/env bash
# =============================================================================
# 步骤3：使用 iSeq 下载宏基因组原始 FASTQ 数据
# 下载后还需要运行 QIIME2 + PICRUSt2 才能得到流程需要的 EC 贡献表
#
# 安装 iSeq：conda create -n iseq -c conda-forge -c bioconda iseq
#
# 用法：bash scripts/03_download_metagenomics_iseq.sh
# =============================================================================

set -euo pipefail

OUTDIR="data/COPD/raw/meta/fastq"
LOG="logs/download_metagenomics.log"
THREADS=8       # iSeq fasterq-dump 线程数
PARALLEL=5      # 并发下载连接数

mkdir -p "$OUTDIR" logs

log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

# ── 两个宏基因组数据集 ────────────────────────────────────────────
# SRP066375: Human Lung Microbiome (Pragman et al.)
# SRP136124: COPD Airway Microbiome (CAMRA study)
STUDIES=(
    "SRP066375"
    "SRP136124"
)

# 检查 iSeq 是否可用
if ! command -v iseq &> /dev/null; then
    log "❌ 未找到 iseq 命令，请先安装："
    log "   conda create -n iseq -c conda-forge -c bioconda iseq"
    log "   conda activate iseq"
    exit 1
fi

log "开始下载宏基因组数据集（共 ${#STUDIES[@]} 个）"
log "输出目录：$OUTDIR"

for study in "${STUDIES[@]}"; do
    study_dir="${OUTDIR}/${study}"
    mkdir -p "$study_dir"

    if [ -f "${study_dir}/.download_complete" ]; then
        log "跳过（已完成）：$study"
        continue
    fi

    log "━━━ 下载 $study ━━━"
    cd "$study_dir"

    # 使用 iSeq 下载 gzip 格式 FASTQ，多线程并行
    iseq -i "$study" \
         -g \
         -t "$THREADS" \
         -p "$PARALLEL"

    if [ $? -eq 0 ]; then
        touch .download_complete
        log "✓ $study 下载完成"
        fastq_count=$(find . -name "*.fastq.gz" | wc -l)
        log "  FASTQ 文件数：$fastq_count"
    else
        log "✗ $study 下载可能有部分失败，请检查 fail.log"
    fi

    cd - > /dev/null
done

log ""
log "✅ 原始数据下载完成"
log ""
log "══════════════════════════════════════════════"
log "  后续步骤：QIIME2 + PICRUSt2 处理"
log "══════════════════════════════════════════════"
log "  下载完成后，还需要运行以下流程才能得到"
log "  我们的流程所需要的 EC 贡献表："
log ""
log "  1. QIIME2 质控 + OTU/ASV 聚类"
log "     → 输出：taxonomy.tsv（分类注释文件）"
log ""
log "  2. PICRUSt2 功能预测"
log "     → 输出：EC_metagenome_out/pred_metagenome_unstrat.tsv"
log "             EC_metagenome_out/pred_metagenome_strat.tsv  ← 流程需要这个"
log ""
log "  3. 将输出文件重命名并放置到："
log "     data/COPD/raw/meta/SRP066375.tsv         (pred_metagenome_strat.tsv)"
log "     data/COPD/raw/meta/SRP066375.taxonomy.tsv (taxonomy.tsv)"
log "     data/COPD/raw/meta/SRP136124.tsv"
log "     data/COPD/raw/meta/SRP136124.taxonomy.tsv"
log ""
log "  参考命令（QIIME2 + PICRUSt2 安装后）："
log "  bash scripts/04_run_qiime2_picrust2.sh"
