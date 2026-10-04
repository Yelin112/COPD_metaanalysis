#!/usr/bin/env bash
# =============================================================================
# 步骤7（OLP专用）：下载口腔扁平苔藓 16S 微生物组原始 FASTQ 数据
# 使用 iSeq 从 SRA/ENA 下载
#
# 数据集：
#   PRJNA542018  唾液 16S V4（MiSeq），20 OLP + 10 RAU + 10 HC
#                Yu et al. BMC Microbiology 2020 (PMID 32245419)
#                注：原始 SRP 编号请在 https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=542018 确认
#
#   PRJNA1097766 口腔微生物组 + 代谢组（16S rRNA）
#                OLP 患者 vs 健康对照（2024年发表）
#
# 依赖：conda create -n iseq -c conda-forge -c bioconda iseq
#
# 用法：bash scripts/07_download_olp_microbiome_sra.sh
# =============================================================================

set -euo pipefail

OUTDIR="data/OLP/raw/meta/fastq"
LOG="logs/download_olp_microbiome.log"
THREADS=8
PARALLEL=5

mkdir -p "$OUTDIR" logs

log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

# ── OLP 16S 数据集（BioProject 编号）────────────────────────────
STUDIES=(
    "PRJNA542018"    # 唾液微生物组，Yu et al. 2020
    "PRJNA1097766"   # 口腔微生物组+代谢组
)

# 检查 iSeq 是否可用
if ! command -v iseq &> /dev/null; then
    log "❌ 未找到 iseq 命令，请先安装："
    log "   conda create -n iseq -c conda-forge -c bioconda iseq"
    log "   conda activate iseq"
    exit 1
fi

log "开始下载 OLP 微生物组数据集（共 ${#STUDIES[@]} 个）"
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

    # iSeq 支持 BioProject（PRJNA）和 SRP 编号
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
log "✅ OLP 微生物组原始数据下载完成"
log ""
log "══════════════════════════════════════════════════════════"
log "  后续步骤：QIIME2 + PICRUSt2 处理（16S扩增子）"
log "══════════════════════════════════════════════════════════"
log ""
log "  下载完成后，需要运行 16S 上游分析流程："
log ""
log "  1. 对每个数据集运行 QIIME2 DADA2 + PICRUSt2："
log "     bash scripts/04_run_qiime2_picrust2.sh PRJNA542018"
log "     bash scripts/04_run_qiime2_picrust2.sh PRJNA1097766"
log ""
log "  2. 输出文件放置到："
log "     data/OLP/raw/meta/PRJNA542018.tsv          (pred_metagenome_strat.tsv)"
log "     data/OLP/raw/meta/PRJNA542018.taxonomy.tsv  (taxonomy.tsv)"
log "     data/OLP/raw/meta/PRJNA1097766.tsv"
log "     data/OLP/raw/meta/PRJNA1097766.taxonomy.tsv"
log ""
log "  ⚠ 注意：PRJNA542018 数据集包含 OLP、RAU（复发性口腔溃疡）、HC 三组。"
log "    请在 metadata 文件中将 OLP 标注为 'OLP'，健康对照标注为 'Control'，"
log "    RAU 组可标注为 'RAU'（流程会自动跳过非 OLP/Control 样本）。"
log ""
log "  若数据为鸟枪法宏基因组，请改用 HUMAnN3 流程："
log "     bash scripts/08_run_biobakery_humann3.sh PRJNA542018"
log "     并在 config/OLP/study_config.yaml 中设置 upstream_tool: humann3"
