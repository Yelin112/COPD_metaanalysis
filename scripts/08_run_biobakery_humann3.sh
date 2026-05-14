#!/usr/bin/env bash
# =============================================================================
# 步骤8（可选）：Biobakery / HUMAnN3 宏基因组功能预测
# 适用于鸟枪法宏基因组（shotgun metagenomics）数据
# 16S 扩增子数据请使用 scripts/04_run_qiime2_picrust2.sh
#
# HUMAnN3 与 PICRUSt2 的区别：
#   PICRUSt2  基于 16S 序列预测功能，间接推断
#   HUMAnN3   直接对鸟枪法reads比对，获得精确 EC 丰度，无需中间预测步骤
#
# HUMAnN3 输出（本流程需要的文件）：
#   *_genefamilies.tsv   UniRef90/EC 丰度，含属级分层（适合 LOGO 分析）
#   *_pathabundance.tsv  MetaCyc 通路丰度（可选）
#
# 安装：
#   conda create -n biobakery -c biobakery biobakery3
#   conda activate biobakery
#   humann_databases --download chocophlan full databases/chocophlan
#   humann_databases --download uniref uniref90_diamond databases/uniref
#
# 用法：bash scripts/08_run_biobakery_humann3.sh SRP066375
# =============================================================================

set -euo pipefail

STUDY=${1:?"用法: bash $0 <STUDY_ID>，如 PRJNA542018"}
FASTQ_DIR="data/OLP/raw/meta/fastq/${STUDY}"
META_OUTDIR="data/OLP/raw/meta"
THREADS=${2:-16}
LOG="logs/humann3_${STUDY}.log"
DB_CHOCO="databases/chocophlan"
DB_UNIREF="databases/uniref"

mkdir -p "$META_OUTDIR" "tmp/humann3_${STUDY}" logs

log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

log "═══ HUMAnN3 处理：$STUDY ═══"

# ── 检查数据库 ─────────────────────────────────────────────────
if [ ! -d "$DB_CHOCO" ] || [ ! -d "$DB_UNIREF" ]; then
    log "⚠ 未找到数据库，请先下载："
    log "  humann_databases --download chocophlan full $DB_CHOCO"
    log "  humann_databases --download uniref uniref90_diamond $DB_UNIREF"
    exit 1
fi

# ── 合并配对 reads（HUMAnN3 可接受合并后的 FASTQ）─────────────
log "1/4 合并配对端 reads..."
merged="${META_OUTDIR}/tmp_${STUDY}_merged.fastq.gz"
cat "${FASTQ_DIR}"/*_R1*.fastq.gz "${FASTQ_DIR}"/*_R2*.fastq.gz > "$merged" 2>/dev/null || \
    cat "${FASTQ_DIR}"/*.fastq.gz > "$merged"

# ── KneadData 质控（可选但推荐）────────────────────────────────
# log "1.5/4 KneadData 去宿主..."
# kneaddata --input "$merged" \
#     --reference-db databases/kneaddata_human \
#     --output "tmp/kneaddata_${STUDY}" \
#     --threads "$THREADS"
# merged="tmp/kneaddata_${STUDY}/$(basename $merged .fastq.gz)_kneaddata_paired_1.fastq"

# ── HUMAnN3 主流程 ─────────────────────────────────────────────
log "2/4 运行 HUMAnN3（此步骤较慢，约 2-8 小时）..."
humann \
    --input "$merged" \
    --output "tmp/humann3_${STUDY}" \
    --nucleotide-database "$DB_CHOCO" \
    --protein-database "$DB_UNIREF" \
    --threads "$THREADS" \
    --output-basename "${STUDY}"

log "3/4 重新分组为 EC 编号..."
humann_regroup_table \
    --input "tmp/humann3_${STUDY}/${STUDY}_genefamilies.tsv" \
    --groups uniref90_level4ec \
    --output "tmp/humann3_${STUDY}/${STUDY}_ec.tsv"

log "4/4 整理输出文件..."

# 保留 EC 分层表（行格式：EC:X.X.X.X|g__Genus.s__Species 或 EC:X.X.X.X|unclassified）
# 此文件供 modules/meta/metagenomics/01b_humann3_to_ec.py 处理
cp "tmp/humann3_${STUDY}/${STUDY}_ec.tsv" \
   "${META_OUTDIR}/${STUDY}.humann3_ec.tsv"
log "  ✓ EC 分层表 → ${META_OUTDIR}/${STUDY}.humann3_ec.tsv"

# 清理合并的临时文件
rm -f "$merged"

log "═══ $STUDY 处理完成 ═══"
log "输出文件："
log "  ${META_OUTDIR}/${STUDY}.humann3_ec.tsv  ← 流程步骤需要"
log ""
log "⚠ 完成后请更新 config/OLP/study_config.yaml："
log "   将 upstream_tool: picrust2 改为 upstream_tool: humann3"
log "   将 feature_file 指向 ${STUDY}.humann3_ec.tsv"
log "   HUMAnN3 数据无需 taxonomy_file（分层信息已内嵌）"
