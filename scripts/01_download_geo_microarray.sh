#!/usr/bin/env bash
# =============================================================================
# 步骤1：下载 GEO 微阵列 series_matrix.txt 文件
# 直接从 GEO FTP 下载，无需 iSeq（iSeq 下载的是原始 CEL，不是 series_matrix）
#
# 用法：bash scripts/01_download_geo_microarray.sh
#        bash scripts/01_download_geo_microarray.sh --threads 4
# =============================================================================

set -euo pipefail

THREADS=${2:-4}
OUTDIR="data/COPD/raw/host"
GEO_FTP="https://ftp.ncbi.nlm.nih.gov/geo/series"
LOG="logs/download_geo_microarray.log"

mkdir -p "$OUTDIR" logs

# ── 19 个 GEO 数据集 ──────────────────────────────────────────────
DATASETS=(
    GSE103174
    GSE106986
    GSE112260
    GSE11784
    GSE119040
    GSE11906
    GSE12472
    GSE13896
    GSE16972
    GSE37147
    GSE37768
    GSE38974
    GSE47460
    GSE56341
    GSE57148
    GSE73395
    GSE76925
    GSE8581
    GSE86064
)

log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

# GEO FTP 路径规则：GSE11906 → GSE11nnn/GSE11906
geo_prefix() {
    local gse=$1
    local num=${gse#GSE}
    local prefix_len=$(( ${#num} - 3 ))
    if [ "$prefix_len" -le 0 ]; then
        echo "${gse%???}nnn"
    else
        echo "GSE${num:0:$prefix_len}nnn"
    fi
}

download_series_matrix() {
    local gse=$1
    local prefix
    prefix=$(geo_prefix "$gse")
    local url="${GEO_FTP}/${prefix}/${gse}/matrix/${gse}_series_matrix.txt.gz"
    local outfile="${OUTDIR}/${gse}_series_matrix.txt"

    if [ -f "$outfile" ]; then
        log "跳过（已存在）：$gse"
        return 0
    fi

    log "下载 $gse series_matrix ..."
    if wget -q --tries=3 --timeout=60 -O "${outfile}.gz" "$url"; then
        gunzip -f "${outfile}.gz"
        log "✓ $gse"
    else
        # 部分数据集有多个平台，尝试带平台后缀的文件名
        log "⚠ 标准路径失败，尝试列目录: $gse"
        local dir_url="${GEO_FTP}/${prefix}/${gse}/matrix/"
        # 获取目录下所有 series_matrix 文件
        wget -q -O /tmp/geo_matrix_list.html "$dir_url" || true
        local files
        files=$(grep -oP 'href="\K[^"]*series_matrix[^"]*' /tmp/geo_matrix_list.html 2>/dev/null || true)
        if [ -z "$files" ]; then
            log "✗ 无法下载 $gse，请手动检查：$dir_url"
            echo "$gse" >> logs/download_failed.log
            return 1
        fi
        # 下载第一个匹配文件（通常只有一个）
        local first_file
        first_file=$(echo "$files" | head -1)
        wget -q --tries=3 -O "${outfile}.gz" "${dir_url}${first_file}"
        gunzip -f "${outfile}.gz"
        log "✓ $gse（备用路径）"
    fi
}

log "开始下载 ${#DATASETS[@]} 个 GEO 微阵列数据集..."
log "输出目录：$OUTDIR"

# 并行下载（使用 xargs 控制并发数）
printf '%s\n' "${DATASETS[@]}" | \
    xargs -P "$THREADS" -I{} bash -c 'source scripts/01_download_geo_microarray.sh; download_series_matrix "$@"' _ {}

log "✅ series_matrix 下载完成"
log "失败列表（如有）：logs/download_failed.log"
