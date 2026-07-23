#!/bin/bash
# PRJNA542018 全量下载脚本（16S V4，单端 merged PE，SRR 编号）
# 断点续传：已存在 .fastq.gz 的 accession 自动跳过
#
# 依赖：prefetch, fasterq-dump（SRA Toolkit，conda activate iseq 环境中已有）
#
# 用法：
#   conda activate iseq
#   cd /home/usr/yel/exp-OLP_meta/pipeline
#   nohup bash scripts/16_download_prjna542018.sh > logs/download_PRJNA542018.log 2>&1 &
#
# 查看进度：
#   tail -f logs/download_PRJNA542018.log

set -euo pipefail

PROJ="PRJNA542018"
OUTDIR="data/OLP/raw/meta/${PROJ}"
SRA_TMP="${OUTDIR}/.sra_tmp"
THREADS=8

mkdir -p "${OUTDIR}" "${SRA_TMP}" logs

echo "[$(date '+%H:%M:%S')] ── ${PROJ} 下载开始 ──"
echo "  输出目录：${OUTDIR}"

# ── 第一步：获取该项目所有 SRR accession ──────────────────────────
ACC_LIST="${OUTDIR}/accession_list.txt"
if [[ ! -f "${ACC_LIST}" ]]; then
    echo "[$(date '+%H:%M:%S')] 查询 ${PROJ} 的 Run 列表..."
    esearch -db sra -query "${PROJ}[BioProject]" \
        | efetch -format runinfo \
        | tail -n +2 \
        | cut -d',' -f1 \
        | grep -v '^$' \
        > "${ACC_LIST}"
    echo "  共找到 $(wc -l < ${ACC_LIST}) 个 run"
else
    echo "[$(date '+%H:%M:%S')] 使用已有 accession list（$(wc -l < ${ACC_LIST}) 个）"
fi

TOTAL=$(grep -c . "${ACC_LIST}" || true)
DONE=0; SKIP=0; FAIL=0

while read -r ACC; do
    [[ -z "${ACC}" ]] && continue

    # ── 断点续传 ────────────────────────────────────────────────────
    if ls "${OUTDIR}/${ACC}"*.fastq.gz 2>/dev/null | grep -q .; then
        echo "[$(date '+%H:%M:%S')] SKIP ${ACC}（已存在）"
        SKIP=$((SKIP + 1))
        continue
    fi

    echo "[$(date '+%H:%M:%S')] ▶ ${ACC} prefetch..."
    if ! prefetch "${ACC}" \
            --output-directory "${SRA_TMP}" \
            --max-size 50G \
            --progress 2>&1; then
        echo "[$(date '+%H:%M:%S')] ❌ ${ACC} prefetch 失败"
        FAIL=$((FAIL + 1))
        continue
    fi

    SRA_FILE="${SRA_TMP}/${ACC}/${ACC}.sra"
    [[ ! -f "${SRA_FILE}" ]] && SRA_FILE="${SRA_TMP}/${ACC}.sra"

    echo "[$(date '+%H:%M:%S')] ▶ ${ACC} fasterq-dump..."
    if ! fasterq-dump "${SRA_FILE}" \
            --split-files \
            --threads "${THREADS}" \
            --outdir "${OUTDIR}" \
            --progress 2>&1; then
        echo "[$(date '+%H:%M:%S')] ❌ ${ACC} fasterq-dump 失败"
        FAIL=$((FAIL + 1))
        rm -rf "${SRA_TMP}/${ACC}" "${SRA_TMP}/${ACC}.sra"
        continue
    fi

    for fq in "${OUTDIR}/${ACC}"*.fastq; do
        [[ -f "${fq}" ]] && gzip "${fq}"
    done
    rm -rf "${SRA_TMP}/${ACC}" "${SRA_TMP}/${ACC}.sra"

    DONE=$((DONE + 1))
    echo "[$(date '+%H:%M:%S')] ✅ ${ACC} 完成 （${DONE}/${TOTAL}）"

done < "${ACC_LIST}"

rmdir "${SRA_TMP}" 2>/dev/null || true

echo ""
echo "[$(date '+%H:%M:%S')] ── ${PROJ} 下载结束 ──"
echo "  ✅ 成功：${DONE}  ⏭ 跳过：${SKIP}  ❌ 失败：${FAIL}"
echo "  下一步：bash scripts/11_qiime2_PRJNA542018.sh"

if [[ ${FAIL} -gt 0 ]]; then
    echo "[提示] 重新运行脚本即可断点续传"
fi
