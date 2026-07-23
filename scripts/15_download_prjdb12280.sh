#!/bin/bash
# PRJDB12280 Presort fraction 下载脚本
# 断点续传：已存在 .fastq.gz 的 accession 自动跳过
#
# 用法：
#   conda activate iseq
#   cd /home/usr/yel/exp-OLP_meta/pipeline
#   nohup bash scripts/15_download_prjdb12280.sh > logs/download_PRJDB12280.log 2>&1 &

set -euo pipefail

PROJ="PRJDB12280"
OUTDIR="data/OLP/raw/meta/${PROJ}"
ACC_LIST="${OUTDIR}/accession_list.txt"
THREADS=8

mkdir -p "${OUTDIR}" logs

if [[ ! -f "${ACC_LIST}" ]]; then
    echo "[ERROR] 找不到 accession list：${ACC_LIST}"
    echo "  请先运行：python3 scripts/14_filter_prjdb12280.py ..."
    exit 1
fi

TOTAL=$(wc -l < "${ACC_LIST}")
echo "[$(date '+%H:%M:%S')] ── ${PROJ} 下载开始，共 ${TOTAL} 个 accession ──"

DONE=0
SKIP=0
FAIL=0

while read -r ACC; do
    [[ -z "${ACC}" ]] && continue

    # 断点续传：双端或单端任意一个 .fastq.gz 存在则跳过
    if ls "${OUTDIR}/${ACC}"*.fastq.gz 2>/dev/null | grep -q .; then
        echo "[$(date '+%H:%M:%S')] SKIP ${ACC}（已存在）"
        SKIP=$((SKIP + 1))
        continue
    fi

    echo "[$(date '+%H:%M:%S')] 下载 ${ACC} ..."

    # ── 优先尝试 iSeq ──────────────────────────────────────────────
    if command -v iseq &>/dev/null; then
        if iseq download --accession "${ACC}" --output "${OUTDIR}" --threads "${THREADS}" 2>&1; then
            DONE=$((DONE + 1))
            echo "[$(date '+%H:%M:%S')] ✅ ${ACC} 完成"
            continue
        fi
        echo "[$(date '+%H:%M:%S')] iSeq 失败，切换 prefetch..."
    fi

    # ── 备用：prefetch + fasterq-dump ─────────────────────────────
    TMP_SRA="${OUTDIR}/sra_tmp/${ACC}"
    mkdir -p "${OUTDIR}/sra_tmp"

    if prefetch "${ACC}" --output-directory "${OUTDIR}/sra_tmp" --max-size 50G 2>&1 && \
       fasterq-dump "${TMP_SRA}/${ACC}.sra" \
           --split-files --threads "${THREADS}" \
           --outdir "${OUTDIR}" 2>&1; then

        # 压缩
        for fq in "${OUTDIR}/${ACC}"*.fastq; do
            [[ -f "${fq}" ]] && gzip "${fq}"
        done
        rm -rf "${TMP_SRA}"
        DONE=$((DONE + 1))
        echo "[$(date '+%H:%M:%S')] ✅ ${ACC} 完成"
    else
        FAIL=$((FAIL + 1))
        echo "[$(date '+%H:%M:%S')] ❌ ${ACC} 下载失败，跳过"
    fi

done < "${ACC_LIST}"

echo ""
echo "[$(date '+%H:%M:%S')] ── ${PROJ} 下载完成 ──"
echo "  成功：${DONE}  跳过（已存在）：${SKIP}  失败：${FAIL}"
echo "  文件目录：${OUTDIR}"

if [[ ${FAIL} -gt 0 ]]; then
    echo ""
    echo "[提示] 有 ${FAIL} 个 accession 下载失败，可重新运行脚本（断点续传会自动跳过已完成的）"
fi
