#!/bin/bash
# PRJDB12280 Presort fraction 下载脚本
# 使用 SRA Toolkit（prefetch + fasterq-dump），原生支持 DRR 编号
# 断点续传：已存在 .fastq.gz 的 accession 自动跳过
#
# 依赖：prefetch, fasterq-dump（SRA Toolkit）
#
# 用法：
#   conda activate iseq   # 或任何含 SRA Toolkit 的环境
#   cd /home/usr/yel/exp-OLP_meta/pipeline
#   nohup bash scripts/15_download_prjdb12280.sh > logs/download_PRJDB12280.log 2>&1 &
#
# 查看进度：
#   tail -f logs/download_PRJDB12280.log

set -euo pipefail

PROJ="PRJDB12280"
OUTDIR="data/OLP/raw/meta/${PROJ}"
ACC_LIST="${OUTDIR}/accession_list.txt"
SRA_TMP="${OUTDIR}/.sra_tmp"
THREADS=8

mkdir -p "${OUTDIR}" "${SRA_TMP}" logs

if [[ ! -f "${ACC_LIST}" ]]; then
    echo "[ERROR] 找不到 accession list：${ACC_LIST}"
    echo "  请先运行：python3 scripts/14_filter_prjdb12280.py ..."
    exit 1
fi

TOTAL=$(grep -c . "${ACC_LIST}" || true)
echo "[$(date '+%H:%M:%S')] ── ${PROJ} 下载开始，共 ${TOTAL} 个 accession ──"
echo "  输出目录：${OUTDIR}"

DONE=0; SKIP=0; FAIL=0

while read -r ACC; do
    [[ -z "${ACC}" ]] && continue

    # ── 断点续传：任一 .fastq.gz 已存在则跳过 ──────────────────────
    if ls "${OUTDIR}/${ACC}"*.fastq.gz 2>/dev/null | grep -q .; then
        echo "[$(date '+%H:%M:%S')] SKIP ${ACC}（已存在）"
        SKIP=$((SKIP + 1))
        continue
    fi

    echo "[$(date '+%H:%M:%S')] ▶ ${ACC} prefetch..."

    # ── prefetch：下载 .sra 文件 ───────────────────────────────────
    if ! prefetch "${ACC}" \
            --output-directory "${SRA_TMP}" \
            --max-size 50G \
            --progress 2>&1; then
        echo "[$(date '+%H:%M:%S')] ❌ ${ACC} prefetch 失败，跳过"
        FAIL=$((FAIL + 1))
        continue
    fi

    # ── fasterq-dump：.sra → .fastq，自动分 _1/_2 ─────────────────
    SRA_FILE="${SRA_TMP}/${ACC}/${ACC}.sra"
    if [[ ! -f "${SRA_FILE}" ]]; then
        # prefetch 有时直接写成 ACC.sra 而不建子目录
        SRA_FILE="${SRA_TMP}/${ACC}.sra"
    fi

    echo "[$(date '+%H:%M:%S')] ▶ ${ACC} fasterq-dump..."
    if ! fasterq-dump "${SRA_FILE}" \
            --split-files \
            --threads "${THREADS}" \
            --outdir "${OUTDIR}" \
            --progress 2>&1; then
        echo "[$(date '+%H:%M:%S')] ❌ ${ACC} fasterq-dump 失败，跳过"
        FAIL=$((FAIL + 1))
        rm -rf "${SRA_TMP}/${ACC}" "${SRA_TMP}/${ACC}.sra"
        continue
    fi

    # ── 压缩 + 清理 SRA 缓存 ──────────────────────────────────────
    for fq in "${OUTDIR}/${ACC}"*.fastq; do
        [[ -f "${fq}" ]] && gzip "${fq}"
    done
    rm -rf "${SRA_TMP}/${ACC}" "${SRA_TMP}/${ACC}.sra"

    DONE=$((DONE + 1))
    echo "[$(date '+%H:%M:%S')] ✅ ${ACC} 完成"

done < "${ACC_LIST}"

# 清理空的临时目录
rmdir "${SRA_TMP}" 2>/dev/null || true

echo ""
echo "[$(date '+%H:%M:%S')] ── ${PROJ} 下载结束 ──"
echo "  ✅ 成功：${DONE}  ⏭ 跳过：${SKIP}  ❌ 失败：${FAIL}"

if [[ ${FAIL} -gt 0 ]]; then
    echo ""
    echo "[提示] ${FAIL} 个失败，重新运行脚本即可断点续传（已完成的自动跳过）"
fi
