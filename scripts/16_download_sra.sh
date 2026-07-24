#!/bin/bash
# 通用 SRA/DDBJ 数据下载脚本（SRR/DRR/ERR 均支持）
# 断点续传：已存在 .fastq.gz 的 accession 自动跳过
#
# 用法：
#   # 下载某项目全部样本
#   bash scripts/16_download_sra.sh PRJNA542018
#   bash scripts/16_download_sra.sh PRJNA1049117
#
#   # 下载指定 accession list（每行一个 Run ID）
#   bash scripts/16_download_sra.sh PRJDB12280 \
#       --accession-list data/OLP/raw/meta/PRJDB12280/accession_list.txt
#
#   # 后台运行
#   nohup bash scripts/16_download_sra.sh PRJNA542018 \
#       > logs/download_PRJNA542018.log 2>&1 &
#
# 依赖：prefetch, fasterq-dump（SRA Toolkit）；esearch/efetch（NCBI E-utils，仅全量下载时需要）

set -euo pipefail

# ── 参数解析 ──────────────────────────────────────────────────────
usage() {
    echo "用法: bash $0 <PROJECT_ID> [--accession-list <file>] [--threads <N>]"
    exit 1
}

PROJ="${1:?$(usage)}"
ACC_LIST_ARG=""
THREADS=8

shift
while [[ $# -gt 0 ]]; do
    case "$1" in
        --accession-list) ACC_LIST_ARG="$2"; shift 2 ;;
        --threads)        THREADS="$2";      shift 2 ;;
        *) echo "[ERROR] 未知参数: $1"; usage ;;
    esac
done

OUTDIR="data/OLP/raw/meta/${PROJ}"
SRA_TMP="${OUTDIR}/.sra_tmp"
ACC_LIST="${OUTDIR}/accession_list.txt"

mkdir -p "${OUTDIR}" "${SRA_TMP}" logs

echo "[$(date '+%H:%M:%S')] ── ${PROJ} 下载开始 ──"
echo "  输出目录：${OUTDIR}  线程：${THREADS}"

# ── 获取 accession list ───────────────────────────────────────────
if [[ -n "${ACC_LIST_ARG}" ]]; then
    # 用户指定的列表
    if [[ ! -f "${ACC_LIST_ARG}" ]]; then
        echo "[ERROR] 找不到 accession list：${ACC_LIST_ARG}"; exit 1
    fi
    cp "${ACC_LIST_ARG}" "${ACC_LIST}"
    echo "  使用指定 accession list：${ACC_LIST_ARG}（$(wc -l < ${ACC_LIST}) 个）"

elif [[ -f "${ACC_LIST}" ]]; then
    # 已有缓存列表（上次运行生成）
    echo "  使用缓存 accession list（$(wc -l < ${ACC_LIST}) 个）"

else
    # 自动查询项目全部 Run：优先 esearch，备用 curl SRA API
    echo "[$(date '+%H:%M:%S')] 查询 ${PROJ} 的 Run 列表..."

    if command -v esearch &>/dev/null; then
        esearch -db sra -query "${PROJ}[BioProject]" \
            | efetch -format runinfo \
            | tail -n +2 \
            | cut -d',' -f1 \
            | grep -v '^$' \
            > "${ACC_LIST}"
    else
        echo "  esearch 不可用，改用 NCBI SRA API（curl）..."
        # NCBI SRA API：按 BioProject 查询，最多返回 10000 条
        curl -fsSLg \
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=${PROJ}[BioProject]&retmax=10000&retmode=json" \
            | python3 -c "
import sys, json
data = json.load(sys.stdin)
ids = data['esearchresult']['idlist']
print(f'  找到 {len(ids)} 个 SRA ID', file=sys.stderr)
# 分批获取 Run accession
import urllib.request, time
batch = 200
for i in range(0, len(ids), batch):
    chunk = ids[i:i+batch]
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=' + ','.join(chunk) + '&rettype=runinfo&retmode=text'
    resp = urllib.request.urlopen(url).read().decode()
    for line in resp.splitlines()[1:]:
        run = line.split(',')[0].strip()
        if run and run != 'Run':
            print(run)
    time.sleep(0.4)
" > "${ACC_LIST}"
    fi

    COUNT=$(grep -c . "${ACC_LIST}" || true)
    if [[ "${COUNT}" -eq 0 ]]; then
        echo "[ERROR] 未查询到任何 Run，请检查项目号或网络"; exit 1
    fi
    echo "  共找到 ${COUNT} 个 run，已保存至 ${ACC_LIST}"
fi

TOTAL=$(grep -c . "${ACC_LIST}" || true)
DONE=0; SKIP=0; FAIL=0

# ── 主下载循环 ────────────────────────────────────────────────────
while read -r ACC; do
    [[ -z "${ACC}" ]] && continue

    # 断点续传：任一 .fastq.gz 已存在则跳过
    if ls "${OUTDIR}/${ACC}"*.fastq.gz 2>/dev/null | grep -q .; then
        echo "[$(date '+%H:%M:%S')] SKIP ${ACC}（已存在）"
        SKIP=$((SKIP + 1))
        continue
    fi

    echo "[$(date '+%H:%M:%S')] ▶ ${ACC} prefetch... [${DONE}+${SKIP}+${FAIL}/${TOTAL}]"
    if ! prefetch "${ACC}" \
            --output-directory "${SRA_TMP}" \
            --max-size 50G \
            --progress 2>&1; then
        echo "[$(date '+%H:%M:%S')] ❌ ${ACC} prefetch 失败"
        FAIL=$((FAIL + 1)); continue
    fi

    # prefetch 输出位置兼容：子目录或直接文件
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
        rm -rf "${SRA_TMP:?}/${ACC}" "${SRA_TMP}/${ACC}.sra"
        continue
    fi

    # 压缩并清理 SRA 缓存
    for fq in "${OUTDIR}/${ACC}"*.fastq; do
        [[ -f "${fq}" ]] && gzip "${fq}"
    done
    rm -rf "${SRA_TMP:?}/${ACC}" "${SRA_TMP}/${ACC}.sra"

    DONE=$((DONE + 1))
    echo "[$(date '+%H:%M:%S')] ✅ ${ACC} 完成 （${DONE}/${TOTAL}）"

done < "${ACC_LIST}"

rmdir "${SRA_TMP}" 2>/dev/null || true

echo ""
echo "[$(date '+%H:%M:%S')] ── ${PROJ} 下载结束 ──"
echo "  ✅ 成功：${DONE}  ⏭ 跳过：${SKIP}  ❌ 失败：${FAIL}"
[[ ${FAIL} -gt 0 ]] && echo "  [提示] 有失败项，重新运行脚本即可断点续传"
