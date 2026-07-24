#!/bin/bash
# 通用 SRA/DDBJ 数据下载脚本（SRR/DRR/ERR 均支持）
# 断点续传：已存在 .fastq.gz 的 accession 自动跳过
#
# 下载策略：
#   全量下载（不指定列表）→ 优先用 iSeq（项目级，网络兼容性最好）
#   指定 accession list  → 逐条 prefetch + fasterq-dump
#   指定 SraRunTable.csv → 提取 Run 列后逐条 prefetch + fasterq-dump
#
# 用法：
#   # 全量下载（推荐，iSeq 自动处理）
#   bash scripts/16_download_sra.sh PRJNA542018
#   bash scripts/16_download_sra.sh PRJNA1049117
#
#   # 下载指定子集（accession list，每行一个 Run ID）
#   bash scripts/16_download_sra.sh PRJDB12280 \
#       --accession-list data/OLP/raw/meta/PRJDB12280/accession_list.txt
#
#   # 从 SraRunTable.csv 提取全部 Run 并下载
#   bash scripts/16_download_sra.sh PRJNA1049117 \
#       --sra-table /path/to/SraRunTable.csv
#
#   # 后台运行
#   nohup bash scripts/16_download_sra.sh PRJNA542018 \
#       > logs/download_PRJNA542018.log 2>&1 &

set -euo pipefail

# ── 参数解析 ──────────────────────────────────────────────────────
usage() {
    echo "用法: bash $0 <PROJECT_ID> [--accession-list <file>] [--sra-table <csv>] [--threads <N>]"
    exit 1
}

PROJ="${1:?$(usage)}"
ACC_LIST_ARG=""
SRA_TABLE_ARG=""
THREADS=8

shift
while [[ $# -gt 0 ]]; do
    case "$1" in
        --accession-list) ACC_LIST_ARG="$2"; shift 2 ;;
        --sra-table)      SRA_TABLE_ARG="$2"; shift 2 ;;
        --threads)        THREADS="$2";      shift 2 ;;
        *) echo "[ERROR] 未知参数: $1"; usage ;;
    esac
done

OUTDIR="data/OLP/raw/meta/${PROJ}"
ACC_LIST="${OUTDIR}/accession_list.txt"

mkdir -p "${OUTDIR}" logs

echo "[$(date '+%H:%M:%S')] ── ${PROJ} 下载开始 ──"
echo "  输出目录：${OUTDIR}  线程：${THREADS}"

# ── 情形1：全量下载且 iSeq 可用 → 直接用 iSeq 项目级下载 ──────────
if [[ -z "${ACC_LIST_ARG}" && -z "${SRA_TABLE_ARG}" ]] && command -v iseq &>/dev/null; then
    echo "[$(date '+%H:%M:%S')] 使用 iSeq 全量下载 ${PROJ}..."
    iseq -i "${PROJ}" -g -p 5 -t "${THREADS}" -o "${OUTDIR}"
    echo ""
    echo "[$(date '+%H:%M:%S')] ── ${PROJ} 下载结束 ──"
    exit 0
fi

# ── 情形2：按列表下载（指定子集）────────────────────────────────────
if [[ -n "${SRA_TABLE_ARG}" ]]; then
    # 从 SraRunTable.csv 提取 Run 列
    if [[ ! -f "${SRA_TABLE_ARG}" ]]; then
        echo "[ERROR] 找不到 SraRunTable.csv：${SRA_TABLE_ARG}"; exit 1
    fi
    python3 - "${SRA_TABLE_ARG}" "${ACC_LIST}" <<'PYEOF'
import csv, sys
src, dst = sys.argv[1], sys.argv[2]
with open(src, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    run_col = next((c for c in reader.fieldnames if c.lower() in ("run", "run_accession")), None)
    if not run_col:
        sys.exit(f"[ERROR] 找不到 Run 列，列名：{reader.fieldnames}")
    runs = [row[run_col].strip() for row in reader if row[run_col].strip()]
with open(dst, "w") as f:
    f.write("\n".join(runs) + "\n")
print(f"  从 CSV 提取到 {len(runs)} 个 run")
PYEOF
    echo "  accession list 已保存：${ACC_LIST}"

elif [[ -n "${ACC_LIST_ARG}" ]]; then
    if [[ ! -f "${ACC_LIST_ARG}" ]]; then
        echo "[ERROR] 找不到 accession list：${ACC_LIST_ARG}"; exit 1
    fi
    cp "${ACC_LIST_ARG}" "${ACC_LIST}"
    echo "  使用指定 accession list（$(wc -l < ${ACC_LIST}) 个）"

elif [[ -f "${ACC_LIST}" && $(grep -c . "${ACC_LIST}" || true) -gt 0 ]]; then
    echo "  使用缓存 accession list（$(wc -l < ${ACC_LIST}) 个）"

else
    echo "[ERROR] 未提供 accession list 或 SraRunTable.csv"
    echo "  请从 NCBI SRA Run Selector 下载 SraRunTable.csv："
    echo "  https://www.ncbi.nlm.nih.gov/Traces/study/?acc=${PROJ}"
    echo "  然后运行：bash $0 ${PROJ} --sra-table /path/to/SraRunTable.csv"
    exit 1
fi

COUNT=$(grep -c . "${ACC_LIST}" || true)
echo "[$(date '+%H:%M:%S')] 使用 iSeq 按列表下载（${COUNT} 个）..."

# iSeq -i 支持传入文件（每行一个 accession），-g 直接输出 .fastq.gz，断点续传由 iSeq 自身处理
iseq -i "${ACC_LIST}" -g -p 5 -t "${THREADS}" -o "${OUTDIR}"

echo ""
echo "[$(date '+%H:%M:%S')] ── ${PROJ} 下载结束 ──"
echo "  [提示] 如有失败项，重新运行脚本即可断点续传（iSeq 自动跳过已存在文件）"
