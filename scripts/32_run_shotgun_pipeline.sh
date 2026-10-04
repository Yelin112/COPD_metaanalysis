#!/bin/bash
# 流程 A：fastp QC → kneaddata 去宿主 → humann3 功能定量 → EC regroup（per-sample）
#
# 用法：
#   bash scripts/32_run_shotgun_pipeline.sh <sample_id> [线程数] [选项]
#   bash scripts/32_run_shotgun_pipeline.sh CJ_PBS 16 \
#       --rawdir /home/exs/data3/seq_backup/upload/Metagenome/Rawdata \
#       --name1 .raw.1.fq.gz --name2 .raw.2.fq.gz \
#       --outbase /home/usr/yel/exp-OLK/mgx
#
# 默认行为（不传选项）= OLP PRJNA1201607 的旧用法（_1/_2.fastq.gz 命名，输出 data/OLP/shotgun/）
# 输出: <outbase>/<sample>/{fastp,kneaddata,humann}，日志 <outbase>/logs/32_<sample>.log
#
# 断点续跑（两级）：
#   - humann/.done_catPE 存在           → 整个样本已按"R1+R2 拼接"新流程完成，直接退出
#   - kneaddata paired_1/2 存在且非空   → 跳过 fastp + kneaddata，只跑 humann
# 旧的 R2-only 产物没有 .done_catPE 标记，会被自动识别为未完成并重跑，无需手动删除。
set -euo pipefail

SAMPLE=${1:?"用法: bash $0 <sample_id> [线程数] [--rawdir DIR] [--name1 S] [--name2 S] [--outbase DIR]"}
THREADS=${2:-8}
shift 2 2>/dev/null || shift

RAWDIR="/home/usr/yel/exp-OLP_meta/pipeline/data/OLP/raw/meta/PRJNA1201607"
NAME1="_1.fastq.gz"
NAME2="_2.fastq.gz"
OUTBASE="/home/usr/yel/exp-OLP_meta/pipeline/data/OLP/shotgun"
HUMAN_DB="/home/exs/data3/Software_Database/kneaddata/human"
# kneaddata 的 trimmomatic 必须走真 jar（conda wrapper 硬编码 -Xmx1g，2G+ 输入必 OOM）：
# 它按 fnmatch "trimmomatic*" 取目录里第一个匹配——share 目录里 wrapper 排在 jar 前会命中 wrapper，
# 所以专建一个只含 jar 软链的目录喂给它，--max-memory 的 -Xmx 才真正生效（2026-09 实踩三次 OOM）
TRIM_DIR="/home/usr/yel/.local/share/kneaddata_trimmomatic"
TRIM_JAR="/home/usr/yel/miniconda3/envs/biobakery/share/trimmomatic/trimmomatic.jar"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --rawdir)  RAWDIR="$2"; shift 2 ;;
        --name1)   NAME1="$2"; shift 2 ;;
        --name2)   NAME2="$2"; shift 2 ;;
        --outbase) OUTBASE="$2"; shift 2 ;;
        --db)      HUMAN_DB="$2"; shift 2 ;;
        *) echo "未知参数: $1"; exit 1 ;;
    esac
done

OUTDIR="${OUTBASE}/${SAMPLE}"
HDIR="${OUTDIR}/humann"
LOGDIR="${OUTBASE}/logs"
LOG="${LOGDIR}/32_${SAMPLE}.log"
DONE_FLAG="${HDIR}/.done_catPE"          # 新流程完成标记
CAT_FQ="${HDIR}/${SAMPLE}_paired_cat.fastq"
mkdir -p "${LOGDIR}"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biobakery

mkdir -p "${OUTDIR}/fastp" "${OUTDIR}/kneaddata" "${HDIR}"
mkdir -p "${TRIM_DIR}"
ln -sf "${TRIM_JAR}" "${TRIM_DIR}/trimmomatic.jar"   # 幂等：保证 TRIM_DIR 里只有 jar

# 关键：切到本样本自己的 humann 目录再跑。metaphlan 把 fifo_map.mapout.txt 按固定相对名
# 写进 cwd——若多个样本共享 cwd 会互相踩踏（2026-09 实踩：12 并发全崩 ValueError）
cd "${HDIR}"

exec > >(tee -a "$LOG") 2>&1
echo "[$(date '+%H:%M:%S')] ==== ${SAMPLE} 开始 ===="
echo "  rawdir=${RAWDIR}  name1=${NAME1}  name2=${NAME2}  out=${OUTDIR}  threads=${THREADS}"

# 0) 整样本已完成（新流程）→ 跳过。放在任何 rm 之前，避免批量重跑时清掉已完成样本
if [[ -f "${DONE_FLAG}" && -s "${HDIR}/${SAMPLE}_ec.tsv" && -s "${HDIR}/${SAMPLE}_genefamilies.tsv" ]]; then
    echo "[$(date '+%H:%M:%S')] 已完成（${DONE_FLAG} 存在），跳过"
    exit 0
fi

# 拼接文件体积大（未压缩），无论成功失败都删掉
trap 'rm -f "${CAT_FQ}"' EXIT

# 断点续跑：kneaddata 配对产物在就直接进 humann（前两步失败过的样本不重算）
PAIRED1="${OUTDIR}/kneaddata/${SAMPLE}_R1_kneaddata_paired_1.fastq"
PAIRED2="${OUTDIR}/kneaddata/${SAMPLE}_R1_kneaddata_paired_2.fastq"
if [[ -s "${PAIRED1}" && -s "${PAIRED2}" ]]; then
    echo "[$(date '+%H:%M:%S')] fastp+kneaddata 产物已存在，直接进 humann"
else

# 1) fastp QC（fastp 线程上限 16）
echo "[$(date '+%H:%M:%S')] fastp..."
FASTP_W=$(( THREADS < 16 ? THREADS : 16 ))
fastp -i "${RAWDIR}/${SAMPLE}${NAME1}" -o "${OUTDIR}/fastp/${SAMPLE}_R1.fastq.gz" \
      -I "${RAWDIR}/${SAMPLE}${NAME2}" -O "${OUTDIR}/fastp/${SAMPLE}_R2.fastq.gz" \
      -w "${FASTP_W}" -j "${OUTDIR}/fastp/${SAMPLE}.json" -h "${OUTDIR}/fastp/${SAMPLE}.html"

# 2) kneaddata 去宿主
# --max-memory 必给：Trimmomatic 默认 java 堆仅 500m，2G+ 的 gz 输入会 OOM（2026-09 实踩）
echo "[$(date '+%H:%M:%S')] kneaddata..."
kneaddata --input1 "${OUTDIR}/fastp/${SAMPLE}_R1.fastq.gz" \
          --input2 "${OUTDIR}/fastp/${SAMPLE}_R2.fastq.gz" \
          -db "${HUMAN_DB}" \
          --output "${OUTDIR}/kneaddata" \
          --threads "${THREADS}" --remove-intermediate-output \
          --max-memory 8g --trimmomatic "${TRIM_DIR}"

fi  # 断点续跑 else 结束

# 3) humann3
# 注意：--taxonomic-profile 是【输入】选项（跳过内部 metaphlan 用现成 profile），不是输出！
# 要物种表就事后从 humann_temp 拷贝（humann 内部 metaphlan 的 bugs_list 留在里面，免费）
# HUMAnN 3 不支持双端输入：--input 是 argparse store（传两次后者覆盖前者），旧版本只分析了 R2。
# 双端需拼接成单文件喂入（humann 只做读段注释，pair 信息无意义）；
# 用 --output-basename 固定输出命名，避免 *_paired_2_* 这种靠 glob 猜的名字。
# 先清残留：上次中断留下的空 humann_temp 会让 humann 假退出 0 却不产任何文件（2026-09 实踩）；
# 旧的 R2-only 产物（*.tsv / *.log）一并清掉，保证重跑干净。
rm -rf "${HDIR}"/*_humann_temp "${HDIR}"/*.tsv "${HDIR}"/*_humann.log "${DONE_FLAG}"

echo "[$(date '+%H:%M:%S')] 拼接 R1+R2 → $(basename "${CAT_FQ}")..."
cat "${PAIRED1}" "${PAIRED2}" > "${CAT_FQ}"

echo "[$(date '+%H:%M:%S')] humann..."
humann --input "${CAT_FQ}" \
       --output "${HDIR}" --output-basename "${SAMPLE}" \
       --threads "${THREADS}"
rm -f "${CAT_FQ}"

# humann 可能"假成功"（退出 0 但无产物），显式检查
GENEFAM="${HDIR}/${SAMPLE}_genefamilies.tsv"
if [[ ! -s "${GENEFAM}" ]]; then
    echo "[error] humann 未产出 ${GENEFAM}，检查 ${HDIR}/${SAMPLE}_humann_temp/ 下的日志"
    exit 1
fi

# 3b) 从 humann_temp 提取 metaphlan 物种表 + humann 日志（含比对率），然后删掉 temp
# temp 里的 bowtie2 .sam / diamond 比对 / 未比对 .fa 每样本可达十几~几十 GB
TEMP="${HDIR}/${SAMPLE}_humann_temp"
if [[ -d "${TEMP}" ]]; then
    BUGS=$(ls "${TEMP}"/*_metaphlan_bugs_list.tsv 2>/dev/null | head -1 || true)
    if [[ -n "${BUGS}" ]]; then
        cp "${BUGS}" "${HDIR}/${SAMPLE}_metaphlan_bugs_list.tsv"
        echo "[$(date '+%H:%M:%S')] metaphlan bugs_list → humann/${SAMPLE}_metaphlan_bugs_list.tsv"
    else
        echo "[warn] humann_temp 里没有 metaphlan bugs_list（metaphlan 可能失败）"
    fi
    HLOG=$(ls "${TEMP}"/*.log 2>/dev/null | head -1 || true)
    [[ -n "${HLOG}" ]] && cp "${HLOG}" "${HDIR}/${SAMPLE}_humann.log"
    rm -rf "${TEMP}"
    echo "[$(date '+%H:%M:%S')] 已删除 humann_temp"
fi

# 4) EC 重分组
echo "[$(date '+%H:%M:%S')] EC regroup..."
humann_regroup_table \
    --input "${GENEFAM}" \
    --groups uniref90_level4ec \
    --output "${HDIR}/${SAMPLE}_ec.tsv"

# 全部成功才写完成标记
date '+%F %T' > "${DONE_FLAG}"
echo "[$(date '+%H:%M:%S')] ==== ${SAMPLE} 完成 ===="
ls -la "${HDIR}/"
