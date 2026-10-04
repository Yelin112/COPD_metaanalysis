#!/bin/bash
# PICRUSt2 功能预测脚本（QIIME2 输出 → EC 丰度表）
#
# 用法：
#   conda activate picrust2
#   cd /home/usr/yel/exp-OLP_meta/pipeline
#   bash scripts/13_picrust2.sh PRJNA542018 [--config config/OLP/16s_datasets.yaml]

set -euo pipefail

PROJ="${1:?请传入数据集名称，如：bash scripts/13_picrust2.sh PRJNA542018}"
shift
CONFIG="config/OLP/16s_datasets.yaml"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) CONFIG="$2"; shift 2 ;;
        *) echo "未知参数: $1"; exit 1 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ "${CONFIG}" = /* ]] || CONFIG="${SCRIPT_DIR}/../${CONFIG}"

# 路径从 config 读（defaults），不再硬编码 OLP
eval "$(/usr/bin/python3 - "${CONFIG}" <<'PYEOF'
import sys, yaml, shlex
d = yaml.safe_load(open(sys.argv[1]))["defaults"]
out = {
    "QIIME_BASE": d.get("out_dir", ""),
    "RAW_DIR": d.get("raw_dir", ""),
    "PICRU2_DIR": d.get("picrust2_dir", ""),
}
for k, v in out.items():
    print(f"{k}={shlex.quote(str(v))}")
PYEOF
)"

QIIME_OUT="${QIIME_BASE}/${PROJ}/exported"
OUTDIR="${PICRU2_DIR}/${PROJ}"
THREADS=8

SEQS="${QIIME_OUT}/rep-seqs/dna-sequences.fasta"
TABLE="${QIIME_OUT}/table/feature-table.biom"

mkdir -p logs
# 注意：不要预建 OUTDIR——picrust2_pipeline.py 要求输出目录不存在，由它自己创建

echo "[$(date '+%H:%M:%S')] ── ${PROJ} PICRUSt2 开始 ──"

# ── PICRUSt2 一步式流程 ───────────────────────────────────────────
# 2.6.x 参数：-s/-i/-o 短形式；无 --threads/--study_fasta/--input_table（老版脚本里的长参数全都不存在）
# --wide_table 必加：不加时 --stratified 输出的是 long 格式 pred_metagenome_contrib.tsv.gz，
# 加了才写宽表 pred_metagenome_strat.tsv.gz（列=function/taxon/样本），下游 01_picrust2_to_ec.py 要的就是这个
picrust2_pipeline.py \
    -s "${SEQS}" \
    -i "${TABLE}" \
    -o "${OUTDIR}" \
    --stratified \
    --wide_table \
    --verbose

# ── 提取 EC 丰度表 ────────────────────────────────────────────────
# 必须用 strat 表：下游 modules/meta/metagenomics/01_picrust2_to_ec.py 期望
# 列1=EC 列2=ASV 列3+=样本（stratified 三层表）；unstrat 表缺 ASV 列会静默丢一个样本
EC_FILE="${OUTDIR}/EC_metagenome_out/pred_metagenome_strat.tsv.gz"
EC_OUT="${RAW_DIR}/${PROJ}.tsv"

echo "[$(date '+%H:%M:%S')] 提取 EC 丰度表..."
gunzip -c "${EC_FILE}" > "${EC_OUT}"

echo "[$(date '+%H:%M:%S')] ✅ ${PROJ} PICRUSt2 完成"
echo "  EC 丰度表 → ${EC_OUT}"
echo "  （可直接用于 Snakemake meta_metagenomics 流程）"
