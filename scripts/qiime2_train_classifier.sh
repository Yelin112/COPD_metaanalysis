#!/bin/bash
# SILVA 138 区域分类器自训（参数化，可复用：换引物即可训其他区域）
# 步骤：extract-reads 截取引物区域（10-20min）→ fit-classifier-naive-bayes（1-2h，内存 16GB+）
#
# 用法：
#   conda activate qiime2-amplicon-2023
#   bash scripts/qiime2_train_classifier.sh \
#       --i-sequences   /path/to/silva-138-99-seqs.qza \
#       --i-taxonomy    /path/to/silva-138-99-tax.qza \
#       --p-f-primer    AGRGTTYGATYMTGGCTCAG \
#       --p-r-primer    TGCTGCCTCCCGTAGGAGT \
#       --o-classifier  /path/to/out-classifier.qza \
#       [--p-min-length 100] [--p-max-length 400]
#
# 训练完成后把 --o-classifier 的路径写入 config/OLP/16s_datasets.yaml 的 classifiers.<区域>
set -euo pipefail
export R_ENVIRON_USER=/dev/null  # 跳过用户 .Renviron（系统 R 4.5.1 包库会污染 qiime2 env 的 R 4.2.2）

SEQ=""; TAX=""; F=""; R=""; OUT=""; MINLEN=100; MAXLEN=400
while [[ $# -gt 0 ]]; do
    case "$1" in
        --i-sequences)  SEQ="$2"; shift 2 ;;
        --i-taxonomy)   TAX="$2"; shift 2 ;;
        --p-f-primer)   F="$2"; shift 2 ;;
        --p-r-primer)   R="$2"; shift 2 ;;
        --o-classifier) OUT="$2"; shift 2 ;;
        --p-min-length) MINLEN="$2"; shift 2 ;;
        --p-max-length) MAXLEN="$2"; shift 2 ;;
        *) echo "未知参数: $1"; exit 1 ;;
    esac
done
[[ -n "${SEQ}" && -n "${TAX}" && -n "${F}" && -n "${R}" && -n "${OUT}" ]] \
    || { echo "用法: $0 --i-sequences SEQ --i-taxonomy TAX --p-f-primer F --p-r-primer R --o-classifier OUT"; exit 1; }

mkdir -p logs
EXTRACTED="${OUT%.qza}-reads.qza"

echo "[$(date '+%H:%M:%S')] 步骤1：extract-reads 截取引物区域..."
qiime feature-classifier extract-reads \
    --i-sequences "${SEQ}" \
    --p-f-primer "${F}" \
    --p-r-primer "${R}" \
    --p-trunc-len 0 \
    --p-min-length "${MINLEN}" \
    --p-max-length "${MAXLEN}" \
    --p-n-jobs 8 \
    --o-reads "${EXTRACTED}"

echo "[$(date '+%H:%M:%S')] 步骤2：fit-classifier-naive-bayes（约 1-2h）..."
nohup qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "${EXTRACTED}" \
    --i-reference-taxonomy "${TAX}" \
    --o-classifier "${OUT}" \
    > "logs/train_classifier_$(basename ${OUT}).log" 2>&1 &

echo "训练已在后台启动（PID: $!），日志: logs/train_classifier_$(basename ${OUT}).log"
echo "完成后把分类器路径写入 config 的 classifiers 对应区域"
