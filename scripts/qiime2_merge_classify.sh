#!/bin/bash
# 16S 合并 + 分区物种分类
# 按 amplicon_region 选分类器对每个 active 数据集的 rep-seqs 分类，
# 导出 taxonomy.tsv 到下游契约路径 data/OLP/raw/meta/{id}.taxonomy.tsv；
# 并合并所有 feature table 为总览表 data/OLP/qiime2/merged/feature-table.tsv
#
# 用法：
#   conda activate qiime2-amplicon-2023
#   bash scripts/qiime2_merge_classify.sh [--config config/OLP/16s_datasets.yaml]
# 幂等：可重跑，缺失 rep-seqs 的数据集自动跳过（其 DADA2 完成后重跑补全）

set -euo pipefail
export R_ENVIRON_USER=/dev/null  # 跳过用户 .Renviron（系统 R 4.5.1 包库会污染 qiime2 env 的 R 4.2.2）

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/../config/OLP/16s_datasets.yaml"
SKIP_CHECK=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) CONFIG="$2"; shift 2 ;;
        --skip-classifier-check) SKIP_CHECK=1; shift ;;
        *) echo "未知参数: $1"; exit 1 ;;
    esac
done

eval "$(/usr/bin/python3 - "${CONFIG}" <<'PYEOF'
import sys, yaml, shlex
cfg = yaml.safe_load(open(sys.argv[1]))
out = {}
for region, path in (cfg.get("classifiers") or {}).items():
    out[f"CLS_{region.replace('-', '_')}"] = path   # bash 变量名不允许连字符
for d in cfg["datasets"]:
    if d.get("status") != "active":
        continue
    d = {**cfg["defaults"], **d}                     # defaults 兜底（与 qiime2_16s.sh 一致）
    out[f"DS_{d['id']}_REGION"] = d.get("amplicon_region", "").replace("-", "_")
    out[f"DS_{d['id']}_OUT"] = d.get("out_dir", "")
out["OUT_BASE"] = cfg["defaults"].get("out_dir", "data/OLP/qiime2")
out["RAW_DIR"] = cfg["defaults"].get("raw_dir", "data/OLP/raw/meta")
for k, v in out.items():
    print(f"{k}={shlex.quote(str(v))}")
PYEOF
)"

# active 数据集清单
mapfile -t DATASETS < <(/usr/bin/python3 - "${CONFIG}" <<'PYEOF'
import sys, yaml
cfg = yaml.safe_load(open(sys.argv[1]))
for d in cfg["datasets"]:
    if d.get("status") == "active":
        print(d["id"])
PYEOF
)

# 预检分类器（缺失时提示自训；--skip-classifier-check 时对应数据集跳过分类）
declare -A NEEDED
for ds in "${DATASETS[@]}"; do
    region_var="DS_${ds}_REGION"
    NEEDED["${!region_var}"]=1
done
declare -A MISSING
for region in "${!NEEDED[@]}"; do
    clsf_var="CLS_${region}"
    clsf="${!clsf_var:-}"
    if [[ -z "${clsf}" || ! -f "${clsf}" ]]; then
        MISSING["${region}"]="${clsf:-未配置}"
    fi
done
if [[ -n "${MISSING[@]+x}" ]]; then
    for region in "${!MISSING[@]}"; do
        echo "警告: 区域 ${region} 分类器缺失: ${MISSING[$region]}"
    done
    if [[ ${SKIP_CHECK} -ne 1 ]]; then
        echo "自训：bash scripts/qiime2_train_classifier.sh --p-f-primer ... --p-r-primer ..."
        echo "或先跑其他区域：加 --skip-classifier-check（缺分类器的数据集将跳过）"
        exit 1
    fi
fi

echo "[$(date '+%H:%M:%S')] ── 分区分类 + 合并开始 ──"

# 逐数据集分类
TABLES=()
for ds in "${DATASETS[@]}"; do
    region_var="DS_${ds}_REGION"
    out_var="DS_${ds}_OUT"
    region="${!region_var}"
    outdir="${!out_var}"
    clsf_var="CLS_${region}"
    clsf="${!clsf_var}"

    if [[ -z "${clsf}" || ! -f "${clsf}" ]]; then
        echo "[skip] ${ds}: 分类器缺失（${region}），自训完成后重跑"
        continue
    fi

    if [[ ! -f "${outdir}/${ds}/rep-seqs.qza" ]]; then
        echo "[skip] ${ds}: rep-seqs.qza 缺失（DADA2 未完成），重跑本脚本补全"
        continue
    fi

    echo "[$(date '+%H:%M:%S')] 分类 ${ds}（${region}）..."
    qiime feature-classifier classify-sklearn \
        --i-classifier "${clsf}" \
        --i-reads "${outdir}/${ds}/rep-seqs.qza" \
        --o-classification "${outdir}/${ds}/taxonomy.qza"

    qiime tools export --input-path "${outdir}/${ds}/taxonomy.qza" \
        --output-path "${outdir}/${ds}/exported/taxonomy"

    # 下游契约路径（study_config.yaml 的 taxonomy_file）
    cp "${outdir}/${ds}/exported/taxonomy/taxonomy.tsv" "${RAW_DIR}/${ds}.taxonomy.tsv"
    echo "  → ${RAW_DIR}/${ds}.taxonomy.tsv"

    [[ -f "${outdir}/${ds}/table.qza" ]] && TABLES+=("${outdir}/${ds}/table.qza")
done

# 合并总表（辅助总览产物；下游主流程仍用 per-dataset 文件）
if [[ ${#TABLES[@]} -gt 0 ]]; then
    echo "[$(date '+%H:%M:%S')] 合并 ${#TABLES[@]} 个 feature table..."
    mkdir -p "${OUT_BASE}/merged"
    qiime feature-table merge \
        --i-tables "${TABLES[@]}" \
        --o-merged-table "${OUT_BASE}/merged/table.qza"
    qiime tools export --input-path "${OUT_BASE}/merged/table.qza" \
        --output-path "${OUT_BASE}/merged/exported"
    biom convert \
        -i "${OUT_BASE}/merged/exported/feature-table.biom" \
        -o "${OUT_BASE}/merged/feature-table.tsv" \
        --to-tsv
    echo "  → ${OUT_BASE}/merged/feature-table.tsv"
else
    echo "警告: 无可用 table.qza，跳过合并"
fi

echo "[$(date '+%H:%M:%S')] ✅ 分区分类 + 合并完成"
echo "  下一步：bash scripts/13_picrust2.sh <PROJ>（各数据集功能预测）"
