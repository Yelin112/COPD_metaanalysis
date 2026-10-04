#!/bin/bash
# 16S QIIME2 通用流程：manifest → import → cutadapt → DADA2 → export
# 参数全部来自 config（默认 config/OLP/16s_datasets.yaml），新数据集只需加 config 条目
#
# 用法：
#   conda activate qiime2-amplicon-2023
#   cd /home/usr/yel/exp-OLP_meta/pipeline
#   bash scripts/qiime2_16s.sh <PROJ> [--config PATH] [--dry-run]
#   nohup bash scripts/qiime2_16s.sh <PROJ> > logs/qiime2_<PROJ>.log 2>&1 &
#
# 输出路径与旧 per-dataset 脚本逐字节一致：data/OLP/qiime2/{PROJ}/...
set -euo pipefail

PROJ="${1:?用法: bash scripts/qiime2_16s.sh <PROJ> [--config PATH] [--dry-run]}"
CONFIG="config/OLP/16s_datasets.yaml"
DRY=0
shift
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) CONFIG="$2"; shift 2 ;;
        --dry-run) DRY=1; shift ;;
        *) echo "未知参数: $1"; exit 1 ;;
    esac
done

# config 路径相对脚本自身目录解析（与 cwd 无关）
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ "${CONFIG}" = /* ]] || CONFIG="${SCRIPT_DIR}/../${CONFIG}"
[[ -f "${CONFIG}" ]] || { echo "错误: config 不存在: ${CONFIG}"; exit 1; }

# R 库路径隔离：跳过用户级 .Renviron（其 R_LIBS_USER 指向系统 R 4.5.1 的包库，
# 会与 qiime2 env 的 R 4.2.2 冲突，导致 Rcpp 等包 dyn.load 失败）
export R_ENVIRON_USER=/dev/null
CONDA_BASE="$(command -v conda >/dev/null && conda info --base 2>/dev/null || true)"
[[ -n "${CONDA_BASE}" && -d "${CONDA_BASE}/envs/qiime2-amplicon-2023/lib/R/library" ]] \
    && export R_LIBS_SITE="${CONDA_BASE}/envs/qiime2-amplicon-2023/lib/R/library"

# 读 config：系统 python3 + pyyaml，扁平化为 shell 变量
eval "$(/usr/bin/python3 - "${PROJ}" "${CONFIG}" <<'PYEOF'
import sys, yaml, shlex
proj, cfg_path = sys.argv[1], sys.argv[2]
cfg = yaml.safe_load(open(cfg_path))
ds = next((d for d in cfg["datasets"] if d["id"] == proj), None)
if ds is None:
    print(f'echo "错误: {proj} 不在 {cfg_path} 中"; exit 1', file=sys.stderr)
    sys.exit(1)
d = {**cfg["defaults"], **ds}          # defaults 兜底，条目覆盖
out = {
    "STATUS": d.get("status", "active"),
    "REASON": d.get("exclude_reason", ""),
    "LAYOUT": d.get("layout", ""),
    "REGION": d.get("amplicon_region", ""),
    "GLOB": d.get("read_glob", ""),
    "SFX1": d.get("sample_suffix", ""),
    "SFX2": d.get("r2_suffix") or "",
    "RAW": d["raw_dir"],
    "FASTQ_DIR": d.get("fastq_dir") or "",
    "OUT": d["out_dir"],
    "THREADS": d["threads"],
    "MIN_LENGTH": d.get("min_length", 100),
    "TRIM_LEFT": d.get("trim_left", 0),
    "CA_ENABLED": "1" if d.get("cutadapt", {}).get("enabled") else "0",
    "CA_MODE": d.get("cutadapt", {}).get("mode", ""),
    "PRIMER_F": (d.get("primers") or {}).get("fwd") or "",
    "PRIMER_R": (d.get("primers") or {}).get("rev") or "",
}
dd = d.get("dada2", {})
if d.get("layout") == "single":
    out["TRUNC"] = dd.get("trunc_len", 0); out["MAXEE"] = dd.get("max_ee", 2.0)
else:
    out["TRUNC_F"] = dd.get("trunc_len_f", 0); out["TRUNC_R"] = dd.get("trunc_len_r", 0)
    out["MAXEE_F"] = dd.get("max_ee_f", 2.0); out["MAXEE_R"] = dd.get("max_ee_r", 2.0)
for k, v in out.items():
    print(f"{k}={shlex.quote(str(v))}")
PYEOF
)"

# excluded 早退
if [[ "${STATUS}" != "active" ]]; then
    echo "[skip] ${PROJ}: ${REASON}"
    exit 0
fi

RAWDIR="${FASTQ_DIR:-${RAW}/${PROJ}}"
OUTDIR="${OUT}/${PROJ}"

# dry-run：只打印参数清单，不建目录不碰数据
if [[ ${DRY} -eq 1 ]]; then
    echo "── ${PROJ} dry-run 参数清单 ──"
    echo "  layout    : ${LAYOUT}"
    echo "  region    : ${REGION}"
    echo "  rawdir    : ${RAWDIR}"
    echo "  glob      : ${GLOB}"
    echo "  suffix    : ${SFX1} / ${SFX2:-无(R2)}"
    if [[ ${CA_ENABLED} -eq 1 ]]; then
        echo "  cutadapt  : 启用, mode=${CA_MODE}"
        echo "  primers   : ${PRIMER_F} / ${PRIMER_R}"
    else
        echo "  cutadapt  : 跳过（引物已截除）"
    fi
    if [[ "${LAYOUT}" == "single" ]]; then
        echo "  dada2     : trunc-len=${TRUNC} max-ee=${MAXEE}"
    else
        echo "  dada2     : trunc ${TRUNC_F}/${TRUNC_R} max-ee ${MAXEE_F}/${MAXEE_R}"
    fi
    echo "  threads   : ${THREADS}"
    exit 0
fi

mkdir -p "${OUTDIR}" logs
echo "[$(date '+%H:%M:%S')] ── ${PROJ} QIIME2 流程开始 ──"

# Step 1: manifest
echo "[$(date '+%H:%M:%S')] 生成 manifest..."
MANIFEST="${OUTDIR}/manifest.tsv"
if [[ "${LAYOUT}" == "single" ]]; then
    echo -e "sample-id\tabsolute-filepath" > "${MANIFEST}"
    for fq in "${RAWDIR}"/${GLOB}; do
        sample=$(basename "${fq}" "${SFX1}")
        echo -e "${sample}\t$(realpath ${fq})"
    done >> "${MANIFEST}"
else
    echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "${MANIFEST}"
    for fwd in "${RAWDIR}"/${GLOB}; do
        sample=$(basename "${fwd}" "${SFX1}")
        rev="${RAWDIR}/${sample}${SFX2}"
        if [[ -f "${rev}" ]]; then
            echo -e "${sample}\t$(realpath ${fwd})\t$(realpath ${rev})"
        fi
    done >> "${MANIFEST}"
fi
echo "  样本数：$(tail -n +2 ${MANIFEST} | wc -l)"

# Step 2: import
echo "[$(date '+%H:%M:%S')] 导入 FASTQ..."
if [[ "${LAYOUT}" == "single" ]]; then
    qiime tools import \
        --type 'SampleData[SequencesWithQuality]' \
        --input-path "${MANIFEST}" \
        --input-format SingleEndFastqManifestPhred33V2 \
        --output-path "${OUTDIR}/demux.qza"
else
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "${MANIFEST}" \
        --input-format PairedEndFastqManifestPhred33V2 \
        --output-path "${OUTDIR}/demux.qza"
fi

# Step 3: cutadapt（可选）
DADA2_INPUT="${OUTDIR}/demux.qza"
if [[ ${CA_ENABLED} -eq 1 ]]; then
    echo "[$(date '+%H:%M:%S')] cutadapt 截除引物 (mode=${CA_MODE})..."
    case "${CA_MODE}" in
        single)
            qiime cutadapt trim-single \
                --i-demultiplexed-sequences "${OUTDIR}/demux.qza" \
                --p-front "${PRIMER_F}" \
                --p-adapter "${PRIMER_R}" \
                --p-discard-untrimmed \
                --p-minimum-length ${MIN_LENGTH} \
                --p-cores ${THREADS} \
                --o-trimmed-sequences "${OUTDIR}/demux_trimmed.qza" \
                --verbose 2>&1 | tail -20 ;;
        paired-front)
            qiime cutadapt trim-paired \
                --i-demultiplexed-sequences "${OUTDIR}/demux.qza" \
                --p-front-f "${PRIMER_F}" \
                --p-front-r "${PRIMER_R}" \
                --p-discard-untrimmed \
                --p-minimum-length ${MIN_LENGTH} \
                --p-cores ${THREADS} \
                --o-trimmed-sequences "${OUTDIR}/demux_trimmed.qza" \
                --verbose 2>&1 | tail -20 ;;
        paired-anywhere-r)
            # R2 可能有未知前缀时用 --p-anywhere-r 而非 --p-front-r（兼容有无前缀）
            qiime cutadapt trim-paired \
                --i-demultiplexed-sequences "${OUTDIR}/demux.qza" \
                --p-front-f "${PRIMER_F}" \
                --p-anywhere-r "${PRIMER_R}" \
                --p-discard-untrimmed \
                --p-minimum-length ${MIN_LENGTH} \
                --p-cores ${THREADS} \
                --o-trimmed-sequences "${OUTDIR}/demux_trimmed.qza" \
                --verbose 2>&1 | tail -20 ;;
        *) echo "错误: 未知 cutadapt mode=${CA_MODE}"; exit 1 ;;
    esac
    DADA2_INPUT="${OUTDIR}/demux_trimmed.qza"
fi

# Step 4: DADA2
echo "[$(date '+%H:%M:%S')] DADA2 去噪（${LAYOUT}）..."
if [[ "${LAYOUT}" == "single" ]]; then
    qiime dada2 denoise-single \
        --i-demultiplexed-seqs "${DADA2_INPUT}" \
        --p-trunc-len ${TRUNC} \
        --p-trim-left ${TRIM_LEFT} \
        --p-max-ee ${MAXEE} \
        --p-n-threads ${THREADS} \
        --o-table "${OUTDIR}/table.qza" \
        --o-representative-sequences "${OUTDIR}/rep-seqs.qza" \
        --o-denoising-stats "${OUTDIR}/dada2-stats.qza"
else
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs "${DADA2_INPUT}" \
        --p-trunc-len-f ${TRUNC_F} \
        --p-trunc-len-r ${TRUNC_R} \
        --p-trim-left-f ${TRIM_LEFT} \
        --p-trim-left-r ${TRIM_LEFT} \
        --p-max-ee-f ${MAXEE_F} \
        --p-max-ee-r ${MAXEE_R} \
        --p-n-threads ${THREADS} \
        --o-table "${OUTDIR}/table.qza" \
        --o-representative-sequences "${OUTDIR}/rep-seqs.qza" \
        --o-denoising-stats "${OUTDIR}/dada2-stats.qza"
fi

# Step 5: export
echo "[$(date '+%H:%M:%S')] 导出结果..."
qiime tools export --input-path "${OUTDIR}/table.qza" \
    --output-path "${OUTDIR}/exported/table"
qiime tools export --input-path "${OUTDIR}/rep-seqs.qza" \
    --output-path "${OUTDIR}/exported/rep-seqs"
biom convert \
    -i "${OUTDIR}/exported/table/feature-table.biom" \
    -o "${OUTDIR}/exported/table/feature-table.tsv" \
    --to-tsv

echo "[$(date '+%H:%M:%S')] ✅ ${PROJ} QIIME2 完成"
echo "  feature table : ${OUTDIR}/exported/table/feature-table.tsv"
echo "  rep seqs      : ${OUTDIR}/exported/rep-seqs/dna-sequences.fasta"
echo "  下一步：bash scripts/qiime2_merge_classify.sh（合并+分区分类）"
