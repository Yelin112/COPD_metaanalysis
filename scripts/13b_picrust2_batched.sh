#!/bin/bash
# PICRUSt2 分批运行（解决大样本量 OOM）
#
# 适用场景：样本数 > ~100 时整批 --stratified --wide_table 内存爆炸
# 策略：拆分 biom 表 → 每批独立跑 picrust2_pipeline → 合并 strat 表
#
# 用法：
#   conda activate picrust2
#   cd /home/usr/yel/exp-OLP_meta/pipeline
#   bash scripts/13b_picrust2_batched.sh CRA008410 [--batch-size 60] [--config config/OLP/16s_datasets.yaml]

set -euo pipefail

PROJ="${1:?请传入数据集名称}"
shift

BATCH_SIZE=60
CONFIG="config/OLP/16s_datasets.yaml"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --batch-size) BATCH_SIZE="$2"; shift 2 ;;
        --config)     CONFIG="$2";     shift 2 ;;
        *) echo "未知参数: $1"; exit 1 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ "${CONFIG}" = /* ]] || CONFIG="${SCRIPT_DIR}/../${CONFIG}"

eval "$(/usr/bin/python3 - "${CONFIG}" <<'PYEOF'
import sys, yaml, shlex
d = yaml.safe_load(open(sys.argv[1]))["defaults"]
for k, v in {"QIIME_BASE": d.get("out_dir",""), "PICRU2_DIR": d.get("picrust2_dir",""), "RAW_DIR": d.get("raw_dir","")}.items():
    print(f"{k}={shlex.quote(str(v))}")
PYEOF
)"

QIIME_OUT="${QIIME_BASE}/${PROJ}/exported"
OUTDIR="${PICRU2_DIR}/${PROJ}"
SEQS="${QIIME_OUT}/rep-seqs/dna-sequences.fasta"
TABLE="${QIIME_OUT}/table/feature-table.biom"
BATCH_DIR="${OUTDIR}/_batches"
EC_OUT="${RAW_DIR}/${PROJ}.tsv"

mkdir -p "${BATCH_DIR}" logs

echo "[$(date '+%H:%M:%S')] ── ${PROJ} 分批 PICRUSt2 开始（batch-size=${BATCH_SIZE}）──"

# ── 1. 获取全部样本 ID ─────────────────────────────────────────────
SAMPLE_LIST="${BATCH_DIR}/all_samples.txt"
python3 - "${TABLE}" "${SAMPLE_LIST}" <<'PYEOF'
import sys
from biom import load_table
t = load_table(sys.argv[1])
with open(sys.argv[2], "w") as f:
    for s in t.ids(axis="sample"):
        f.write(s + "\n")
PYEOF

N_SAMPLES=$(wc -l < "${SAMPLE_LIST}")
echo "  总样本数：${N_SAMPLES}，每批：${BATCH_SIZE}"

# ── 2. 拆分 biom 并逐批运行 ───────────────────────────────────────
python3 - "${SAMPLE_LIST}" "${BATCH_SIZE}" "${BATCH_DIR}" <<'PYEOF'
import sys, math
samples = open(sys.argv[1]).read().splitlines()
bs      = int(sys.argv[2])
outdir  = sys.argv[3]
n_batch = math.ceil(len(samples) / bs)
for i in range(n_batch):
    chunk = samples[i*bs:(i+1)*bs]
    with open(f"{outdir}/batch_{i:02d}_samples.txt", "w") as f:
        f.write("\n".join(chunk) + "\n")
    print(f"batch_{i:02d}: {len(chunk)} 样本")
PYEOF

BATCH_EC_FILES=()
for batch_samples in "${BATCH_DIR}"/batch_*_samples.txt; do
    BNAME=$(basename "${batch_samples}" _samples.txt)
    BIOM_BATCH="${BATCH_DIR}/${BNAME}.biom"
    BATCH_OUT="${BATCH_DIR}/${BNAME}_picrust2"

    # 子集 biom
    if [[ ! -f "${BIOM_BATCH}" ]]; then
        echo "[$(date '+%H:%M:%S')] 拆分 ${BNAME}..."
        biom subset-table \
            -i "${TABLE}" \
            -a sample \
            -s "${batch_samples}" \
            -o "${BIOM_BATCH}"
    fi

    # 跑 picrust2（跳过已完成的批次）
    EC_STRAT="${BATCH_OUT}/EC_metagenome_out/pred_metagenome_strat.tsv.gz"
    if [[ ! -f "${EC_STRAT}" ]]; then
        echo "[$(date '+%H:%M:%S')] 运行 ${BNAME}..."
        picrust2_pipeline.py \
            -s "${SEQS}" \
            -i "${BIOM_BATCH}" \
            -o "${BATCH_OUT}" \
            --stratified \
            --wide_table \
            --verbose
    else
        echo "[$(date '+%H:%M:%S')] ${BNAME} 已完成，跳过"
    fi

    BATCH_EC_FILES+=("${EC_STRAT}")
done

# ── 3. 合并各批次 strat 表 ────────────────────────────────────────
echo "[$(date '+%H:%M:%S')] 合并 ${#BATCH_EC_FILES[@]} 批次 strat 表..."

MERGED_GZ="${OUTDIR}/EC_metagenome_out/pred_metagenome_strat.tsv.gz"
mkdir -p "${OUTDIR}/EC_metagenome_out"

python3 - "${MERGED_GZ}" "${BATCH_EC_FILES[@]}" <<'PYEOF'
import sys, gzip, csv
from collections import defaultdict

out_gz   = sys.argv[1]
in_files = sys.argv[2:]

# 读所有批次，以 (function, taxon) 为 key 聚合
data   = defaultdict(dict)   # (func, taxon) → {sample: value}
header = None
all_samples = []

for fpath in in_files:
    with gzip.open(fpath, "rt") as f:
        reader = csv.reader(f, delimiter="\t")
        h = next(reader)   # function, taxon, samp1, samp2, ...
        batch_samples = h[2:]
        if header is None:
            header = h[:2]
        all_samples.extend(batch_samples)
        for row in reader:
            key = (row[0], row[1])
            for i, s in enumerate(batch_samples):
                v = row[2 + i] if (2 + i) < len(row) else "0"
                data[key][s] = v

# 写合并表
with gzip.open(out_gz, "wt") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(header + all_samples)
    for key in sorted(data.keys()):
        row = list(key) + [data[key].get(s, "0") for s in all_samples]
        writer.writerow(row)

print(f"  合并完成：{len(data)} 行 × {len(all_samples)+2} 列 → {out_gz}")
PYEOF

# ── 4. 提取 EC 表（与 13_picrust2.sh 相同逻辑）─────────────────────
echo "[$(date '+%H:%M:%S')] 提取 EC 丰度表..."
gunzip -c "${MERGED_GZ}" > "${EC_OUT}"

echo "[$(date '+%H:%M:%S')] ✅ ${PROJ} 分批 PICRUSt2 完成"
echo "  EC 丰度表 → ${EC_OUT}"
echo "  中间批次结果保留于 ${BATCH_DIR}（可手动清理）"
