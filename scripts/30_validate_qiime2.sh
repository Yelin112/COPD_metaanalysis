#!/bin/bash
# QIIME2 结果一键验证脚本
# 输出各数据集：样本数、输入reads、过滤率、合并率、去嵌合率、最终保留率、中位reads数
#
# 用法：
#   bash scripts/30_validate_qiime2.sh [--config config/OLP/16s_datasets.yaml]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/../config/OLP/16s_datasets.yaml"
[[ $# -gt 0 && "$1" == "--config" ]] && CONFIG="$2"

# 从 config 读 active 数据集（单一事实来源）
mapfile -t DATASETS < <(python3 - "${CONFIG}" <<'PYEOF'
import sys, yaml
cfg = yaml.safe_load(open(sys.argv[1]))
for d in cfg["datasets"]:
    if d.get("status") == "active":
        print(d["id"])
PYEOF
)
OUTBASE="$(python3 -c "
import yaml
print(yaml.safe_load(open('${CONFIG}'))['defaults']['out_dir'])")"

printf "\n%-16s %5s %9s %8s %8s %9s %7s %10s\n" \
    "Dataset" "Samp" "Med_in" "Filt%" "Merge%" "Chimera%" "Final%" "Med_reads"
printf '%s\n' "$(printf '─%.0s' {1..85})"

for proj in "${DATASETS[@]}"; do
    dir="${OUTBASE}/${proj}"
    ft="${dir}/exported/table/feature-table.tsv"
    stats_qza="${dir}/dada2-stats.qza"

    # 未完成
    if [[ ! -f "${ft}" ]]; then
        printf "%-16s  %s\n" "${proj}" "⏳ 未完成"
        continue
    fi

    # 导出 dada2-stats（如果还没导出）
    stats_tsv="${dir}/exported/dada2-stats/stats.tsv"
    if [[ ! -f "${stats_tsv}" ]] && [[ -f "${stats_qza}" ]]; then
        qiime tools export --input-path "${stats_qza}" \
            --output-path "${dir}/exported/dada2-stats" 2>/dev/null
    fi

    # 解析并输出
    if [[ -f "${stats_tsv}" ]]; then
        python3 - "${proj}" "${ft}" "${stats_tsv}" <<'PYEOF'
import sys, csv, statistics

proj, ft_path, stats_path = sys.argv[1], sys.argv[2], sys.argv[3]

# 按表头名取列（双端表有 merged 列，单端表没有；此前按列位置解析在两种表上都错位）
rows = []
with open(stats_path) as f:
    header = None
    for line in f:
        if line.startswith('#q2:types'):
            continue
        parts = line.strip().split('\t')
        if parts and parts[0] == 'sample-id':
            header = parts
            continue
        if header is None or not parts:
            continue
        d = dict(zip(header, parts))
        try:
            inp  = float(d['input'])
            filt = float(d['filtered'])
            deno = float(d['denoised'])
            nc   = float(d['non-chimeric'])
            mrg  = float(d['merged']) if 'merged' in d else None   # 单端无 merged
            rows.append((inp, filt, deno, mrg, nc))
        except (ValueError, KeyError):
            pass

if not rows:
    print(f"{proj:<16} {'?':>5} {'?':>9} {'?':>8} {'?':>8} {'?':>9} {'?':>7} {'?':>10}")
    sys.exit()

n          = len(rows)
med_in     = statistics.median(r[0] for r in rows)
med_filt   = statistics.median(r[1]/r[0]*100 if r[0] else 0 for r in rows)
# 双端: merge% = merged/denoised, chimera% = non-chimeric/merged；单端: merge 不适用, chimera% = non-chimeric/denoised
med_merge  = statistics.median(r[3]/r[2]*100 if (r[3] is not None and r[2]) else 100 for r in rows)
_den       = lambda r: r[3] if (r[3] is not None and r[3] > 0) else r[2]
med_chim   = statistics.median(r[4]/_den(r)*100 if _den(r) else 0 for r in rows)
med_final  = statistics.median(r[4]/r[0]*100 if r[0] else 0 for r in rows)
med_nc     = statistics.median(r[4] for r in rows)

# 简单判断是否正常
ok = med_final >= 30 and med_nc >= 500 and n >= 5
flag = "✅" if ok else "⚠️ "

print(f"{flag} {proj:<14} {n:>5} {med_in:>9.0f} {med_filt:>7.1f}% "
      f"{med_merge:>7.1f}% {med_chim:>8.1f}% {med_final:>6.1f}% {med_nc:>10.0f}")
PYEOF

    else
        # stats 不可用，只看 feature table
        n_samp=$(awk 'NR==2{print NF-1}' "${ft}" 2>/dev/null || echo "?")
        n_feat=$(awk 'NR>2{c++}END{print c}' "${ft}" 2>/dev/null || echo "?")
        printf "⚠️  %-14s %5s  (stats未导出；ASV数=%s)\n" "${proj}" "${n_samp}" "${n_feat}"
    fi
done

printf '%s\n' "$(printf '─%.0s' {1..85})"
printf "判断标准：Final%% ≥ 30%% 且 中位reads ≥ 500 且 样本数 ≥ 5\n\n"
