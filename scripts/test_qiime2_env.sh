#!/bin/bash
# QIIME2 环境验证脚本（使用 EasyAmplicon 测试数据）
#
# 用法：
#   conda activate qiime2-amplicon-2023
#   cd /home/usr/yel/exp-OLP_meta/pipeline
#   bash scripts/test_qiime2_env.sh

set -euo pipefail

# ── R 路径修复：确保 conda 环境的 R 优先于系统 R ──────────────────
# QIIME2 DADA2 插件调用 run_dada.R，需要 conda 环境中的 R（含 optparse）
CONDA_ENV_BIN="$(conda info --base 2>/dev/null)/envs/qiime2-amplicon-2023/bin"
if [[ -f "${CONDA_ENV_BIN}/Rscript" ]]; then
    export PATH="${CONDA_ENV_BIN}:${PATH}"
fi
unset R_LIBS R_LIBS_USER R_LIBS_SITE  # 清除可能覆盖 conda R 库路径的变量

SEQDIR="/home/usr/yel/repos/EasyAmplicon/seq"
OUTDIR="results/test_qiime2"
SILVA_SEQS="/home/usr/yel/exp-16S_workflow/Qiime2/database/silva/silva-138-99-seqs.qza"
SILVA_TAX="/home/usr/yel/exp-16S_workflow/Qiime2/database/silva/silva-138-99-tax.qza"
THREADS=8

mkdir -p "${OUTDIR}" logs

echo "════════════════════════════════════════"
echo " QIIME2 环境验证  $(date '+%Y-%m-%d %H:%M')"
echo "════════════════════════════════════════"

# ── 0. 环境检查 ────────────────────────────────────────────────────
echo ""
echo "【0】环境检查"
echo "  QIIME2 版本：$(qiime --version 2>&1 | head -1)"
echo "  测试数据目录：${SEQDIR}"
echo ""

# 列出测试数据
echo "  测试数据文件："
ls "${SEQDIR}" | sed 's/^/    /'

# 判断单端还是双端
PE_COUNT=$(ls "${SEQDIR}"/*_2.f*q* 2>/dev/null | wc -l || true)
if [[ ${PE_COUNT} -gt 0 ]]; then
    SEQ_TYPE="paired"
    echo "  → 检测到双端测序数据（paired-end）"
else
    SEQ_TYPE="single"
    echo "  → 检测到单端测序数据（single-end）"
fi

# ── 1. 生成 manifest ───────────────────────────────────────────────
echo ""
echo "【1】生成 QIIME2 manifest"
MANIFEST="${OUTDIR}/manifest.tsv"

if [[ "${SEQ_TYPE}" == "paired" ]]; then
    echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "${MANIFEST}"
    for fwd in "${SEQDIR}"/*_1.f*q*; do
        [[ -f "${fwd}" ]] || continue
        sample=$(basename "${fwd}" | sed 's/_1\..*//')
        rev=$(echo "${fwd}" | sed 's/_1\./_2./')
        if [[ -f "${rev}" ]]; then
            echo -e "${sample}\t$(realpath ${fwd})\t$(realpath ${rev})"
        fi
    done >> "${MANIFEST}"
else
    echo -e "sample-id\tabsolute-filepath" > "${MANIFEST}"
    for fq in "${SEQDIR}"/*.f*q*; do
        [[ -f "${fq}" ]] || continue
        sample=$(basename "${fq}" | sed 's/\..*//')
        echo -e "${sample}\t$(realpath ${fq})"
    done >> "${MANIFEST}"
fi

SAMPLE_COUNT=$(tail -n +2 "${MANIFEST}" | wc -l)
echo "  样本数：${SAMPLE_COUNT}"
cat "${MANIFEST}"

# ── 2. 导入数据 ────────────────────────────────────────────────────
echo ""
echo "【2】导入 FASTQ 数据"
if [[ "${SEQ_TYPE}" == "paired" ]]; then
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "${MANIFEST}" \
        --input-format PairedEndFastqManifestPhred33V2 \
        --output-path "${OUTDIR}/demux.qza"
else
    qiime tools import \
        --type 'SampleData[SequencesWithQuality]' \
        --input-path "${MANIFEST}" \
        --input-format SingleEndFastqManifestPhred33V2 \
        --output-path "${OUTDIR}/demux.qza"
fi
echo "  ✅ demux.qza 生成成功"

# ── 3. 质量可视化 ─────────────────────────────────────────────────
echo ""
echo "【3】生成质量报告"
qiime demux summarize \
    --i-data "${OUTDIR}/demux.qza" \
    --o-visualization "${OUTDIR}/demux-summary.qzv"
echo "  ✅ demux-summary.qzv 生成成功"

# ── 4. DADA2 去噪 ──────────────────────────────────────────────────
echo ""
echo "【4】DADA2 去噪"
if [[ "${SEQ_TYPE}" == "paired" ]]; then
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs "${OUTDIR}/demux.qza" \
        --p-trunc-len-f 230 \
        --p-trunc-len-r 200 \
        --p-trim-left-f 0 \
        --p-trim-left-r 0 \
        --p-max-ee-f 2.0 \
        --p-max-ee-r 2.0 \
        --p-n-threads ${THREADS} \
        --o-table "${OUTDIR}/table.qza" \
        --o-representative-sequences "${OUTDIR}/rep-seqs.qza" \
        --o-denoising-stats "${OUTDIR}/dada2-stats.qza"
else
    qiime dada2 denoise-single \
        --i-demultiplexed-seqs "${OUTDIR}/demux.qza" \
        --p-trunc-len 0 \
        --p-trim-left 0 \
        --p-max-ee 2.0 \
        --p-n-threads ${THREADS} \
        --o-table "${OUTDIR}/table.qza" \
        --o-representative-sequences "${OUTDIR}/rep-seqs.qza" \
        --o-denoising-stats "${OUTDIR}/dada2-stats.qza"
fi
echo "  ✅ DADA2 完成"

# 输出去噪统计
qiime metadata tabulate \
    --m-input-file "${OUTDIR}/dada2-stats.qza" \
    --o-visualization "${OUTDIR}/dada2-stats.qzv"

# ── 5. 分类注释（训练好的分类器存在才运行）──────────────────────────
echo ""
echo "【5】分类注释（SILVA）"

CLASSIFIER=""
# 优先找已训练的区域裁剪版分类器
for f in \
    "/home/usr/yel/exp-16S_workflow/Qiime2/database/silva/silva-138-99-classifier-515-806.qza" \
    "/home/usr/yel/exp-16S_workflow/Qiime2/database/silva/silva-138-99-classifier.qza"; do
    if [[ -f "$f" ]]; then CLASSIFIER="$f"; break; fi
done

if [[ -n "${CLASSIFIER}" ]]; then
    echo "  使用分类器：${CLASSIFIER}"
    qiime feature-classifier classify-sklearn \
        --i-classifier "${CLASSIFIER}" \
        --i-reads "${OUTDIR}/rep-seqs.qza" \
        --p-n-jobs ${THREADS} \
        --o-classification "${OUTDIR}/taxonomy.qza"
    echo "  ✅ 分类注释完成"
else
    echo "  ⚠️ 未找到已训练分类器，跳过分类步骤"
    echo "  （请先训练分类器，参考 scripts/train_silva_classifier.sh）"
fi

# ── 6. 导出结果 ────────────────────────────────────────────────────
echo ""
echo "【6】导出结果"
qiime tools export --input-path "${OUTDIR}/table.qza" \
    --output-path "${OUTDIR}/exported/table"
qiime tools export --input-path "${OUTDIR}/rep-seqs.qza" \
    --output-path "${OUTDIR}/exported/rep-seqs"
biom convert \
    -i "${OUTDIR}/exported/table/feature-table.biom" \
    -o "${OUTDIR}/exported/table/feature-table.tsv" --to-tsv

ASV_COUNT=$(grep -c "^>" "${OUTDIR}/exported/rep-seqs/dna-sequences.fasta" || true)
echo "  ASV 数量：${ASV_COUNT}"
echo "  feature table：${OUTDIR}/exported/table/feature-table.tsv"
echo "  rep seqs：${OUTDIR}/exported/rep-seqs/dna-sequences.fasta"

# ── 7. 汇总 ───────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════"
echo " ✅ 环境验证完成！"
echo "  输出目录：${OUTDIR}"
echo "  可视化文件（需 QIIME2 View 查看）："
echo "    ${OUTDIR}/demux-summary.qzv"
echo "    ${OUTDIR}/dada2-stats.qzv"
echo "════════════════════════════════════════"
