#!/bin/bash
# 修复 QIIME2 conda 环境中 R 路径问题
# 问题：QIIME2 的 DADA2 插件调用 run_dada.R 时，使用的是系统 R 而非 conda 环境 R
# 导致找不到 optparse 等 R 包
#
# 诊断：
#   conda activate qiime2-amplicon-2023
#   bash scripts/fix_qiime2_r_path.sh
#
# 此脚本会：
#   1. 诊断当前 R 路径问题
#   2. 提供修复方案

set -euo pipefail

CONDA_ENV_NAME="qiime2-amplicon-2023"
CONDA_BASE="${CONDA_BASE:-${HOME}/miniconda3}"
ENV_BIN="${CONDA_BASE}/envs/${CONDA_ENV_NAME}/bin"

echo "════════════════════════════════════════"
echo " QIIME2 R 路径诊断"
echo "════════════════════════════════════════"
echo ""

# ── 1. 当前 R 路径 ─────────────────────────────────────────────────
echo "【1】当前 PATH 中的 Rscript"
echo "  which Rscript: $(which Rscript 2>/dev/null || echo '未找到')"
echo "  Rscript 版本：$(Rscript --version 2>&1 | head -1)"
echo ""

# ── 2. conda 环境 R ───────────────────────────────────────────────
echo "【2】conda 环境 R 路径"
CONDA_RSCRIPT="${ENV_BIN}/Rscript"
if [[ -f "${CONDA_RSCRIPT}" ]]; then
    echo "  ✅ 找到：${CONDA_RSCRIPT}"
    echo "  版本：$(${CONDA_RSCRIPT} --version 2>&1 | head -1)"
    echo ""
    echo "  .libPaths()："
    ${CONDA_RSCRIPT} -e 'cat(paste(.libPaths(), collapse="\n"))' 2>/dev/null | sed 's/^/    /'
    echo ""
    echo "  optparse 检查："
    if ${CONDA_RSCRIPT} -e 'library(optparse)' 2>/dev/null; then
        echo "  ✅ optparse 可用"
    else
        echo "  ❌ optparse 不可用"
    fi
else
    echo "  ❌ 未找到 conda 环境 R：${CONDA_RSCRIPT}"
    echo "  尝试其他路径..."
    find "${CONDA_BASE}/envs/${CONDA_ENV_NAME}" -name "Rscript" 2>/dev/null | head -5 | sed 's/^/    /'
fi
echo ""

# ── 3. 问题诊断 ───────────────────────────────────────────────────
echo "【3】问题原因"
CURRENT_R=$(which Rscript 2>/dev/null || echo "")
if [[ "${CURRENT_R}" != "${CONDA_RSCRIPT}" ]]; then
    echo "  ⚠️  当前使用的 Rscript 不是 conda 环境 R"
    echo "  当前：${CURRENT_R}"
    echo "  期望：${CONDA_RSCRIPT}"
    echo ""
    echo "  可能原因："
    echo "  - ~/.bashrc 或 ~/.bash_profile 中的 PATH 设置覆盖了 conda 环境"
    echo "  - 系统 R 在 conda 环境 R 之前"

    # 检查是否真的是 conda env 被激活
    echo ""
    echo "  当前 CONDA_DEFAULT_ENV：${CONDA_DEFAULT_ENV:-未设置}"
    echo "  PATH 前几项："
    echo "${PATH}" | tr ':' '\n' | head -5 | sed 's/^/    /'
else
    echo "  ✅ Rscript 路径正确，问题可能是 R_LIBS 或 R_LIBS_USER 覆盖"
fi
echo ""

# ── 4. 修复方案 ───────────────────────────────────────────────────
echo "【4】修复方案"
echo ""
echo "  方案A（推荐）：清除 R_LIBS 环境变量后重新激活 conda 环境"
echo ""
echo "    unset R_LIBS R_LIBS_USER R_LIBS_SITE"
echo "    conda deactivate"
echo "    conda activate ${CONDA_ENV_NAME}"
echo "    which Rscript   # 应显示 conda env 路径"
echo ""
echo "  方案B：临时修正 PATH（在当前 shell 中）"
echo ""
echo "    export PATH=\"${ENV_BIN}:\${PATH}\""
echo "    which Rscript   # 验证"
echo ""
echo "  方案C：创建符号链接（持久修复）"
echo "    # 仅在 conda env 未正确设置 R 时使用"
echo "    # ln -sf ${CONDA_RSCRIPT} ${ENV_BIN}/Rscript"
echo ""

# ── 5. 自动尝试修复 ───────────────────────────────────────────────
echo "【5】尝试修复（临时设置 PATH）"
export PATH="${ENV_BIN}:${PATH}"
echo "  已将 ${ENV_BIN} 置于 PATH 最前"
echo "  当前 Rscript：$(which Rscript)"
echo ""

if [[ -f "${CONDA_RSCRIPT}" ]]; then
    echo "  重新检查 optparse："
    if ${CONDA_RSCRIPT} -e 'library(optparse); cat("optparse OK\n")' 2>/dev/null; then
        echo ""
        echo "  ✅ 修复成功！现在可以运行 QIIME2 DADA2"
        echo ""
        echo "  ─────────────────────────────────────────"
        echo "  永久修复：在 ~/.bashrc 或 ~/.bash_profile 中添加："
        echo "  export PATH=\"${ENV_BIN}:\${PATH}\""
        echo "  （注意：只在 ${CONDA_ENV_NAME} 环境中有效）"
        echo "  更好的方案：查找为什么 conda activate 没有正确设置 PATH"
        echo "  ─────────────────────────────────────────"
    else
        echo ""
        echo "  ❌ 即使使用 conda R，optparse 仍不可用"
        echo "  需要在 conda 环境中安装 optparse："
        echo "    conda install -n ${CONDA_ENV_NAME} r-optparse"
    fi
fi

echo ""
echo "════════════════════════════════════════"
