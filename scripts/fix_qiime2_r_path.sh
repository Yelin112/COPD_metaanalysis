#!/bin/bash
# QIIME2 DADA2 R 包诊断脚本
#
# 问题：系统 R 被优先使用，但 optparse 只在 conda 环境 R 包库中存在
# 解决：通过 R_LIBS_SITE 让系统 R 也能找到 conda 环境的 R 包
#       不改变 PATH，系统 R 仍然优先，RStudio 不受影响
#
# 用法：
#   conda activate qiime2-amplicon-2023
#   bash scripts/fix_qiime2_r_path.sh

set -euo pipefail

CONDA_ENV_NAME="qiime2-amplicon-2023"
CONDA_BASE="$(conda info --base 2>/dev/null)"
CONDA_R_LIB="${CONDA_BASE}/envs/${CONDA_ENV_NAME}/lib/R/library"

echo "════════════════════════════════════════"
echo " QIIME2 DADA2 R 包诊断"
echo "════════════════════════════════════════"
echo ""

echo "【当前状态】"
echo "  系统 Rscript：$(which Rscript)"
echo "  系统 R 版本：$(Rscript --version 2>&1 | head -1)"
echo ""
echo "  系统 R .libPaths()："
Rscript -e 'cat(paste(.libPaths(), collapse="\n"))' 2>/dev/null | sed 's/^/    /'
echo ""

echo "【conda 环境 R 包库】"
if [[ -d "${CONDA_R_LIB}" ]]; then
    echo "  路径：${CONDA_R_LIB}"
    if [[ -d "${CONDA_R_LIB}/optparse" ]]; then
        echo "  ✅ optparse 存在于 conda 包库"
    else
        echo "  ❌ optparse 不在 conda 包库中"
    fi
else
    echo "  ❌ 找不到 conda 环境 R 包库：${CONDA_R_LIB}"
fi
echo ""

echo "【修复方案：R_LIBS_SITE】"
echo "  原理：R_LIBS_SITE 追加额外的库搜索路径，不改变 PATH"
echo "  系统 R 仍然被使用，RStudio 完全不受影响"
echo ""
echo "  临时测试（仅当前 shell）："
echo "    export R_LIBS_SITE=\"${CONDA_R_LIB}\""
echo "    Rscript -e 'library(optparse); cat(\"OK\\n\")'"
echo ""

export R_LIBS_SITE="${CONDA_R_LIB}"
echo "  正在测试..."
if Rscript -e 'library(optparse); cat("optparse 加载成功\n")' 2>/dev/null; then
    echo ""
    echo "  ✅ 验证通过！DADA2 可以正常运行"
    echo ""
    echo "  所有 QIIME2 脚本（test/11/12）已自动设置 R_LIBS_SITE，"
    echo "  运行前无需手动操作。"
else
    echo ""
    echo "  ❌ 仍然无法加载 optparse"
    echo "  请确认 conda 环境安装正确："
    echo "    conda install -n ${CONDA_ENV_NAME} r-optparse"
fi

echo ""
echo "════════════════════════════════════════"
