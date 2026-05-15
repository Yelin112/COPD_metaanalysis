"""
RNA-seq 数据集合并模块
合并逻辑与微阵列完全相同（均为 gene × sample z-score 矩阵）
Snakemake 规则 host_rnaseq.smk 直接调用 microarray/03_merge_datasets.py，
本文件不被流程调用，仅供直接执行或测试时使用。
"""

import os
import sys

# 将调用转发给微阵列合并模块（逻辑完全相同）
_this_dir = os.path.dirname(os.path.abspath(__file__))
_merge_module = os.path.join(_this_dir, "..", "microarray", "03_merge_datasets.py")

with open(_merge_module) as _f:
    exec(compile(_f.read(), _merge_module, "exec"))
