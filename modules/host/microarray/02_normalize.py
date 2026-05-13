"""
步骤2：log2 变换（可选）+ 跨样本 z-score 标准化
对应原 Perl 脚本：2_transcriptome_normalization.pl

Snakemake 注入变量：
  snakemake.input[0]            — 上一步输出的基因级别矩阵
  snakemake.output[0]           — 标准化后矩阵
  snakemake.params.log2_transform — bool，是否需要 log2 变换
"""

import numpy as np
import pandas as pd


def main():
    input_file     = snakemake.input[0]
    output_file    = snakemake.output[0]
    log2_transform = snakemake.params.log2_transform

    df = pd.read_csv(input_file, sep="\t", index_col=0)

    if log2_transform:
        df = np.log2(df + 1)

    # 对每个样本（列）做 z-score 标准化
    df = (df - df.mean()) / df.std()

    df.to_csv(output_file, sep="\t")


main()
