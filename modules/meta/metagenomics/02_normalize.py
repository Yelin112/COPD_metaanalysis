"""
宏基因组 EC 丰度矩阵标准化：相对丰度 + log2(x+1) 变换

Snakemake 注入变量：
  snakemake.input[0]  — EC 丰度矩阵（行=EC，列=样本）
  snakemake.output[0] — 标准化后矩阵
"""

import numpy as np
import pandas as pd


def main():
    df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)

    # 按列（样本）转为相对丰度
    df = df.div(df.sum(axis=0), axis=1)

    # log2(x+1) 变换
    df = np.log2(df + 1)

    df.to_csv(snakemake.output[0], sep="\t")


main()
