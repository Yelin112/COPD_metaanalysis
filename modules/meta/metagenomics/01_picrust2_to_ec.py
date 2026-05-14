"""
将 PICRUSt2 输出的 ASV-EC 分层贡献表汇总为样本级 EC 丰度矩阵
（仅汇总，不过滤属——属过滤在 LOGO 分析中进行）
对应原 Perl 脚本 10 的数据加载部分

Snakemake 注入变量：
  snakemake.input.feature_file  — PICRUSt2 EC 贡献表（列1=EC, 列2=ASV, 其余=样本）
  snakemake.input.taxonomy_file — QIIME2 taxonomy.tsv（asv_id\ttaxonomy_string）
  snakemake.output[0]           — 汇总后 EC 丰度矩阵（行=EC，列=样本）
"""

import pandas as pd


def main():
    feature_file  = snakemake.input.feature_file
    output_file   = snakemake.output[0]

    df = pd.read_csv(feature_file, sep="\t")

    # 前两列是 EC 和 ASV，其余是样本
    ec_col  = df.columns[0]
    asv_col = df.columns[1]
    sample_cols = df.columns[2:].tolist()

    # 按 EC 汇总所有 ASV 的贡献
    agg = df.groupby(ec_col)[sample_cols].sum()
    agg.index.name = "feature_id"

    agg.to_csv(output_file, sep="\t")


main()
