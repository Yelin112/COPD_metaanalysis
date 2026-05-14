"""
将 HUMAnN3 输出的 EC 分层表转换为样本级 EC 丰度矩阵
（鸟枪法宏基因组替代 PICRUSt2 的处理模块）

输入文件格式（HUMAnN3 EC 分层表）：
  # Gene Family    Sample1    Sample2    ...
  UNMAPPED         100.0      200.0
  UNGROUPED        50.0       60.0
  EC:1.1.1.1       25.0       30.0       ← 汇总行（合计）
  EC:1.1.1.1|g__Streptococcus.s__mutans  15.0  20.0  ← 分层行
  EC:1.1.1.1|unclassified                10.0  10.0
  EC:1.1.1.2       8.0        5.0

处理逻辑：
  - 若存在不含 | 的 EC 汇总行，直接使用（精确度更高）
  - 若只有分层行，则按 EC 编号汇总所有分层的丰度
  - 过滤掉 UNMAPPED / UNGROUPED 行

LOGO 分析兼容性：
  HUMAnN3 的分层行（EC|taxonomy）可直接用于 LOGO 分析，
  无需像 PICRUSt2 那样额外提供 taxonomy 文件。
  LOGO 模块读取分层行时会从行名中解析属级信息。

Snakemake 注入变量：
  snakemake.input[0]   — HUMAnN3 EC 分层表（.humann3_ec.tsv）
  snakemake.output[0]  — 汇总后 EC 丰度矩阵（行=EC，列=样本）
"""

import re
import pandas as pd


def main():
    input_file  = snakemake.input[0]
    output_file = snakemake.output[0]

    rows = []
    header = None

    with open(input_file, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                # 表头行：# Gene Family\tSample1\t...
                header = line.lstrip("# ").split("\t")
                continue
            if not line:
                continue

            parts = line.split("\t")
            feature = parts[0]

            # 跳过 UNMAPPED / UNGROUPED
            if feature in ("UNMAPPED", "UNGROUPED"):
                continue

            # 只保留 EC: 开头的行（跳过 UniRef90 等其他命名空间）
            if not feature.startswith("EC:"):
                continue

            # 跳过分层行（含 |），只保留汇总行
            if "|" in feature:
                continue

            rows.append(parts)

    if not rows:
        raise ValueError(
            f"在 {input_file} 中未找到 EC 汇总行。"
            "请确认已运行 humann_regroup_table --groups uniref90_level4ec"
        )

    # 如果没有汇总行（极少数情况），从分层行重新汇总
    if not rows:
        rows_stratified = []
        with open(input_file, encoding="utf-8") as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("#") or not line:
                    continue
                parts = line.split("\t")
                feature = parts[0]
                if feature.startswith("EC:") and "|" in feature:
                    rows_stratified.append(parts)

        df = pd.DataFrame(rows_stratified)
        df.columns = ["feature"] + (header[1:] if header else list(range(1, df.shape[1])))
        df["ec"] = df["feature"].str.split("|").str[0]
        sample_cols = df.columns[1:-1]
        df[sample_cols] = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        result = df.groupby("ec")[list(sample_cols)].sum()
    else:
        df = pd.DataFrame(rows)
        n_samples = df.shape[1] - 1
        sample_names = header[1:] if header and len(header) > 1 else [f"s{i}" for i in range(n_samples)]
        df.columns = ["feature_id"] + sample_names
        df = df.set_index("feature_id")
        df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
        result = df

    result.index.name = "feature_id"
    result.to_csv(output_file, sep="\t")


main()
