"""
步骤1：微阵列探针级别表达矩阵 → 基因级别矩阵
对于一个基因对应多个探针的情况，选取 IQR 最大的探针代表该基因
对应原 Perl 脚本：1_generate_transcriptome_table.pl

Snakemake 注入变量：
  snakemake.input.matrix    — 原始序列矩阵文件（GEO series_matrix.txt）
  snakemake.input.metadata  — 样本元数据文件（sample_id + condition 等）
  snakemake.input.mapping   — 探针-基因映射文件（probe_id\tgene_symbol）
  snakemake.output[0]       — 输出：基因级别表达矩阵 TSV
  snakemake.params.dataset_id
"""

import pandas as pd
import numpy as np


def load_probe_gene_mapping(mapping_file: str) -> tuple[dict, dict]:
    """读取探针→基因和基因→探针映射，跳过多基因映射行（含 ///）"""
    probe2gene = {}
    gene2probes = {}
    with open(mapping_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if "///" in line:
                continue
            parts = line.split("\t")
            if len(parts) < 2 or not parts[1].strip():
                continue
            probe_id, gene = parts[0], parts[1].strip()
            probe2gene[probe_id] = gene
            gene2probes.setdefault(gene, []).append(probe_id)
    return probe2gene, gene2probes


def load_series_matrix(matrix_file: str, probe2gene: dict) -> pd.DataFrame:
    """
    解析 GEO series_matrix.txt
    跳过 '!' 开头的元数据行，返回 DataFrame（index=probe_id，columns=sample_id）
    """
    rows = []
    headers = None
    with open(matrix_file) as f:
        for line in f:
            line = line.rstrip("\n").replace('"', "")
            if line.startswith("!"):
                continue
            if line.startswith("ID_REF"):
                headers = line.split("\t")
                continue
            if headers is None:
                continue
            parts = line.split("\t")
            probe_id = parts[0]
            if probe_id not in probe2gene:
                continue
            values = parts[1:]
            rows.append([probe_id] + values)

    df = pd.DataFrame(rows, columns=headers)
    df = df.set_index("ID_REF")
    df = df.apply(pd.to_numeric, errors="coerce")
    return df


def select_best_probe(gene2probes: dict, expr_df: pd.DataFrame) -> dict:
    """对每个基因，选取非缺失样本中 IQR 最大的探针"""
    best = {}
    for gene, probes in gene2probes.items():
        valid = [p for p in probes if p in expr_df.index]
        if not valid:
            continue
        best_probe = max(
            valid,
            key=lambda p: float(np.nanpercentile(expr_df.loc[p].dropna(), 75)
                                - np.nanpercentile(expr_df.loc[p].dropna(), 25))
        )
        best[gene] = best_probe
    return best


def main():
    matrix_file  = snakemake.input.matrix
    mapping_file = snakemake.input.mapping
    output_file  = snakemake.output[0]

    probe2gene, gene2probes = load_probe_gene_mapping(mapping_file)
    expr_df = load_series_matrix(matrix_file, probe2gene)

    best_probe = select_best_probe(gene2probes, expr_df)

    selected_probes = list(best_probe.values())
    gene_names      = list(best_probe.keys())

    result = expr_df.loc[selected_probes].copy()
    result.index = gene_names
    result.index.name = "gene_id"

    result.to_csv(output_file, sep="\t")


main()
