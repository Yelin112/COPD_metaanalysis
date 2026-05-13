"""
步骤3：合并所有数据集的标准化矩阵，生成 all_normalized.tsv 和 all_metadata.tsv
all_metadata.tsv 包含 sample_id / study / batch / condition 四列

Snakemake 注入变量：
  snakemake.input.matrices — 所有 *_normalized.tsv 文件路径列表
  snakemake.output.matrix  — 合并后矩阵
  snakemake.output.metadata — 样本元数据
  snakemake.params.datasets — 来自 YAML 的数据集配置列表
"""

import os
import pandas as pd


def main():
    datasets     = snakemake.params.datasets
    input_files  = snakemake.input.matrices
    out_matrix   = snakemake.output.matrix
    out_metadata = snakemake.output.metadata

    # 建立 dataset_id → 文件路径的映射
    id_to_file = {}
    for f in input_files:
        basename = os.path.basename(f)
        dataset_id = basename.replace("_normalized.tsv", "")
        id_to_file[dataset_id] = f

    frames = []
    meta_rows = []

    for ds in datasets:
        ds_id = ds["id"]
        if ds_id not in id_to_file:
            continue
        df = pd.read_csv(id_to_file[ds_id], sep="\t", index_col=0)

        # 读取对应的元数据文件，获取每个样本的 condition
        metadata_path = ds["metadata_file"]
        meta_df = pd.read_csv(metadata_path, sep="\t", index_col=0)

        for sample_id in df.columns:
            condition = meta_df.loc[sample_id, "condition"] if sample_id in meta_df.index else "Unknown"
            meta_rows.append({
                "sample_id": sample_id,
                "study":     ds_id,
                "batch":     ds["batch"],
                "condition": condition,
            })

        frames.append(df)

    combined = pd.concat(frames, axis=1, join="outer")
    combined.index.name = "gene_id"
    combined.to_csv(out_matrix, sep="\t")

    meta = pd.DataFrame(meta_rows).set_index("sample_id")
    meta.to_csv(out_metadata, sep="\t")


main()
