"""
合并所有微生物数据集的标准化矩阵，生成 all_normalized.tsv 和 all_metadata.tsv

Snakemake 注入变量：
  snakemake.input.matrices  — 所有 *_normalized.tsv 路径列表
  snakemake.output.matrix   — 合并矩阵
  snakemake.output.metadata — 样本元数据（sample_id / study / batch / condition）
  snakemake.params.datasets — YAML 中 meta_omics.datasets 列表
"""

import os
import pandas as pd


def main():
    datasets     = snakemake.params.datasets
    input_files  = snakemake.input.matrices
    out_matrix   = snakemake.output.matrix
    out_metadata = snakemake.output.metadata

    id_to_file = {}
    for f in input_files:
        ds_id = os.path.basename(f).replace("_normalized.tsv", "")
        id_to_file[ds_id] = f

    frames = []
    meta_rows = []

    for ds in datasets:
        ds_id = ds["id"]
        if ds_id not in id_to_file:
            continue
        df = pd.read_csv(id_to_file[ds_id], sep="\t", index_col=0)
        for sample_id in df.columns:
            meta_rows.append({
                "sample_id": sample_id,
                "study":     ds_id,
                "batch":     ds["batch"],
                "condition": "Unknown",   # 微生物侧 condition 由宿主元数据决定
            })
        frames.append(df)

    combined = pd.concat(frames, axis=1, join="outer").fillna(0)
    combined.index.name = "feature_id"
    combined.to_csv(out_matrix, sep="\t")

    meta = pd.DataFrame(meta_rows).set_index("sample_id")
    meta.to_csv(out_metadata, sep="\t")


main()
