"""
宿主组学 — RNA-seq 数据处理规则（待实现）
输入：featureCounts / HTSeq 输出的 counts 矩阵
输出格式与 host_microarray.smk 完全相同，下游规则无感知
"""

def _ds(dataset_id, field):
    for d in config["host_omics"]["datasets"]:
        if d["id"] == dataset_id:
            return d[field]
    raise KeyError(f"Dataset {dataset_id} not found")


rule rnaseq_normalize:
    """TMM/DESeq2 VST 标准化，输出与微阵列相同的 z-score 特征矩阵"""
    input:
        counts   = lambda wc: _ds(wc.dataset, "counts_file"),
        metadata = lambda wc: _ds(wc.dataset, "metadata_file")
    output:
        f"{OUT}/processed/host/{{dataset}}_normalized.tsv"
    conda:
        "../envs/r_metaanalysis.yaml"
    script:
        "../../modules/host/rnaseq/01_normalize.R"


rule merge_host_datasets:
    """合并所有 RNA-seq 数据集"""
    input:
        matrices = expand(
            f"{OUT}/processed/host/{{dataset}}_normalized.tsv",
            dataset=DATASETS_HOST
        )
    output:
        matrix   = f"{OUT}/processed/host/all_normalized.tsv",
        metadata = f"{OUT}/processed/host/all_metadata.tsv"
    params:
        datasets = config["host_omics"]["datasets"],
        study    = STUDY
    conda:
        "../envs/python.yaml"
    script:
        "../../modules/host/rnaseq/02_merge_datasets.py"
