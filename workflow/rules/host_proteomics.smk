"""
宿主组学 — 蛋白质组数据处理规则（待实现）
输入：MaxQuant LFQ 强度矩阵 / DIA-NN 输出
输出格式与 host_microarray.smk 完全相同
"""

def _ds(dataset_id, field):
    for d in config["host_omics"]["datasets"]:
        if d["id"] == dataset_id:
            return d[field]
    raise KeyError(f"Dataset {dataset_id} not found")


rule proteomics_normalize:
    """Median centering + log2 变换"""
    input:
        intensity = lambda wc: _ds(wc.dataset, "intensity_file"),
        metadata  = lambda wc: _ds(wc.dataset, "metadata_file")
    output:
        f"{OUT}/processed/host/{{dataset}}_normalized.tsv"
    conda:
        "../envs/r_metaanalysis.yaml"
    script:
        "../../modules/host/proteomics/01_normalize.R"


rule merge_host_datasets:
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
        "../../modules/host/proteomics/02_merge_datasets.py"
