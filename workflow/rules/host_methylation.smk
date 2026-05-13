"""
宿主组学 — 甲基化组数据处理规则（待实现）
输入：Illumina EPIC/450K 芯片 Beta 值矩阵
输出格式与 host_microarray.smk 完全相同
"""

def _ds(dataset_id, field):
    for d in config["host_omics"]["datasets"]:
        if d["id"] == dataset_id:
            return d[field]
    raise KeyError(f"Dataset {dataset_id} not found")


rule methylation_beta_to_m:
    """Beta 值 → M 值转换，过滤高变异 CpG 位点"""
    input:
        beta     = lambda wc: _ds(wc.dataset, "beta_file"),
        metadata = lambda wc: _ds(wc.dataset, "metadata_file")
    output:
        f"{OUT}/processed/host/{{dataset}}_normalized.tsv"
    conda:
        "../envs/r_metaanalysis.yaml"
    script:
        "../../modules/host/methylation/01_beta_to_m.R"


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
        "../../modules/host/methylation/02_merge_datasets.py"
