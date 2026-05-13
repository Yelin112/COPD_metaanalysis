"""
微生物侧组学 — 宏转录组数据处理规则（待实现）
输入：PICRUSt2 宏转录组 EC 功能预测表（基于 mRNA reads）
输出格式与 meta_metagenomics.smk 完全相同
"""

def _meta_ds(dataset_id, field):
    for d in config["meta_omics"]["datasets"]:
        if d["id"] == dataset_id:
            return d[field]
    raise KeyError(f"Meta dataset {dataset_id} not found")


rule metatranscriptomics_normalize:
    input:
        feature_file = lambda wc: _meta_ds(wc.dataset, "feature_file"),
        metadata     = lambda wc: _meta_ds(wc.dataset, "metadata_file")
    output:
        f"{OUT}/processed/meta/{{dataset}}_normalized.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../../modules/meta/metatranscriptomics/01_normalize.py"


rule merge_meta_datasets:
    input:
        matrices = expand(
            f"{OUT}/processed/meta/{{dataset}}_normalized.tsv",
            dataset=DATASETS_META
        )
    output:
        matrix   = f"{OUT}/processed/meta/all_normalized.tsv",
        metadata = f"{OUT}/processed/meta/all_metadata.tsv"
    params:
        datasets = config["meta_omics"]["datasets"],
        study    = STUDY
    conda:
        "../envs/python.yaml"
    script:
        "../../modules/meta/metatranscriptomics/02_merge_datasets.py"
