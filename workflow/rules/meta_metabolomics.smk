"""
微生物/环境侧 — 直接测量代谢组数据处理规则（待实现）
输入：代谢物强度矩阵（LC-MS/GC-MS）
当使用此规则时，bridge 应配置为 direct_metabolomics
"""

def _meta_ds(dataset_id, field):
    for d in config["meta_omics"]["datasets"]:
        if d["id"] == dataset_id:
            return d[field]
    raise KeyError(f"Meta dataset {dataset_id} not found")


rule metabolomics_normalize:
    """log 变换 + median scaling"""
    input:
        intensity = lambda wc: _meta_ds(wc.dataset, "intensity_file"),
        metadata  = lambda wc: _meta_ds(wc.dataset, "metadata_file")
    output:
        f"{OUT}/processed/meta/{{dataset}}_normalized.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../../modules/meta/metabolomics/01_normalize.py"


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
        "../../modules/meta/metabolomics/02_merge_datasets.py"
