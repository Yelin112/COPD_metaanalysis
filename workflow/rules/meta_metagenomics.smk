"""
微生物侧组学 — 宏基因组（PICRUSt2 功能预测）处理规则
"""

def _meta_ds(dataset_id, field):
    for d in config["meta_omics"]["datasets"]:
        if d["id"] == dataset_id:
            return d[field]
    raise KeyError(f"Meta dataset {dataset_id} not found")


rule metagenomics_picrust2_to_ec:
    """
    将 PICRUSt2 输出的 ASV-EC 贡献表汇总为样本级 EC 丰度矩阵
    输入：stratified ASV 贡献表 + QIIME2 taxonomy 文件
    输出：{dataset}_processed.tsv（行=EC号，列=样本）
    """
    input:
        feature_file  = lambda wc: _meta_ds(wc.dataset, "feature_file"),
        taxonomy_file = lambda wc: _meta_ds(wc.dataset, "taxonomy_file")
    output:
        f"{OUT}/processed/meta/{{dataset}}_processed.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../../modules/meta/metagenomics/01_picrust2_to_ec.py"


rule metagenomics_normalize:
    """对宏基因组 EC 丰度矩阵做相对丰度归一化和 log 变换"""
    input:
        f"{OUT}/processed/meta/{{dataset}}_processed.tsv"
    output:
        f"{OUT}/processed/meta/{{dataset}}_normalized.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../../modules/meta/metagenomics/02_normalize.py"


rule merge_meta_datasets:
    """合并所有微生物数据集，生成 all_normalized.tsv + all_metadata.tsv"""
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
        "../../modules/meta/metagenomics/03_merge_datasets.py"
