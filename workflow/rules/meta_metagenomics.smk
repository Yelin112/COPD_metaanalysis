"""
微生物侧组学 — 宏基因组处理规则
支持两种上游工具：
  picrust2  基于16S扩增子+PICRUSt2功能预测（默认，适合16S数据）
  humann3   基于鸟枪法宏基因组+HUMAnN3功能分析（适合 shotgun 数据）

在 study_config.yaml 中通过 meta_omics.upstream_tool 控制
"""

_UPSTREAM = config["meta_omics"].get("upstream_tool", "picrust2")


def _meta_ds(dataset_id, field):
    for d in config["meta_omics"]["datasets"]:
        if d["id"] == dataset_id:
            return d[field]
    raise KeyError(f"Meta dataset {dataset_id} not found")


if _UPSTREAM == "humann3":

    rule metagenomics_humann3_to_ec:
        """
        将 HUMAnN3 EC 分层表转换为样本级 EC 丰度矩阵
        输入：HUMAnN3 运行后经 humann_regroup_table 生成的 EC 表
        输出：{dataset}_processed.tsv（行=EC号，列=样本）
        """
        input:
            lambda wc: _meta_ds(wc.dataset, "feature_file")
        output:
            f"{OUT}/processed/meta/{{dataset}}_processed.tsv"
        conda:
            "../envs/python.yaml"
        script:
            "../../modules/meta/metagenomics/01b_humann3_to_ec.py"

else:  # picrust2（默认）

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
