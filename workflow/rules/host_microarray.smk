"""
宿主组学 — 微阵列数据处理规则
对应原 Perl 脚本 1（探针→基因）和 2（标准化）
"""

# 辅助函数：从配置中查找某个数据集的字段值
def _ds(dataset_id, field):
    for d in config["host_omics"]["datasets"]:
        if d["id"] == dataset_id:
            return d[field]
    raise KeyError(f"Dataset {dataset_id} not found")


rule microarray_probe_to_gene:
    """步骤1：将探针级别表达矩阵折叠为基因级别（取 IQR 最大探针代表基因）"""
    input:
        matrix   = lambda wc: _ds(wc.dataset, "matrix_file"),
        metadata = lambda wc: _ds(wc.dataset, "metadata_file"),
        mapping  = lambda wc: os.path.join(
            config["probe_mappings"]["base_dir"],
            _ds(wc.dataset, "platform") + "_mapping.txt"
        )
    output:
        f"{OUT}/processed/host/{{dataset}}_processed.tsv"
    params:
        dataset_id = lambda wc: wc.dataset
    conda:
        "../envs/python.yaml"
    script:
        "../../modules/host/microarray/01_probe_to_gene.py"


rule microarray_normalize:
    """步骤2：log2 变换（可选）+ z-score 标准化"""
    input:
        f"{OUT}/processed/host/{{dataset}}_processed.tsv"
    output:
        f"{OUT}/processed/host/{{dataset}}_normalized.tsv"
    params:
        log2_transform = lambda wc: _ds(wc.dataset, "log2_transform")
    conda:
        "../envs/python.yaml"
    script:
        "../../modules/host/microarray/02_normalize.py"


rule merge_host_datasets:
    """步骤3：合并所有数据集的标准化矩阵，并生成统一的样本元数据"""
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
        "../../modules/host/microarray/03_merge_datasets.py"
