"""
宿主组学 — Bulk RNA-seq 数据处理规则

输入类型（通过 datasets[].input_type 控制）：
  counts  原始计数矩阵（featureCounts / HTSeq / STARsolo）→ DESeq2 VST → z-score
  tpm     TPM 矩阵（Salmon / kallisto / StringTie）        → log2(x+1)  → z-score
  fpkm    FPKM 矩阵                                        → log2(x+1)  → z-score

study_config.yaml 中 RNA-seq 数据集必填字段：
  id            数据集唯一标识（如 GSE211630）
  batch         批次编号（整数）
  input_type    counts | tpm | fpkm
  counts_file   表达矩阵路径（首列=gene_id，其余列=样本）
  metadata_file 样本元数据路径（含 sample_id 和 condition 列）
"""


def _ds(dataset_id, field):
    for d in config["host_omics"]["datasets"]:
        if d["id"] == dataset_id:
            return d[field]
    raise KeyError(f"Dataset {dataset_id} not found")


rule rnaseq_normalize:
    """
    RNA-seq 标准化
    counts   → DESeq2 VST → 列 z-score
    tpm/fpkm → log2(x+1) → 列 z-score
    输出格式与微阵列标准化结果一致，下游 ComBat + MetaDE 无需修改
    """
    input:
        counts_file   = lambda wc: _ds(wc.dataset, "counts_file"),
        metadata_file = lambda wc: _ds(wc.dataset, "metadata_file"),
    output:
        f"{OUT}/processed/host/{{dataset}}_normalized.tsv"
    params:
        input_type = lambda wc: _ds(wc.dataset, "input_type"),
        dataset_id = lambda wc: wc.dataset
    conda:
        "../envs/r_rnaseq.yaml"
    script:
        "../../modules/host/rnaseq/01_normalize.R"


rule merge_host_datasets:
    """合并所有 RNA-seq 标准化矩阵，生成 all_normalized.tsv + all_metadata.tsv"""
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
