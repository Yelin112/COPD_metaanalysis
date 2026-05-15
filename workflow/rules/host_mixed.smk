"""
宿主组学 — 混合输入类型处理规则（microarray + RNA-seq 同一研究）

每个数据集通过 input_type 字段指定数据类型：
  microarray  已折叠探针矩阵（matrix_file） → 探针→基因 → z-score
  counts      原始计数矩阵（counts_file）   → DESeq2 VST → z-score
  tpm / fpkm  TPM/FPKM 矩阵（counts_file）  → log2(x+1)  → z-score

study_config.yaml 中 mixed 数据集必填字段（所有类型）：
  id            数据集唯一标识
  batch         批次编号（整数）
  input_type    microarray | counts | tpm | fpkm
  metadata_file 样本元数据路径（含 sample_id 和 condition 列）

microarray 类型额外必填：
  matrix_file   series_matrix.txt 路径
  platform      平台 ID（如 GPL10558），用于查找探针映射文件
  log2_transform bool，是否需要 log2 变换

counts / tpm / fpkm 类型额外必填：
  counts_file   表达矩阵路径（首列=gene_id，其余列=样本）
"""

import os


def _ds(dataset_id, field, default=None):
    for d in config["host_omics"]["datasets"]:
        if d["id"] == dataset_id:
            return d.get(field, default) if default is not None else d[field]
    raise KeyError(f"Dataset {dataset_id} not found")


def _input_type(dataset_id):
    return _ds(dataset_id, "input_type")


# ── 步骤1（仅 microarray）：探针→基因折叠 ─────────────────────────
rule microarray_probe_to_gene:
    """将 Illumina/Affymetrix 探针矩阵折叠为基因矩阵（IQR 最大探针代表基因）"""
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


# ── 步骤2：统一标准化（所有类型均产生同格式 z-score 矩阵）──────────
rule normalize_host_dataset:
    """
    微阵列：读取 _processed.tsv → log2（可选）→ z-score
    RNA-seq counts：读取 counts_file → DESeq2 VST → z-score
    RNA-seq tpm/fpkm：读取 counts_file → log2(x+1) → z-score
    """
    input:
        data_file = lambda wc: (
            f"{OUT}/processed/host/{wc.dataset}_processed.tsv"
            if _input_type(wc.dataset) == "microarray"
            else _ds(wc.dataset, "counts_file")
        ),
        metadata_file = lambda wc: _ds(wc.dataset, "metadata_file"),
    output:
        f"{OUT}/processed/host/{{dataset}}_normalized.tsv"
    params:
        input_type     = lambda wc: _input_type(wc.dataset),
        log2_transform = lambda wc: _ds(wc.dataset, "log2_transform", False),
        dataset_id     = lambda wc: wc.dataset
    conda:
        "../envs/r_rnaseq.yaml"
    script:
        "../../modules/host/mixed/01_normalize.R"


# ── 步骤3：合并所有数据集 ─────────────────────────────────────────
rule merge_host_datasets:
    """合并全部标准化矩阵，生成 all_normalized.tsv + all_metadata.tsv"""
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
