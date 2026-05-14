"""
荟萃分析规则 — 宿主侧和微生物侧通用
依赖上游各 host_*.smk / meta_*.smk 生成标准格式的
all_normalized.tsv + all_metadata.tsv
"""


# ── 宿主侧荟萃分析 ────────────────────────────────────────────────

rule host_impute:
    """缺失值填补（knn）"""
    input:
        f"{OUT}/processed/host/all_normalized.tsv"
    output:
        f"{OUT}/processed/host/all_imputed.tsv"
    params:
        method = config["host_omics"]["metaanalysis"]["imputation"]["method"],
        k      = config["host_omics"]["metaanalysis"]["imputation"]["k"]
    conda:
        "../envs/r_metaanalysis.yaml"
    script:
        "../../modules/metaanalysis/impute.R"


rule host_combat:
    """ComBat 批次效应校正"""
    input:
        matrix   = f"{OUT}/processed/host/all_imputed.tsv",
        metadata = f"{OUT}/processed/host/all_metadata.tsv"
    output:
        f"{OUT}/processed/host/all_combat.tsv"
    params:
        batch_column = "batch",
        run_combat   = config["host_omics"]["metaanalysis"]["batch_correction"]
    conda:
        "../envs/r_metaanalysis.yaml"
    script:
        "../../modules/metaanalysis/combat_correction.R"


rule host_metaDE:
    """随机效应荟萃分析（MetaDE）"""
    input:
        matrix   = f"{OUT}/processed/host/all_combat.tsv",
        metadata = f"{OUT}/processed/host/all_metadata.tsv"
    output:
        fdr      = f"{OUT}/host_metaanalysis/metaFDR.tsv",
        zval     = f"{OUT}/host_metaanalysis/metaZval.tsv",
        pval     = f"{OUT}/host_metaanalysis/metaPval.tsv",
        combined = f"{OUT}/host_metaanalysis/metaanalysis_results.tsv"
    params:
        method           = config["host_omics"]["metaanalysis"]["method"],
        individual_method = config["host_omics"]["metaanalysis"]["individual_method"],
        nperm            = config["host_omics"]["metaanalysis"]["nperm"],
        study_column     = "study"
    conda:
        "../envs/r_metaanalysis.yaml"
    script:
        "../../modules/metaanalysis/metaDE_analysis.R"


# ── 微生物侧荟萃分析（不做 ComBat，直接 MetaDE）──────────────────

rule meta_metaDE:
    """微生物侧随机效应荟萃分析"""
    input:
        matrix   = f"{OUT}/processed/meta/all_normalized.tsv",
        metadata = f"{OUT}/processed/meta/all_metadata.tsv"
    output:
        fdr      = f"{OUT}/meta_metaanalysis/metaFDR.tsv",
        zval     = f"{OUT}/meta_metaanalysis/metaZval.tsv",
        pval     = f"{OUT}/meta_metaanalysis/metaPval.tsv",
        combined = f"{OUT}/meta_metaanalysis/metaanalysis_results.tsv"
    params:
        method            = config["meta_omics"]["metaanalysis"]["method"],
        individual_method = config["meta_omics"]["metaanalysis"]["individual_method"],
        nperm             = config["meta_omics"]["metaanalysis"]["nperm"],
        study_column      = "study"
    conda:
        "../envs/r_metaanalysis.yaml"
    script:
        "../../modules/metaanalysis/metaDE_analysis.R"
