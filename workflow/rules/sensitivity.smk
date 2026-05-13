"""
敏感性分析规则 — LOGO（Leave-One-Genus-Out）分析
对应原 Perl 脚本 10
"""

_logo_cfg  = config["sensitivity"]["logo_analysis"]
_logo_dsets = _logo_cfg["datasets"]


rule logo_analysis:
    """
    步骤10：逐一排除每个属，重算宏基因组 EC 功能谱
    输出：每个属对应一个 {study}_{genus}_logo.tsv
    再将所有结果合并为 LOGO_summary.tsv
    """
    input:
        feature_files  = [
            next(d["feature_file"]  for d in config["meta_omics"]["datasets"] if d["id"] == ds)
            for ds in _logo_dsets
        ],
        taxonomy_files = [
            next(d["taxonomy_file"] for d in config["meta_omics"]["datasets"] if d["id"] == ds)
            for ds in _logo_dsets
        ],
        genus_list = _logo_cfg["genus_list_file"],
        ec_list    = _logo_cfg["ec_list_file"],
        emm        = f"{OUT}/integration/EMM.tsv",
        meta_results = f"{OUT}/meta_metaanalysis/metaanalysis_results.tsv"
    output:
        summary = f"{OUT}/sensitivity/LOGO_summary.tsv",
        outdir  = directory(f"{OUT}/sensitivity/per_genus")
    params:
        datasets = _logo_dsets
    conda:
        "../envs/python.yaml"
    script:
        "../../modules/integration/10_logo_sensitivity.py"
