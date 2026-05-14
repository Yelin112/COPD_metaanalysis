"""
桥接规则 — 直接代谢组桥接（跳过 EMM，代谢物差异直接映射宿主基因）
当 meta_omics.type: metabolomics 时使用此规则替代 bridge_metacyc_stitch.smk
"""

_bridge_cfg = config["bridge"]["direct_metabolomics"]


rule map_metabolite_to_gene:
    """
    利用 HMDB 或自定义代谢物-基因互作数据库，
    将代谢物荟萃分析结果直接映射到宿主基因靶点
    输出格式与 generate_emm 完全相同（行=代谢物，列=宿主基因，值=互作权重）
    """
    input:
        meta_results  = f"{OUT}/meta_metaanalysis/metaanalysis_results.tsv",
        gene_db       = _bridge_cfg["compound_gene_database"]
    output:
        emm = f"{OUT}/integration/EMM.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../../modules/bridge/direct_metabolomics/01_map_metabolite_gene.py"
