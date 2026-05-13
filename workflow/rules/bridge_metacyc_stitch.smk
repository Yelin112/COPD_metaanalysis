"""
桥接规则 — MetaCyc + STITCH 数据库整合
对应原 Perl 脚本 4-7
"""

_bridge_cfg = config["bridge"]["metacyc_stitch"]
_stitch_taxid = _bridge_cfg["host_species_id"]


rule extract_metacyc:
    """步骤4：从 MetaCyc reactions.tab 中提取酶-代谢反应信息"""
    input:
        config["metacyc"]["reactions_file"]
    output:
        f"{OUT}/bridge/microbial_metabolic_reactions.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../../modules/bridge/metacyc_stitch/04_extract_metacyc.py"


rule convert_metacyc_to_stitch:
    """步骤5：将 MetaCyc 化合物 ID 转换为 STITCH ID（基于 PubChem/ChEBI）"""
    input:
        reactions      = f"{OUT}/bridge/microbial_metabolic_reactions.tsv",
        chemical_sources = config["stitch"]["chemical_sources"]
    output:
        compounds = f"{OUT}/bridge/cmpd_with_stitch_ids.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../../modules/bridge/metacyc_stitch/05_convert_to_stitch.py"


rule parse_stitch:
    """步骤6：从 STITCH 数据库提取化合物-宿主基因相互作用"""
    input:
        compounds   = f"{OUT}/bridge/cmpd_with_stitch_ids.tsv",
        prot_chem   = config["stitch"]["protein_chemical_links"].format(taxid=_stitch_taxid),
        actions     = config["stitch"]["actions"].format(taxid=_stitch_taxid)
    output:
        f"{OUT}/bridge/compound_target_match.tsv"
    params:
        confidence_threshold = _bridge_cfg["stitch_confidence_threshold"]
    conda:
        "../envs/python.yaml"
    script:
        "../../modules/bridge/metacyc_stitch/06_parse_stitch.py"


rule link_enzyme_gene:
    """步骤7：整合酶-化合物-宿主基因链接，过滤宿主源化合物"""
    input:
        reactions    = f"{OUT}/bridge/microbial_metabolic_reactions.tsv",
        compounds    = f"{OUT}/bridge/cmpd_with_stitch_ids.tsv",
        targets      = f"{OUT}/bridge/compound_target_match.tsv",
        host_meta    = f"{OUT}/host_metaanalysis/metaanalysis_results.tsv",
        human_filter = _bridge_cfg["human_compound_filter"]
    output:
        links_dir = directory(f"{OUT}/bridge/enzyme_compound_gene_links")
    conda:
        "../envs/python.yaml"
    script:
        "../../modules/bridge/metacyc_stitch/07_link_enzyme_gene.py"
