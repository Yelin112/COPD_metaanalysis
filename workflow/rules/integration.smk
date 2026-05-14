"""
整合规则 — EMM 生成和 EPCS（原 PRMT）计算
对应原 Perl 脚本 8-9
"""


rule generate_emm:
    """
    步骤8：生成归一化酶-化合物贡献矩阵（EMM）
    每个化合物对每个酶的正/负贡献按反应方向归一化到 [-1, 1]
    """
    input:
        reactions = f"{OUT}/bridge/microbial_metabolic_reactions.tsv"
    output:
        emm = f"{OUT}/integration/EMM.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../../modules/integration/08_generate_EMM.py"


rule calculate_epcs:
    """
    步骤9：计算 EPCS 分数（原 PRMT）
    EPCS(compound) = Σ EMM(compound, enzyme) × EffectSize(enzyme)
    """
    input:
        emm        = f"{OUT}/integration/EMM.tsv",
        effect_size = f"{OUT}/meta_metaanalysis/metaanalysis_results.tsv"
    output:
        epcs = f"{OUT}/integration/EPCS_scores.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../../modules/integration/09_calculate_EPCS.py"
