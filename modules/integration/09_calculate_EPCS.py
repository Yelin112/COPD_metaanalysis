"""
步骤9：计算 EPCS 分数（Environmental-to-Phenotype Contribution Score，原 PRMT）
  EPCS(compound_j) = Σ_i [ EMM(compound_j, enzyme_i) × EffectSize(enzyme_i) ]
对应原 Perl 脚本：9_calculate_PRMT.pl

Snakemake 注入变量：
  snakemake.input.emm         — EMM.tsv（行=化合物，列=EC酶）
  snakemake.input.effect_size — 微生物荟萃分析结果（含 feature_id + zval 列）
  snakemake.output.epcs       — EPCS_scores.tsv
"""

import pandas as pd


def main():
    emm = pd.read_csv(snakemake.input.emm, sep="\t", index_col=0)
    meta = pd.read_csv(snakemake.input.effect_size, sep="\t", index_col=0)

    # 取 zval 列作为效应量
    if "zval" in meta.columns:
        effect = meta["zval"]
    elif "effect_size" in meta.columns:
        effect = meta["effect_size"]
    else:
        raise ValueError("荟萃分析结果文件中未找到 'zval' 或 'effect_size' 列")

    # 对齐 EMM 列与效应量索引
    common = emm.columns.intersection(effect.index)
    emm_aligned    = emm[common]
    effect_aligned = effect.loc[common]

    # 矩阵乘法：每行（化合物）的加权求和
    epcs = emm_aligned.dot(effect_aligned)
    epcs = epcs.sort_values(ascending=False)
    epcs.name = "epcs_score"
    epcs.index.name = "compound_id"

    result = epcs.reset_index()
    result["rank"] = range(1, len(result) + 1)
    result.to_csv(snakemake.output.epcs, sep="\t", index=False)


main()
