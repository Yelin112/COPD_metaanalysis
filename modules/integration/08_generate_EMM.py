"""
步骤8：生成归一化酶-化合物贡献矩阵（EMM）
  - 行 = 化合物，列 = EC 酶
  - 值 = 该酶对该化合物的净贡献方向（正=产生，负=消耗），归一化到 [-1, 1]
对应原 Perl 脚本：8_generate_EMM_for_PRMT.pl

Snakemake 注入变量：
  snakemake.input.reactions — microbial_metabolic_reactions.tsv
  snakemake.output.emm      — EMM.tsv
"""

import pandas as pd
from collections import defaultdict


def main():
    reactions = pd.read_csv(snakemake.input.reactions, sep="\t")

    # contrib[ec][compound] = 净贡献计数
    contrib: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))

    for _, row in reactions.iterrows():
        ec        = row["ec"]
        left      = [c.strip() for c in str(row.get("left_compounds", "")).split(";") if c.strip()]
        right     = [c.strip() for c in str(row.get("right_compounds", "")).split(";") if c.strip()]
        direction = str(row.get("direction", "")).upper()

        if "RIGHT-TO-LEFT" in direction:
            # 反向：left 为产物，right 为底物
            for c in left:
                contrib[ec][c] += 1
            for c in right:
                contrib[ec][c] -= 1
        else:
            # 正向或未知：left 为底物，right 为产物
            for c in left:
                contrib[ec][c] -= 1
            for c in right:
                contrib[ec][c] += 1

    # 收集所有化合物和酶
    all_ecs   = sorted(contrib.keys())
    all_cmpds = sorted({c for ec_dict in contrib.values() for c in ec_dict})

    # 构建原始贡献矩阵
    emm = pd.DataFrame(0, index=all_cmpds, columns=all_ecs, dtype=float)
    for ec, cmpd_dict in contrib.items():
        for cmpd, val in cmpd_dict.items():
            emm.loc[cmpd, ec] = val

    # 按行（化合物）归一化：正贡献之和归一，负贡献之和归一
    def normalize_row(row):
        pos = row[row > 0].sum()
        neg = row[row < 0].abs().sum()
        new = row.copy()
        if pos > 0:
            new[new > 0] = new[new > 0] / pos
        if neg > 0:
            new[new < 0] = new[new < 0] / neg
        return new

    emm = emm.apply(normalize_row, axis=1)
    emm.index.name = "compound_id"
    emm.to_csv(snakemake.output.emm, sep="\t")


main()
