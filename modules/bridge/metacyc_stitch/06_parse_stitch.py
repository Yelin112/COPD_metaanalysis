"""
步骤6：从 STITCH 数据库提取化合物→宿主蛋白靶点的相互作用
过滤条件：combined_score >= confidence_threshold
对应原 Perl 脚本：6_parse_stitch_database.pl

Snakemake 注入变量：
  snakemake.input.compounds  — 步骤5 输出的带 STITCH ID 化合物表
  snakemake.input.prot_chem  — STITCH protein_chemical.links.detailed 文件
  snakemake.input.actions    — STITCH actions 文件
  snakemake.output[0]        — 化合物→靶点匹配表
  snakemake.params.confidence_threshold — 置信度阈值（默认 700）
"""

import pandas as pd


def main():
    threshold = snakemake.params.confidence_threshold

    cmpd_df = pd.read_csv(snakemake.input.compounds, sep="\t")
    stitch_ids = set(cmpd_df["stitch_id"].dropna())

    # STITCH protein_chemical.links 格式：
    # chemical | protein | similarity | ... | combined_score
    links = pd.read_csv(snakemake.input.prot_chem, sep="\t")
    links.columns = links.columns.str.lower().str.replace(".", "_", regex=False)

    chem_col  = links.columns[0]   # chemical
    prot_col  = links.columns[1]   # protein
    score_col = "combined_score"

    links = links[links[chem_col].isin(stitch_ids)]
    links = links[links[score_col] >= threshold][[chem_col, prot_col, score_col]]
    links.columns = ["stitch_id", "protein_id", "combined_score"]

    # 关联 action 类型（activation / inhibition / binding ...）
    actions = pd.read_csv(snakemake.input.actions, sep="\t")
    actions.columns = actions.columns.str.lower()
    act_cols = [c for c in actions.columns if "item_id" in c or "action" in c or "mode" in c]
    if len(act_cols) >= 2:
        actions = actions[act_cols[:2]]
        actions.columns = ["stitch_id", "action"]
        links = links.merge(actions, on="stitch_id", how="left")

    links.to_csv(snakemake.output[0], sep="\t", index=False)


main()
