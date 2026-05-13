"""
步骤7：整合酶(EC)-化合物-宿主基因链接
  - 过滤宿主源化合物（human_compound_filter 黑名单）
  - 对每个 EC 号生成一个独立的关联结果文件（写入 links_dir/）
  - 同时输出汇总文件 all_enzyme_gene_links.tsv
对应原 Perl 脚本：7_get_enzyme_compound_gene_links.pl

Snakemake 注入变量：
  snakemake.input.reactions    — microbial_metabolic_reactions.tsv
  snakemake.input.compounds    — cmpd_with_stitch_ids.tsv
  snakemake.input.targets      — compound_target_match.tsv
  snakemake.input.host_meta    — 宿主荟萃分析结果（含 effect_size / fdr 等）
  snakemake.input.human_filter — 宿主源化合物 ID 列表（黑名单，每行一个）
  snakemake.output.links_dir   — 输出目录
"""

import os
import pandas as pd


def load_human_filter(path: str) -> set[str]:
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def main():
    reactions  = pd.read_csv(snakemake.input.reactions, sep="\t")
    compounds  = pd.read_csv(snakemake.input.compounds, sep="\t")
    targets    = pd.read_csv(snakemake.input.targets, sep="\t")
    host_meta  = pd.read_csv(snakemake.input.host_meta, sep="\t", index_col=0)
    human_cmpds = load_human_filter(snakemake.input.human_filter)

    links_dir = snakemake.output.links_dir
    os.makedirs(links_dir, exist_ok=True)

    # metacyc_id → stitch_id 映射
    cmpd_map = dict(zip(compounds["metacyc_id"], compounds["stitch_id"]))

    # stitch_id → [(protein_id, score)] 映射
    target_map = targets.groupby("stitch_id").apply(
        lambda g: list(zip(g["protein_id"], g["combined_score"]))
    ).to_dict()

    all_rows = []

    for _, rxn in reactions.iterrows():
        ec = rxn["ec"]
        all_cmpds = set()
        for col in ("left_compounds", "right_compounds"):
            if pd.notna(rxn[col]):
                all_cmpds.update(c.strip() for c in rxn[col].split(";") if c.strip())

        # 过滤宿主源化合物
        all_cmpds -= human_cmpds

        ec_rows = []
        for metacyc_id in all_cmpds:
            stitch_id = cmpd_map.get(metacyc_id)
            if not stitch_id:
                continue
            for protein_id, score in target_map.get(stitch_id, []):
                gene = protein_id.split(".")[-1] if "." in protein_id else protein_id
                if gene not in host_meta.index:
                    continue
                row = {
                    "ec":          ec,
                    "metacyc_id":  metacyc_id,
                    "stitch_id":   stitch_id,
                    "protein_id":  protein_id,
                    "gene":        gene,
                    "stitch_score": score,
                }
                row.update(host_meta.loc[gene].to_dict())
                ec_rows.append(row)

        if ec_rows:
            ec_df = pd.DataFrame(ec_rows)
            safe_ec = ec.replace(":", "_").replace("/", "_")
            ec_df.to_csv(os.path.join(links_dir, f"{safe_ec}.tsv"), sep="\t", index=False)
            all_rows.extend(ec_rows)

    summary = pd.DataFrame(all_rows)
    summary.to_csv(os.path.join(links_dir, "all_enzyme_gene_links.tsv"), sep="\t", index=False)


main()
