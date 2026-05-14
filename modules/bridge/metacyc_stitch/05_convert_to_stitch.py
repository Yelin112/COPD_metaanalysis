"""
步骤5：将 MetaCyc 化合物 ID 转换为 STITCH CID（基于 PubChem/ChEBI 匹配）
对应原 Perl 脚本：5_convert_metacyc_to_stitch.pl

Snakemake 注入变量：
  snakemake.input.reactions       — 步骤4 输出的反应表
  snakemake.input.chemical_sources — STITCH chemical.sources.v5.0.tsv
  snakemake.output.compounds      — 输出：带 STITCH ID 的化合物表
"""

import pandas as pd


def extract_all_compounds(reactions_df: pd.DataFrame) -> set[str]:
    """从反应表中提取所有唯一化合物 ID"""
    cmpds = set()
    for col in ("left_compounds", "right_compounds"):
        for val in reactions_df[col].dropna():
            cmpds.update(c.strip() for c in val.split(";") if c.strip())
    return cmpds


def main():
    reactions_df = pd.read_csv(snakemake.input.reactions, sep="\t")
    all_cmpds = extract_all_compounds(reactions_df)

    # STITCH chemical.sources 格式：
    # flat_chemical_id | source_db | source_id | ...
    sources = pd.read_csv(
        snakemake.input.chemical_sources,
        sep="\t",
        usecols=[0, 1, 2],
        names=["stitch_id", "source_db", "source_id"],
        skiprows=1,
    )

    # MetaCyc 化合物 ID 通常对应 PubChem 或 ChEBI
    pubchem = sources[sources["source_db"].isin(["PubChem", "ChEBI"])]
    mapping = pubchem[pubchem["source_id"].isin(all_cmpds)][["source_id", "stitch_id"]]
    mapping.columns = ["metacyc_id", "stitch_id"]
    mapping = mapping.drop_duplicates(subset="metacyc_id")

    mapping.to_csv(snakemake.output.compounds, sep="\t", index=False)


main()
