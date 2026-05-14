"""
步骤10：LOGO（Leave-One-Genus-Out）敏感性分析
逐一排除每个属，重新汇总 EC 丰度，再用 EMM 和荟萃分析效应量
重新计算 EPCS，评估该属对整体结果的贡献
对应原 Perl 脚本：10_LOGO.pl

Snakemake 注入变量：
  snakemake.input.feature_files  — 各微生物数据集的原始 EC 贡献表列表
  snakemake.input.taxonomy_files — 对应的 QIIME2 taxonomy.tsv 列表
  snakemake.input.genus_list     — 待逐一排除的属列表（每行一个）
  snakemake.input.ec_list        — 关注的 EC 基因列表（每行一个）
  snakemake.input.emm            — EMM.tsv（步骤8 输出）
  snakemake.input.meta_results   — 微生物荟萃分析结果
  snakemake.output.summary       — 汇总 EPCS 结果（所有属的对比）
  snakemake.output.outdir        — 每个属的详细 EPCS 输出目录
  snakemake.params.datasets      — 数据集 ID 列表（与 feature/taxonomy 文件对应）
"""

import os
import pandas as pd
import numpy as np


def load_taxonomy(taxonomy_file: str) -> dict[str, str]:
    """返回 {asv_id: taxonomy_string} 映射"""
    taxa = {}
    with open(taxonomy_file) as f:
        next(f)  # 跳过表头
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                taxa[parts[0]] = parts[1]
    return taxa


def load_ec_list(ec_list_file: str) -> set[str]:
    with open(ec_list_file) as f:
        return {line.strip() for line in f if line.strip()}


def load_genus_list(genus_list_file: str) -> list[str]:
    with open(genus_list_file) as f:
        return [line.strip() for line in f if line.strip()]


def compute_ec_profile_without_genus(
    feature_file: str,
    taxonomy: dict[str, str],
    ec_whitelist: set[str],
    exclude_genus: str,
) -> pd.DataFrame:
    """
    从 PICRUSt2 分层贡献表中排除指定属的 ASV，
    再对 EC 进行汇总，返回 DataFrame（行=EC，列=样本）
    """
    df = pd.read_csv(feature_file, sep="\t")
    ec_col  = df.columns[0]
    asv_col = df.columns[1]
    sample_cols = df.columns[2:].tolist()

    # 过滤属于 exclude_genus 的 ASV
    def is_excluded(asv):
        return f"g__{exclude_genus}" in taxonomy.get(asv, "")

    mask = ~df[asv_col].apply(is_excluded)
    df = df[mask]

    # 只保留关注的 EC
    if ec_whitelist:
        df = df[df[ec_col].isin(ec_whitelist)]

    agg = df.groupby(ec_col)[sample_cols].sum()
    return agg


def calculate_epcs(ec_profile: pd.DataFrame, emm: pd.DataFrame, effect: pd.Series) -> pd.Series:
    """基于给定 EC profile 的效应量（直接用全局效应量），重新计算 EPCS"""
    # 对齐
    common_ec = emm.columns.intersection(effect.index).intersection(ec_profile.index)
    emm_aligned = emm[common_ec]
    effect_aligned = effect.loc[common_ec]
    return emm_aligned.dot(effect_aligned)


def main():
    feature_files  = snakemake.input.feature_files
    taxonomy_files = snakemake.input.taxonomy_files
    datasets       = snakemake.params.datasets
    genus_list     = load_genus_list(snakemake.input.genus_list)
    ec_whitelist   = load_ec_list(snakemake.input.ec_list)

    emm  = pd.read_csv(snakemake.input.emm, sep="\t", index_col=0)
    meta = pd.read_csv(snakemake.input.meta_results, sep="\t", index_col=0)
    effect = meta["zval"] if "zval" in meta.columns else meta["effect_size"]

    outdir = snakemake.output.outdir
    os.makedirs(outdir, exist_ok=True)

    # taxonomy 合并（多数据集共享 ASV 命名空间）
    combined_taxonomy: dict[str, str] = {}
    for tf in taxonomy_files:
        combined_taxonomy.update(load_taxonomy(tf))

    summary_rows = []

    for genus in genus_list:
        # 对每个数据集，排除该属后重新汇总 EC 丰度
        profiles = []
        for ds_id, ff in zip(datasets, feature_files):
            profile = compute_ec_profile_without_genus(
                ff, combined_taxonomy, ec_whitelist, genus
            )
            profiles.append(profile)

        # 合并多数据集（取均值或求和）
        if profiles:
            merged = pd.concat(profiles, axis=1).fillna(0)
            # 取每个 EC 在所有样本的均值作为代表丰度（简化）
            ec_mean = merged.mean(axis=1)

            # 重新计算 EPCS（仅用 EMM 权重，效应量不变）
            common = emm.columns.intersection(effect.index).intersection(ec_mean.index)
            epcs = emm[common].dot(effect.loc[common])

            genus_df = epcs.sort_values(ascending=False).reset_index()
            genus_df.columns = ["compound_id", "epcs_score"]
            genus_df["excluded_genus"] = genus
            genus_df.to_csv(os.path.join(outdir, f"{genus}_EPCS.tsv"), sep="\t", index=False)

            top = genus_df.head(10)
            for _, row in top.iterrows():
                summary_rows.append({
                    "excluded_genus": genus,
                    "compound_id":    row["compound_id"],
                    "epcs_score":     row["epcs_score"],
                })

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(snakemake.output.summary, sep="\t", index=False)


main()
