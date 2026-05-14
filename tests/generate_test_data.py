"""
生成完整的模拟测试数据，覆盖流程所有输入文件
运行：python3 tests/generate_test_data.py
"""

import os
import numpy as np
import pandas as pd

np.random.seed(42)
BASE = os.path.join(os.path.dirname(__file__), "data")

# ── 基本参数 ──────────────────────────────────────────────────────
N_PROBES   = 80
N_GENES    = 50
N_EC       = 10
N_ASVS     = 15
N_CMPDS    = 8

DATASETS_HOST = [
    {"id": "GSE_TEST1", "n_copd": 6, "n_ctrl": 6, "platform": "GPL_TEST"},
    {"id": "GSE_TEST2", "n_copd": 5, "n_ctrl": 5, "platform": "GPL_TEST"},
]
DATASETS_META = [
    {"id": "SRP_TEST1", "n_samples": 12},
]

GENES  = [f"GENE{i:03d}" for i in range(1, N_GENES + 1)]
PROBES = [f"probe_{i:04d}" for i in range(1, N_PROBES + 1)]
ECS    = [f"{a}.{b}.{c}.{d}" for a, b, c, d in [
    (1,1,1,1),(1,2,1,3),(2,3,1,5),(3,1,3,2),(4,2,1,1),
    (1,3,1,27),(2,1,1,1),(3,2,1,1),(4,1,1,1),(5,1,1,1),
]]
ASVS   = [f"ASV{i:03d}" for i in range(1, N_ASVS + 1)]
CMPDS  = [f"CPD{i:04d}" for i in range(1, N_CMPDS + 1)]
GENERA = ["Streptococcus", "Haemophilus", "Pseudomonas"]

# ── 探针-基因映射 ─────────────────────────────────────────────────
def make_probe_mapping():
    rows = []
    # 每个基因分配 1-2 个探针
    probe_idx = 0
    for gene in GENES:
        n = np.random.choice([1, 2])
        for _ in range(n):
            if probe_idx >= len(PROBES):
                break
            rows.append(f"{PROBES[probe_idx]}\t{gene}")
            probe_idx += 1
    path = os.path.join(BASE, "platform_mappings", "GPL_TEST_mapping.txt")
    with open(path, "w") as f:
        f.write("probe_id\tgene_symbol\n")
        f.write("\n".join(rows) + "\n")
    print(f"  ✓ {path}")


# ── 宿主微阵列数据 ────────────────────────────────────────────────
def make_host_dataset(ds_id, n_copd, n_ctrl):
    samples_copd = [f"{ds_id}_COPD{i:02d}" for i in range(1, n_copd + 1)]
    samples_ctrl = [f"{ds_id}_CTRL{i:02d}" for i in range(1, n_ctrl + 1)]
    all_samples  = samples_copd + samples_ctrl

    # 元数据
    meta = pd.DataFrame({
        "sample_id": all_samples,
        "condition": ["COPD"] * n_copd + ["Control"] * n_ctrl,
        "age":       np.random.randint(50, 80, len(all_samples)),
        "sex":       np.random.choice(["M", "F"], len(all_samples)),
    })
    meta_path = os.path.join(BASE, "host", f"{ds_id}_metadata.txt")
    meta.to_csv(meta_path, sep="\t", index=False)
    print(f"  ✓ {meta_path}")

    # series_matrix：COPD 组中特定基因表达偏高
    signal = np.random.normal(7, 1.5, (N_PROBES, len(all_samples)))
    for i in range(10):   # 前10个探针在 COPD 中偏高
        signal[i, :n_copd] += np.random.uniform(1.5, 3.0, n_copd)

    matrix_path = os.path.join(BASE, "host", f"{ds_id}_series_matrix.txt")
    with open(matrix_path, "w") as f:
        f.write(f'!Series_title\t"{ds_id} test data"\n')
        f.write("!series_matrix_table_begin\n")
        header = "ID_REF\t" + "\t".join(all_samples)
        f.write(header + "\n")
        for j, probe in enumerate(PROBES):
            vals = "\t".join(f"{v:.4f}" for v in signal[j])
            f.write(f"{probe}\t{vals}\n")
        f.write("!series_matrix_table_end\n")
    print(f"  ✓ {matrix_path}")


# ── 微生物宏基因组数据 ────────────────────────────────────────────
def make_meta_dataset(ds_id, n_samples):
    samples = [f"{ds_id}_S{i:02d}" for i in range(1, n_samples + 1)]

    # taxonomy：每个 ASV 分配一个属
    tax_rows = []
    asv_genus = {}
    for asv in ASVS:
        genus = np.random.choice(GENERA)
        asv_genus[asv] = genus
        tax_rows.append({
            "Feature ID": asv,
            "Taxon": f"d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;"
                     f"f__Streptococcaceae;g__{genus};s__sp"
        })
    tax_df = pd.DataFrame(tax_rows)
    tax_path = os.path.join(BASE, "meta", f"{ds_id}.taxonomy.tsv")
    tax_df.to_csv(tax_path, sep="\t", index=False)
    print(f"  ✓ {tax_path}")

    # PICRUSt2 stratified EC 贡献表
    rows = []
    for ec in ECS:
        for asv in ASVS:
            contrib = np.random.exponential(0.01, n_samples)
            contrib = np.round(contrib, 6)
            rows.append([ec, asv] + list(contrib))
    feat_df = pd.DataFrame(rows, columns=["function", "taxon"] + samples)
    feat_path = os.path.join(BASE, "meta", f"{ds_id}.tsv")
    feat_df.to_csv(feat_path, sep="\t", index=False)
    print(f"  ✓ {feat_path}")


# ── MetaCyc reactions.tab（模拟块格式）────────────────────────────
def make_metacyc_reactions():
    lines = []
    for i, ec in enumerate(ECS):
        # 每个酶参与一个反应
        left_cmpds  = [CMPDS[i % N_CMPDS], CMPDS[(i + 1) % N_CMPDS]]
        right_cmpds = [CMPDS[(i + 2) % N_CMPDS]]
        direction   = np.random.choice(["LEFT-TO-RIGHT", "RIGHT-TO-LEFT", ""])
        lines += [
            f"UNIQUE-ID - RXN-TEST{i:04d}",
            f"EC-NUMBER - {ec}",
        ]
        for c in left_cmpds:
            lines.append(f"LEFT - {c}")
        for c in right_cmpds:
            lines.append(f"RIGHT - {c}")
        if direction:
            lines.append(f"REACTION-DIRECTION - {direction}")
        lines.append("//")

    path = os.path.join(BASE, "bridge", "reactions.tab")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"  ✓ {path}")


# ── STITCH chemical.sources（化合物 ID 映射）──────────────────────
def make_stitch_sources():
    rows = []
    for cmpd in CMPDS:
        stitch_id = f"CIDm{abs(hash(cmpd)) % 100000:08d}"
        rows.append({"stitch_id": stitch_id, "source_db": "PubChem", "source_id": cmpd})
    df = pd.DataFrame(rows)
    path = os.path.join(BASE, "bridge", "chemical.sources.v5.0.tsv")
    df.to_csv(path, sep="\t", index=False)
    print(f"  ✓ {path}")
    return {r["source_id"]: r["stitch_id"] for _, r in df.iterrows()}


# ── STITCH protein_chemical.links ─────────────────────────────────
def make_stitch_links(cmpd_to_stitch):
    proteins = [f"9606.ENSP{i:011d}" for i in range(1, N_GENES + 1)]
    rows = []
    for cmpd_id, stitch_id in cmpd_to_stitch.items():
        for gene, prot in zip(GENES[:10], proteins[:10]):
            score = np.random.randint(400, 1000)
            rows.append({
                "chemical":       stitch_id,
                "protein":        prot,
                "combined_score": score,
            })
    df = pd.DataFrame(rows)
    path = os.path.join(BASE, "bridge", "9606.protein_chemical.links.detailed.v5.0.tsv")
    df.to_csv(path, sep="\t", index=False)
    print(f"  ✓ {path}")

    # actions 文件（简化版）
    act_rows = [{"item_id_a": r["chemical"], "action": "binding"}
                for _, r in df.drop_duplicates("chemical").iterrows()]
    act_df = pd.DataFrame(act_rows)
    act_path = os.path.join(BASE, "bridge", "9606.actions.v5.0.tsv")
    act_df.to_csv(act_path, sep="\t", index=False)
    print(f"  ✓ {act_path}")


# ── 宿主源化合物黑名单 ────────────────────────────────────────────
def make_human_cmpd_list():
    path = os.path.join(BASE, "bridge", "human_cmpd_list.txt")
    with open(path, "w") as f:
        # 把最后两个化合物标记为宿主源，测试过滤逻辑
        f.write("\n".join(CMPDS[-2:]) + "\n")
    print(f"  ✓ {path}")


# ── LOGO 分析辅助文件 ─────────────────────────────────────────────
def make_sensitivity_files():
    genus_path = os.path.join(BASE, "sensitivity", "genus_list.txt")
    with open(genus_path, "w") as f:
        f.write("\n".join(GENERA) + "\n")
    print(f"  ✓ {genus_path}")

    ec_path = os.path.join(BASE, "sensitivity", "ec_list.txt")
    with open(ec_path, "w") as f:
        f.write("\n".join(ECS) + "\n")
    print(f"  ✓ {ec_path}")


# ── 模拟荟萃分析结果（替代 R 步骤输出，用于测试下游 Python 步骤）──
def make_mock_meta_results():
    """
    跳过真实的 ComBat/MetaDE R 步骤，
    直接生成格式正确的荟萃分析结果供下游测试使用
    """
    os.makedirs(os.path.join(BASE, "mock_results"), exist_ok=True)

    # 宿主荟萃分析结果（行=基因）
    host_meta = pd.DataFrame({
        "feature_id": GENES,
        "zval":       np.random.normal(0, 2, N_GENES),
        "pval":       np.random.uniform(0, 1, N_GENES),
        "fdr":        np.random.uniform(0, 1, N_GENES),
    })
    host_path = os.path.join(BASE, "mock_results", "host_metaanalysis_results.tsv")
    host_meta.to_csv(host_path, sep="\t", index=False)
    print(f"  ✓ {host_path}")

    # 微生物荟萃分析结果（行=EC号）
    meta_meta = pd.DataFrame({
        "feature_id": ECS,
        "zval":       np.random.normal(0, 1.5, N_EC),
        "pval":       np.random.uniform(0, 1, N_EC),
        "fdr":        np.random.uniform(0, 1, N_EC),
    })
    meta_path = os.path.join(BASE, "mock_results", "meta_metaanalysis_results.tsv")
    meta_meta.to_csv(meta_path, sep="\t", index=False)
    print(f"  ✓ {meta_path}")


if __name__ == "__main__":
    print("生成模拟测试数据...\n")

    print("[ 探针映射 ]")
    make_probe_mapping()

    print("\n[ 宿主微阵列数据 ]")
    for ds in DATASETS_HOST:
        make_host_dataset(ds["id"], ds["n_copd"], ds["n_ctrl"])

    print("\n[ 微生物宏基因组数据 ]")
    for ds in DATASETS_META:
        make_meta_dataset(ds["id"], ds["n_samples"])

    print("\n[ 桥接数据库 ]")
    make_metacyc_reactions()
    cmpd_map = make_stitch_sources()
    make_stitch_links(cmpd_map)
    make_human_cmpd_list()

    print("\n[ 敏感性分析辅助文件 ]")
    make_sensitivity_files()

    print("\n[ 模拟荟萃分析结果（跳过R步骤用）]")
    make_mock_meta_results()

    print("\n✅ 全部测试数据生成完成")
