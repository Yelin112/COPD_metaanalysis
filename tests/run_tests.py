"""
端到端流程测试（跳过 R 步骤，用模拟荟萃分析结果替代）
运行：python3 tests/run_tests.py
"""

import os
import sys
import traceback
import importlib.util
import pandas as pd

# 路径设置
ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA  = os.path.join(ROOT, "tests", "data")
OUT   = os.path.join(ROOT, "tests", "outputs")
MODS  = os.path.join(ROOT, "modules")

sys.path.insert(0, os.path.join(ROOT, "tests"))
from mock_snakemake import MockSnakemake

os.makedirs(OUT, exist_ok=True)

# ── 工具函数 ──────────────────────────────────────────────────────
passed = []
failed = []

def run_step(name, script_path, snakemake_obj):
    """执行一个模块脚本，注入 mock snakemake 对象"""
    print(f"\n{'─'*55}")
    print(f"  步骤：{name}")
    print(f"  脚本：{os.path.relpath(script_path, ROOT)}")

    # 创建输出目录
    for val in vars(snakemake_obj.output)["_d"].values():
        if isinstance(val, list):
            for v in val:
                os.makedirs(os.path.dirname(v) or ".", exist_ok=True)
        elif val.endswith("/") or "." not in os.path.basename(val):
            os.makedirs(val, exist_ok=True)
        else:
            os.makedirs(os.path.dirname(val) or ".", exist_ok=True)

    try:
        with open(script_path) as f:
            code = f.read()
        ns = {"snakemake": snakemake_obj, "__file__": script_path}
        exec(compile(code, script_path, "exec"), ns)

        # 验证输出文件存在
        for key, path in vars(snakemake_obj.output)["_d"].items():
            if os.path.isdir(path):
                n = len(os.listdir(path))
                print(f"  ✓ 输出目录 [{key}]: {os.path.relpath(path, ROOT)} ({n} 个文件)")
            elif os.path.exists(path):
                df = pd.read_csv(path, sep="\t", nrows=2)
                print(f"  ✓ 输出文件 [{key}]: {os.path.relpath(path, ROOT)} "
                      f"({os.path.getsize(path)//1024+1}KB, {df.shape[1]} 列)")
            else:
                raise FileNotFoundError(f"输出文件不存在: {path}")

        passed.append(name)
    except Exception as e:
        print(f"  ✗ 失败: {e}")
        traceback.print_exc()
        failed.append(name)


# ═══════════════════════════════════════════════════════════════
print("=" * 55)
print("  多组学荟萃分析流程 — 模拟数据测试")
print("=" * 55)

# ── 步骤 1：探针 → 基因 ───────────────────────────────────────────
for ds_id in ["GSE_TEST1", "GSE_TEST2"]:
    run_step(
        f"[宿主] 探针→基因 ({ds_id})",
        os.path.join(MODS, "host/microarray/01_probe_to_gene.py"),
        MockSnakemake(
            input={
                "matrix":   os.path.join(DATA, f"host/{ds_id}_series_matrix.txt"),
                "metadata": os.path.join(DATA, f"host/{ds_id}_metadata.txt"),
                "mapping":  os.path.join(DATA, "platform_mappings/GPL_TEST_mapping.txt"),
            },
            output={"0": os.path.join(OUT, f"host/{ds_id}_processed.tsv")},
            params={"dataset_id": ds_id},
        )
    )

# ── 步骤 2：标准化 ─────────────────────────────────────────────────
for ds_id in ["GSE_TEST1", "GSE_TEST2"]:
    run_step(
        f"[宿主] 标准化 ({ds_id})",
        os.path.join(MODS, "host/microarray/02_normalize.py"),
        MockSnakemake(
            input={"0": os.path.join(OUT, f"host/{ds_id}_processed.tsv")},
            output={"0": os.path.join(OUT, f"host/{ds_id}_normalized.tsv")},
            params={"log2_transform": False},
        )
    )

# ── 步骤 3：合并数据集 ──────────────────────────────────────────────
run_step(
    "[宿主] 合并数据集",
    os.path.join(MODS, "host/microarray/03_merge_datasets.py"),
    MockSnakemake(
        input={
            "matrices": [
                os.path.join(OUT, "host/GSE_TEST1_normalized.tsv"),
                os.path.join(OUT, "host/GSE_TEST2_normalized.tsv"),
            ]
        },
        output={
            "matrix":   os.path.join(OUT, "host/all_normalized.tsv"),
            "metadata": os.path.join(OUT, "host/all_metadata.tsv"),
        },
        params={
            "datasets": [
                {"id": "GSE_TEST1", "batch": 1,
                 "metadata_file": os.path.join(DATA, "host/GSE_TEST1_metadata.txt")},
                {"id": "GSE_TEST2", "batch": 2,
                 "metadata_file": os.path.join(DATA, "host/GSE_TEST2_metadata.txt")},
            ],
            "study": "TEST",
        }
    )
)

# ── 步骤 4（微生物）：PICRUSt2 → EC ──────────────────────────────
run_step(
    "[微生物] PICRUSt2 → EC 丰度",
    os.path.join(MODS, "meta/metagenomics/01_picrust2_to_ec.py"),
    MockSnakemake(
        input={
            "feature_file":  os.path.join(DATA, "meta/SRP_TEST1.tsv"),
            "taxonomy_file": os.path.join(DATA, "meta/SRP_TEST1.taxonomy.tsv"),
        },
        output={"0": os.path.join(OUT, "meta/SRP_TEST1_processed.tsv")},
    )
)

# ── 步骤 5（微生物）：标准化 ───────────────────────────────────────
run_step(
    "[微生物] 标准化",
    os.path.join(MODS, "meta/metagenomics/02_normalize.py"),
    MockSnakemake(
        input={"0": os.path.join(OUT, "meta/SRP_TEST1_processed.tsv")},
        output={"0": os.path.join(OUT, "meta/SRP_TEST1_normalized.tsv")},
    )
)

# ── 步骤 6（微生物）：合并 ────────────────────────────────────────
run_step(
    "[微生物] 合并数据集",
    os.path.join(MODS, "meta/metagenomics/03_merge_datasets.py"),
    MockSnakemake(
        input={"matrices": [os.path.join(OUT, "meta/SRP_TEST1_normalized.tsv")]},
        output={
            "matrix":   os.path.join(OUT, "meta/all_normalized.tsv"),
            "metadata": os.path.join(OUT, "meta/all_metadata.tsv"),
        },
        params={
            "datasets": [{"id": "SRP_TEST1", "batch": 1}],
            "study": "TEST",
        }
    )
)

# ── 步骤 7：MetaCyc 解析 ───────────────────────────────────────────
run_step(
    "[桥接] MetaCyc 解析",
    os.path.join(MODS, "bridge/metacyc_stitch/04_extract_metacyc.py"),
    MockSnakemake(
        input={"0": os.path.join(DATA, "bridge/reactions.tab")},
        output={"0": os.path.join(OUT, "bridge/microbial_metabolic_reactions.tsv")},
    )
)

# ── 步骤 8：MetaCyc → STITCH ID ────────────────────────────────────
run_step(
    "[桥接] 化合物 ID 转换",
    os.path.join(MODS, "bridge/metacyc_stitch/05_convert_to_stitch.py"),
    MockSnakemake(
        input={
            "reactions":       os.path.join(OUT, "bridge/microbial_metabolic_reactions.tsv"),
            "chemical_sources": os.path.join(DATA, "bridge/chemical.sources.v5.0.tsv"),
        },
        output={"compounds": os.path.join(OUT, "bridge/cmpd_with_stitch_ids.tsv")},
    )
)

# ── 步骤 9：STITCH 解析 ────────────────────────────────────────────
run_step(
    "[桥接] STITCH 化合物-靶点解析",
    os.path.join(MODS, "bridge/metacyc_stitch/06_parse_stitch.py"),
    MockSnakemake(
        input={
            "compounds": os.path.join(OUT, "bridge/cmpd_with_stitch_ids.tsv"),
            "prot_chem": os.path.join(DATA, "bridge/9606.protein_chemical.links.detailed.v5.0.tsv"),
            "actions":   os.path.join(DATA, "bridge/9606.actions.v5.0.tsv"),
        },
        output={"0": os.path.join(OUT, "bridge/compound_target_match.tsv")},
        params={"confidence_threshold": 400},
    )
)

# ── 步骤 10：酶-化合物-宿主基因整合 ──────────────────────────────
run_step(
    "[桥接] 酶-化合物-基因整合",
    os.path.join(MODS, "bridge/metacyc_stitch/07_link_enzyme_gene.py"),
    MockSnakemake(
        input={
            "reactions":    os.path.join(OUT, "bridge/microbial_metabolic_reactions.tsv"),
            "compounds":    os.path.join(OUT, "bridge/cmpd_with_stitch_ids.tsv"),
            "targets":      os.path.join(OUT, "bridge/compound_target_match.tsv"),
            "host_meta":    os.path.join(DATA, "mock_results/host_metaanalysis_results.tsv"),
            "human_filter": os.path.join(DATA, "bridge/human_cmpd_list.txt"),
        },
        output={"links_dir": os.path.join(OUT, "bridge/enzyme_compound_gene_links/")},
    )
)

# ── 步骤 11：EMM 生成 ──────────────────────────────────────────────
run_step(
    "[整合] EMM 矩阵生成",
    os.path.join(MODS, "integration/08_generate_EMM.py"),
    MockSnakemake(
        input={"reactions": os.path.join(OUT, "bridge/microbial_metabolic_reactions.tsv")},
        output={"emm": os.path.join(OUT, "integration/EMM.tsv")},
    )
)

# ── 步骤 12：EPCS 计算 ─────────────────────────────────────────────
run_step(
    "[整合] EPCS 分数计算",
    os.path.join(MODS, "integration/09_calculate_EPCS.py"),
    MockSnakemake(
        input={
            "emm":         os.path.join(OUT, "integration/EMM.tsv"),
            "effect_size": os.path.join(DATA, "mock_results/meta_metaanalysis_results.tsv"),
        },
        output={"epcs": os.path.join(OUT, "integration/EPCS_scores.tsv")},
    )
)

# ── 步骤 13：LOGO 敏感性分析 ───────────────────────────────────────
run_step(
    "[整合] LOGO 敏感性分析",
    os.path.join(MODS, "integration/10_logo_sensitivity.py"),
    MockSnakemake(
        input={
            "feature_files":  [os.path.join(DATA, "meta/SRP_TEST1.tsv")],
            "taxonomy_files": [os.path.join(DATA, "meta/SRP_TEST1.taxonomy.tsv")],
            "genus_list":     os.path.join(DATA, "sensitivity/genus_list.txt"),
            "ec_list":        os.path.join(DATA, "sensitivity/ec_list.txt"),
            "emm":            os.path.join(OUT, "integration/EMM.tsv"),
            "meta_results":   os.path.join(DATA, "mock_results/meta_metaanalysis_results.tsv"),
        },
        output={
            "summary": os.path.join(OUT, "sensitivity/LOGO_summary.tsv"),
            "outdir":  os.path.join(OUT, "sensitivity/per_genus/"),
        },
        params={"datasets": ["SRP_TEST1"]},
    )
)

# ── 结果汇总 ──────────────────────────────────────────────────────
print(f"\n{'='*55}")
print(f"  测试结果")
print(f"{'='*55}")
print(f"  通过：{len(passed)}/{len(passed)+len(failed)} 步")
for name in passed:
    print(f"    ✓ {name}")
if failed:
    print(f"\n  失败：")
    for name in failed:
        print(f"    ✗ {name}")
    sys.exit(1)
else:
    print("\n  🎉 所有步骤通过！")
