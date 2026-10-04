# OLP 多组学荟萃分析项目

## 项目概述

本项目将 Wang et al. COPD 多组学荟萃分析框架适配用于**口腔扁平苔藓（Oral Lichen Planus, OLP）**研究。

- **疾病**：OLP，T细胞介导的慢性口腔黏膜炎症性疾病
- **分析框架**：ComBat 批次校正 + MetaDE 随机效应荟萃分析 + EPCS 整合评分
- **宿主组学**：口腔黏膜基因表达微阵列（GEO）
- **微生物组学**：口腔16S rRNA测序（SRA/GSA）→ PICRUSt2 功能预测
- **桥接**：MetaCyc 代谢通路 + STITCH 蛋白-化合物互作

## 服务器环境

- **项目根目录**：`/home/usr/yel/exp-OLP_meta/pipeline`
- **系统 R**：`/usr/lib/R`（RStudio 使用，包库在 `/home/usr/yel/R/library-system`）
- **Conda R**：`/home/usr/yel/miniconda3/bin/Rscript`（与系统 R 不同，不要混用）
- **运行 R 脚本**：`/usr/lib/R/bin/Rscript scripts/xxx.R`

### Conda 环境

| 环境名 | 用途 |
|--------|------|
| `iseq` | 下载 SRA/GSA/ENA 数据 |
| `qiime2-amplicon-2023` | 16S 数据质控和去噪 |
| `picrust2` | 功能预测（ASV → EC） |
| `biobakery` | 鸟枪法：HUMAnN 3.9 + MetaPhlAn 4.2.5 + kneaddata 0.12.4 + fastp 1.3.6（python 钉 3.12；三处兼容性补丁见 docs/biobakery_env配置与补丁.md，env 重建必须重打）|

数据库位置：`/home/exs/data3/Software_Database/`（chocophlan/uniref/utility_mapping/kneaddata T2T/GTDB-tk/CARD/CAZy；kraken2+bracken 55G；eggnog 未下）

网络要点：conda 用 USTC http 镜像 + fetch_threads=1；data3 已 fstab 自动挂载（nofail）。

## 数据状态

### 宿主转录组（已完成处理）

| 数据集 | 平台 | 样本 | 状态 |
|--------|------|------|------|
| GSE52130 | GPL10558 Illumina HumanHT-12 V4.0 | 7 OLP + 7 Control（已排除9个生殖器样本）| ✅ 标准化完成 |
| GSE38616 | GPL6244 Affymetrix HuGene-1_0-st | 7 OLP + 7 Control | ✅ 标准化完成 |

输出文件：
- `results/OLP/processed/host/GSE52130_normalized.tsv`
- `results/OLP/processed/host/GSE38616_normalized.tsv`
- `data/OLP/raw/host/GSE52130_metadata.txt`
- `data/OLP/raw/host/GSE38616_metadata.txt`

处理脚本：`scripts/10_process_olp_microarray.R`（tinyarray + AnnoProbe 方案）
运行方式：`/usr/lib/R/bin/Rscript scripts/10_process_olp_microarray.R`

**GSE52130 元数据注意事项**：混有口腔 OLP 和生殖器 LP 样本，9 个生殖器样本 condition = Unknown，流程自动过滤。

**探针注释缓存**：AnnoProbe 下载的 `.rda` 文件缓存在运行目录（`GPL10558_bioc.rda`、`GPL6244_bioc.rda`），首次运行后无需重新下载。

### ComBat + MetaDE（scripts/31，2026-10-04 完成 ✅）

运行：`/usr/lib/R/bin/Rscript scripts/31_combat_metaDE_host.R`

**结果**：
- 公共基因 18,156 个；FDR < 0.05 差异基因 **220 个**（前列：C9ORF24、EIF1、NME2、TUFM、EEF2）
- ComBat 验证：校正前 PC1 R²=1.00（批次主导）→ 校正后批次 R²=0.00、disease PC1 R²=0.38（✅ 有效）
- 输出：`results/OLP/processed/host/metaDE/`（Meta.RDS / metaFDR/Zval/Pval.txt / pca_combat.pdf）

**已知注意事项**：
- 样本列名带研究前缀（`GSE52130.GSM1260095`），下游读 metaDE 输出时需处理
- FDR 有并列值（nperm=300 + 2 个研究，置换粒度粗，属正常）
- **MetaDE 已从 CRAN 下架**，安装命令：`remotes::install_github("cran/MetaDE")`；bnstruct、sva 用 `pak::pak()` 安装

### 单细胞（GSE211630，2026-09-22 完成 ✅）

用 Seurat 5 + Harmony + SingleR 分析（脚本 39-44，`config/OLP/scrnaseq_gse211630.yaml` 驱动）：

- 57,439 细胞（92% 保留）、26 簇、9 类型（7 文献类型 + **Mast 肥大细胞** + 骨骼肌活检带入）
- T 亚群 8 类；**核心结论：EOLP 的 CD8_Trm 中 IFNG(+2.2, p=6e-14)/IL17A(+3.95)/GZMB(+3.02)/CCL20(+3.30) 上调，IL17 模块 EOLP>NEOLP p=0.02** —— 复现 eLife 2023 的 CD8 Trm 驱动 IFN-γ/IL-17 轴结论
- 产物：`results/OLP/processed/host/GSE211630/`（rds checkpoint、图+表、NOTES.md；含 CellChat 1.x 环状图 16 张，脚本 47）
- 数据坑见 NOTES.md：5/6 文件末行截断；ITGAE 缺失于 3 个样本（GEO 小面板）；EOLP/NEOLP 性别不平衡（Y 基因全部下调）
- **局限**：Normal n=1、EOLP vs NEOLP 2v3，结论标注为探索性

### 微生物组（QIIME2 全部完成 ✅，2026-10-04 验证 8/9 通过；框架总结见 docs/项目框架与流程总结_20261004.md）

| 数据集 | 来源 | 类型 | 状态 |
|--------|------|------|------|
| PRJNA542018 | SRA | 16S V4，单端（merged PE），40样本 | ✅ Final 63.3%，中位 11,280 reads |
| CRA008410 | GSA（国内） | 16S V3-V4，双端 2×250bp | ✅ Final 43.3%，中位 34,493 reads |
| PRJNA306560 | SRA | 16S V4，双端 2×251bp，唾液 | ❌ 排除（原始中位 62 reads，终表 14 reads，实质无数据）|
| PRJNA555458 | SRA | 16S V3-V4，双端 2×300bp，组织 | ✅ Final 85.1%，中位 13,878 reads |
| PRJNA556311 | SRA | 16S V3-V4，双端 2×300bp，唾液 | ✅ Final 40.6%，中位 29,742 reads |
| PRJNA598825 | SRA | 16S V3-V4，双端 2×300bp | ✅ Final 62.9%，中位 60,764 reads |
| PRJNA690677 | SRA | 16S V3-V4，双端 2×241bp | ✅ Final 49.0%，中位 52,578 reads |
| PRJNA1049117 | SRA | 16S V3-V4，双端 2×250bp | ✅ Final 51.4%，中位 16,042 reads |
| PRJDB12280 | SRA（DDBJ） | 16S V1-V2，双端 2×251bp | ✅ Final 85.7%，中位 31,149 reads |
| PRJNA1201607 | SRA | 鸟枪法宏基因组 | ❌ 排除（需独立流程）→ shotgun 流程 A 试跑中 |
| PRJEB90477 | ENA | 16S | ❌ 排除（无临床分组信息） |
| PRJNA512581 | SRA | 16S，单端+双端混合，两种平台 | ❌ 排除（无法统一处理） |
| PRJNA1043432 | SRA | 16S V3-V4，双端 2×250bp | ❌ 排除（引物预先去除+质量值均一，DADA2 无法处理） |

数据位置：`data/OLP/raw/meta/{项目ID}/`

**PRJNA542018 特殊说明**：reads 为 merged paired-end 以单端格式提交，长度可变（43-558bp），引物（515F/806R）仍在，需 cutadapt 截除。

**CRA008410 特殊说明**：引物已在上游截除，reads 固定 250bp。

**PRJDB12280 特殊说明**：含多分组，保留 OLP/Cont，剔除 HypoT 与月经周期组（`scripts/14_filter_prjdb12280.py`）。

## QIIME2 脚本说明（config 驱动，2026-09 收敛）

16S 上游已收敛为 config 驱动：**参数唯一事实来源 = `config/OLP/16s_datasets.yaml`**（9 active + 4 excluded 条目，新数据集加条目即可）。

```bash
conda activate qiime2-amplicon-2023
bash scripts/qiime2_16s.sh <PROJ> [--config PATH] [--dry-run]   # 去噪
bash scripts/qiime2_merge_classify.sh [--skip-classifier-check]  # 合并+分区分类
bash scripts/qiime2_train_classifier.sh ...                      # 区域分类器自训（V1-V2 训练中）
bash scripts/13_picrust2.sh <PROJ>                               # 功能预测（已修 strat 表 bug）
bash scripts/30_validate_qiime2.sh                               # 结果验证（读 config）
```

旧 per-dataset 脚本已归档 `legacy/qiime2/`（参数溯源原文）。

各数据集参数速览（以 config 为准）：

| 数据集 | 区域 | 关键参数 |
|--------|------|----------|
| PRJNA542018 | V4 | 单端 merged，515F/806R，trunc 0，max-ee 2.0 |
| CRA008410 | V3-V4 | 引物已截除跳过 cutadapt，trunc 240/240，max-ee 2.0 |
| PRJNA306560 | V4 | R2 有 7bp 前缀用 **--p-anywhere-r**，trunc 200/180，max-ee 5.0 |
| PRJNA555458 | V3-V4 | 组织样本，341F/806R，trunc 260/210，max-ee 5.0 |
| PRJNA556311 | V3-V4 | 341F/806R，trunc 280/230，max-ee 2.0 |
| PRJNA598825 | V3-V4 | 341F/**805R**（少一个 G），trunc 280/230，max-ee 2.0 |
| PRJNA690677 | V3-V4 | 2×241bp 短读，338F/806R，trunc 0/0，max-ee 2.0 |
| PRJNA1049117 | V3-V4 | 343F/798R（798R 待核对原文），trunc 240/200，max-ee 2.0 |
| PRJDB12280 | V1-V2 | 27Fmod/338R，trunc 230/160，max-ee 2.0（分类器自训中） |

排除：PRJNA1043432（引物预去除+质量值均一）、PRJEB90477（无临床分组）、PRJNA1201607（WGS 走 shotgun）、PRJNA512581（单双端混合双平台）。

## 当前进度与下一步（2026-10-04）

### 阶段一：QIIME2 ✅（8/9 完成）

```bash
# 验证结果（已修复列解析 bug）：
bash scripts/30_validate_qiime2.sh
```

PRJNA306560 ⚠️：中位 14 reads，建议排除（低于 500 reads 门槛）。

### 阶段二：PICRUSt2（QIIME2 完成后）

```bash
conda activate picrust2
cd /home/usr/yel/exp-OLP_meta/pipeline

# 7 个数据集后台运行中（542018/CRA008410/556311/598825/690677/1049117/PRJDB12280）
# PRJNA555458 已完成

# 批量提交（验证通过的数据集）：
for proj in PRJNA542018 CRA008410 PRJNA555458 PRJNA556311 PRJNA598825 PRJNA690677 PRJNA1049117 PRJDB12280; do
    nohup bash scripts/13_picrust2.sh ${proj} > logs/picrust2_${proj}.log 2>&1 &
done
```

### 阶段三：下游整合（数据齐后）

- V1-V2 分类器自训完成 → `qiime2_merge_classify.sh` 全套
- PICRUSt2 EC 表汇总
- shotgun 剩余 3 样本（SRR31818480-82，PRJNA1201607，流程 A）
- **legacy 1-10 桥接链路适配 OLP（核心工作）**：PICRUSt2 EC 表 + 宿主差异基因 → EPCS 整合评分
- Snakemake 主流程适配：`snakemake --configfile config/OLP/study_config.yaml --config dbconfig=config/databases.yaml`

## 项目结构

```
pipeline/
├── config/OLP/
│   ├── 16s_datasets.yaml       # 16S 参数唯一事实来源（9 active + 4 excluded）
│   ├── scrnaseq_gse211630.yaml # scRNA 分析参数
│   ├── study_config.yaml       # 主流程配置
│   └── ../databases.yaml       # 数据库路径
├── workflow/
│   ├── Snakefile               # 主流程入口
│   └── rules/                  # host_microarray/rnaseq/mixed/meta_metagenomics/metaanalysis/integration/sensitivity
├── scripts/
│   ├── 10_process_olp_microarray.R     # 宿主微阵列（tinyarray + AnnoProbe）
│   ├── 13_picrust2.sh                  # PICRUSt2 功能预测
│   ├── 14_filter_prjdb12280.py         # PRJDB12280 样本筛选
│   ├── 21_apply_dataset_rules.py       # metadata condition 标注
│   ├── 30_validate_qiime2.sh           # QIIME2 结果验证
│   ├── 39-44_scrnaseq_*.R              # scRNA 分析（Seurat5+Harmony+SingleR）
│   ├── 47_cellchat_*.R                 # CellChat 细胞通讯
│   ├── qiime2_16s.sh                   # config 驱动的 QIIME2 去噪（新）
│   ├── qiime2_merge_classify.sh        # 合并+分区分类（新）
│   └── qiime2_train_classifier.sh      # 分类器自训（新）
├── legacy/
│   ├── 1-10/                   # COPD 框架 EPCS/STITCH 桥接链路（待适配 OLP）
│   └── qiime2/                 # 旧 per-dataset 脚本归档（参数溯源）
├── data/OLP/
│   ├── raw/host/               # 微阵列原始文件 + metadata
│   ├── raw/meta/               # 16S FASTQ
│   ├── raw/shotgun/            # shotgun FASTQ（PRJNA1201607）
│   ├── qiime2/<项目>/          # QIIME2 输出
│   └── picrust2/               # PICRUSt2 EC 表
├── results/OLP/processed/host/
│   ├── GSE52130_normalized.tsv
│   ├── GSE38616_normalized.tsv
│   └── GSE211630/              # scRNA 产物（rds/图/表/NOTES.md）
└── docs/                       # 系统地图、biobakery 补丁记录、QIIME2 指南、周报
```

## Git 分支

当前开发分支：`claude/define-project-purpose-gv8GK`
远程仓库：`git@github.com:Yelin112/COPD_metaanalysis.git`
