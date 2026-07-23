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
| `biobakery` | HUMAnN3（鸟枪法，本项目暂未用） |

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

### 微生物组（处理中）

| 数据集 | 来源 | 类型 | 状态 |
|--------|------|------|------|
| PRJNA542018 | SRA | 16S V4，单端（merged PE），40样本 | ✅ 已下载，待 QIIME2 |
| CRA008410 | GSA（国内） | 16S V3-V4，双端 2×250bp | ✅ 已下载，待 QIIME2 |
| PRJNA1201607 | SRA | 16S，待确认 | 🔄 下载中 |

数据位置：`data/OLP/raw/meta/{项目ID}/`

**PRJNA542018 特殊说明**：reads 为 merged paired-end 以单端格式提交，长度可变（43-558bp），引物（515F/806R）仍在，QIIME2 需要 cutadapt 截除。

**CRA008410 特殊说明**：引物已在上游截除，reads 固定 250bp。

## 下一步待运行

```bash
# 1. QIIME2 处理（两个可并行）
conda activate qiime2-amplicon-2023
cd /home/usr/yel/exp-OLP_meta/pipeline
nohup bash scripts/11_qiime2_PRJNA542018.sh > logs/qiime2_PRJNA542018.log 2>&1 &
nohup bash scripts/12_qiime2_CRA008410.sh   > logs/qiime2_CRA008410.log 2>&1 &

# 2. PICRUSt2 功能预测（QIIME2 完成后）
conda activate picrust2
bash scripts/13_picrust2.sh PRJNA542018
bash scripts/13_picrust2.sh CRA008410

# 3. Snakemake 主流程（所有数据准备好后）
conda activate qiime2-amplicon-2023  # 或新建 snakemake 环境
snakemake --configfile config/OLP/study_config.yaml \
          --config dbconfig=config/databases.yaml \
          --cores 8 --use-conda
```

## 项目结构

```
pipeline/
├── config/
│   ├── databases.yaml          # 数据库路径（MetaCyc、STITCH、探针注释）
│   └── OLP/
│       └── study_config.yaml   # OLP 研究配置（数据集、参数、关键词）
├── workflow/
│   ├── Snakefile               # 主流程入口（支持 microarray/rnaseq/mixed）
│   └── rules/
│       ├── host_microarray.smk # 微阵列规则
│       ├── host_rnaseq.smk     # RNA-seq 规则
│       ├── host_mixed.smk      # 混合类型规则（已实现，本项目未用）
│       ├── meta_metagenomics.smk
│       ├── metaanalysis.smk
│       ├── integration.smk
│       └── sensitivity.smk
├── modules/
│   └── host/
│       ├── microarray/         # 探针折叠、标准化、合并
│       ├── rnaseq/             # DESeq2 VST 标准化
│       └── mixed/              # 混合类型统一标准化
├── scripts/
│   ├── 06_download_olp_geo_microarray.py   # 下载 GEO 微阵列 + GSE213346 补充文件
│   ├── 07_download_olp_microbiome_sra.sh   # 微生物组下载（已被 iSeq 替代）
│   ├── 10_process_olp_microarray.R         # 宿主转录组处理（tinyarray + AnnoProbe）
│   ├── 11_qiime2_PRJNA542018.sh            # QIIME2：单端 V4
│   ├── 12_qiime2_CRA008410.sh             # QIIME2：双端 V3-V4
│   └── 13_picrust2.sh                     # PICRUSt2 功能预测
├── data/OLP/
│   ├── raw/host/               # 微阵列原始文件 + metadata
│   ├── raw/meta/               # 16S FASTQ 文件
│   └── raw/sensitivity/
│       └── genus_list.txt      # LOGO 敏感性分析菌属列表
└── results/OLP/
    └── processed/host/         # 标准化后矩阵（Snakemake 输入）
```

## 重要配置说明

### study_config.yaml 关键字段

```yaml
host_omics:
  type: "microarray"            # 当前：microarray（框架支持 rnaseq/mixed）
  datasets:
    - id: "GSE52130"
      metadata_file: "data/OLP/raw/host/GSE52130_metadata.txt"
      matrix_file:   "data/OLP/raw/host/GSE52130_series_matrix.txt"
    # condition 标签：disease_label="OLP", control_label="Control"

meta_omics:
  type: "metagenomics"
  upstream_tool: "picrust2"     # 16S → picrust2；鸟枪法 → humann3
```

### GSE52130 元数据注意事项

该数据集混有口腔 OLP 和生殖器 LP 样本，metadata 中：
- 14 个口腔样本（7 OLP + 7 Control）→ condition 已正确标注
- 9 个生殖器样本 → condition = Unknown，流程自动过滤

### 探针注释缓存

AnnoProbe 下载的 `.rda` 文件会缓存在运行目录，首次运行后无需重新下载：
- `GPL10558_bioc.rda`（Illumina HumanHT-12 V4.0）
- `GPL6244_bioc.rda`（Affymetrix HuGene-1_0-st）

## Git 分支

当前开发分支：`claude/define-project-purpose-gv8GK`
远程仓库：`git@github.com:Yelin112/COPD_metaanalysis.git`
