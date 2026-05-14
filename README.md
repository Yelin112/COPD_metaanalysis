# 多组学荟萃分析框架

配置驱动的多组学荟萃分析流程，支持**换疾病**和**换组学类型**。
原始项目：Wang et al. "Multi-omic meta-analysis identifies functional signatures of airway microbiome in COPD"

---

## 目录结构

```
├── config/
│   ├── databases.yaml          # 数据库路径（MetaCyc、STITCH 等，与疾病无关）
│   └── COPD/
│       └── study_config.yaml   # COPD 研究配置（换疾病只需新建同格式文件）
├── workflow/
│   ├── Snakefile               # 主工作流入口
│   ├── rules/                  # 各组学类型的 Snakemake 规则模块
│   └── envs/                   # Conda 环境配置
├── modules/
│   ├── host/                   # 宿主组学处理（microarray/rnaseq/proteomics/methylation）
│   ├── meta/                   # 微生物/环境侧处理（metagenomics/metatranscriptomics/metabolomics）
│   ├── metaanalysis/           # 荟萃分析（ComBat + MetaDE，R）
│   ├── bridge/                 # 跨组学桥接数据库整合
│   └── integration/            # EPCS 分数计算 + LOGO 敏感性分析
├── tests/
│   ├── generate_test_data.py   # 生成模拟测试数据
│   ├── run_tests.py            # 端到端测试（跳过 R 步骤）
│   └── data/                   # 模拟输入数据
└── legacy/                     # 原始 Perl/R 脚本（仅供参考）
```

---

## 环境要求

### Python 依赖
```bash
pip install -r requirements.txt
```

### R 依赖（仅荟萃分析步骤需要）
```r
install.packages("BiocManager")
BiocManager::install(c("sva", "MetaDE", "bnstruct"))
```

### Snakemake（运行完整流程）
```bash
pip install snakemake
# 或使用 conda
conda install -c bioconda snakemake
```

---

## 快速开始

### 1. 克隆仓库
```bash
git clone https://github.com/Yelin112/COPD_metaanalysis.git
cd COPD_metaanalysis
pip install -r requirements.txt
```

### 2. 运行模拟数据测试（验证安装）
```bash
python3 tests/generate_test_data.py   # 生成模拟数据
python3 tests/run_tests.py            # 运行所有 Python 步骤测试
```
预期输出：`🎉 所有步骤通过！`（15/15）

### 3. 准备真实数据
按以下结构放置数据文件（路径与 `study_config.yaml` 中一致）：
```
data/COPD/raw/
├── host/
│   ├── GSE11906_metadata.txt        # 样本元数据（见格式说明）
│   ├── GSE11906_series_matrix.txt   # GEO 直接下载
│   └── ...（19 个数据集）
├── meta/
│   ├── SRP066375.tsv                # PICRUSt2 EC 贡献表
│   ├── SRP066375.taxonomy.tsv       # QIIME2 分类注释
│   └── SRP136124.*
└── bridge/
    └── human_cmpd_list.txt          # 宿主源化合物黑名单
```

探针映射文件（放在 `databases.yaml` 中 `probe_mappings.base_dir` 指定的目录）：
```
/你的路径/platform_mappings/
├── GPL96_mapping.txt     # Affymetrix HG-U133A
├── GPL570_mapping.txt    # Affymetrix HG-U133_Plus_2
└── GPL14550_mapping.txt  # Agilent 4x44K
```

### 4. 下载并配置数据库
编辑 `config/databases.yaml`，将路径指向本地安装的数据库：
- **MetaCyc**：https://metacyc.org/download.shtml （需注册，免费学术使用）
- **STITCH v5**：http://stitch.embl.de/cgi/download.pl

### 5. 运行完整流程
```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile config/COPD/study_config.yaml \
  --config dbconfig=config/databases.yaml \
  --cores 8 \
  --use-conda
```

结果输出在 `results/COPD/`。

---

## 输入文件格式说明

### 样本元数据文件（`*_metadata.txt`）
制表符分隔，**必须包含 `sample_id` 和 `condition` 列**：
```
sample_id       condition    age    sex
GSM123456       COPD         65     M
GSM123457       Control      62     F
```
`condition` 列的值必须与 `study_config.yaml` 中的 `disease_label`/`control_label` 完全一致。

### 微阵列原始矩阵（`*_series_matrix.txt`）
从 NCBI GEO 直接下载，无需修改。

### 探针-基因映射文件（`{平台号}_mapping.txt`）
制表符分隔，两列，**不含多基因映射行（含 `///` 的行）**：
```
probe_id      gene_symbol
1007_s_at     DDR1
1053_at       RFC2
```

### PICRUSt2 EC 贡献表
PICRUSt2 的 stratified 分层输出，前两列为 EC 号和 ASV ID：
```
function    taxon      Sample1    Sample2
1.1.1.1     ASV_001    0.023      0.041
```

### QIIME2 分类注释文件
```
Feature ID    Taxon
ASV_001       d__Bacteria;...;g__Streptococcus;s__sp
```

---

## 换疾病

只需创建新的配置文件，代码无需改动：
```bash
cp -r config/COPD config/新疾病名
# 编辑 config/新疾病名/study_config.yaml
#   修改：study.name, disease_label, control_label
#   替换：host_omics.datasets 列表（数据集 ID、路径、批次）
#   替换：meta_omics.datasets 列表
```

运行时指定新配置：
```bash
snakemake --configfile config/新疾病名/study_config.yaml ...
```

## 换组学类型

修改 `study_config.yaml` 中的 `type` 字段：

| 字段 | 当前支持 | 待实现（有占位脚本） |
|------|---------|-------------------|
| `host_omics.type` | `microarray` ✅ | `rnaseq` / `proteomics` / `methylation` |
| `meta_omics.type` | `metagenomics` ✅ | `metatranscriptomics` / `metabolomics` |
| `bridge.type` | `metacyc_stitch` ✅ | `direct_metabolomics` |

---

## 流程说明

```
宿主组学数据（多数据集）          微生物/环境侧数据（多数据集）
  ↓ 步骤1-3：处理+标准化+合并        ↓ 步骤1-3：处理+标准化+合并
  ↓ 步骤4：缺失值填补（R）            ↓
  ↓ 步骤5：ComBat 批次校正（R）        ↓ 步骤4：MetaDE 荟萃分析（R）
  ↓ 步骤6：MetaDE 荟萃分析（R）        ↓ → meta_metaanalysis_results.tsv
  ↓ → host_metaanalysis_results.tsv   ↓
        ↓                              ↓
        ├── MetaCyc 解析 ──────────────┤
        ├── STITCH ID 转换              │
        ├── STITCH 靶点解析             │
        └── 酶-化合物-基因整合          │
                ↓                      ↓
           EMM 矩阵生成 ← 代谢反应方向权重
                ↓
           EPCS 分数计算（= EMM × meta 效应量）
                ↓
           LOGO 敏感性分析（逐属排除验证）
                ↓
           results/疾病名/integration/EPCS_scores.tsv
```

**EPCS**（Environmental-to-Phenotype Contribution Score，原 PRMT）：
预测微生物产生的代谢物对宿主基因表达的贡献分数，正值表示在疾病状态下该代谢物产生增加。

---

## 输出文件说明

| 路径 | 内容 |
|------|------|
| `results/{疾病}/host_metaanalysis/metaanalysis_results.tsv` | 宿主基因差异表达荟萃分析结果（zval + FDR） |
| `results/{疾病}/meta_metaanalysis/metaanalysis_results.tsv` | 微生物 EC 基因差异荟萃分析结果 |
| `results/{疾病}/integration/EMM.tsv` | 酶-化合物贡献矩阵（归一化） |
| `results/{疾病}/integration/EPCS_scores.tsv` | 代谢物 EPCS 分数排名 |
| `results/{疾病}/sensitivity/LOGO_summary.tsv` | LOGO 敏感性分析汇总 |

---

## 引用

如使用本流程，请引用原始论文：
> Wang et al. "Multi-omic meta-analysis identifies functional signatures of airway microbiome in COPD." *[期刊]* (年份).
