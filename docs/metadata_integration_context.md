# OLP 微生物组荟萃分析 — Metadata 整合上下文

本文档供 Codex/AI 助手理解并继续 Metadata 整合任务。

---

## 项目背景

**疾病**：口腔扁平苔藓（Oral Lichen Planus, OLP）  
**目标**：将 13 个公开微生物组数据集的样本元数据整合为统一表格，
标注每个样本的分组（OLP 患者 / 健康对照 / 排除），
用于后续多数据集荟萃分析。

---

## 目标输出文件

`results/OLP/metadata/combined_metadata_annotated.tsv`

### 列定义（16 列，制表符分隔）

| 列名 | 说明 | 示例值 |
|------|------|--------|
| `sample_id` | SRR/DRR/ERR/CRR run accession | SRR8952025 |
| `study_id` | 项目 accession | PRJNA542018 |
| `biosample_id` | SAMN/SAMD/ERS accession | SAMN11606556 |
| `condition` | **核心字段**，见下方规则 | OLP / Control / FILL_ME |
| `sample_type` | 样本类型（文献固定值） | saliva / swab / tissue / plaque_supra / plaque_sub |
| `amplicon_region` | 扩增区域（文献固定值） | V3-V4 / V4 / V1-V2 / WGS |
| `primers` | 引物（文献固定值） | 341F/806R |
| `library_layout` | PAIRED / SINGLE | PAIRED |
| `platform` | 测序平台 | ILLUMINA |
| `instrument_model` | 仪器型号 | Illumina MiSeq |
| `read_count` | reads 数量 | 85234 |
| `sample_title` | 原始标题（辅助判断） | OLP patient saliva |
| `library_name` | 文库名（辅助判断） | OLP.12 |
| `age` | 年龄 | 45 |
| `sex` | 性别 | M / F |
| `notes` | 备注 / 无法推断时的原始属性提示 | disease_state=oral lichen planus |

---

## condition 字段赋值规则

优先级从高到低：

### 1. 排除规则（Exclude_*）
以下关键词匹配 → 标为对应 Exclude 标签，**不参与主分析**

| 关键词 | 标签 |
|--------|------|
| genital / vulvar / penile / vaginal / esophageal | `Exclude_nonoral` |
| treated / treatment / after therapy / post-treatment | `Exclude_treated` |
| hypothyroid / hypot / hypo-t | `Exclude_other_disease` |
| menstrual / menstruation / cycle / follicular / luteal | `Exclude_other_disease` |
| periodontitis / gingivitis / caries / cancer / tumor | `Exclude_other_disease` |

### 2. OLP 关键词
`olp` / `oral lichen planus` / `lichen planus` / `patient` / `case` / `diseased` / `affected` → `OLP`

### 3. 对照关键词
`control` / `ctrl` / `healthy` / `normal` / `hc` / `unaffected` / `disease free` → `Control`

### 4. 无法推断
→ `FILL_ME`（需手工标注）

**重要**：样本**保留**在表格中，下游分析用 `condition in ["OLP", "Control"]` 过滤。

---

## 各数据集状态

### ✅ 已完成标注（无需处理）

| 数据集 | 来源 | condition 来源 | 样本数（估计） |
|--------|------|----------------|----------------|
| **CRA008410** | 中国 GSA | LibraryName 前缀：`OLP.N`→OLP，`C.N`/`H.N`→Control | ~120 |
| **PRJDB12280** | DDBJ/日本 | 论文 Supp Table 2，已筛选 Presort fraction：SAMD01578002–022（OLP n=21）+ SAMD01578047–102（Control n=56）| 77 |

### ❌ 排除（无临床分组信息）

| 数据集 | 原因 |
|--------|------|
| **PRJEB90477** | ENA 元数据无任何临床字段，联系作者成本高 |

### ⏳ 需要从 NCBI 下载 SraRunTable

以下数据集需在 **NCBI SRA Run Selector** 下载 SraRunTable：

| 数据集 | 下载链接 | 重点查找字段 | 备注 |
|--------|----------|--------------|------|
| PRJNA542018 | https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA542018 | LibraryName / disease_state | V4, 唾液, n≈40 |
| PRJNA306560 | https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA306560 | source_name / Sample Name | V4, 唾液 |
| PRJNA555458 | https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA555458 | disease_state / condition | V3-V4, 唾液 |
| PRJNA556311 | https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA556311 | source_name / tissue | V3-V4, 组织 |
| PRJNA598825 | https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA598825 | disease_state | V3-V4, 拭子 |
| PRJNA690677 | https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA690677 | Sample Name / condition | V3-V4, 龈下菌斑 |
| PRJNA1049117 | https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1049117 | source_name / LibraryName | V3-V4, 混合 |
| PRJNA1043432 | https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1043432 | disease_state / condition | V3-V4, 拭子 |
| PRJNA1201607 | https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1201607 | source_name / disease | WGS（鸟枪法），龈上菌斑 |
| PRJNA512581 | https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA512581 | Sample Name / tissue | 混合扩增区域，慎用 |

**下载方式**：
1. 打开链接
2. 右上角点 **"Metadata"** 按钮
3. 下载得到 `SraRunTable.txt`（实际为 CSV 格式，逗号分隔）

---

## 文献固定字段（已确认，不需要从 SraRunTable 提取）

```python
DATASET_FIXED = {
    "CRA008410":    {"amplicon_region": "V3-V4", "sample_type": "swab",         "primers": "341F/806R"},
    "PRJDB12280":   {"amplicon_region": "V1-V2", "sample_type": "saliva",        "primers": "27Fmod/338R"},
    "PRJEB90477":   {"amplicon_region": "V1-V2", "sample_type": "plaque_supra",  "primers": "27F/338R"},
    "PRJNA1043432": {"amplicon_region": "V3-V4", "sample_type": "swab",          "primers": "341F/806R"},
    "PRJNA1049117": {"amplicon_region": "V3-V4", "sample_type": "mixed",         "primers": "343F/798R"},
    "PRJNA1201607": {"amplicon_region": "WGS",   "sample_type": "plaque_supra",  "primers": "NA"},
    "PRJNA306560":  {"amplicon_region": "V4",    "sample_type": "saliva",        "primers": "515F/806R"},
    "PRJNA512581":  {"amplicon_region": "V1-V3_or_V3-V4", "sample_type": "swab_or_tissue", "primers": "mixed"},
    "PRJNA542018":  {"amplicon_region": "V4",    "sample_type": "saliva",        "primers": "515F/806R"},
    "PRJNA555458":  {"amplicon_region": "V3-V4", "sample_type": "saliva",        "primers": "341F/806R"},
    "PRJNA556311":  {"amplicon_region": "V3-V4", "sample_type": "tissue",        "primers": "341F/806R"},
    "PRJNA598825":  {"amplicon_region": "V3-V4", "sample_type": "swab",          "primers": "V3-V4_Illumina"},
    "PRJNA690677":  {"amplicon_region": "V3-V4", "sample_type": "plaque_sub",    "primers": "338F/806R"},
}
```

---

## 任务：从 SraRunTable 提取 condition

### 输入
- 每个数据集目录下的 `SraRunTable.txt`（NCBI 下载，CSV 格式）
- `combined_metadata.tsv`（已有技术字段，condition=FILL_ME 待填）

### 关键逻辑

SraRunTable 中 Run accession（`Run` 列，如 SRR8952025）
与 `combined_metadata.tsv` 的 `sample_id` 列一一对应。

需要在 SraRunTable 中找含分组信息的列，常见列名：
- `disease_state`、`disease`、`health_state`、`host_health_state`
- `condition`、`status`、`diagnosis`、`phenotype`
- `source_name`、`isolation_source`、`tissue`
- `LibraryName`、`Sample Name`、`sample_name`
- `subject_is_affected`、`affected`

对找到的值按上述关键词规则推断 OLP / Control / Exclude_* / FILL_ME。

同时提取 `age`（列名：age / host_age / Age）
和 `sex`（列名：sex / gender / Sex / host_sex / Gender），
标准化：male/m → M，female/f → F。

### 输出

更新后的 `combined_metadata_annotated.tsv`，
以及一份 `sra_merge_report.tsv` 记录每行的推断依据：

```
sample_id  study_id  cond_before  cond_after  evidence           age  sex
SRR8952025 PRJNA542018 FILL_ME    OLP         LibraryName='OLP_1' 45  F
SRR8952026 PRJNA542018 FILL_ME    Control     disease_state='healthy' 38 M
```

---

## 注意事项

1. **不要删行**：无法判断的保留 `FILL_ME`，在 `notes` 列写原始属性值
2. **PRJNA512581 慎用**：扩增区域混杂（V1-V3 + V3-V4），可能需要按样本分层分析或排除
3. **PRJNA1201607 单独处理**：WGS 鸟枪法，不走 QIIME2，需 KneadData + MetaPhlAn4 流程
4. **CRA008410 已完成**：LibraryName `OLP.N`/`C.N`/`H.N` 已处理，**不要覆盖**
5. **PRJDB12280 已完成**：OLP n=21，Control n=56，**不要覆盖**
6. **去重**：同一 biosample_id 可能对应多个 run（技术重复），每个 run 单独一行但共享同一 condition
7. 分隔符：输出文件用 **制表符**（`\t`），SraRunTable 输入为逗号分隔
