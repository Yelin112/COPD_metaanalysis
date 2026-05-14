"""
宏转录组标准化模块（待实现）
功能：PICRUSt2 宏转录组 EC 贡献表（基于 mRNA reads）→ TPM 标准化 → log2 矩阵
输出格式与 meta/metagenomics/02_normalize.py 完全一致

输入字段（来自 study_config.yaml）：
  datasets[].feature_file  — PICRUSt2 宏转录组 stratified EC 贡献表
  datasets[].metadata_file — 样本元数据
"""
raise NotImplementedError(
    "宏转录组标准化模块尚未实现。\n"
    "请实现 modules/meta/metatranscriptomics/01_normalize.py 后再使用 meta_omics.type: metatranscriptomics。"
)
