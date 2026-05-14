"""
直接测量代谢组标准化模块（待实现）
功能：LC-MS/GC-MS 代谢物强度矩阵 → log 变换 → median scaling → z-score 矩阵
输出格式与 meta/metagenomics/02_normalize.py 完全一致

注意：使用此模块时，bridge.type 应配置为 direct_metabolomics

输入字段（来自 study_config.yaml）：
  datasets[].intensity_file — 代谢物峰面积/强度矩阵（行=代谢物，列=样本）
  datasets[].metadata_file  — 样本元数据
"""
raise NotImplementedError(
    "代谢组标准化模块尚未实现。\n"
    "请实现 modules/meta/metabolomics/01_normalize.py 后再使用 meta_omics.type: metabolomics。"
)
