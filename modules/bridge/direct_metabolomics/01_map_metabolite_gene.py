"""
直接代谢组桥接模块（待实现）
功能：利用 HMDB 或自定义数据库，将代谢物荟萃分析结果映射到宿主基因靶点
输出格式与 08_generate_EMM.py 完全一致（行=代谢物，列=宿主基因，值=互作权重）
下游 09_calculate_EPCS.py 无需任何修改

输入字段（来自 study_config.yaml bridge.direct_metabolomics）：
  compound_gene_database — 代谢物-基因互作数据库（TSV 格式：metabolite_id | gene_id | score）
"""
raise NotImplementedError(
    "直接代谢组桥接模块尚未实现。\n"
    "请实现 modules/bridge/direct_metabolomics/01_map_metabolite_gene.py 后再使用 bridge.type: direct_metabolomics。"
)
