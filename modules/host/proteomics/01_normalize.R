# 蛋白质组标准化模块（待实现）
# 功能：MaxQuant LFQ 强度 / DIA-NN 输出 → median centering + log2 → z-score 矩阵
# 输出格式与 host/microarray/02_normalize.py 完全一致
#
# 实现时需要的 R 包：limma / NormalyzerDE
# 输入字段（来自 study_config.yaml）：
#   datasets[].intensity_file — MaxQuant proteinGroups.txt 或 DIA-NN report.tsv
#   datasets[].metadata_file  — 样本元数据

stop("蛋白质组标准化模块尚未实现。",
     "请实现 modules/host/proteomics/01_normalize.R 后再使用 host_omics.type: proteomics。")
