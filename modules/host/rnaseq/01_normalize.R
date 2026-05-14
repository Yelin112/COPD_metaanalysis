# RNA-seq 标准化模块（待实现）
# 功能：featureCounts/HTSeq counts 矩阵 → TMM/DESeq2 VST 标准化 → z-score 矩阵
# 输出格式与 host/microarray/02_normalize.py 完全一致，下游无需改动
#
# 实现时需要的 R 包：DESeq2 / edgeR
# 输入字段（来自 study_config.yaml）：
#   datasets[].counts_file   — featureCounts 输出的 counts 矩阵
#   datasets[].metadata_file — 样本元数据

stop("RNA-seq 标准化模块尚未实现。",
     "请实现 modules/host/rnaseq/01_normalize.R 后再使用 host_omics.type: rnaseq。")
