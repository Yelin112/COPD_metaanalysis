# 甲基化组预处理模块（待实现）
# 功能：Illumina EPIC/450K Beta 值矩阵 → M 值转换 → 过滤低变异 CpG → z-score 矩阵
# 输出格式与 host/microarray/02_normalize.py 完全一致
#
# 实现时需要的 R 包：minfi / ChAMP
# 输入字段（来自 study_config.yaml）：
#   datasets[].beta_file     — idat 处理后的 Beta 值矩阵
#   datasets[].metadata_file — 样本元数据

stop("甲基化组预处理模块尚未实现。",
     "请实现 modules/host/methylation/01_beta_to_m.R 后再使用 host_omics.type: methylation。")
