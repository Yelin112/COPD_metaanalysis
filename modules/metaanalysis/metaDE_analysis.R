# MetaDE 随机效应荟萃分析
# 对应原脚本3 第三段（MetaDE.rawdata + 输出 FDR/Zval/Pval）
# 数据集列表从元数据文件动态读取，不再硬编码
#
# Snakemake 注入：
#   snakemake@input$matrix    — ComBat 校正后矩阵（行=基因，列=样本）
#   snakemake@input$metadata  — 样本元数据（含 study / condition 列）
#   snakemake@output$fdr      — metaFDR.tsv
#   snakemake@output$zval     — metaZval.tsv
#   snakemake@output$pval     — metaPval.tsv
#   snakemake@output$combined — metaanalysis_results.tsv（FDR + zval 合并）
#   snakemake@params$method           — "REM" | "FEM"
#   snakemake@params$individual_method — "modt" | "t" | "Wilcoxon"
#   snakemake@params$nperm            — 置换次数
#   snakemake@params$study_column     — 元数据中表示研究来源的列名

library(MetaDE)

matrix_file   <- snakemake@input$matrix
metadata_file <- snakemake@input$metadata
out_fdr       <- snakemake@output$fdr
out_zval      <- snakemake@output$zval
out_pval      <- snakemake@output$pval
out_combined  <- snakemake@output$combined

meta_method  <- snakemake@params$method
ind_method   <- snakemake@params$individual_method
nperm        <- as.integer(snakemake@params$nperm)
study_col    <- snakemake@params$study_column

# 读取数据
data <- read.table(matrix_file, header = TRUE, row.names = 1,
                   check.names = FALSE, sep = "\t")
meta <- read.table(metadata_file, header = TRUE, row.names = 1,
                   check.names = FALSE, sep = "\t")

# 从元数据动态获取数据集列表（不再硬编码）
study_names <- unique(meta[[study_col]])

# 为每个研究单独写出子矩阵文件，供 MetaDE.Read 读取
tmpdir <- tempdir()
for (study in study_names) {
  samples <- rownames(meta)[meta[[study_col]] == study]
  samples <- intersect(samples, colnames(data))
  sub_data <- data[, samples, drop = FALSE]
  write.table(sub_data, file = file.path(tmpdir, paste0(study, ".txt")),
              quote = FALSE, sep = "\t")
}

# 切换工作目录到临时目录后读取
old_wd <- getwd()
setwd(tmpdir)

n_studies <- length(study_names)
data_raw    <- MetaDE.Read(study_names, skip = rep(1, n_studies),
                            via = "txt", matched = TRUE, log = FALSE)
data_merged <- MetaDE.merge(data_raw)

result <- MetaDE.rawdata(
  data_merged,
  ind.method  = rep(ind_method, n_studies),
  paired      = rep(FALSE, n_studies),
  meta.method = meta_method,
  nperm       = nperm
)

setwd(old_wd)

# 输出各统计量
write.table(result$meta.analysis$FDR,  file = out_fdr,  sep = "\t", quote = FALSE)
write.table(result$meta.analysis$zval, file = out_zval, sep = "\t", quote = FALSE)
write.table(result$meta.analysis$pval, file = out_pval, sep = "\t", quote = FALSE)

# 合并为单一结果文件：feature_id + zval + fdr
fdr_vec  <- result$meta.analysis$FDR
zval_vec <- result$meta.analysis$zval
pval_vec <- result$meta.analysis$pval

combined <- data.frame(
  feature_id   = names(fdr_vec),
  zval         = zval_vec,
  pval         = pval_vec,
  fdr          = fdr_vec,
  row.names    = NULL,
  check.names  = FALSE
)
write.table(combined, file = out_combined, sep = "\t", quote = FALSE, row.names = FALSE)
