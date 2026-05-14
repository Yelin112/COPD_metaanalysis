# ComBat 批次效应校正
# 对应原脚本3 第二段（ComBat）
# 批次列名从参数读取，不再硬编码
#
# Snakemake 注入：
#   snakemake@input$matrix   — 填补后的矩阵
#   snakemake@input$metadata — 样本元数据（含批次列）
#   snakemake@output[[1]]    — 校正后矩阵
#   snakemake@params$batch_column — 批次列名（默认 "batch"）
#   snakemake@params$run_combat   — TRUE/FALSE

library(sva)

matrix_file   <- snakemake@input$matrix
metadata_file <- snakemake@input$metadata
output_file   <- snakemake@output[[1]]
batch_col     <- snakemake@params$batch_column
run_combat    <- as.logical(snakemake@params$run_combat)

data <- read.table(matrix_file, header = TRUE, row.names = 1,
                   check.names = FALSE, sep = "\t")
meta <- read.table(metadata_file, header = TRUE, row.names = 1,
                   check.names = FALSE, sep = "\t")

if (run_combat) {
  batch <- meta[[batch_col]]
  data_corrected <- ComBat(as.matrix(data), batch,
                            mod         = NULL,
                            par.prior   = TRUE,
                            prior.plots = FALSE)
} else {
  data_corrected <- as.matrix(data)
}

write.table(data_corrected, file = output_file,
            append = FALSE, quote = FALSE, sep = "\t")
