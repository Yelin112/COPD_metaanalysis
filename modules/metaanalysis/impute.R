# 缺失值填补
# 对应原脚本3 第一段（knn.impute）
# Snakemake 注入：snakemake@input[[1]], snakemake@output[[1]]
# snakemake@params$method ("knn" | "none"), snakemake@params$k

library(bnstruct)

input_file  <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
method      <- snakemake@params$method
k           <- as.integer(snakemake@params$k)

data <- read.table(input_file, header = TRUE, row.names = 1,
                   check.names = FALSE, sep = "\t")

if (method == "knn") {
  data_imputed <- knn.impute(as.matrix(data), k = k,
                              cat.var   = seq_len(ncol(data)),
                              to.impute = seq_len(nrow(data)),
                              using     = seq_len(nrow(data)))
} else {
  data_imputed <- as.matrix(data)
}

write.table(data_imputed, file = output_file,
            append = FALSE, quote = FALSE, sep = "\t")
