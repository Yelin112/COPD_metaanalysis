# 混合输入类型标准化模块（microarray + RNA-seq 统一输出）
#
# input_type: microarray  — 读取已折叠的基因矩阵 (_processed.tsv)
#                           → log2（可选）→ 列 z-score
# input_type: counts      — 原始计数矩阵
#                           → 低表达过滤 → DESeq2 VST → 列 z-score
# input_type: tpm | fpkm  — TPM/FPKM 矩阵
#                           → log2(x+1) → 列 z-score
#
# Snakemake 注入变量：
#   snakemake@input[["data_file"]]     — 表达矩阵（microarray: _processed.tsv；其他: counts_file）
#   snakemake@input[["metadata_file"]] — 样本元数据（sample_id + condition）
#   snakemake@output[[1]]              — 标准化后矩阵
#   snakemake@params[["input_type"]]   — "microarray" | "counts" | "tpm" | "fpkm"
#   snakemake@params[["log2_transform"]]— bool，仅 microarray 类型使用
#   snakemake@params[["dataset_id"]]   — 数据集 ID（日志）

suppressPackageStartupMessages({
    library(data.table)
})

input_type    <- snakemake@params[["input_type"]]
data_file     <- snakemake@input[["data_file"]]
meta_file     <- snakemake@input[["metadata_file"]]
output_file   <- snakemake@output[[1]]
dataset_id    <- snakemake@params[["dataset_id"]]

# ── 读取表达矩阵 ───────────────────────────────────────────────────
mat <- as.matrix(data.frame(fread(data_file), row.names = 1))
storage.mode(mat) <- "numeric"

# ── 对齐样本：只保留元数据和矩阵共有的样本 ────────────────────────
meta <- read.table(meta_file, sep = "\t", header = TRUE, row.names = 1,
                   check.names = FALSE)
common <- intersect(colnames(mat), rownames(meta))
if (length(common) == 0) {
    stop(sprintf("[%s] 表达矩阵与元数据无公共样本，请检查列名与 sample_id 是否一致",
                 dataset_id))
}
mat  <- mat[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]

message(sprintf("[%s] 输入类型: %s | 基因: %d | 样本: %d",
                dataset_id, input_type, nrow(mat), ncol(mat)))

# ── 按 input_type 标准化 ──────────────────────────────────────────
if (input_type == "microarray") {

    log2_transform <- snakemake@params[["log2_transform"]]

    # 过滤全 NA 基因
    keep <- rowSums(is.na(mat)) < ncol(mat)
    mat  <- mat[keep, , drop = FALSE]

    if (isTRUE(log2_transform)) {
        mat <- log2(mat + 1)
        message("  已应用 log2(x+1) 变换")
    }

    expr_mat <- mat

} else if (input_type == "counts") {

    suppressPackageStartupMessages(library(DESeq2))

    mat_int <- round(mat)
    storage.mode(mat_int) <- "integer"

    # 过滤低表达基因：至少 25% 样本中 count > 0
    min_samples <- max(2L, floor(ncol(mat_int) / 4L))
    keep        <- rowSums(mat_int > 0L) >= min_samples
    mat_int     <- mat_int[keep, , drop = FALSE]
    message(sprintf("  低表达过滤后保留 %d / %d 基因", sum(keep), length(keep)))

    dds <- DESeqDataSetFromMatrix(
        countData = mat_int,
        colData   = meta,
        design    = ~ 1
    )
    dds      <- estimateSizeFactors(dds)
    vsd      <- vst(dds, blind = TRUE)
    expr_mat <- assay(vsd)

} else if (input_type %in% c("tpm", "fpkm")) {

    keep     <- rowSums(mat > 0) >= 1
    mat      <- mat[keep, , drop = FALSE]
    expr_mat <- log2(mat + 1)
    message(sprintf("  全零过滤后保留 %d / %d 基因", sum(keep), length(keep)))

} else {
    stop(sprintf("未知 input_type: '%s'，应为 microarray | counts | tpm | fpkm",
                 input_type))
}

# ── 列 z-score ────────────────────────────────────────────────────
expr_mat <- scale(expr_mat, center = TRUE, scale = TRUE)

# 移除常数基因（SD=0 → NaN）
na_genes <- rowSums(is.nan(expr_mat)) > 0
if (any(na_genes)) {
    message(sprintf("  移除 %d 个常数基因（z-score 无法计算）", sum(na_genes)))
    expr_mat <- expr_mat[!na_genes, , drop = FALSE]
}

# ── 写出 ──────────────────────────────────────────────────────────
result           <- as.data.frame(expr_mat)
result           <- cbind(gene_id = rownames(result), result)
rownames(result) <- NULL

write.table(result, file = output_file, sep = "\t",
            row.names = FALSE, quote = FALSE)

message(sprintf("[%s] 完成：输出 %d 基因 × %d 样本 → %s",
                dataset_id, nrow(result), ncol(result) - 1L, output_file))
