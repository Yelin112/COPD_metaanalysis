# OLP 微阵列数据处理（tinyarray + AnnoProbe 方案）
#
# 替代路径：比原版 Python 脚本更简洁，读取已下载的 series_matrix.txt，
# 输出与 Snakemake 流程完全兼容的文件。
#
# 产出：
#   data/OLP/raw/host/{GSE}_metadata.txt
#   results/OLP/processed/host/{GSE}_normalized.tsv
#
# 用法：
#   /usr/lib/R/bin/Rscript scripts/10_process_olp_microarray.R   # 终端
#   source("scripts/10_process_olp_microarray.R")                 # RStudio

library(stringr)
library(AnnoProbe)
library(tinyarray)

# ── 1. 数据集配置 ──────────────────────────────────────────────────
DATASETS <- list(

    GSE52130 = list(
        matrix_file    = "data/OLP/raw/host/GSE52130_series_matrix.txt",
        log2_transform = FALSE,
        # disease 信息在 Sample_title 中
        infer_condition = function(pd) {
            title <- tolower(pd$title)
            cond  <- rep(NA_character_, nrow(pd))
            cond[str_detect(title, "oral lichen planus")] <- "OLP"
            cond[str_detect(title, "control oral")]       <- "Control"
            cond   # 生殖器样本（genital）保持 NA，后续自动排除
        }
    ),

    GSE38616 = list(
        matrix_file    = "data/OLP/raw/host/GSE38616_series_matrix.txt",
        log2_transform = FALSE,
        # disease 信息在 characteristics_ch1 中
        infer_condition = function(pd) {
            # 合并所有 characteristics 列为一个字符串
            char_cols <- grep("characteristics_ch1", colnames(pd), value = TRUE)
            txt  <- tolower(apply(pd[, char_cols, drop = FALSE], 1,
                                  paste, collapse = " "))
            cond <- rep(NA_character_, nrow(pd))
            cond[str_detect(txt, "oral lichen planus")] <- "OLP"
            cond[str_detect(txt, "healthy")]            <- "Control"
            cond
        }
    )
)

# ── 2. 探针→基因折叠（IQR 最大探针代表基因）─────────────────────────
collapse_probes <- function(expr_mat, ids) {
    # ids: data.frame，列 probe_id / symbol
    common <- intersect(rownames(expr_mat), ids$probe_id)
    if (length(common) == 0) stop("探针 ID 与注释无交集，请检查平台")

    expr_mat <- expr_mat[common, , drop = FALSE]
    ids      <- ids[ids$probe_id %in% common, ]
    ids      <- ids[!is.na(ids$symbol) & ids$symbol != "", ]

    iqr_vals <- apply(expr_mat[ids$probe_id, , drop = FALSE], 1,
                      IQR, na.rm = TRUE)
    df <- data.frame(probe  = ids$probe_id,
                     gene   = ids$symbol,
                     iqr    = iqr_vals,
                     stringsAsFactors = FALSE)
    df <- df[order(df$iqr, decreasing = TRUE), ]
    best <- df[!duplicated(df$gene), ]

    result <- expr_mat[best$probe, , drop = FALSE]
    rownames(result) <- best$gene
    result
}

# ── 3. 列 z-score ─────────────────────────────────────────────────
col_zscore <- function(mat) {
    mat <- scale(mat, center = TRUE, scale = TRUE)
    mat[rowSums(is.nan(mat)) == 0, , drop = FALSE]
}

# ── 4. 主循环 ─────────────────────────────────────────────────────
dir.create("results/OLP/processed/host", recursive = TRUE, showWarnings = FALSE)

for (gse_id in names(DATASETS)) {
    cfg <- DATASETS[[gse_id]]
    message("\n══ ", gse_id, " ══")

    # ── 4a. 读取数据（GEOquery 读本地文件）─────────────────────────
    gse_obj  <- GEOquery::getGEO(filename = cfg$matrix_file,
                                  GSEMatrix = TRUE, getGPL = FALSE)
    expr_mat <- Biobase::exprs(gse_obj)
    pd       <- Biobase::pData(gse_obj)
    gpl      <- as.character(gse_obj@annotation)   # 平台号

    message("  平台: ", gpl,
            "  探针: ", nrow(expr_mat),
            "  样本: ", ncol(expr_mat))

    # ── 4b. 推断 condition，过滤无关样本 ───────────────────────────
    conditions <- cfg$infer_condition(pd)
    keep       <- !is.na(conditions)
    expr_mat   <- expr_mat[, keep, drop = FALSE]
    pd         <- pd[keep, , drop = FALSE]
    conditions <- conditions[keep]

    message("  OLP=",     sum(conditions == "OLP"),
            "  Control=", sum(conditions == "Control"),
            "  排除=",    sum(!keep))

    # ── 4c. log2 变换（可选）──────────────────────────────────────
    if (isTRUE(cfg$log2_transform)) {
        expr_mat <- log2(expr_mat + 1)
        message("  已应用 log2(x+1)")
    }

    # ── 4d. 获取探针注释并折叠为基因 ──────────────────────────────
    ids      <- AnnoProbe::idmap(gpl)
    gene_mat <- collapse_probes(expr_mat, ids)
    message("  基因数（折叠后）: ", nrow(gene_mat))

    # ── 4e. 列 z-score ────────────────────────────────────────────
    norm_mat <- col_zscore(gene_mat)
    message("  基因数（z-score 后）: ", nrow(norm_mat))

    # ── 4f. 写出文件 ───────────────────────────────────────────────
    meta_out <- file.path("data/OLP/raw/host",
                          paste0(gse_id, "_metadata.txt"))
    norm_out <- file.path("results/OLP/processed/host",
                          paste0(gse_id, "_normalized.tsv"))

    meta <- data.frame(sample_id = rownames(pd),
                       condition = conditions,
                       stringsAsFactors = FALSE)
    write.table(meta, meta_out, sep = "\t", row.names = FALSE, quote = FALSE)
    message("  元数据 → ", meta_out)

    result <- cbind(gene_id = rownames(norm_mat), as.data.frame(norm_mat))
    rownames(result) <- NULL
    write.table(result, norm_out, sep = "\t", row.names = FALSE, quote = FALSE)
    message("  标准化矩阵 → ", norm_out)
}

message("\n✅ 全部完成。")
