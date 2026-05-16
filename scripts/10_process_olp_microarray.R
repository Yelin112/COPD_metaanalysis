# OLP 微阵列数据处理（GEOquery + limma 方案）
#
# 替代路径：Python 脚本（02 + Snakemake microarray 规则）出现问题时使用
# 输出文件路径和格式与 Snakemake 流程完全兼容，可直接衔接后续步骤
#
# 产出：
#   data/OLP/raw/host/{GSE}_metadata.txt        （sample_id / condition）
#   results/OLP/processed/host/{GSE}_normalized.tsv  （gene × sample z-score）
#
# 用法：
#   Rscript scripts/10_process_olp_microarray.R          # 终端
#   source("scripts/10_process_olp_microarray.R")         # RStudio

# ── 0. 包检查与加载 ────────────────────────────────────────────────
required_pkgs <- c(
    "GEOquery", "limma", "Biobase", "AnnotationDbi",
    "illuminaHumanv4.db",           # GPL10558
    "hugene10sttranscriptcluster.db" # GPL6244
)

missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
    message("安装缺失包: ", paste(missing, collapse = ", "))
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(missing, ask = FALSE)
}

suppressPackageStartupMessages({
    library(GEOquery)
    library(limma)
    library(Biobase)
    library(AnnotationDbi)
    library(illuminaHumanv4.db)
    library(hugene10sttranscriptcluster.db)
})

# ── 1. 数据集配置 ──────────────────────────────────────────────────
DATASETS <- list(

    GSE52130 = list(
        platform      = "GPL10558",
        log2_transform = FALSE,   # Illumina BeadChip 通常已 log 变换
        matrix_file   = "data/OLP/raw/host/GSE52130_series_matrix.txt",
        # condition 推断：disease 信息在 Sample_title 中
        infer_condition = function(title, chars) {
            txt <- tolower(paste(title, paste(chars, collapse = " ")))
            if (grepl("oral lichen planus", txt))  return("OLP")
            if (grepl("control oral",       txt))  return("Control")
            return(NA_character_)   # 生殖器样本等直接排除
        }
    ),

    GSE38616 = list(
        platform       = "GPL6244",
        log2_transform = FALSE,
        matrix_file    = "data/OLP/raw/host/GSE38616_series_matrix.txt",
        # condition 推断：disease 信息在 Sample_characteristics_ch1 中
        infer_condition = function(title, chars) {
            txt <- tolower(paste(title, paste(chars, collapse = " ")))
            if (grepl("oral lichen planus", txt)) return("OLP")
            if (grepl("healthy",            txt)) return("Control")
            return(NA_character_)
        }
    )
)

# ── 2. 探针→基因映射函数 ────────────────────────────────────────────
#   每个基因保留 IQR 最大的探针（与 Python 版 01_probe_to_gene.py 逻辑一致）
probe_to_gene <- function(expr_mat, platform) {

    if (platform == "GPL10558") {
        db  <- illuminaHumanv4.db
        key <- "SYMBOL"
    } else if (platform == "GPL6244") {
        db  <- hugene10sttranscriptcluster.db
        key <- "SYMBOL"
    } else {
        stop("未支持的平台: ", platform)
    }

    symbols <- mapIds(db, keys = rownames(expr_mat),
                      column = key, keytype = "PROBEID",
                      multiVals = "first")

    # 去掉无法映射的探针
    valid <- !is.na(symbols) & symbols != ""
    expr_mat <- expr_mat[valid, , drop = FALSE]
    symbols  <- symbols[valid]

    # 每个基因取 IQR 最大探针
    iqr_vals <- apply(expr_mat, 1, IQR, na.rm = TRUE)
    df <- data.frame(probe = rownames(expr_mat),
                     gene  = symbols,
                     iqr   = iqr_vals,
                     stringsAsFactors = FALSE)
    df <- df[order(df$iqr, decreasing = TRUE), ]
    best <- df[!duplicated(df$gene), ]

    expr_mat[best$probe, , drop = FALSE] |>
        (\(m) { rownames(m) <- best$gene; m })()
}

# ── 3. 列 z-score 函数 ─────────────────────────────────────────────
col_zscore <- function(mat) {
    mat <- scale(mat, center = TRUE, scale = TRUE)
    # 移除常数基因（SD=0 → NaN）
    mat[rowSums(is.nan(mat)) == 0, , drop = FALSE]
}

# ── 4. 主循环 ─────────────────────────────────────────────────────
dir.create("results/OLP/processed/host", recursive = TRUE, showWarnings = FALSE)
dir.create("data/OLP/raw/host",          recursive = TRUE, showWarnings = FALSE)

for (gse_id in names(DATASETS)) {
    cfg <- DATASETS[[gse_id]]
    message("\n══ ", gse_id, " (", cfg$platform, ") ══")

    # ── 4a. 读取表达数据 ────────────────────────────────────────────
    gse <- tryCatch(
        getGEO(filename = cfg$matrix_file, GSEMatrix = TRUE, getGPL = FALSE),
        error = function(e) {
            message("  本地文件读取失败，尝试在线下载...")
            getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
        }
    )
    # getGEO(filename=) 返回单个对象；getGEO(GEO=) 返回列表
    if (is.list(gse)) gse <- gse[[1]]

    expr_mat <- exprs(gse)          # 探针 × 样本
    pd       <- pData(gse)          # 样本元数据 data.frame

    message("  原始矩阵: ", nrow(expr_mat), " 探针 × ", ncol(expr_mat), " 样本")

    # ── 4b. 提取样本 condition ──────────────────────────────────────
    # characteristics_ch1 可能有多列，统一取所有 "characteristics_ch1.*" 列
    char_cols <- grep("^characteristics_ch1", colnames(pd), value = TRUE)

    conditions <- mapply(function(title, ...) {
        chars <- unlist(list(...))
        cfg$infer_condition(title, chars)
    }, pd$title, pd[, char_cols, drop = FALSE],
    SIMPLIFY = TRUE)

    meta <- data.frame(
        sample_id = rownames(pd),
        condition = conditions,
        stringsAsFactors = FALSE
    )

    # 只保留 OLP / Control，排除 NA（生殖器样本等）
    keep_samples <- meta$sample_id[!is.na(meta$condition)]
    meta         <- meta[!is.na(meta$condition), ]
    expr_mat     <- expr_mat[, keep_samples, drop = FALSE]

    message("  保留样本: OLP=", sum(meta$condition == "OLP"),
            "  Control=", sum(meta$condition == "Control"),
            "  排除=", sum(is.na(conditions)))

    # ── 4c. log2 变换（可选）──────────────────────────────────────
    if (isTRUE(cfg$log2_transform)) {
        expr_mat <- log2(expr_mat + 1)
        message("  已应用 log2(x+1)")
    }

    # ── 4d. 探针→基因折叠 ─────────────────────────────────────────
    gene_mat <- probe_to_gene(expr_mat, cfg$platform)
    message("  折叠后: ", nrow(gene_mat), " 基因")

    # ── 4e. 列 z-score ────────────────────────────────────────────
    norm_mat <- col_zscore(gene_mat)
    message("  z-score 后: ", nrow(norm_mat), " 基因（移除常数基因）")

    # ── 4f. 写出文件 ───────────────────────────────────────────────
    meta_out <- file.path("data/OLP/raw/host",
                          paste0(gse_id, "_metadata.txt"))
    norm_out <- file.path("results/OLP/processed/host",
                          paste0(gse_id, "_normalized.tsv"))

    write.table(meta, meta_out, sep = "\t", row.names = FALSE, quote = FALSE)
    message("  元数据 → ", meta_out)

    result <- cbind(gene_id = rownames(norm_mat), as.data.frame(norm_mat))
    rownames(result) <- NULL
    write.table(result, norm_out, sep = "\t", row.names = FALSE, quote = FALSE)
    message("  标准化矩阵 → ", norm_out)
}

message("\n✅ 全部完成。后续直接运行 Snakemake 的 merge / metaanalysis 步骤即可。")
