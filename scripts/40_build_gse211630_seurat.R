# 步骤40（GSE211630 scRNA-seq）：读入 + 合并 + QC → Seurat 对象
#
# 方案：6 个基因×细胞 count 矩阵（symbol、表头偏移、基因集不一致）→
#       union 基因集合并（26,599）→ 细胞 ID 加样本前缀去冲突 →
#       QC 过滤（normal-1 基因集小，nFeature 上界自适应）→ rds checkpoint
#
# 产出：
#   results/OLP/processed/host/GSE211630/checkpoints/02_qc_filtered.rds
#   results/OLP/processed/host/GSE211630/tables/qc_summary.tsv
#   results/OLP/processed/host/GSE211630/figures/01_qc_violin.png
#   results/OLP/metadata/GSE211630_samples.tsv
#
# 用法：
#   nohup /usr/lib/R/bin/Rscript scripts/40_build_gse211630_seurat.R \
#     > logs/40_gse211630_build.log 2>&1 &

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(future)
})
plan("multicore", workers = 8)
options(future.globals.maxSize = 4 * 1024^3)

source("modules/host/scrnaseq/utils.R")
CFG <- load_sc_config("config/OLP/scrnaseq_gse211630.yaml")
OUT <- CFG$out_dir

ids <- sapply(CFG$samples, `[[`, "id")

# ── 1. 读入 6 个矩阵，加样本前缀，union 合并 ──────────────────────────
# 检查点：合并矩阵已存在则直接复用（读入+转换 ~7 min，重跑不重复烧）
raw_rds <- file.path(OUT, "checkpoints/01_merged_raw.rds")
if (file.exists(raw_rds)) {
  message("══ 复用已有合并矩阵 ══ ", raw_rds)
  merged <- readRDS(raw_rds)
} else {
mats <- list()
gene_all <- character(0)
for (s in CFG$samples) {
  f <- file.path(CFG$data_dir, s$file)
  message("══ ", s$id, " ══ 读取 ", f)
  m <- read_gz_count_matrix(f)
  message("  原始: ", nrow(m), " 基因 x ", ncol(m), " 细胞")

  # barcode 去 _N 全局序号后缀 + 样本前缀（实测跨文件有 42 个冲突）
  colnames(m) <- paste0(s$id, "_", sub("_[0-9]+$", "", colnames(m)))
  if (anyDuplicated(colnames(m))) stop("barcode 去重后仍有冲突: ", s$id)

  mats[[s$id]] <- m
  gene_all <- union(gene_all, rownames(m))
}
message("union 基因数: ", length(gene_all))   # 预期 26,599

merged <- do.call(cbind, lapply(mats, function(m) {
  add <- setdiff(gene_all, rownames(m))
  if (length(add) > 0) {
    m <- rbind(m, Matrix::Matrix(0, nrow = length(add), ncol = ncol(m),
                                 dimnames = list(add, NULL), sparse = TRUE))
  }
  m[gene_all, , drop = FALSE]
}))
message("合并后: ", nrow(merged), " 基因 x ", ncol(merged), " 细胞")  # 预期 62,357
  dir.create(file.path(OUT, "checkpoints"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(merged, raw_rds)
  message("  合并矩阵缓存 → ", raw_rds)
}

seu <- CreateSeuratObject(counts = merged, project = "GSE211630",
                          min.cells = 3, min.features = 200)
message("CreateSeuratObject 后细胞数: ", ncol(seu))

# ── 2. 元数据（cell ID 首个 "_" 前即样本 id）─────────────────────────
sample_of <- sub("_.*", "", colnames(seu))
stopifnot(all(sample_of %in% ids))
cond_map <- setNames(sapply(CFG$samples, `[[`, "condition"), ids)
grp_map  <- setNames(sapply(CFG$samples, `[[`, "group"), ids)
seu$sample    <- sample_of
seu$condition <- unname(cond_map[sample_of])
seu$group     <- unname(grp_map[sample_of])

samp_tab <- data.frame(
  sample_id = ids,
  file      = sapply(CFG$samples, `[[`, "file"),
  condition = unname(cond_map[ids]),
  group     = unname(grp_map[ids]),
  n_cells   = as.numeric(table(sample_of)[ids]))
save_tsv(samp_tab, "results/OLP/metadata/GSE211630_samples.tsv")

# ── 3. QC（normal-1 无 MT- 基因，percent.mito 恒 0 属预期）────────────
seu[["percent.mito"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
meta <- seu@meta.data
p1 <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
              group.by = "sample", pt.size = 0, ncol = 1) & NoLegend()
save_fig(p1, file.path(OUT, "figures/01_qc_violin.png"), 7, 10)

qc <- CFG$qc
qc_summary <- data.frame()
keep_cells <- character(0)
for (s in ids) {
  m <- meta[meta$sample == s, ]
  max_genes <- qc$max_genes
  # 自适应：基因集小的样本（normal-1）用 3×自身 nFeature 中位数作上界
  if (3 * median(m$nFeature_RNA) < max_genes) {
    max_genes <- 3 * median(m$nFeature_RNA)
    message(sprintf("  %-9s nFeature 上界自适应为 %.0f（3×中位数 %.0f）",
                    s, max_genes, median(m$nFeature_RNA)))
  }
  keep <- m$nFeature_RNA > qc$min_genes & m$nFeature_RNA < max_genes &
          m$nCount_RNA  > qc$min_counts & m$nCount_RNA < qc$max_counts &
          m$percent.mito < qc$max_mito
  keep_cells <- c(keep_cells, rownames(m)[keep])
  qc_summary <- rbind(qc_summary, data.frame(
    sample = s, before = nrow(m), after = sum(keep),
    pct_kept = round(100 * mean(keep), 1),
    med_nFeature = round(median(m$nFeature_RNA), 1),
    med_nCount = round(median(m$nCount_RNA), 1),
    med_mito = round(median(m$percent.mito), 2),
    max_genes_used = max_genes))
  message(sprintf("  %-9s %5d → %5d (%.1f%%)  med: %.0f genes, %.0f counts, %.2f%% mito",
                  s, nrow(m), sum(keep), 100 * mean(keep),
                  median(m$nFeature_RNA), median(m$nCount_RNA),
                  median(m$percent.mito)))
}
save_tsv(qc_summary, file.path(OUT, "tables/qc_summary.tsv"))

seu <- subset(seu, cells = keep_cells)
message("QC 后总细胞数: ", ncol(seu))

# ── 4. checkpoint ─────────────────────────────────────────────────────
saveRDS(seu, file.path(OUT, "checkpoints/02_qc_filtered.rds"))
message("✅ 40 完成 → ", file.path(OUT, "checkpoints/02_qc_filtered.rds"))
