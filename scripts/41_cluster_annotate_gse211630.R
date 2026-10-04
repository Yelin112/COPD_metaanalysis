# 步骤41（GSE211630 scRNA-seq）：标准化 → 聚类（Harmony）→ SingleR 注释
#
# 方案：Normalize → HVG(2000) → PCA(50) → Harmony(仅按 sample 整合) →
#       聚类(0.8) → UMAP → SingleR(BlueprintEncodeData, 簇级) →
#       marker dotplot 验证 → curated 7 大类型
#
# 产出：
#   results/OLP/processed/host/GSE211630/checkpoints/03_clustered_annotated.rds
#   results/OLP/processed/host/GSE211630/tables/cluster_labels.tsv
#   results/OLP/processed/host/GSE211630/tables/cluster_markers.tsv
#   results/OLP/processed/host/GSE211630/figures/02_elbow.png
#   .../figures/03_umap_by_sample.png（批次混合检查）
#   .../figures/04_umap_by_cluster.png
#   .../figures/05_umap_singler_label.png
#   .../figures/06_dotplot_major_markers.png
#
# 用法：
#   nohup /usr/lib/R/bin/Rscript scripts/41_cluster_annotate_gse211630.R \
#     > logs/41_gse211630_cluster.log 2>&1 &

suppressMessages({
  library(Seurat)
  library(harmony)
  library(SingleR)
  library(celldex)
  library(dplyr)
  library(ggplot2)
  library(future)
})
plan("multicore", workers = 8)
options(future.globals.maxSize = 4 * 1024^3)
options(future.seed = TRUE)

source("modules/host/scrnaseq/utils.R")
CFG <- load_sc_config("config/OLP/scrnaseq_gse211630.yaml")
OUT <- CFG$out_dir
ig <- CFG$integration

seu <- readRDS(file.path(OUT, "checkpoints/02_qc_filtered.rds"))

# 检查点：聚类完成的未注释对象已存在则复用（Normalize→UMAP ~30 min）
pre_rds <- file.path(OUT, "checkpoints/03b_preannotation.rds")
if (file.exists(pre_rds)) {
  message("══ 复用已有聚类结果 ══ ", pre_rds)
  seu <- readRDS(pre_rds)
} else {

# ── 1. 标准化 + HVG + PCA ─────────────────────────────────────────────
seu <- NormalizeData(seu, scale.factor = 10000)
seu <- FindVariableFeatures(seu, nfeatures = CFG$hvg_n, selection.method = "vst")
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = ig$n_pcs)
save_fig(ElbowPlot(seu, ndims = ig$n_pcs), file.path(OUT, "figures/02_elbow.png"), 7, 5)

# ── 2. Harmony 整合（只按 sample，绝不按 condition）────────────────────
seu <- RunHarmony(seu, group.by.vars = ig$group, reduction.use = "pca",
                  reduction.save = "harmony")

# ── 3. 聚类 + UMAP ────────────────────────────────────────────────────
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:ig$dims_use)
seu <- FindClusters(seu, resolution = ig$cluster_res)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:ig$dims_use)
message("聚类数: ", length(unique(seu$seurat_clusters)))

save_fig(DimPlot(seu, group.by = "sample", shuffle = TRUE),
         file.path(OUT, "figures/03_umap_by_sample.png"), 8, 6)
save_fig(DimPlot(seu, group.by = "seurat_clusters", label = TRUE, repel = TRUE),
         file.path(OUT, "figures/04_umap_by_cluster.png"), 8, 6)

  saveRDS(seu, pre_rds)
  message("  预注释检查点 → ", pre_rds)
}

# ── 4. SingleR 簇级注释（新版 celldex rownames 即 gene symbol，无需映射）─
message("══ SingleR ══ 加载 BlueprintEncodeData")
ref <- celldex::BlueprintEncodeData()
message("  参考: ", nrow(ref), " symbol x ", ncol(ref), " 样本")

avg <- AverageExpression(seu, group.by = "seurat_clusters", layer = "data")
avg <- as.matrix(avg[[1]])
common <- intersect(rownames(avg), rownames(ref))
message("  基因交集: ", length(common), "/", nrow(avg))

pred_main <- SingleR::SingleR(test = avg[common, ], ref = ref[common, ],
                              labels = ref$label.main, de.method = "classic")
pred_fine <- SingleR::SingleR(test = avg[common, ], ref = ref[common, ],
                              labels = ref$label.fine, de.method = "classic")

# ── 5. SingleR 标签 → curated 7 大类型 ────────────────────────────────
map_curated <- function(x) {
  x <- tolower(x)
  dplyr::case_when(
    grepl("t-cell|t cell|nk", x) ~ "TNK",
    grepl("b-cell|b cell|plasma", x) ~ "B_Plasma",
    grepl("mono|macro|dendritic|dc|neutro|granulo", x) ~ "Myeloid",
    grepl("keratino|epithel", x) ~ "Epithelial",
    grepl("fibroblast", x) ~ "Fibroblast",
    grepl("endothel", x) ~ "Endothelial",
    grepl("smooth muscle", x) ~ "SmoothMuscle",
    TRUE ~ paste0("Other_", x))
}
labels_tab <- data.frame(
  cluster = sub("^g", "", rownames(pred_main)),   # Seurat 给数字簇名加的 g 前缀
  singler_main = pred_main$labels,
  singler_score = round(pred_main$scores, 2),     # classic 模式 scores 为数值向量
  singler_fine = pred_fine$labels,
  curated = map_curated(pred_main$labels),
  stringsAsFactors = FALSE)

# 手工定稿：SingleR 与 marker 矛盾时以 marker 为准（2026-09-22 按 dotplot+top markers 核对）
overrides <- c(
  "11" = "SmoothMuscle",     # CASQ2/ACTG2 血管平滑肌样（SingleR 误判 Adipocytes）
  "13" = "Mast",             # TPSAB1/TPSB2/CPA3/CMA1 肥大细胞（SingleR 误判 HSC）
  "22" = "SkeletalMuscle")   # TNNI2/TNNT3/DES 活检带入骨骼肌（SingleR 误判 Chondrocytes）
labels_tab$curated[labels_tab$cluster %in% names(overrides)] <-
  unname(overrides[labels_tab$cluster[labels_tab$cluster %in% names(overrides)]])
save_tsv(labels_tab, file.path(OUT, "tables/cluster_labels.tsv"))

# 簇 → curated 映射回每个细胞
cl2lab <- setNames(labels_tab$curated, labels_tab$cluster)
seu$curated_label <- unname(cl2lab[as.character(seu$seurat_clusters)])
message("curated 类型: ", paste(sort(unique(seu$curated_label)), collapse = ", "))
save_fig(DimPlot(seu, group.by = "curated_label", label = TRUE, repel = TRUE),
         file.path(OUT, "figures/05_umap_singler_label.png"), 8, 6)

# ── 6. marker dotplot 验证（以 marker 为准，人工核对 cluster_labels.tsv）──
mk <- CFG$major_markers
mk_vec <- unlist(mk)
names(mk_vec) <- rep(names(mk), lengths(mk))
DotPlot(seu, features = mk_vec, group.by = "seurat_clusters") +
  RotatedAxis() + ggplot2::ggtitle("Major cell-type markers (see cluster_labels.tsv)")
save_fig(last_plot(), file.path(OUT, "figures/06_dotplot_major_markers.png"), 14, 8)

# ── 7. 各簇 marker（验证表）───────────────────────────────────────────
all_markers <- FindAllMarkers(seu, only.pos = TRUE, logfc.threshold = 0.25,
                              min.pct = 0.1)
save_tsv(all_markers, file.path(OUT, "tables/cluster_markers.tsv"))

# ── 8. checkpoint ─────────────────────────────────────────────────────
saveRDS(seu, file.path(OUT, "checkpoints/03_clustered_annotated.rds"))
message("✅ 41 完成 → ", file.path(OUT, "checkpoints/03_clustered_annotated.rds"))
