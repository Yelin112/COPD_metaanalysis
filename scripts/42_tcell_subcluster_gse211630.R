# 步骤42（GSE211630 scRNA-seq）：T/NK 亚群细分 + 自动命名
#
# 方案：TNK subset → 专属 HVG/PCA/Harmony 重跑 → 聚类(0.8) →
#       11 组 marker AddModuleScore → 簇级规则命名
#       （NK / Proliferating / Treg / CD4_Naive|Th1|Th17 / CD8_Naive|Effector|Trm|Exhausted）
#       → dotplot 供人工核对（以 marker 为准，可改 tcell_cluster_scores.tsv 后手动定稿）
#
# 产出：
#   results/OLP/processed/host/GSE211630/checkpoints/04_tcell.rds
#   results/OLP/processed/host/GSE211630/tables/tcell_cluster_scores.tsv
#   results/OLP/processed/host/GSE211630/tables/tcell_subtype_proportions.tsv
#   .../figures/09_tcell_umap_subtypes.png
#   .../figures/10_dotplot_tcell_markers.png
#   .../figures/12_featureplot_key_genes.png
#
# 用法：
#   nohup /usr/lib/R/bin/Rscript scripts/42_tcell_subcluster_gse211630.R \
#     > logs/42_gse211630_tcell.log 2>&1 &

suppressMessages({
  library(Seurat)
  library(harmony)
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

seu <- readRDS(file.path(OUT, "checkpoints/03_clustered_annotated.rds"))

# ── 1. TNK subset ─────────────────────────────────────────────────────
seu_t <- subset(seu, curated_label == "TNK")
message("TNK 细胞数: ", ncol(seu_t))

# ── 2. subset 专属重跑（HVG 不继承全数据）──────────────────────────────
seu_t <- NormalizeData(seu_t, scale.factor = 10000)
seu_t <- FindVariableFeatures(seu_t, nfeatures = 2000, selection.method = "vst")
seu_t <- ScaleData(seu_t)
seu_t <- RunPCA(seu_t, npcs = 30)
seu_t <- RunHarmony(seu_t, group.by.vars = "sample", reduction.use = "pca")
seu_t <- FindNeighbors(seu_t, reduction = "harmony", dims = 1:20)
seu_t <- FindClusters(seu_t, resolution = 0.8)
seu_t <- RunUMAP(seu_t, reduction = "harmony", dims = 1:20)
message("T 亚簇数: ", length(unique(seu_t$seurat_clusters)))

# ── 3. marker 模块评分 + 簇级规则命名 ─────────────────────────────────
tk <- CFG$tcell_markers
seu_t <- AddModuleScore(seu_t, features = unname(tk), name = "tmod",
                        ctrl = 100, search = FALSE)
# 注意：AddModuleScore 的 features 顺序与 YAML tcell_markers 顺序一致
mod_cols <- paste0("tmod", seq_along(tk))
names(mod_cols) <- names(tk)

cl_scores <- seu_t@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(mod_cols), mean), .groups = "drop") %>%
  as.data.frame()
rownames(cl_scores) <- cl_scores$seurat_clusters
colnames(cl_scores)[-1] <- names(tk)
# ITGAE+ 比例（Trm 定义验证，eLife：ITGAE+CD69+）
ita_expr <- GetAssayData(seu_t, layer = "counts")["ITGAE", ]
ita_df <- data.frame(cl = seu_t$seurat_clusters,
                     ita_pos = as.numeric(ita_expr) > 0,
                     row.names = NULL)
ita <- ita_df %>%
  group_by(cl) %>%
  summarise(pct_ITGAE_pos = 100 * mean(ita_pos), .groups = "drop")
cl_scores$pct_ITGAE_pos <- ita$pct_ITGAE_pos[match(cl_scores$seurat_clusters, ita$cl)]

assign_subtype <- function(r) {
  states <- c("Naive", "Treg", "Th1", "Th17", "Cytotoxic", "Trm",
              "Exhausted", "Proliferating")
  # 1) 增殖优先（谱系标志常全阴性）
  if (which.max(r[states]) == which(states == "Proliferating")) return("Proliferating")
  # 2) NK 需分数明确为正且是谱系最大（避免全负簇误判）
  if (r[["NK"]] > 0.2 & r[["NK"]] >= r[["CD4"]] & r[["NK"]] >= r[["CD8"]]) return("NK")
  # 3) 谱系
  lineage <- if (r[["CD4"]] >= r[["CD8"]]) "CD4" else "CD8"
  # 4) Trm 覆盖：ITGAE+ >=20% 的簇定为组织驻留（eLife 定义为 ITGAE+CD69+；
  #    Trm 常与耗竭共表达，模块 argmax 会漏判，如簇15 ITGAE+ 52.9% 但 Exhausted 分更高）
  if (r[["pct_ITGAE_pos"]] >= 20) return(paste0(lineage, "_Trm"))
  # 5) 其余按谱系内状态 argmax
  if (lineage == "CD4") {
    st <- c("Naive", "Treg", "Th1", "Th17")[which.max(r[c("Naive", "Treg", "Th1", "Th17")])]
    return(if (st == "Treg") "Treg" else paste0("CD4_", st))
  }
  st <- c("Naive", "Cytotoxic", "Trm", "Exhausted")[
    which.max(r[c("Naive", "Cytotoxic", "Trm", "Exhausted")])]
  paste0("CD8_", if (st == "Cytotoxic") "Effector" else st)
}
cl_scores$subtype <- apply(cl_scores, 1, assign_subtype)
save_tsv(cl_scores, file.path(OUT, "tables/tcell_cluster_scores.tsv"))

cl2st <- setNames(cl_scores$subtype, cl_scores$seurat_clusters)
seu_t$subtype <- unname(cl2st[as.character(seu_t$seurat_clusters)])
message("亚群: ", paste(sort(unique(seu_t$subtype)), collapse = ", "))

# ── 4. 图：UMAP + dotplot 验证 + 关键基因 ──────────────────────────────
save_fig(DimPlot(seu_t, group.by = "subtype", label = TRUE, repel = TRUE),
         file.path(OUT, "figures/09_tcell_umap_subtypes.png"), 8, 6)

tk_vec <- unlist(tk)
tk_vec <- tk_vec[!duplicated(tk_vec)]   # GNLY 在 NK 和 Cytotoxic 中重复，DotPlot 会因重复 gene 报错
names(tk_vec) <- rep(names(tk), lengths(tk))[!duplicated(unlist(tk))]
DotPlot(seu_t, features = tk_vec, group.by = "subtype") + RotatedAxis() +
  ggplot2::ggtitle("T-cell subtype markers (see tcell_cluster_scores.tsv)")
save_fig(last_plot(), file.path(OUT, "figures/10_dotplot_tcell_markers.png"), 16, 8)

key <- c("CD8A", "ITGAE", "CXCR6", "GZMB", "PDCD1", "IFNG", "IL17A", "FOXP3")
FeaturePlot(seu_t, features = key, ncol = 4)
save_fig(last_plot(), file.path(OUT, "figures/12_featureplot_key_genes.png"), 16, 10)

# ── 5. 亚群比例（每样本）──────────────────────────────────────────────
prop <- table(seu_t$sample, seu_t$subtype)
prop <- prop / rowSums(prop)
save_tsv(data.frame(sample = rownames(prop),
                    as.data.frame.matrix(round(prop, 4)),
                    check.names = FALSE),
         file.path(OUT, "tables/tcell_subtype_proportions.tsv"))

# ── 6. checkpoint ─────────────────────────────────────────────────────
saveRDS(seu_t, file.path(OUT, "checkpoints/04_tcell.rds"))
message("✅ 42 完成 → ", file.path(OUT, "checkpoints/04_tcell.rds"))
