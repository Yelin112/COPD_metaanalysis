# 步骤44（GSE211630 scRNA-seq）：T 细胞定制通路活性（AddModuleScore）
#
# 方案：8 组 curated 基因集（TCR / IFN-γ / TNF-α / IL-17 / 细胞毒 / 耗竭 /
#       Trm / 增殖）→ 每细胞 AddModuleScore → 样本均值（n=6 视图）+
#       亚群均值（亚群×通路热图，核心结论图）。
#       检验：EOLP vs NEOLP 样本级 Welch t（2v3 探索性）；OLP vs Normal 只出
#       描述（均值+差值，n=1 不做检验）。
#
# 产出：
#   results/OLP/processed/host/GSE211630/tables/module_scores_sample.tsv
#   .../tables/module_scores_subtype.tsv
#   .../tables/module_stats.tsv
#   .../figures/14_module_violin_by_condition.png
#   .../figures/15_module_heatmap_subtype.png
#   .../figures/16_featureplot_trm_score.png
#
# 用法：
#   nohup /usr/lib/R/bin/Rscript scripts/44_module_scores_gse211630.R \
#     > logs/44_gse211630_module.log 2>&1 &

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(future)
})
plan("multicore", workers = 8)
options(future.globals.maxSize = 4 * 1024^3)

source("modules/host/scrnaseq/utils.R")
CFG <- load_sc_config("config/OLP/scrnaseq_gse211630.yaml")
OUT <- CFG$out_dir

seu_t <- readRDS(file.path(OUT, "checkpoints/04_tcell.rds"))

# ── 1. 模块评分（缺失基因自动忽略，normal-1 基因集小交集少几个）───────
sets <- CFG$module_sets
seu_t <- AddModuleScore(seu_t, features = unname(sets), name = "mod",
                        ctrl = 100, search = FALSE)
mod_cols <- paste0("mod", seq_along(sets))
names(mod_cols) <- names(sets)

meta <- seu_t@meta.data

# ── 2. 样本均值 + 统计 ────────────────────────────────────────────────
ids <- sapply(CFG$samples, `[[`, "id")
sample_means <- meta %>%
  group_by(sample) %>%
  summarise(across(all_of(mod_cols), mean), .groups = "drop") %>%
  as.data.frame()
colnames(sample_means)[-1] <- names(sets)
cond_map <- setNames(sapply(CFG$samples, `[[`, "condition"), ids)
grp_map  <- setNames(sapply(CFG$samples, `[[`, "group"), ids)
sample_means$condition <- unname(cond_map[sample_means$sample])
sample_means$group     <- unname(grp_map[sample_means$sample])
save_tsv(sample_means, file.path(OUT, "tables/module_scores_sample.tsv"))

stats <- data.frame(module = names(sets))
stats$welch_p_EOLP_vs_NEOLP <- sapply(names(sets), function(mod) {
  x <- sample_means[[mod]][sample_means$group == "EOLP"]
  y <- sample_means[[mod]][sample_means$group == "NEOLP"]
  if (length(x) < 2 | length(y) < 2) NA else t.test(x, y)$p.value
})
# 模块分是中心化分数，组间比较用差值而非比值
stats$diff_OLP_vs_Normal <- sapply(names(sets), function(mod) {
  mean(sample_means[[mod]][sample_means$condition == "OLP"]) -
    mean(sample_means[[mod]][sample_means$condition == "Normal"])
})
save_tsv(stats, file.path(OUT, "tables/module_stats.tsv"))

# ── 3. 图 14：8 通路 × 3 组（每样本一点）──────────────────────────────
long <- tidyr::pivot_longer(sample_means, cols = all_of(names(sets)),
                            names_to = "module", values_to = "score")
long$group <- factor(long$group, levels = c("Normal", "NEOLP", "EOLP"))
long$module <- factor(long$module, levels = names(sets))
p_viol <- ggplot(long, aes(group, score, color = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.12, size = 1.5) +
  facet_wrap(~module, scales = "free_y", ncol = 4) +
  labs(x = NULL, y = "Module score (sample mean)",
       caption = "Normal n=1 descriptive only; EOLP vs NEOLP in module_stats.tsv") +
  theme(legend.position = "none")
save_fig(p_viol, file.path(OUT, "figures/14_module_violin_by_condition.png"), 12, 7)

# ── 4. 图 15：亚群 × 通路 z-score 热图（核心结论图）───────────────────
subtype_means <- meta %>%
  group_by(subtype) %>%
  summarise(across(all_of(mod_cols), mean), .groups = "drop") %>%
  as.data.frame()
rownames(subtype_means) <- subtype_means$subtype
subtype_means <- subtype_means[, -1]
colnames(subtype_means) <- names(sets)
save_tsv(data.frame(subtype = rownames(subtype_means),
                    round(subtype_means, 4), check.names = FALSE),
         file.path(OUT, "tables/module_scores_subtype.tsv"))

# 注释条：每亚群在 EOLP/NEOLP/Normal 的细胞占比（样本均值未加权）
subtype_means <- subtype_means[order(-subtype_means$Trm), ]
ann <- meta %>%
  group_by(subtype, group) %>% summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = "group", values_from = "n", values_fill = 0)
ann_m <- as.data.frame(ann[, -1]); rownames(ann_m) <- ann$subtype
ann_m <- ann_m[rownames(subtype_means), ]
ann_m[] <- lapply(ann_m, function(x) sprintf("%d", x))

mat_z <- scale(subtype_means)   # 每通路（列）z-score
p_hm <- pheatmap(mat_z, cluster_cols = FALSE, cluster_rows = FALSE,
                 angle_col = 45, fontsize = 9,
                 annotation_row = ann_m,
                 color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
                 main = "T-cell subtype x pathway module scores (column z-score, rows sorted by Trm)",
                 filename = file.path(OUT, "figures/15_module_heatmap_subtype.pdf"),
                 width = 9, height = 7)
message("  热图 → figures/15_module_heatmap_subtype.png")

# ── 5. 图 16：Trm 评分 UMAP ───────────────────────────────────────────
trm_col <- mod_cols[["Trm"]]
seu_t$Trm_score <- meta[[trm_col]]
FeaturePlot(seu_t, features = "Trm_score") +
  scale_color_viridis_c() +
  ggplot2::ggtitle("Trm module score (ITGAE/CD69/CXCR6/ZNF683...)")
save_fig(last_plot(), file.path(OUT, "figures/16_featureplot_trm_score.png"), 8, 6)

message("✅ 44 完成")
