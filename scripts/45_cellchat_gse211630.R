# 步骤45（GSE211630 scRNA-seq）：CellChat 三组细胞通讯比较
#
# 基础流程参考官方教程（~/Inference and analysis of cell-cell communication
# using CellChat.md：Seurat 输入 → CellChatDB.human → triMean → 通路级聚合），
# 三组比较用官方多数据集流程（mergeCellChat + 差异可视化）。
#
# 组：Normal（n=1，13.5k 细胞）/ NEOLP（n=3，25.7k → 抽样 20k）/ EOLP（n=2，18.2k）
# 分组注释用 curated_label（本数据含上皮/基质，T-上皮互作是 OLP 病理重点）。
#
# 产出（results/OLP/processed/host/GSE211630/）：
#   tables/45_cellchat_LR_delta_*.csv（三对两两比较的 L-R 差异表）
#   figures/pdf|png/45_cellchat_*.pdf/.png（官方可视化 + 差异热图）
#   checkpoints/cp_45_cellchat_list.rds（复用 4-6h 的 computeCommunProb）
#
# 用法：
#   nohup /usr/lib/R/bin/Rscript scripts/45_cellchat_gse211630.R \
#     > logs/45_cellchat_gse211630.log 2>&1 &

suppressMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
  library(pheatmap)
  library(future)
})
plan("multicore", workers = 8)
options(future.globals.maxSize = 32 * 1024^3)

source("modules/host/scrnaseq/utils.R")
OUT <- "results/OLP/processed/host/GSE211630"
TB_DIR <- file.path(OUT, "tables")
FG_DIR <- file.path(OUT, "figures")
CP_DIR <- file.path(OUT, "checkpoints")

seu <- readRDS(file.path(CP_DIR, "03_clustered_annotated.rds"))
seu <- JoinLayers(seu)
seu$group <- factor(seu$group, levels = c("Normal", "NEOLP", "EOLP"))
message("各组细胞数:")
print(table(seu$group))

# ── 1. 抽样（每组 ≤20k，按样本等比）────────────────────────────────────
set.seed(42)
CAP <- 20000
sample_cells <- function(seu_obj, grp) {
  s <- subset(seu_obj, group == grp)
  n <- ncol(s)
  if (n <= CAP) return(s)
  frac <- CAP / n
  keep <- unlist(lapply(unique(s$sample), function(sm) {
    idx <- which(s$sample == sm)
    idx[sample.int(length(idx), max(1, round(length(idx) * frac)))]
  }))
  subset(s, cells = keep)
}
seu_normal <- sample_cells(seu, "Normal")
seu_neolp  <- sample_cells(seu, "NEOLP")
seu_eolp   <- sample_cells(seu, "EOLP")
message("抽样后: Normal=", ncol(seu_normal), ", NEOLP=", ncol(seu_neolp),
        ", EOLP=", ncol(seu_eolp))

# ── 2. 组级 CellChat（官方教程基础流程）────────────────────────────────
run_cellchat <- function(seu_obj) {
  seu_obj$samples <- seu_obj$sample    # CellChat 要求 samples（复数）列
  cc <- createCellChat(seu_obj, group.by = "curated_label")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc, do.fast = FALSE)   # 未装 presto
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

cc_cp <- file.path(CP_DIR, "cp_45_cellchat_list.rds")
if (file.exists(cc_cp)) {
  message("══ 复用已有 CellChat 计算 ══ ", cc_cp)
  cc_list <- readRDS(cc_cp)
} else {
  cc_list <- list(Normal = run_cellchat(seu_normal),
                  NEOLP  = run_cellchat(seu_neolp),
                  EOLP   = run_cellchat(seu_eolp))
  saveRDS(cc_list, cc_cp)
}

merged <- mergeCellChat(cc_list, add.names = c("Normal", "NEOLP", "EOLP"))

# ── 3. 官方多组比较可视化（tryCatch：本机 CellChat 版本部分图函数有 bug）──
safe_plot <- function(expr, outname, w = 10, h = 8) {
  p <- tryCatch(expr, error = function(e) {
    message("  ", outname, " 跳过（", conditionMessage(e), "）"); NULL })
  if (is.null(p)) return(invisible(NULL))
  if (inherits(p, "ggplot")) {
    ggplot2::ggsave(file.path(FG_DIR, "pdf", paste0(outname, ".pdf")), p,
                    width = w, height = h)
    ggplot2::ggsave(file.path(FG_DIR, "png", paste0(outname, ".png")), p,
                    width = w, height = h, dpi = 300)
  } else if (inherits(p, "Heatmap")) {
    # ComplexHeatmap 对象（新版 CellChat 的 netVisual_heatmap 等）
    grDevices::pdf(file.path(FG_DIR, "pdf", paste0(outname, ".pdf")),
                   width = w, height = h)
    ComplexHeatmap::draw(p); grDevices::dev.off()
    grDevices::png(file.path(FG_DIR, "png", paste0(outname, ".png")),
                   width = w, height = h, units = "in", res = 300)
    ComplexHeatmap::draw(p); grDevices::dev.off()
  } else {
    message("  ", outname, " 跳过（未知图形对象类型）")
  }
  message("  图 → ", outname, "(.pdf/.png)")
}

# 互作总数比较
safe_plot(compareInteractions(merged, show.legend = FALSE, group = c(1, 2, 3)),
          "45_cellchat_compare_interactions", 7, 5)
# 总体差异热图（EOLP vs Normal / EOLP vs NEOLP）
safe_plot(netVisual_heatmap(merged, comparison = c(3, 1)), "45_cellchat_diffheat_EOLP_vs_Normal", 10, 8)
safe_plot(netVisual_heatmap(merged, comparison = c(3, 2)), "45_cellchat_diffheat_EOLP_vs_NEOLP", 10, 8)
safe_plot(netVisual_heatmap(merged, comparison = c(2, 1)), "45_cellchat_diffheat_NEOLP_vs_Normal", 10, 8)
# 信号强度排名
safe_plot(rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = TRUE),
          "45_cellchat_rankNet", 12, 8)

# ── 4. 稳健 L-R 差异表（直接从概率数组提取；两两比较）──────────────────
prob_of <- function(grp) merged@net[[grp]]$prob   # sources x targets x interactions

pair_delta <- function(g1, g2) {
  p1 <- prob_of(g1); p2 <- prob_of(g2)
  lr <- intersect(dimnames(p1)[[3]], dimnames(p2)[[3]])
  src <- intersect(dimnames(p1)[[1]], dimnames(p2)[[1]])
  tgt <- intersect(dimnames(p1)[[2]], dimnames(p2)[[2]])
  m1 <- apply(p1[src, tgt, lr], 3, mean, na.rm = TRUE)
  m2 <- apply(p2[src, tgt, lr], 3, mean, na.rm = TRUE)
  tab <- data.frame(interaction = lr, mean_g1 = m1, mean_g2 = m2,
                    delta = m1 - m2)
  tab[order(-abs(tab$delta)), ]
}

pairs <- list(c("EOLP", "Normal"), c("EOLP", "NEOLP"), c("NEOLP", "Normal"))
pair_heat <- list()
for (pr in pairs) {
  tab <- pair_delta(pr[1], pr[2])
  save_tsv(tab, file.path(TB_DIR, paste0("45_cellchat_LR_delta_",
                                         pr[1], "_vs_", pr[2], ".csv")))
  message("── ", pr[1], " vs ", pr[2], " Top 10 ──")
  print(head(tab, 10), row.names = FALSE)
  # 差异热图数据（亚群对 × top L-R）
  top20 <- head(tab$interaction, 20)
  p1 <- prob_of(pr[1]); p2 <- prob_of(pr[2])
  src <- intersect(dimnames(p1)[[1]], dimnames(p2)[[1]])
  tgt <- intersect(dimnames(p1)[[2]], dimnames(p2)[[2]])
  mat <- apply(p1[src, tgt, top20] - p2[src, tgt, top20], c(1, 2), mean,
               na.rm = TRUE)
  if (nrow(mat) > 1) {
    hm <- pheatmap(mat, cluster_cols = TRUE, cluster_rows = TRUE,
                   color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
                   main = paste(pr[1], "-", pr[2], "mean prob, top 20 L-R"))
    ggplot2::ggsave(file.path(FG_DIR, "pdf", paste0("45_cellchat_LRheat_",
          pr[1], "_vs_", pr[2], ".pdf")), hm, width = 11, height = 8)
    ggplot2::ggsave(file.path(FG_DIR, "png", paste0("45_cellchat_LRheat_",
          pr[1], "_vs_", pr[2], ".png")), hm, width = 11, height = 8, dpi = 300)
  }
}

# ── 5. TNK 参与的细胞对差异（细胞类型在 source/target 维度，不在互作名里）──
message("── TNK 参与的细胞对差异（平均概率差）──")
for (pr in pairs) {
  p1 <- prob_of(pr[1]); p2 <- prob_of(pr[2])
  src <- intersect(dimnames(p1)[[1]], dimnames(p2)[[1]])
  tgt <- intersect(dimnames(p1)[[2]], dimnames(p2)[[2]])
  lr  <- intersect(dimnames(p1)[[3]], dimnames(p2)[[3]])
  d <- apply(p1[src, tgt, lr] - p2[src, tgt, lr], c(1, 2), mean, na.rm = TRUE)
  message(pr[1], " vs ", pr[2], " —— TNK 作为 source（→各类型）:")
  print(round(sort(d["TNK", ], decreasing = TRUE), 4))
  message("—— TNK 作为 target（各类型→TNK）:")
  print(round(sort(d[, "TNK"], decreasing = TRUE), 4))
}

message("✅ 45 完成")
