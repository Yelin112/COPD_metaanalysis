# 步骤46（GSE211630 scRNA-seq）：CellChat 1.5.0 三组比较（验证版）
#
# 背景：CellChat 2.x（GitHub main）的 netVisual_bubble 等绘图函数有 bug。
# 本脚本用 1.5.0（sqjin/CellChat tag v1.5.0，装于 ~/R/library-cellchat1）
# 按官方教程（Inference and analysis of cell-cell communication using
# CellChat.md）A 路线（matrix + meta）重跑三组（Normal/NEOLP/EOLP），
# 验证官方全套可视化是否可用。验证通过后将以 1.x 取代 2.x。
#
# 产出（results/OLP/processed/host/GSE211630/cellchat_v1/）：
#   tables/（显著 L-R 表）、figures/pdf|png/（官方图全套）
#   checkpoints/cp_46_cellchat_v1_list.rds（复用 computeCommunProb）
#
# 用法（必须用隔离库优先）：
#   nohup /usr/lib/R/bin/Rscript scripts/46_cellchat_v1_gse211630.R \
#     > logs/46_cellchat_v1_gse211630.log 2>&1 &

# 隔离库优先：CellChat 1.5.0 与 2.x 同名不能共存
.libPaths(c("/home/usr/yel/R/library-cellchat1", .libPaths()))
suppressMessages({
  library(CellChat)   # 应为 1.5.0
  library(Seurat)
  library(patchwork)
  library(future)
})
stopifnot(packageVersion("CellChat") == "1.5.0")
message("CellChat 版本: ", packageVersion("CellChat"))
plan("multicore", workers = 8)
options(future.globals.maxSize = 32 * 1024^3)

OUT <- "results/OLP/processed/host/GSE211630/cellchat_v1"
TB_DIR <- file.path(OUT, "tables")
FG_DIR <- file.path(OUT, "figures")
CP_DIR <- file.path(OUT, "checkpoints")
for (d in c(TB_DIR, file.path(FG_DIR, "pdf"), file.path(FG_DIR, "png"), CP_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
save_fig <- function(p, name, w = 10, h = 8) {   # 按对象类型分流保存
  if (inherits(p, "ggplot")) {
    ggplot2::ggsave(file.path(FG_DIR, "pdf", paste0(name, ".pdf")), p, width = w, height = h)
    ggplot2::ggsave(file.path(FG_DIR, "png", paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  } else if (inherits(p, "Heatmap")) {   # netVisual_heatmap 返回 ComplexHeatmap
    grDevices::pdf(file.path(FG_DIR, "pdf", paste0(name, ".pdf")), width = w, height = h)
    ComplexHeatmap::draw(p); grDevices::dev.off()
    grDevices::png(file.path(FG_DIR, "png", paste0(name, ".png")),
                   width = w, height = h, units = "in", res = 300)
    ComplexHeatmap::draw(p); grDevices::dev.off()
  } else {
    message("  ", name, " 跳过（未知对象类型）")
  }
  message("  图 → ", name)
}
save_tsv <- function(x, path) { write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE); message("  表 → ", path) }

seu <- readRDS("results/OLP/processed/host/GSE211630/checkpoints/03_clustered_annotated.rds")
seu$group <- factor(seu$group, levels = c("Normal", "NEOLP", "EOLP"))

# ── 1. 教程 A 路线：matrix + meta ──────────────────────────────────────
run_v1 <- function(seu_obj, grp) {
  s <- subset(seu_obj, group == grp)
  data.input <- as.matrix(GetAssayData(s, layer = "data"))   # 归一化数据
  meta <- data.frame(labels = s$curated_label, row.names = colnames(s))
  cc <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

cp <- file.path(CP_DIR, "cp_46_cellchat_v1_list.rds")
if (file.exists(cp)) {
  message("══ 复用 1.x 计算 ══")
  cc_list <- readRDS(cp)
} else {
  cc_list <- list(Normal = run_v1(seu, "Normal"),
                  NEOLP  = run_v1(seu, "NEOLP"),
                  EOLP   = run_v1(seu, "EOLP"))
  saveRDS(cc_list, cp)
}
# 角色分析前置：1.x 要求先算网络中心性
cc_list <- lapply(cc_list, function(cc)
  tryCatch(netAnalysis_computeCentrality(cc, slot.name = "netP"),
           error = function(e) { message("  centrality 跳过: ", conditionMessage(e)); cc }))
merged <- mergeCellChat(cc_list, add.names = c("Normal", "NEOLP", "EOLP"))
message("组: ", paste(names(merged@net), collapse = ", "))

# ── 2. 官方多组比较可视化（1.x 应全部可用）─────────────────────────────
gg1 <- compareInteractions(merged, show.legend = FALSE, group = c(1, 2, 3))
save_fig(gg1, "46_compare_interactions", 7, 5)

gg2 <- rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
save_fig(gg2, "46_rankNet", 12, 8)

gg3 <- netVisual_heatmap(merged, comparison = c(3, 1))
save_fig(gg3, "46_diffheat_EOLP_vs_Normal", 10, 8)
gg3b <- netVisual_heatmap(merged, comparison = c(3, 2))
save_fig(gg3b, "46_diffheat_EOLP_vs_NEOLP", 10, 8)

# 官方气泡图（2.x 中报 wrong sign 的那个；T 细胞相关对）
t_types <- "TNK"
gg4 <- netVisual_bubble(merged, sources.use = t_types, targets.use = 1:9,
                        comparison = c(1, 3), angle.x = 45)
save_fig(gg4, "46_bubble_TNK_source_Normal_vs_EOLP", 12, 10)
gg4b <- netVisual_bubble(merged, sources.use = 1:9, targets.use = t_types,
                         comparison = c(1, 3), angle.x = 45)
save_fig(gg4b, "46_bubble_TNK_target_Normal_vs_EOLP", 12, 10)

# 信号角色热图：1.x 的 merged 对象 centr 槽有结构问题 → 每组单独出
# （官方教程同款展示方式）；1.x 要求真实通路名（"ALL" 无效）
for (grp in names(cc_list)) {
  if (length(slot(cc_list[[grp]], "netP")$centr) == 0) {
    message("  ", grp, " 无 centr，跳过角色热图")
    next
  }
  tp <- head(cc_list[[grp]]@netP$pathways, 15)
  gg5 <- tryCatch(netAnalysis_signalingRole_heatmap(cc_list[[grp]],
                                                    pattern = "outgoing",
                                                    signaling = tp,
                                                    width = 10, height = 12),
                  error = function(e) { message("  outgoing 跳过: ", conditionMessage(e)); NULL })
  if (!is.null(gg5)) save_fig(gg5, paste0("46_role_heatmap_outgoing_", grp), 10, 12)
  gg5b <- tryCatch(netAnalysis_signalingRole_heatmap(cc_list[[grp]],
                                                     pattern = "incoming",
                                                     signaling = tp,
                                                     width = 10, height = 12),
                   error = function(e) { message("  incoming 跳过: ", conditionMessage(e)); NULL })
  if (!is.null(gg5b)) save_fig(gg5b, paste0("46_role_heatmap_incoming_", grp), 10, 12)
}

# 全局通讯模式（NMF，outgoing）
cc_eolp <- cc_list$EOLP
patterns <- tryCatch(
  identifyCommunicationPatterns(cc_eolp, pattern = "outgoing", k = 4, width = 10, height = 8),
  error = function(e) { message("  模式分析跳过: ", conditionMessage(e)); NULL })
if (!is.null(patterns)) {
  # 1.5.0 签名：pattern 是字符串（"outgoing"），模式槽已由 identify 存入对象
  gg6 <- tryCatch(netAnalysis_river(cc_eolp, pattern = "outgoing", font.size = 3),
                  error = function(e) { message("  river 跳过: ", conditionMessage(e)); NULL })
  if (!is.null(gg6)) save_fig(gg6, "46_river_outgoing_EOLP", 12, 10)
  gg7 <- tryCatch(netAnalysis_dot(cc_eolp, pattern = "outgoing"),
                  error = function(e) { message("  dot 跳过: ", conditionMessage(e)); NULL })
  if (!is.null(gg7)) save_fig(gg7, "46_dot_outgoing_EOLP", 12, 8)
}

# ── 3. 显著 L-R 表（subsetCommunication 对三数据集合并对象有问题 → 保护；
#    两两差异表才是科学产出）──────────────────────────────────────────
df_net <- tryCatch(subsetCommunication(merged, slot.name = "net"),
                   error = function(e) { message("  subsetCommunication 跳过: ",
                                                 conditionMessage(e)); NULL })
if (!is.null(df_net)) {
  tryCatch({
    save_tsv(df_net, file.path(TB_DIR, "46_cellchat_net_interactions.csv"))
    message("net 表行数: ", nrow(df_net))
  }, error = function(e) message("  net 表保存跳过: ", conditionMessage(e)))
}

# 差异 L-R（两两比较，用 1.x 惯例）
for (pr in list(c("EOLP", "Normal"), c("EOLP", "NEOLP"), c("NEOLP", "Normal"))) {
  tab <- tryCatch({
    p1 <- merged@net[[pr[1]]]$prob; p2 <- merged@net[[pr[2]]]$prob
    lr <- unique(intersect(dimnames(p1)[[3]], dimnames(p2)[[3]]))  # 去重防索引扩展
    src <- intersect(dimnames(p1)[[1]], dimnames(p2)[[1]])
    tgt <- intersect(dimnames(p1)[[2]], dimnames(p2)[[2]])
    d <- data.frame(interaction = lr,
                    g1 = apply(p1[src, tgt, lr], 3, mean, na.rm = TRUE),
                    g2 = apply(p2[src, tgt, lr], 3, mean, na.rm = TRUE))
    d$delta <- d$g1 - d$g2
    d[order(-abs(d$delta)), ]
  }, error = function(e) { message("  ", pr[1], " vs ", pr[2], " 跳过: ", conditionMessage(e)); NULL })
  if (!is.null(tab)) save_tsv(tab, file.path(TB_DIR, paste0("46_cellchat_LR_delta_",
                                                            pr[1], "_vs_", pr[2], ".csv")))
}

message("✅ 46（1.x 验证）完成")
