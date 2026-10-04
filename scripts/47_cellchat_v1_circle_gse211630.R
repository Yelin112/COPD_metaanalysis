# 47：GSE211630 CellChat 1.x 环状图（基于 46 的 checkpoint，不重跑计算）
# 与 21 样本项目的 13 脚本同款修补：igraph 2.x 兼容 + attach 环境同步 + fun() 包装
.libPaths(c("/home/usr/yel/R/library-cellchat1", .libPaths()))
suppressMessages(library(CellChat))
stopifnot(packageVersion("CellChat") == "1.5.0")

patch_circle <- function() {
  ns <- asNamespace("CellChat")
  src <- paste(deparse(get("netVisual_circle", envir = ns)), collapse = "\n")
  src <- gsub("vertex\\.label\\.degree = label\\.locs,", "", src)
  src <- gsub("(if \\(sum\\(edge\\.start\\[, 2\\] == edge\\.start\\[, 1\\]\\) != 0\\) \\{)",
              "igraph::E(g)$loop.angle <- 0\n    \\1", src)
  nf <- eval(parse(text = src))
  environment(nf) <- ns
  unlockBinding("netVisual_circle", ns)
  assign("netVisual_circle", nf, envir = ns)
  if ("package:CellChat" %in% search()) {
    ae <- as.environment("package:CellChat")
    unlockBinding("netVisual_circle", ae)
    assign("netVisual_circle", nf, envir = ae)
  }
  invisible(TRUE)
}
patch_circle()

OUT <- "results/OLP/processed/host/GSE211630/cellchat_v1"
FG_PDF <- file.path(OUT, "figures", "pdf")
FG_PNG <- file.path(OUT, "figures", "png")
cc_list <- readRDS(file.path(OUT, "checkpoints", "cp_46_cellchat_v1_list.rds"))

# 注：expr 若为 promise，第二次求值不重绘 → cairo 空设备不落盘；必须用 function() 包装
save_base <- function(fun, name, w = 9, h = 9) {
  fp <- file.path(FG_PDF, paste0(name, ".pdf"))
  pn <- file.path(FG_PNG, paste0(name, ".png"))
  ok_pdf <- tryCatch({
    grDevices::pdf(fp, w, h); fun(); grDevices::dev.off(); TRUE
  }, error = function(e) {
    try(grDevices::dev.off(), silent = TRUE)
    if (file.exists(fp)) file.remove(fp)
    message("  ", name, " 跳过: ", conditionMessage(e)); FALSE
  })
  if (!ok_pdf) return(invisible(FALSE))
  ok_png <- tryCatch({
    grDevices::png(pn, w, h, units = "in", res = 300); fun(); grDevices::dev.off(); TRUE
  }, error = function(e) {
    try(grDevices::dev.off(), silent = TRUE)
    if (file.exists(pn)) file.remove(pn)
    message("  ", name, " png 跳过: ", conditionMessage(e)); FALSE
  })
  if (ok_png) message("  图 → ", name)
  invisible(ok_png)
}

# 1) 每组总网络环状图（数量 + 强度）
for (grp in names(cc_list)) {
  cc <- cc_list[[grp]]
  vw <- as.numeric(table(cc@idents))
  save_base(function() netVisual_circle(cc@net$count, vertex.weight = vw, weight.scale = TRUE,
                             label.edge = FALSE, title.name = paste0(grp, " (number of interactions)")),
            paste0("47_circle_count_", grp))
  save_base(function() netVisual_circle(cc@net$weight, vertex.weight = vw, weight.scale = TRUE,
                             label.edge = FALSE, title.name = paste0(grp, " (interaction strength)")),
            paste0("47_circle_weight_", grp))
}

# 2) 重点通路环状图（与 46 结果呼应：IFN-II/IL17/MIF/MHC/CCL，取与 netP 的交集）
avail <- unique(unlist(lapply(cc_list, function(cc) cc@netP$pathways)))
pws <- intersect(c("MIF", "MHC-I", "MHC-II", "CCL", "CLEC", "CD22"), avail)
message("通路: ", paste(pws, collapse = ", "))
for (pw in pws) {
  for (grp in names(cc_list)) {
    cc <- cc_list[[grp]]
    vw <- as.numeric(table(cc@idents))
    tryCatch({
      save_base(function() netVisual_aggregate(cc, signaling = pw, layout = "circle",
                                    vertex.weight = vw, label.edge = FALSE),
                paste0("47_circle_", pw, "_", grp))
    }, error = function(e) message("  ", pw, "/", grp, " 跳过: ", conditionMessage(e)))
  }
}
message("✅ 47（GSE211630 环状图）完成")
