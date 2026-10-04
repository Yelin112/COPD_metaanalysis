# 步骤43（GSE211630 scRNA-seq）：细胞组成变化 + T 亚群 pseudobulk DESeq2
#
# 方案：
#   组成：每样本比例表 + 分层 bootstrap 95% CI（B=1000，按样本整群重抽样；
#         Normal n=1 时 CI 来自单样本重抽样，图注如实标注）。
#         OLP vs Normal 为描述性（n=1 局限），EOLP vs NEOLP 为探索性趋势。
#   差异表达：每 T 亚群 pseudobulk → DESeq2（OLP_vs_Normal 5v1 只作排序参考；
#         EOLP_vs_NEOLP 2v3 探索性）。亚群在某组 <2 样本或每样本 <30 细胞则跳过。
#
# 产出：
#   results/OLP/processed/host/GSE211630/tables/celltype_proportions.tsv
#   .../tables/de_{subtype}_{comparison}.tsv
#   .../figures/07_proportion_major_stacked.png
#   .../figures/08_proportion_olp_vs_normal.png
#   .../figures/11_tcell_proportion_by_condition.png
#   .../figures/13_de_volcano_cd8_effector_eolp_vs_neolp.png
#
# 用法：
#   nohup /usr/lib/R/bin/Rscript scripts/43_composition_deseq2_gse211630.R \
#     > logs/43_gse211630_composition.log 2>&1 &

suppressMessages({
  library(Seurat)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(future)
})
plan("multicore", workers = 8)
options(future.globals.maxSize = 4 * 1024^3)

source("modules/host/scrnaseq/utils.R")
CFG <- load_sc_config("config/OLP/scrnaseq_gse211630.yaml")
OUT <- CFG$out_dir

seu   <- readRDS(file.path(OUT, "checkpoints/03_clustered_annotated.rds"))
seu_t <- readRDS(file.path(OUT, "checkpoints/04_tcell.rds"))

ids <- sapply(CFG$samples, `[[`, "id")
grp_map <- setNames(sapply(CFG$samples, `[[`, "group"), ids)

# ── 1. 分层 bootstrap：组内细胞类型比例均值 ± 95% CI ────────────────────
# 按样本整群重抽样（样本→样本内细胞），Normal n=1 时退化为单样本内细胞重抽样
boot_props <- function(meta, label_col, group_col, B = 1000) {
  labels <- sort(unique(meta[[label_col]]))
  out <- list()
  for (g in unique(meta[[group_col]])) {
    m <- meta[meta[[group_col]] == g, ]
    samples <- unique(m$sample)
    mat <- matrix(NA_real_, B, length(labels), dimnames = list(NULL, labels))
    for (b in seq_len(B)) {
      ss <- sample(samples, length(samples), replace = TRUE)
      cells <- unlist(lapply(ss, function(s) {
        idx <- which(m$sample == s)
        idx[sample.int(length(idx), length(idx), replace = TRUE)]
      }))
      tb <- table(factor(m[[label_col]][cells], levels = labels))
      mat[b, ] <- tb / sum(tb)
    }
    ci <- t(apply(mat, 2, quantile, probs = c(0.025, 0.975)))
    out[[g]] <- data.frame(label = labels, group = g,
                           mean = colMeans(mat), lo = ci[, 1], hi = ci[, 2],
                           row.names = NULL)
  }
  do.call(rbind, out)
}

prop_table <- function(seu_obj, label_col) {
  p <- table(seu_obj@meta.data$sample, seu_obj@meta.data[[label_col]])
  p <- p / rowSums(p)
  data.frame(sample = rownames(p), as.data.frame.matrix(round(p, 4)), check.names = FALSE)
}

message("══ 组成分析：大类型 ══")
save_tsv(prop_table(seu, "curated_label"),
         file.path(OUT, "tables/celltype_proportions.tsv"))

p_stack <- ggplot(seu@meta.data, aes(x = sample, fill = curated_label)) +
  geom_bar(position = "fill") +
  labs(y = "Proportion", x = NULL, fill = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_stack, file.path(OUT, "figures/07_proportion_major_stacked.png"), 9, 5)

bp <- boot_props(seu@meta.data, "curated_label", "group")
bp$group <- factor(bp$group, levels = c("Normal", "NEOLP", "EOLP"))
bp_w <- tidyr::pivot_wider(bp[, c("label", "group", "mean")],
                           names_from = "group", values_from = "mean")
bp_w$log2_ratio_OLP_vs_Normal <-
  log2((bp_w$EOLP + bp_w$NEOLP + 1e-3) / (bp_w$Normal + 1e-3))
save_tsv(bp_w, file.path(OUT, "tables/celltype_group_proportions.tsv"))
bp$mean <- pmax(bp$mean, 1e-4)
p_mean <- ggplot(bp, aes(x = label, y = mean, color = group)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(0.5)) +
  scale_y_log10() + labs(x = NULL, y = "Proportion (mean \u00b1 95% CI, log10)",
                         caption = "Normal n=1: CI from within-sample resampling, reference only") +
  coord_flip()
save_fig(p_mean, file.path(OUT, "figures/08_proportion_olp_vs_normal.png"), 9, 7)

# ── 2. T 亚群比例（三组，含 CD8_Trm 正对照面板）───────────────────────
message("══ 组成分析：T 亚群 ══")
bp_t <- boot_props(seu_t@meta.data, "subtype", "group")
bp_t$group <- factor(bp_t$group, levels = c("Normal", "NEOLP", "EOLP"))
bp_t$mean <- pmax(bp_t$mean, 1e-4)
p_t <- ggplot(bp_t, aes(x = label, y = mean, color = group)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(0.5)) +
  scale_y_log10() + labs(x = NULL, y = "Proportion (mean \u00b1 95% CI, log10)",
                         caption = "Normal n=1: CI reference only; CD8_Trm = eLife positive control") +
  coord_flip()
save_fig(p_t, file.path(OUT, "figures/11_tcell_proportion_by_condition.png"), 10, 8)
message("  CD8_Trm 比例（%）: ",
        if ("CD8_Trm" %in% bp_t$label) {
          trm_rows <- bp_t[bp_t$label == "CD8_Trm", ]
          paste(trm_rows$group, round(100 * trm_rows$mean, 2),
                sep = "=", collapse = ", ")
        } else "（未产生 CD8_Trm 亚群，请检查 42 的 cluster scores）")

# ── 3. pseudobulk DESeq2（T 亚群 × 两个比较）────────────────────────────
run_deseq2 <- function(seu_t, st, cmp) {
  # cmp: list(id, col, groups=[[g1], [g2]])；col 指定用哪列元数据做分组
  # （OLP_vs_Normal 用 condition 列，EOLP_vs_NEOLP 用 group 列）
  g1 <- unlist(cmp$groups[[1]]); g2 <- unlist(cmp$groups[[2]])
  cmp_col <- seu_t@meta.data[[cmp$col]]
  m <- seu_t@meta.data[seu_t$subtype == st & cmp_col %in% c(g1, g2), ]
  if (nrow(m) == 0) {
    message(sprintf("  跳过 %s %s：该比较下无细胞", st, cmp$id))
    return(NULL)
  }
  m$cmp_grp <- ifelse(m[[cmp$col]] %in% g1, "g1", "g2")
  cnt <- table(m$sample)
  samp_cmp <- tapply(m$cmp_grp, m$sample, function(x) x[1])
  ok <- all(c("g1", "g2") %in% samp_cmp) && all(cnt >= 30)
  if (!ok) {
    message(sprintf("  跳过 %s %s：某组无样本或存在 <30 细胞/样本", st, cmp$id))
    return(NULL)
  }
  if (g2 == "Normal" & sum(samp_cmp == "g2") == 1) {
    message(sprintf("  ⚠ %s %s：Normal n=1，结果仅作排序参考", st, cmp$id))
  }
  sub <- subset(seu_t, cells = rownames(m))
  sub$cmp_grp <- ifelse(sub@meta.data[[cmp$col]] %in% g1, "g1", "g2")
  agg <- AggregateExpression(sub, group.by = "sample", return.seurat = FALSE)
  # 注意：v5 下显式传 layer 会与 PseudobulkExpression 形参冲突，默认层即 counts
  cts <- as.matrix(agg[[1]])
  cd <- data.frame(sample = colnames(cts), row.names = colnames(cts))
  cd$cmp_grp <- factor(tapply(sub$cmp_grp, sub$sample, function(x) x[1])[cd$sample],
                       levels = c("g2", "g1"))   # 基线=对照
  dds <- DESeqDataSetFromMatrix(cts, cd, design = ~cmp_grp)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("cmp_grp", "g1", "g2"))
  out <- as.data.frame(res)
  out$gene <- rownames(out)
  out <- out[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat",
                 "pvalue", "padj")]
  save_tsv(out, file.path(OUT, "tables",
                          paste0("de_", st, "_", cmp$id, ".tsv")))
  message(sprintf("  ✓ %s %s：%d 基因，padj<0.05 有 %d 个",
                  st, cmp$id, nrow(out), sum(out$padj < 0.05, na.rm = TRUE)))
  out
}

message("══ pseudobulk DESeq2 ══")
subtypes <- sort(unique(seu_t$subtype))
de_list <- list()
for (st in subtypes) {
  for (cmp in CFG$comparisons) {
    res <- run_deseq2(seu_t, st, cmp)
    de_list[[paste(st, cmp$id)]] <- res
  }
}

# ── 4. 火山图（CD8_Effector, EOLP vs NEOLP；缺失则换第一个可用亚群）────
vol <- de_list[["CD8_Effector EOLP_vs_NEOLP"]]
st_use <- "CD8_Effector"
if (is.null(vol)) {
  hits <- names(de_list)[vapply(de_list, Negate(is.null), logical(1))]
  hits <- hits[grepl("EOLP_vs_NEOLP", hits)]
  if (length(hits) == 0) {
    message("无可用 EOLP_vs_NEOLP DE 表，跳过火山图")
  } else {
    st_use <- strsplit(hits[1], " ")[[1]][1]
    vol <- de_list[[hits[1]]]
  }
}
if (!is.null(vol)) {
  vol <- vol %>%
    mutate(sig = ifelse(!is.na(padj) & padj < 0.05,
                        ifelse(log2FoldChange > 0, "up", "down"), "ns"))
  top <- vol %>% filter(padj < 0.05) %>%
    slice_max(abs(log2FoldChange), n = 10)
  p_vol <- ggplot(vol, aes(log2FoldChange, -log10(pvalue), color = sig)) +
    geom_point(size = 0.6) +
    scale_color_manual(values = c(up = "#d62728", down = "#1f77b4", ns = "grey80")) +
    ggrepel::geom_text_repel(data = top, aes(label = gene), size = 3,
                             max.overlaps = 20, show.legend = FALSE) +
    labs(title = paste(st_use, "EOLP vs NEOLP (exploratory 2v3)"),
         x = "log2FC", y = "-log10(p)") + theme(legend.position = "none")
  save_fig(p_vol,
           file.path(OUT, "figures/13_de_volcano_cd8_effector_eolp_vs_neolp.png"), 8, 6)
}

message("✅ 43 完成")
