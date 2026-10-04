# OLP 宿主转录组 ComBat 批次校正 + MetaDE 荟萃分析
#
# 适配自 legacy/3_combat_metaDE.r（COPD 框架，25 数据集）
# 主要改动：
#   1. 直接读 scripts/10 产出的 normalized.tsv（跳过 legacy perl 1-2 步）
#   2. 参数化数据集列表（n=2，可扩展）
#   3. 修复 PCA ggplot 代码（原版 pca$tPC1 不存在）
#   4. MetaDE 参数全部由 n_studies 驱动，不再硬编码
#
# 用法：
#   /usr/lib/R/bin/Rscript scripts/31_combat_metaDE_host.R
#
# 输出（results/OLP/processed/host/metaDE/）：
#   all_imputed.txt  all_combat.txt  pca_combat.pdf
#   GSE52130.txt  GSE38616.txt      （MetaDE 输入）
#   metaFDR.txt  metaZval.txt  metaPval.txt  Meta.RDS

library(bnstruct)
library(sva)
library(MetaDE)
library(ggplot2)

# ── 参数（唯一需要改的地方）────────────────────────────────────────────
STUDY_IDS <- c("GSE52130", "GSE38616")
NORM_DIR  <- "results/OLP/processed/host"
META_DIR  <- "data/OLP/raw/host"
OUT_DIR   <- "results/OLP/processed/host/metaDE"
N         <- length(STUDY_IDS)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 1. 读取标准化矩阵 ──────────────────────────────────────────────────
# 格式：第一列 gene_id，其余列为样本（row.names=FALSE 写出的）
read_norm <- function(id) {
    f <- file.path(NORM_DIR, paste0(id, "_normalized.tsv"))
    d <- read.table(f, header=TRUE, check.names=FALSE, sep="\t",
                    stringsAsFactors=FALSE)
    rownames(d) <- d$gene_id
    d[, -1, drop=FALSE]
}

mats <- lapply(STUDY_IDS, read_norm)
names(mats) <- STUDY_IDS
message("读入矩阵：", paste(sapply(mats, ncol), "样本", collapse=" / "))

# ── 2. 读取元数据 ──────────────────────────────────────────────────────
# 格式：sample_id  condition（OLP/Control）
read_meta <- function(id) {
    f <- file.path(META_DIR, paste0(id, "_metadata.txt"))
    read.table(f, header=TRUE, check.names=FALSE, sep="\t",
               stringsAsFactors=FALSE)
}

metas <- lapply(STUDY_IDS, read_meta)
names(metas) <- STUDY_IDS

# ── 3. 对齐公共基因 ────────────────────────────────────────────────────
common_genes <- Reduce(intersect, lapply(mats, rownames))
message("公共基因数：", length(common_genes))
if (length(common_genes) < 100)
    stop("公共基因数过少（<100），请检查基因注释平台是否一致")

mats_aligned <- lapply(mats, function(m) m[common_genes, , drop=FALSE])

# ── 4. 拼合矩阵 + 元数据 ──────────────────────────────────────────────
all_mat <- do.call(cbind, mats_aligned)

all_meta <- do.call(rbind, lapply(STUDY_IDS, function(id) {
    samps <- colnames(mats_aligned[[id]])
    cond  <- metas[[id]]$condition[match(samps, metas[[id]]$sample_id)]
    data.frame(
        sample  = samps,
        batch   = id,
        disease = ifelse(cond == "OLP", 1, 0),
        row.names = samps,
        stringsAsFactors = FALSE
    )
}))
message("总样本数：", nrow(all_meta),
        "（OLP=", sum(all_meta$disease==1),
        " Control=", sum(all_meta$disease==0), "）")

# ── 5. knn 填补缺失值 ──────────────────────────────────────────────────
data_impute <- knn.impute(as.matrix(all_mat), k=10,
                          cat.var  = 1:ncol(all_mat),
                          to.impute= 1:nrow(all_mat),
                          using    = 1:nrow(all_mat))
write.table(data_impute,
            file  = file.path(OUT_DIR, "all_imputed.txt"),
            quote = FALSE, sep="\t")

# ── 6. ComBat 批次校正 ─────────────────────────────────────────────────
all_combat <- ComBat(data_impute,
                     batch      = all_meta$batch,
                     mod        = NULL,
                     par.prior  = TRUE,
                     prior.plots= FALSE)
write.table(all_combat,
            file  = file.path(OUT_DIR, "all_combat.txt"),
            quote = FALSE, sep="\t")
message("ComBat 完成")

# ── 7. PCA 图（校正前后对比）──────────────────────────────────────────
make_pca_df <- function(mat, label) {
    pca <- prcomp(t(mat), scale.=FALSE)
    var_explained <- round(summary(pca)$importance[2, 1:2] * 100, 1)
    df <- data.frame(
        PC1     = pca$x[, 1],
        PC2     = pca$x[, 2],
        batch   = all_meta[rownames(pca$x), "batch"],
        disease = ifelse(all_meta[rownames(pca$x), "disease"]==1, "OLP", "Control"),
        stage   = label,
        xlab    = paste0("PC1 (", var_explained[1], "%)"),
        ylab    = paste0("PC2 (", var_explained[2], "%)")
    )
    df
}

df_before <- make_pca_df(as.matrix(all_mat), "Before ComBat")
df_after  <- make_pca_df(all_combat,          "After ComBat")
df_plot   <- rbind(df_before, df_after)
df_plot$stage <- factor(df_plot$stage, levels=c("Before ComBat","After ComBat"))

p <- ggplot(df_plot, aes(PC1, PC2, color=batch, shape=disease)) +
    geom_point(size=2.5, alpha=0.8) +
    facet_wrap(~stage, scales="free") +
    scale_color_brewer(palette="Set1") +
    theme_bw(base_size=12) +
    labs(color="Study", shape="Condition")
ggsave(file.path(OUT_DIR, "pca_combat.pdf"), p, width=10, height=5)
message("PCA 图 → ", file.path(OUT_DIR, "pca_combat.pdf"))

# ── 8. 写 MetaDE 输入文件 ─────────────────────────────────────────────
# MetaDE.Read 格式（skip=1）：
#   第 1 行：空列名 + 各样本的分组标签（1=disease, 0=control）
#   后续行：基因 ID + 表达值
for (id in STUDY_IDS) {
    samps  <- colnames(mats_aligned[[id]])
    labels <- all_meta[samps, "disease"]
    label_row <- matrix(labels, nrow=1,
                        dimnames=list("", samps))
    mat_out <- rbind(label_row, mats_aligned[[id]][common_genes, samps])
    write.table(mat_out,
                file  = file.path(OUT_DIR, paste0(id, ".txt")),
                quote = FALSE, sep="\t", col.names=NA)
}
message("MetaDE 输入文件已写出（", N, " 个数据集）")

# ── 9. MetaDE 荟萃分析 ────────────────────────────────────────────────
old_wd <- getwd()
setwd(OUT_DIR)   # MetaDE.Read 从 CWD 读 {id}.txt

data.raw    <- MetaDE.Read(STUDY_IDS,
                           skip    = rep(1, N),
                           via     = "txt",
                           matched = TRUE,
                           log     = FALSE)
data.merged <- MetaDE.merge(data.raw)

MetaDE.res  <- MetaDE.rawdata(data.merged,
                              ind.method  = rep("modt", N),
                              paired      = rep(FALSE,  N),
                              meta.method = "REM",
                              nperm       = 300)

write.table(MetaDE.res$meta.analysis$FDR,
            file="metaFDR.txt",  quote=FALSE, sep="\t")
write.table(MetaDE.res$meta.analysis$zval,
            file="metaZval.txt", quote=FALSE, sep="\t")
write.table(MetaDE.res$meta.analysis$pval,
            file="metaPval.txt", quote=FALSE, sep="\t")
saveRDS(MetaDE.res, file="Meta.RDS")

setwd(old_wd)

# ── 10. 简单汇总 ──────────────────────────────────────────────────────
fdr <- read.table(file.path(OUT_DIR, "metaFDR.txt"),
                  header=TRUE, sep="\t")
n_sig <- sum(fdr[, 1] < 0.05, na.rm=TRUE)
message("\n✅ MetaDE 完成")
message("  显著基因（FDR < 0.05）：", n_sig, " / ", nrow(fdr))
message("  结果目录：", OUT_DIR)
