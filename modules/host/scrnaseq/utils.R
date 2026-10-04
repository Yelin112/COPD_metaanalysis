# GSE211630 scRNA-seq 分析共享函数（40-44 脚本 source）
# 只收敛 4 个函数：读矩阵（格式坑全在此）、读配置、存表、存图

# 读 GSE211630 count 矩阵（格式坑实测：表头无前导空字段、barcode 带引号、
# 基因列为 symbol、末行无换行且被截断——GEO 源文件末尾缺 ~510 字段，
# fread 会丢弃该行，此处显式记录）→ dgCMatrix（基因×细胞）
read_gz_count_matrix <- function(file) {
  hdr  <- readLines(file, n = 1)
  barcodes <- gsub('"', "", strsplit(hdr, "\t")[[1]])
  dat  <- data.table::fread(file, skip = 1, header = FALSE, showProgress = FALSE,
                            fill = FALSE)
  genes <- gsub('"', "", dat[[1]])
  # 末行截断检查：正常行 = 1 基因 + N 细胞字段；fread 已把截断行当 footer 丢弃
  dropped <- length(readLines(file)) - 1 - nrow(dat)
  if (dropped > 0) message("  ⚠ 末行截断被丢弃 ", dropped, " 行（源文件缺尾，GEO 源头问题）")
  m <- Matrix::Matrix(as.matrix(dat[, -1]), sparse = TRUE)
  rownames(m) <- genes
  colnames(m) <- barcodes
  m
}

load_sc_config <- function(path) {
  yaml::read_yaml(path)
}

save_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("  表 → ", path)
}

sc_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = 8))
}

save_fig <- function(p, path, width = 8, height = 6) {
  # 用户偏好：研究绘图优先 PDF（矢量）；.png 路径自动换 .pdf，图内勿用中文。
  # 输出分目录（规范第3条）：figures/pdf/（矢量）+ figures/png/（Markdown 嵌入）
  stem <- sub("\\.png$", "", path)
  dir.create(file.path(dirname(stem), "pdf"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dirname(stem), "png"), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(dirname(stem), "pdf", paste0(basename(stem), ".pdf")),
                  p + sc_theme(), width = width, height = height, bg = "white")
  ggplot2::ggsave(file.path(dirname(stem), "png", paste0(basename(stem), ".png")),
                  p + sc_theme(), width = width, height = height, dpi = 300, bg = "white")
  message("  图 → ", basename(stem), "(.pdf/.png)")
}
