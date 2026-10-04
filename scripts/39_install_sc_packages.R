# GSE211630 scRNA-seq：安装 Seurat 全家桶（系统 R 4.5.1 → ~/R/library-system）
#
# 产出：
#   装好 Seurat / SeuratObject / harmony / future / SingleR / celldex
#   参考数据 BlueprintEncodeData 预热到 ~/.cache/R/ExperimentHub
#
# 用法（走 USTC 国内镜像，不设代理）：
#   /usr/lib/R/bin/Rscript scripts/39_install_sc_packages.R

options(repos = c(CRAN = "https://mirrors.ustc.edu.cn/CRAN"))
options(BioC_mirror = "https://mirrors.ustc.edu.cn/bioc")
options(Ncpus = 8)
options(timeout = 900)

message("== 1/3 CRAN：Seurat + harmony + future ==")
install.packages(c("Seurat", "harmony", "future", "future.apply"),
                 dependencies = TRUE, ask = FALSE)

message("== 2/3 Bioconductor：SingleR + celldex ==")
BiocManager::install(c("SingleR", "celldex"), ask = FALSE, update = FALSE)

message("== 3/3 预热 celldex 参考数据（几百 MB，经 ExperimentHub）==")
BlueprintEncodeData <- NULL  # 触发类名，避免 R CMD check 提示
ref <- celldex::BlueprintEncodeData()
message("  参考维度：", nrow(ref), " x ", ncol(ref))

for (p in c("Seurat", "harmony", "future", "SingleR", "celldex")) {
  message(sprintf("  %-10s %s", p,
                  if (requireNamespace(p, quietly = TRUE)) "OK" else "FAILED"))
}
