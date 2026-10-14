# 00b_GSE333892_recover_counts.R —— GSE333892 的 counts 层与 data 层相同（log1p(CP10K) 小数），检验能否还原整数 counts
# 还原公式：counts = expm1(data) * nCount_RNA / 1e4；若还原值与最近整数之差 < 0.01 的比例接近 100%，视为可还原。
# 可还原时：逐样本（isolate）写出 10x 风格的 matrix.mtx.gz / features.tsv.gz / barcodes.tsv.gz 与 metadata.csv。
# 用法: Rscript 00b_GSE333892_recover_counts.R <rds> <输出目录> <结果目录>
suppressMessages({library(Seurat); library(SeuratObject); library(Matrix)})
args <- commandArgs(trailingOnly = TRUE); outd <- args[2]; resd <- args[3]
obj <- readRDS(args[1]); md <- obj@meta.data
X <- LayerData(obj, assay = "RNA", layer = "data"); rm(obj); gc()
set.seed(1); idx <- sample(ncol(X), 3000)
sub <- X[, idx]; nc <- md$nCount_RNA[idx]
rec <- t(t(expm1(sub)) * nc / 1e4)
v <- rec@x; frac <- mean(abs(v - round(v)) < 0.01)
csum <- Matrix::colSums(round(rec)); dsum <- abs(csum - nc) / nc
lines <- c(sprintf("cells tested=%d; nonzero entries=%d; fraction within 0.01 of integer=%.4f", length(idx), length(v), frac),
           sprintf("min recovered nonzero=%.4f; median=%.3f", min(v), median(v)),
           sprintf("colSums(round(recovered)) vs nCount_RNA: median rel diff=%.2e, max=%.2e", median(dsum), max(dsum)),
           sprintf("nCount_RNA integer=%s; nFeature_RNA matches nnz=%s", all(abs(md$nCount_RNA - round(md$nCount_RNA)) < 1e-6),
                   isTRUE(all.equal(as.numeric(diff(sub@p)), md$nFeature_RNA[idx]))))
cat(lines, sep = "\n"); writeLines(lines, file.path(resd, "GSE333892_count_recovery.txt"))
if (frac > 0.99) {
  dir.create(outd, showWarnings = FALSE, recursive = TRUE)
  for (s in sort(unique(md$isolate))) {
    j <- which(md$isolate == s); d <- file.path(outd, s); dir.create(d, showWarnings = FALSE)
    if (file.exists(file.path(d, "matrix.mtx.gz"))) next
    m <- X[, j]; m <- t(t(expm1(m)) * md$nCount_RNA[j] / 1e4); m@x <- round(m@x); m <- drop0(m)
    writeMM(m, file.path(d, "matrix.mtx")); R.utils::gzip(file.path(d, "matrix.mtx"), overwrite = TRUE)
    write.table(data.frame(rownames(m), rownames(m), "Gene Expression"), gzfile(file.path(d, "features.tsv.gz")),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    writeLines(colnames(m), gzfile(file.path(d, "barcodes.tsv.gz")))
    write.csv(md[j, ], file.path(d, "metadata.csv"))
    cat("wrote", s, length(j), "cells\n")
  }
}
cat("DONE\n")
