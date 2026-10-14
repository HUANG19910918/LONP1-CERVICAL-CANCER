# 00_inspect_GSE333892.R —— 读取 GSE333892 整体 Seurat 对象，检查是否含原始 counts、样本数与分组
# 用法: Rscript 00_inspect_GSE333892.R <rds 路径> <输出目录>
suppressMessages({library(Seurat); library(SeuratObject); library(Matrix)})
args <- commandArgs(trailingOnly = TRUE)
t0 <- Sys.time(); obj <- readRDS(args[1]); cat("read in", format(Sys.time() - t0), "\n")
cat("class:", class(obj), "\n")
o <- capture.output({
  print(obj)
  cat("assays:", Assays(obj), "\n")
  for (a in Assays(obj)) {
    cat("== assay", a, "class", class(obj[[a]]), "\n")
    ly <- tryCatch(Layers(obj[[a]]), error = function(e) slotNames(obj[[a]]))
    cat("layers:", ly, "\n")
    for (l in intersect(ly, c("counts", "data", grep("^counts", ly, value = TRUE)))) {
      m <- tryCatch(LayerData(obj, assay = a, layer = l), error = function(e) slot(obj[[a]], l))
      if (length(m) == 0 || nrow(m) == 0) { cat(l, ": empty\n"); next }
      sub <- m[, seq_len(min(2000, ncol(m))), drop = FALSE]
      v <- if (inherits(sub, "dgCMatrix")) sub@x else as.vector(sub)
      cat(sprintf("%s: %d x %d; nonzero sample min=%.3f max=%.3f; all integer=%s\n", l, nrow(m), ncol(m),
                  min(v), max(v), all(abs(v - round(v)) < 1e-8)))
    }
  }
  md <- obj@meta.data
  cat("meta columns:", colnames(md), "\n")
  for (cn in colnames(md)) if (!is.numeric(md[[cn]]) || length(unique(md[[cn]])) <= 60) {
    u <- table(md[[cn]], useNA = "ifany"); if (length(u) <= 80) { cat("--", cn, "--\n"); print(u) }
  }
  grep_hpv <- grep("HPV|hpv", rownames(obj), value = TRUE); cat("HPV features:", head(grep_hpv, 30), "\n")
  cat("LONP1 present:", "LONP1" %in% rownames(obj), "\n")
})
writeLines(o, file.path(args[2], "GSE333892_inspect.txt")); cat(o, sep = "\n")
saveRDS(obj@meta.data, file.path(args[2], "GSE333892_metadata.rds"))
