# 02_fgsea.R —— LONP1 共表达排序的 fgsea 富集（GO:BP、Reactome），原始 vs 偏相关对照（补充分析 2）
# 输入: results/LONP1_partial_spearman_genomewide.csv（01_partial_spearman.py 输出）
# 排序统计量: rho_raw、rho_partial_A（LONP1 CN + g CN）、rho_partial_B（A + 19p 臂 + 纯度）
# 基因集: msigdbr（Homo sapiens；C5 GO:BP，C2 CP:REACTOME），基因集大小 15–500；fgseaMultilevel，seed=1
# 用法: Rscript 02_fgsea.R <02_partial_correlation 目录>
suppressMessages({library(fgsea); library(msigdbr); library(data.table)})
args <- commandArgs(trailingOnly = TRUE); base <- args[1]
res <- file.path(base, "results")
d <- fread(file.path(res, "LONP1_partial_spearman_genomewide.csv"))
d <- d[!is.na(gene) & gene != ""]
get_sets <- function(cat, sub) {
  m <- tryCatch(msigdbr(species = "Homo sapiens", collection = cat, subcollection = sub),
                error = function(e) msigdbr(species = "Homo sapiens", category = cat, subcategory = sub))
  split(m$gene_symbol, m$gs_name)
}
sets <- list(GOBP = get_sets("C5", "GO:BP"), REACTOME = get_sets("C2", "CP:REACTOME"))
cat("gene sets:", sapply(sets, length), "\n")
out <- list(); lines <- c()
for (stat in c("rho_raw", "rho_partial_A", "rho_partial_B")) {
  x <- d[!is.na(get(stat))]
  r <- setNames(x[[stat]], x$gene); r <- r[!duplicated(names(r))]
  set.seed(1)
  for (lib in names(sets)) {
    f <- fgseaMultilevel(sets[[lib]], r, minSize = 15, maxSize = 500, eps = 0, nproc = 4)
    f$stat <- stat; f$lib <- lib
    f$leadingEdge <- vapply(f$leadingEdge, paste, "", collapse = ";")
    out[[paste(stat, lib)]] <- f
    top <- f[order(padj, -NES)][NES > 0][1:10]
    lines <- c(lines, sprintf("[%s | %s] n_sig(padj<0.05, NES>0)=%d, NES<0=%d", stat, lib,
                              sum(f$padj < 0.05 & f$NES > 0), sum(f$padj < 0.05 & f$NES < 0)))
    lines <- c(lines, sprintf("   %-70s NES=%.2f padj=%.1e size=%d", substr(top$pathway, 1, 70), top$NES, top$padj, top$size))
  }
}
all <- rbindlist(out)
fwrite(all[order(stat, lib, padj)], file.path(res, "fgsea_GOBP_Reactome_by_layer.csv"))
writeLines(lines, file.path(res, "02_fgsea_summary.txt")); cat(lines, sep = "\n")
