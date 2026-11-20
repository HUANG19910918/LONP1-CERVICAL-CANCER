# 00_dump_genesets_all.R -- export GO:BP, Reactome and Hallmark gene sets (15-500 genes) for the generalization analysis
# usage: Rscript 00_dump_genesets_all.R genesets_all.gmt
suppressMessages(library(msigdbr))
args <- commandArgs(trailingOnly = TRUE); out <- if (length(args)) args[1] else "genesets_all.gmt"
get <- function(cat, sub) {
  m <- tryCatch(msigdbr(species = "Homo sapiens", collection = cat, subcollection = sub),
    error = function(e) {
      if (is.null(sub)) msigdbr(species = "Homo sapiens", category = cat)
      else msigdbr(species = "Homo sapiens", category = cat, subcategory = sub)
    })
  split(m$gene_symbol, m$gs_name)
}
hall <- tryCatch(get("H", NULL), error = function(e) list())
s <- c(get("C5", "GO:BP"), get("C2", "CP:REACTOME"), hall)
s <- lapply(s, unique)
n <- vapply(s, length, 0L)
s <- s[n >= 15 & n <= 500]
writeLines(vapply(names(s), function(k) paste(c(k, "msigdbr", s[[k]]), collapse = "\t"), ""), out)
cat("wrote", length(s), "gene sets to", out, "\n")
