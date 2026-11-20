# 00_dump_genesets.R -- export the MSigDB gene sets used for the control-gene analysis (same msigdbr source as 02_fgsea.R)
suppressMessages(library(msigdbr))
args <- commandArgs(trailingOnly = TRUE); out <- if (length(args)) args[1] else "genesets_selected.gmt"
get <- function(cat, sub) { m <- tryCatch(msigdbr(species = "Homo sapiens", collection = cat, subcollection = sub),
  error = function(e) msigdbr(species = "Homo sapiens", category = cat, subcategory = sub)); split(m$gene_symbol, m$gs_name) }
s <- c(get("C5", "GO:BP"), get("C2", "CP:REACTOME"), get("C5", "GO:CC"))
keep <- c("GOBP_MITOCHONDRIAL_GENE_EXPRESSION", "GOBP_MITOCHONDRIAL_TRANSLATION", "GOBP_OXIDATIVE_PHOSPHORYLATION",
  "GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY", "GOBP_PROTEIN_FOLDING", "REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION",
  "REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT", "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT", "REACTOME_TRANSLATION",
  "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY", "GOBP_CYTOPLASMIC_TRANSLATION", "GOBP_RIBOSOME_BIOGENESIS", "GOCC_MITOCHONDRION")
keep <- keep[keep %in% names(s)]
writeLines(vapply(keep, function(k) paste(c(k, "msigdbr", unique(s[[k]])), collapse = "\t"), ""), out)
cat("wrote", length(keep), "gene sets to", out, "\n")
