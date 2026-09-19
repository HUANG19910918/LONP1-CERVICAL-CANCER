# 01_verify_stats.R —— 用R复算实验记录中R1/R2/R3/R6/R7/R8全部关键统计量
# 输出 raw_data/verification_report.txt，与实验记录.md中数值逐项核对
suppressMessages({library(survival)})
args <- commandArgs(trailingOnly=TRUE); base <- args[1]
rd <- file.path(base, "raw_data")
sink(file.path(rd, "verification_report.txt"))
cat("=== R verification report,", format(Sys.time()), R.version.string, "===\n\n")

# R1: TCGA+GTEx DE
d1 <- read.csv(file.path(rd, "01_TCGA_GTEx_LONP1_expression_tumor_vs_normal.csv"))
t1 <- d1$LONP1_log2TPM[d1$group=="Tumor"]; n1 <- d1$LONP1_log2TPM[d1$group=="Normal"]
tt <- t.test(t1, n1); wt <- wilcox.test(t1, n1)
cat(sprintf("R1 TCGA+GTEx: nT=%d nN=%d meanT=%.3f meanN=%.3f log2FC=%.3f\n   Welch t=%.3f df=%.1f p=%.3e | MWU W=%.0f p=%.3e\n\n",
    length(t1), length(n1), mean(t1), mean(n1), mean(t1)-mean(n1), tt$statistic, tt$parameter, tt$p.value, wt$statistic, wt$p.value))

# R2: survival
d2 <- read.csv(file.path(rd, "02_TCGA_CESC_clinical_LONP1_survival.csv"))
for(ep in list(c("OS_months","OS_event","OS"), c("DSS_months","DSS_event","DSS"), c("PFS_months","PFS_event","PFS"))){
  x <- d2[!is.na(d2[[ep[1]]]) & !is.na(d2[[ep[2]]]) & d2[[ep[1]]]>0,]
  cx <- coxph(Surv(x[[ep[1]]], x[[ep[2]]]) ~ LONP1_log2RSEM, data=x)
  s <- summary(cx)
  cat(sprintf("R2 %s: n=%d ev=%d HR=%.3f (%.3f-%.3f) p=%.4f\n", ep[3], s$n, s$nevent,
      s$conf.int[1], s$conf.int[3], s$conf.int[4], s$coefficients[5]))
}
x <- d2[!is.na(d2$OS_months) & !is.na(d2$OS_event) & d2$OS_months>0,]
x$grp <- ifelse(x$LONP1_log2RSEM > median(x$LONP1_log2RSEM), "High", "Low")
sd1 <- survdiff(Surv(OS_months, OS_event) ~ grp, data=x)
cat(sprintf("R2 OS median-split logrank: chi2=%.3f p=%.4f (cutoff=%.3f, High=%d Low=%d)\n\n",
    sd1$chisq, 1-pchisq(sd1$chisq,1), median(x$LONP1_log2RSEM), sum(x$grp=="High"), sum(x$grp=="Low")))

# R3: mito panel Spearman
d8 <- read.csv(file.path(rd, "08_CESC_gene_panel_expression_RSEM.csv"))
lg <- log2(d8$LONP1+1)
mito <- c("CLPP","HSPD1","FIS1","HSPA9","NRF1","AFG3L2","MFF","MIEF2","MIEF1","PRKN","SPG7","DNM1L","PPARGC1A","PINK1","TFAM","YME1L1","MFN2","MFN1","OPA1")
cat("R3 mito panel Spearman (log2(RSEM+1)):\n")
for(g in mito){
  v <- log2(d8[[g]]+1); ct <- suppressWarnings(cor.test(lg, v, method="spearman"))
  cat(sprintf("  %s rho=%.3f p=%.2e\n", g, ct$estimate, ct$p.value))
}
cat("\n")

# R6: CNA
d7 <- read.csv(file.path(rd, "07_CESC_LONP1_CNA.csv"))
m <- merge(d7, data.frame(sample=d8$sample, expr=lg))
ct <- suppressWarnings(cor.test(m$log2CNA, m$expr, method="spearman"))
cat(sprintf("R6 CNA~expr: n=%d Spearman rho=%.3f p=%.2e\n", sum(complete.cases(m[,c("log2CNA","expr")])), ct$estimate, ct$p.value))
print(aggregate(expr~gistic, m, function(v) c(n=length(v), mean=round(mean(v),3))))
cat("\n")

# R7: subtype
d10 <- read.csv(file.path(rd, "10_CESC_subtype.csv"))
d8$patient <- substr(d8$sample, 1, 12)
ms <- merge(d8, d10)
a <- log2(ms$LONP1[ms$subtype=="CESC_AdenoCarcinoma"]+1); s2 <- log2(ms$LONP1[ms$subtype=="CESC_SquamousCarcinoma"]+1)
tt2 <- t.test(a, s2)
cat(sprintf("R7 Adeno(n=%d, %.3f±%.3f) vs Squamous(n=%d, %.3f±%.3f): Welch p=%.4f\n\n",
    length(a), mean(a), sd(a), length(s2), mean(s2), sd(s2), tt2$p.value))

# R8: GEO
g63 <- read.csv(file.path(rd, "03_GSE63514_LONP1_209017_s_at.csv"))
w63 <- wilcox.test(g63$value[g63$group=="Cancer"], g63$value[g63$group=="Normal"])
cat(sprintf("R8 GSE63514 Cancer(n=%d) vs Normal(n=%d): MWU p=%.4f\n", sum(g63$group=="Cancer"), sum(g63$group=="Normal"), w63$p.value))
g63$grpF <- factor(g63$group, levels=c("Normal","CIN1","CIN2","CIN3","Cancer"))
kw <- kruskal.test(g63$value_log2norm, g63$grpF)
cat(sprintf("   GSE63514 Kruskal-Wallis across 5 groups: p=%.4f\n", kw$p.value))
g97 <- read.csv(file.path(rd, "04_GSE9750_LONP1_209017_s_at.csv"))
w97 <- wilcox.test(g97$value_MAS5[g97$group=="cancer"], g97$value_MAS5[g97$group=="normal"])
cat(sprintf("R8 GSE9750 cancer(n=%d, mean=%.1f) vs normal(n=%d, mean=%.1f): MWU p=%.3e\n",
    sum(g97$group=="cancer"), mean(g97$value_MAS5[g97$group=="cancer"]), sum(g97$group=="normal"), mean(g97$value_MAS5[g97$group=="normal"]), w97$p.value))
g78 <- read.csv(file.path(rd, "05_GSE7803_LONP1_209017_s_at.csv"))
w78 <- wilcox.test(g78$value_log2[g78$group=="scc"], g78$value_log2[g78$group=="normal"])
cat(sprintf("R8 GSE7803 SCC(n=%d) vs normal(n=%d): MWU p=%.4f\n", sum(g78$group=="scc"), sum(g78$group=="normal"), w78$p.value))
sink()
cat("VERIFY_DONE\n")
