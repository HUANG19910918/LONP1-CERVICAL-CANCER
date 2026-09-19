# 02_fig_expression.R —— Figure 1 (表达) + Figure 2 (调控/亚型)
# 输入: raw_data/01,03,04,05,07,08,10  输出: figures/Fig1A-D, Fig2A-C (pdf+png)
suppressMessages({library(ggplot2); library(scales)})
args <- commandArgs(trailingOnly=TRUE); base <- args[1]
rd <- file.path(base,"raw_data"); fig <- file.path(base,"figures","子图原文件"); dir.create(fig, showWarnings=FALSE, recursive=TRUE)
dir.create(fig, showWarnings=FALSE)
COL <- c(Normal="#4DBBD5", Tumor="#E64B35", Cancer="#E64B35",
         CIN1="#91D1C2", CIN2="#F39B7F", CIN3="#8491B4",
         normal="#4DBBD5", cancer="#E64B35", cellline="#3C5488",
         hsil="#F39B7F", scc="#E64B35")
th <- theme_classic(base_size=13) + theme(
  axis.text=element_text(color="black", size=11),
  axis.title=element_text(size=12, face="bold"),
  plot.title=element_text(size=12, face="bold", hjust=0.5),
  legend.position="none")
sv <- function(p, name, w, h){
  ggsave(file.path(fig, paste0(name,".pdf")), p, width=w, height=h)
  ggsave(file.path(fig, paste0(name,".png")), p, width=w, height=h, dpi=300)
}
pfmt <- function(p) ifelse(p<0.001, sprintf("p = %.2e", p), sprintf("p = %.3f", p))
box_jitter <- function(d, xcol, ycol, ylab, title, order, pv, plab_y){
  d$x <- factor(d[[xcol]], levels=order)
  ggplot(d, aes(x=x, y=.data[[ycol]], fill=x)) +
    geom_boxplot(width=0.55, outlier.shape=NA, alpha=0.85, linewidth=0.5) +
    geom_jitter(width=0.15, size=1.1, alpha=0.45, shape=16) +
    scale_fill_manual(values=COL) +
    labs(x=NULL, y=ylab, title=title) + th +
    annotate("text", x=1.5, y=plab_y, label=pv, size=3.8)
}

# Fig1A TCGA+GTEx
d1 <- read.csv(file.path(rd,"01_TCGA_GTEx_LONP1_expression_tumor_vs_normal.csv"))
w <- wilcox.test(d1$LONP1_log2TPM[d1$group=="Tumor"], d1$LONP1_log2TPM[d1$group=="Normal"])
p1a <- ggplot(d1, aes(x=factor(group, levels=c("Normal","Tumor")), y=LONP1_log2TPM, fill=group)) +
  geom_violin(trim=FALSE, alpha=0.6, linewidth=0.4) +
  geom_boxplot(width=0.18, outlier.shape=NA, alpha=0.9, linewidth=0.5) +
  scale_fill_manual(values=COL) +
  labs(x=NULL, y=expression(bold("LONP1 expression, log"[2]*"(TPM+0.001)")),
       title="TCGA-CESC + GTEx") + th +
  scale_x_discrete(labels=c(sprintf("Normal\n(n=%d)", sum(d1$group=="Normal")),
                            sprintf("Tumor\n(n=%d)", sum(d1$group=="Tumor")))) +
  annotate("text", x=1.5, y=max(d1$LONP1_log2TPM)+0.25, label=pfmt(w$p.value), size=3.8)
sv(p1a, "Fig1A_TCGA_GTEx_expression", 3.2, 3.8)

# Fig1B GSE9750 (tissue only for stats; MAS5值取log2显示)
d4 <- read.csv(file.path(rd,"04_GSE9750_LONP1_209017_s_at.csv"))
d4t <- d4[d4$group!="cellline",]; d4t$l2 <- log2(d4t$value_MAS5)
w4 <- wilcox.test(d4t$l2[d4t$group=="cancer"], d4t$l2[d4t$group=="normal"])
p1b <- box_jitter(d4t, "group", "l2", expression(bold("LONP1, log"[2]*"(MAS5 signal)")),
  "GSE9750", c("normal","cancer"), pfmt(w4$p.value), max(d4t$l2)+0.3)
p1b <- p1b + scale_x_discrete(labels=c(sprintf("Normal\n(n=%d)",sum(d4t$group=="normal")),
                                       sprintf("Cancer\n(n=%d)",sum(d4t$group=="cancer"))))
sv(p1b, "Fig1B_GSE9750_expression", 3.0, 3.8)

# Fig1C GSE7803
d5 <- read.csv(file.path(rd,"05_GSE7803_LONP1_209017_s_at.csv"))
w5 <- wilcox.test(d5$value_log2[d5$group=="scc"], d5$value_log2[d5$group=="normal"])
d5$x <- factor(d5$group, levels=c("normal","hsil","scc"))
p1c <- ggplot(d5, aes(x=x, y=value_log2, fill=x)) +
  geom_boxplot(width=0.55, outlier.shape=NA, alpha=0.85, linewidth=0.5) +
  geom_jitter(width=0.15, size=1.1, alpha=0.45, shape=16) +
  scale_fill_manual(values=COL) +
  scale_x_discrete(labels=c(sprintf("Normal\n(n=%d)",sum(d5$group=="normal")),
                            sprintf("HSIL\n(n=%d)",sum(d5$group=="hsil")),
                            sprintf("SCC\n(n=%d)",sum(d5$group=="scc")))) +
  labs(x=NULL, y=expression(bold("LONP1, log"[2]*" expression")), title="GSE7803") + th +
  annotate("text", x=2, y=max(d5$value_log2)+0.15,
           label=paste0("SCC vs Normal ", pfmt(w5$p.value)), size=3.5)
sv(p1c, "Fig1C_GSE7803_expression", 3.4, 3.8)

# Fig1D GSE63514 (阴性结果，如实展示)
d3 <- read.csv(file.path(rd,"03_GSE63514_LONP1_209017_s_at.csv"))
w3 <- wilcox.test(d3$value_log2norm[d3$group=="Cancer"], d3$value_log2norm[d3$group=="Normal"])
d3$x <- factor(d3$group, levels=c("Normal","CIN1","CIN2","CIN3","Cancer"))
p1d <- ggplot(d3, aes(x=x, y=value_log2norm, fill=x)) +
  geom_boxplot(width=0.6, outlier.shape=NA, alpha=0.85, linewidth=0.5) +
  geom_jitter(width=0.15, size=1.0, alpha=0.45, shape=16) +
  scale_fill_manual(values=COL) +
  labs(x=NULL, y=expression(bold("LONP1, log"[2]*" expression")),
       title="GSE63514 (LCM epithelium)") + th +
  annotate("text", x=3, y=max(d3$value_log2norm)+0.4,
           label=paste0("Cancer vs Normal ", pfmt(w3$p.value)), size=3.5)
sv(p1d, "Fig1D_GSE63514_expression", 4.6, 3.8)

# Fig2A CNA scatter
d7 <- read.csv(file.path(rd,"07_CESC_LONP1_CNA.csv"))
d8 <- read.csv(file.path(rd,"08_CESC_gene_panel_expression_RSEM.csv"))
m <- merge(d7, data.frame(sample=d8$sample, expr=log2(d8$LONP1+1)))
m <- m[complete.cases(m[,c("log2CNA","expr")]),]
ct <- suppressWarnings(cor.test(m$log2CNA, m$expr, method="spearman"))
p2a <- ggplot(m, aes(x=log2CNA, y=expr)) +
  geom_point(size=1.4, alpha=0.55, color="#3C5488", shape=16) +
  geom_smooth(method="lm", formula=y~x, color="#E64B35", fill="#F39B7F", linewidth=0.8) +
  labs(x=expression(bold("LONP1 copy number, log"[2]*"(CNA ratio)")),
       y=expression(bold("LONP1 mRNA, log"[2]*"(RSEM+1)")),
       title=sprintf("TCGA-CESC (n=%d)", nrow(m))) + th +
  annotate("text", x=min(m$log2CNA)+0.28, y=max(m$expr)-0.1,
           label=sprintf("Spearman rho = %.2f\np = %.2e", ct$estimate, ct$p.value), size=3.6, hjust=0)
sv(p2a, "Fig2A_CNA_expression_scatter", 3.8, 3.8)

# Fig2B expr by GISTIC
m$g <- factor(m$gistic, levels=c(-2,-1,0,1,2),
              labels=c("Deep del","Shallow del","Diploid","Gain","Amp"))
kw <- kruskal.test(m$expr, m$g)
p2b <- ggplot(m[!is.na(m$g),], aes(x=g, y=expr, fill=g)) +
  geom_boxplot(width=0.6, outlier.shape=NA, alpha=0.85, linewidth=0.5) +
  geom_jitter(width=0.15, size=1.0, alpha=0.4, shape=16) +
  scale_fill_manual(values=c("#4DBBD5","#91D1C2","grey85","#F39B7F","#E64B35")) +
  labs(x=NULL, y=expression(bold("LONP1 mRNA, log"[2]*"(RSEM+1)")),
       title="Expression by copy-number class") + th +
  theme(axis.text.x=element_text(angle=30, hjust=1)) +
  annotate("text", x=3, y=max(m$expr)+0.15, label=paste0("Kruskal-Wallis ", pfmt(kw$p.value)), size=3.6)
sv(p2b, "Fig2B_expression_by_GISTIC", 3.8, 3.9)

# Fig2C subtype
d10 <- read.csv(file.path(rd,"10_CESC_subtype.csv"))
d8$patient <- substr(d8$sample,1,12)
ms <- merge(data.frame(patient=d8$patient, expr=log2(d8$LONP1+1)), d10)
ms <- ms[ms$subtype %in% c("CESC_AdenoCarcinoma","CESC_SquamousCarcinoma"),]
ms$x <- factor(ifelse(ms$subtype=="CESC_AdenoCarcinoma","Adenocarcinoma","Squamous"),
               levels=c("Squamous","Adenocarcinoma"))
t2 <- t.test(expr~x, data=ms)
p2c <- ggplot(ms, aes(x=x, y=expr, fill=x)) +
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.85, linewidth=0.5) +
  geom_jitter(width=0.12, size=1.1, alpha=0.45, shape=16) +
  scale_fill_manual(values=c("#8491B4","#F39B7F")) +
  scale_x_discrete(labels=c(sprintf("SCC\n(n=%d)",sum(ms$x=="Squamous")),
                            sprintf("Adeno\n(n=%d)",sum(ms$x=="Adenocarcinoma")))) +
  labs(x=NULL, y=expression(bold("LONP1 mRNA, log"[2]*"(RSEM+1)")), title="Histological subtype") + th +
  annotate("text", x=1.5, y=max(ms$expr)+0.15, label=pfmt(t2$p.value), size=3.8)
sv(p2c, "Fig2C_subtype", 3.0, 3.8)
cat("FIG12_DONE\n")
