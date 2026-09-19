# 04_fig_survival_immune.R —— Figure 4：生存（阴性结果，如实呈现）；Figure 5：免疫
suppressMessages({library(ggplot2); library(survival)})
args <- commandArgs(trailingOnly=TRUE); base <- args[1]
rd <- file.path(base,"raw_data"); fig <- file.path(base,"figures","子图原文件"); dir.create(fig, showWarnings=FALSE, recursive=TRUE)
th <- theme_classic(base_size=13) + theme(
  axis.text=element_text(color="black", size=11),
  axis.title=element_text(size=12, face="bold"),
  plot.title=element_text(size=12, face="bold", hjust=0.5))
sv <- function(p,name,w,h){ ggsave(file.path(fig,paste0(name,".pdf")),p,width=w,height=h)
  ggsave(file.path(fig,paste0(name,".png")),p,width=w,height=h,dpi=300) }
d2 <- read.csv(file.path(rd,"02_TCGA_CESC_clinical_LONP1_survival.csv"))

# Fig4A-C: KM 曲线（中位数分组），手绘KM（ggplot2+survfit）
km_panel <- function(tk, ek, label){
  x <- d2[!is.na(d2[[tk]]) & !is.na(d2[[ek]]) & d2[[tk]]>0,]
  x$grp <- factor(ifelse(x$LONP1_log2RSEM > median(x$LONP1_log2RSEM), "High", "Low"), levels=c("Low","High"))
  sf <- survfit(Surv(x[[tk]], x[[ek]]) ~ grp, data=x)
  sd <- survdiff(Surv(x[[tk]], x[[ek]]) ~ grp, data=x)
  pv <- 1-pchisq(sd$chisq, 1)
  cx <- summary(coxph(Surv(x[[tk]], x[[ek]]) ~ LONP1_log2RSEM, data=x))
  stnames <- sub("^grp=", "", names(sf$strata))
  dd <- data.frame(t=sf$time, s=sf$surv, c=sf$n.censor,
                   g=rep(stnames, sf$strata))
  dd0 <- rbind(data.frame(t=0, s=1, c=0, g="Low"), data.frame(t=0, s=1, c=0, g="High"), dd)
  dd0$g <- factor(dd0$g, levels=c("Low","High"))
  p <- ggplot(dd0, aes(x=t, y=s, color=g)) +
    geom_step(linewidth=0.9) +
    geom_point(data=dd0[dd0$c>0,], shape=3, size=1.6, show.legend=FALSE) +
    scale_color_manual(values=c(Low="#4DBBD5", High="#E64B35"), breaks=c("Low","High"),
                       labels=c(sprintf("Low (n=%d)", sum(x$grp=="Low")),
                                sprintf("High (n=%d)", sum(x$grp=="High")))) +
    scale_y_continuous(limits=c(0,1), labels=scales::percent) +
    labs(x="Time (months)", y=sprintf("%s probability", label),
         title=sprintf("%s by LONP1 (median split)", label)) + th +
    theme(legend.position=c(0.78,0.9), legend.title=element_blank(),
          legend.text=element_text(size=9), legend.background=element_blank()) +
    annotate("text", x=max(dd0$t)*0.55, y=0.14,
             label=sprintf("Log-rank p = %.3f\nCox (continuous) HR = %.2f\n(95%% CI %.2f-%.2f), p = %.3f",
                           pv, cx$conf.int[1], cx$conf.int[3], cx$conf.int[4], cx$coefficients[5]),
             size=3.1, hjust=0)
  p
}
sv(km_panel("OS_months","OS_event","Overall survival"), "Fig4A_KM_OS", 4.3, 4.0)
sv(km_panel("DSS_months","DSS_event","Disease-specific survival"), "Fig4B_KM_DSS", 4.3, 4.0)
sv(km_panel("PFS_months","PFS_event","Progression-free survival"), "Fig4C_KM_PFS", 4.3, 4.0)

# Fig4D: 森林图（三终点连续变量Cox）
fr <- do.call(rbind, lapply(list(c("OS_months","OS_event","OS"),
                                 c("DSS_months","DSS_event","DSS"),
                                 c("PFS_months","PFS_event","PFS")), function(ep){
  x <- d2[!is.na(d2[[ep[1]]]) & !is.na(d2[[ep[2]]]) & d2[[ep[1]]]>0,]
  s <- summary(coxph(Surv(x[[ep[1]]], x[[ep[2]]]) ~ LONP1_log2RSEM, data=x))
  data.frame(ep=ep[3], n=s$n, ev=s$nevent, hr=s$conf.int[1], lo=s$conf.int[3], hi=s$conf.int[4], p=s$coefficients[5])
}))
write.csv(fr, file.path(rd,"13_cox_forest_data.csv"), row.names=FALSE)
fr$lab <- sprintf("%s (n=%d, events=%d)", fr$ep, fr$n, fr$ev)
fr$lab <- factor(fr$lab, levels=rev(fr$lab))
p4d <- ggplot(fr, aes(x=hr, y=lab)) +
  geom_vline(xintercept=1, linetype=2, color="grey50") +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.15, linewidth=0.7, color="#3C5488") +
  geom_point(size=3, shape=15, color="#E64B35") +
  geom_text(aes(label=sprintf("HR %.2f (%.2f-%.2f), p=%.2f", hr, lo, hi, p)),
            vjust=-1.1, size=3.2) +
  scale_x_log10(limits=c(0.4, 2.2)) +
  labs(x="Hazard ratio (per log2 unit LONP1)", y=NULL,
       title="Univariate Cox regression (TCGA-CESC)") + th
sv(p4d, "Fig4D_forest", 5.4, 3.4)

# Fig5A: 免疫签名相关性条形图；Fig5B: 代表性散点(M1)
d8 <- read.csv(file.path(rd,"08_CESC_gene_panel_expression_RSEM.csv"))
lg <- log2(d8$LONP1+1)
sets <- list("CD8 T cell"=c("CD8A","CD8B"), "CD4/Th"=c("CD4","IL7R"), "B cell"=c("CD19","MS4A1","CD79A"),
  "NK cell"=c("NCAM1","KLRD1","NKG7"), "Macrophage M1"=c("NOS2","IRF5","PTGS2"),
  "Macrophage M2"=c("CD163","MRC1","MS4A4A"), "Monocyte"=c("CD86","CSF1R"),
  "Neutrophil"=c("CEACAM8","ITGAM","CCR7"), "Dendritic cell"=c("ITGAX","CD1C","NRP1"),
  "Treg"=c("FOXP3","CCR8","IL2RA"), "Exhaustion"=c("PDCD1","CTLA4","LAG3","HAVCR2","TIGIT"),
  "Checkpoint ligand"=c("CD274","PDCD1LG2"), "Cytotoxicity"=c("GZMA","GZMB","PRF1","IFNG"))
zs <- function(v){ (v-mean(v))/sd(v) }
ir <- do.call(rbind, lapply(names(sets), function(nm){
  m <- sapply(sets[[nm]], function(g) zs(log2(d8[[g]]+1)))
  sig <- rowMeans(m)
  ct <- suppressWarnings(cor.test(lg, sig, method="spearman"))
  data.frame(sig=nm, rho=unname(ct$estimate), p=ct$p.value)}))
ir$padj <- p.adjust(ir$p, "BH")
write.csv(ir, file.path(rd,"14_immune_signature_spearman.csv"), row.names=FALSE)
ir <- ir[order(ir$rho),]; ir$sig <- factor(ir$sig, levels=ir$sig)
ir$col <- ifelse(ir$padj<0.05, ifelse(ir$rho<0, "negS","posS"), "ns")
p5a <- ggplot(ir, aes(x=rho, y=sig, fill=col)) +
  geom_col(width=0.68, alpha=0.9) +
  geom_vline(xintercept=0, color="grey40", linewidth=0.4) +
  geom_text(aes(label=ifelse(padj<0.05, sprintf("FDR=%.3f", padj), "")),
            x=ifelse(ir$rho<0, 0.012, -0.012), hjust=ifelse(ir$rho<0, 0, 1), size=2.9) +
  scale_fill_manual(values=c(negS="#4DBBD5", posS="#E64B35", ns="grey80")) +
  labs(x="Spearman correlation with LONP1", y=NULL,
       title="Immune marker signatures (TCGA-CESC)") + th + theme(legend.position="none")
sv(p5a, "Fig5A_immune_bar", 5.4, 4.4)
m1 <- rowMeans(sapply(sets[["Macrophage M1"]], function(g) zs(log2(d8[[g]]+1))))
ct <- suppressWarnings(cor.test(lg, m1, method="spearman"))
p5b <- ggplot(data.frame(x=lg, y=m1), aes(x,y)) +
  geom_point(size=1.4, alpha=0.55, color="#3C5488", shape=16) +
  geom_smooth(method="lm", formula=y~x, color="#E64B35", fill="#F39B7F", linewidth=0.8) +
  labs(x=expression(bold("LONP1 mRNA, log"[2]*"(RSEM+1)")), y="M1 macrophage signature (z-score)",
       title=sprintf("rho = %.2f, p = %.1e", ct$estimate, ct$p.value)) + th
sv(p5b, "Fig5B_M1_scatter", 3.8, 3.8)
cat("FIG45_DONE\n")
