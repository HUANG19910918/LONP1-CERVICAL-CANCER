# 03_fig_mito_enrich.R —— Figure 3：线粒体动态基因相关性 + 共表达富集
# 富集口径：CSV06中rho>=0.4的基因(n=165)提交Enrichr(GO_BP_2023/GO_CC_2023/KEGG_2021_Human)
suppressMessages({library(ggplot2); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE); base <- args[1]
rd <- file.path(base,"raw_data"); fig <- file.path(base,"figures","子图原文件"); dir.create(fig, showWarnings=FALSE, recursive=TRUE)
th <- theme_classic(base_size=13) + theme(
  axis.text=element_text(color="black", size=11),
  axis.title=element_text(size=12, face="bold"),
  plot.title=element_text(size=12, face="bold", hjust=0.5),
  legend.position="none")
sv <- function(p,name,w,h){ ggsave(file.path(fig,paste0(name,".pdf")),p,width=w,height=h)
  ggsave(file.path(fig,paste0(name,".png")),p,width=w,height=h,dpi=300) }

# Fig3A: 线粒体动态/质控基因 Spearman 棒棒糖图
d8 <- read.csv(file.path(rd,"08_CESC_gene_panel_expression_RSEM.csv"))
lg <- log2(d8$LONP1+1)
genes <- c("CLPP","HSPD1","HSPA9","AFG3L2","SPG7","YME1L1",       # 质控蛋白酶/伴侣
           "FIS1","MFF","DNM1L","MIEF1","MIEF2",                   # 分裂
           "MFN1","MFN2","OPA1",                                   # 融合
           "PINK1","PRKN","PPARGC1A","NRF1","TFAM")                # 自噬/生物发生
cls <- c(rep("Protease/chaperone",6), rep("Fission",5), rep("Fusion",3), rep("Mitophagy/biogenesis",5))
res <- do.call(rbind, lapply(seq_along(genes), function(i){
  ct <- suppressWarnings(cor.test(lg, log2(d8[[genes[i]]]+1), method="spearman"))
  data.frame(gene=genes[i], class=cls[i], rho=unname(ct$estimate), p=ct$p.value)}))
res$sig <- cut(res$p, c(-Inf,0.001,0.01,0.05,Inf), labels=c("***","**","*",""))
res <- res[order(res$rho),]; res$gene <- factor(res$gene, levels=res$gene)
write.csv(res, file.path(rd,"12_mito_panel_spearman.csv"), row.names=FALSE)
p3a <- ggplot(res, aes(x=rho, y=gene, color=class)) +
  geom_segment(aes(x=0, xend=rho, yend=gene), linewidth=0.7) +
  geom_point(size=2.6) +
  geom_vline(xintercept=0, linetype=2, color="grey50", linewidth=0.4) +
  geom_text(aes(label=sig, x=rho+ifelse(rho>0,0.045,-0.045)), size=3.6, color="black") +
  scale_color_manual(values=c("Protease/chaperone"="#E64B35","Fission"="#F39B7F",
                              "Fusion"="#4DBBD5","Mitophagy/biogenesis"="#8491B4")) +
  labs(x="Spearman correlation with LONP1", y=NULL,
       title="Mitochondrial dynamics & quality control") +
  th + theme(legend.position="right", legend.title=element_blank(),
             legend.text=element_text(size=9))
sv(p3a, "Fig3A_mito_panel_lollipop", 6.2, 4.6)

# Fig3B/C: Enrichr富集（重新提交，保证可复现）
co <- read.csv(file.path(rd,"06_coexpression_LONP1_CESC.csv"))
pos <- co$gene[co$spearman_rho>=0.4 & !is.na(co$gene) & co$gene!=""]
gl <- tempfile(); writeLines(paste(pos, collapse="\n"), gl)
tf <- tempfile()
system(paste0("curl -s -F 'list=<", gl, "' -F 'description=LONP1posCESC' https://maayanlab.cloud/Enrichr/addList -o ", tf))
uid <- fromJSON(tf)$userListId
enr <- function(libname){
  tf2 <- tempfile()
  system(paste0("curl -s 'https://maayanlab.cloud/Enrichr/enrich?userListId=", uid,
                "&backgroundType=", libname, "' -o ", tf2))
  j <- fromJSON(tf2)[[libname]]
  data.frame(lib=libname, term=sapply(j, `[[`, 2), p=as.numeric(sapply(j, `[[`, 3)),
             padj=as.numeric(sapply(j, `[[`, 7)),
             n_overlap=sapply(j, function(x) length(x[[6]])))[1:12,]
}
Sys.sleep(1)
e1 <- enr("GO_Biological_Process_2023"); e2 <- enr("GO_Cellular_Component_2023"); e3 <- enr("KEGG_2021_Human")
ea <- rbind(e1,e2,e3)
write.csv(ea, file.path(rd,"11_enrichment_enrichr_pos.csv"), row.names=FALSE)
mk_bar <- function(e, title, color, n=8){
  e <- e[1:n,]; e$term <- sub(" \\(GO:\\d+\\)","",e$term)
  e$term <- factor(e$term, levels=rev(e$term))
  ggplot(e, aes(x=-log10(padj), y=term)) +
    geom_col(fill=color, alpha=0.85, width=0.7) +
    geom_text(aes(label=n_overlap), hjust=-0.3, size=3.2) +
    labs(x=expression(bold("-log"[10]*" adjusted p")), y=NULL, title=title) +
    th + theme(axis.text.y=element_text(size=9.5)) +
    scale_x_continuous(expand=expansion(mult=c(0,0.12)))
}
sv(mk_bar(rbind(e2[1:4,],e1[c(1,6,8,10),]), "GO enrichment (LONP1-positive genes)", "#E64B35"),
   "Fig3B_GO_enrichment", 5.6, 3.6)
sv(mk_bar(e3, "KEGG enrichment (LONP1-positive genes)", "#3C5488"),
   "Fig3C_KEGG_enrichment", 5.6, 3.6)
cat("FIG3_DONE uid=", uid, "\n")
