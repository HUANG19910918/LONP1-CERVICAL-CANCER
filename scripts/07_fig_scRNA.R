# 07_fig_scRNA.R —— Figure 6：GSE208653单细胞LONP1定位（5 panel）
# 输入: raw_data/17_GSE208653_LONP1_per_cell.csv.gz（80,435细胞，由06脚本生成）
suppressMessages({library(ggplot2); library(scales)})
args <- commandArgs(trailingOnly=TRUE); base <- args[1]
rd <- file.path(base,"raw_data"); fig <- file.path(base,"figures","子图原文件"); dir.create(fig, showWarnings=FALSE, recursive=TRUE)
th <- theme_classic(base_size=13) + theme(
  axis.text=element_text(color="black", size=11),
  axis.title=element_text(size=12, face="bold"),
  plot.title=element_text(size=12, face="bold", hjust=0.5))
sv <- function(p,name,w,h){ ggsave(file.path(fig,paste0(name,".pdf")),p,width=w,height=h)
  ggsave(file.path(fig,paste0(name,".png")),p,width=w,height=h,dpi=300) }
d <- read.csv(file.path(rd,"17_GSE208653_LONP1_per_cell.csv.gz"))
d$stage <- factor(d$stage, levels=c("Normal_HPVneg","Normal_HPVpos","HSIL","Cancer"))
stageLab <- c(Normal_HPVneg="Normal HPV-", Normal_HPVpos="Normal HPV+", HSIL="HSIL", Cancer="Cancer")
ctCols <- c(Epithelial="#E64B35", T_NK="#4DBBD5", Myeloid="#00A087", Plasma="#3C5488",
  B="#F39B7F", Fibroblast="#8491B4", Endothelial="#91D1C2", Mast="#DC0000", SmoothMuscle="#7E6148")

# 6A UMAP by cell type
p6a <- ggplot(d, aes(umap1, umap2, color=cell_type)) +
  geom_point(size=0.12, alpha=0.5, shape=16) +
  scale_color_manual(values=ctCols) +
  guides(color=guide_legend(override.aes=list(size=3, alpha=1))) +
  labs(x="UMAP1", y="UMAP2", title="GSE208653: 80,435 cells, 9 samples", color=NULL) +
  th + theme(legend.text=element_text(size=9))
sv(p6a, "Fig6A_UMAP_celltype", 6.0, 4.4)

# 6B UMAP by LONP1 (order low->high so positives on top)
d2 <- d[order(d$LONP1),]
p6b <- ggplot(d2, aes(umap1, umap2, color=LONP1)) +
  geom_point(size=0.12, alpha=0.6, shape=16) +
  scale_color_gradient(low="grey88", high="#B2182B", name="LONP1") +
  labs(x="UMAP1", y="UMAP2", title="LONP1 expression") + th
sv(p6b, "Fig6B_UMAP_LONP1", 5.2, 4.4)

# 6C violin by cell type (ordered by mean)
ord <- names(sort(tapply(d$LONP1, d$cell_type, mean), decreasing=TRUE))
d$ctF <- factor(d$cell_type, levels=ord)
p6c <- ggplot(d, aes(ctF, LONP1, fill=ctF)) +
  geom_violin(scale="width", linewidth=0.3, adjust=2) +
  stat_summary(fun=mean, geom="point", size=1.6, color="black") +
  scale_fill_manual(values=ctCols) +
  labs(x=NULL, y="LONP1, log1p(CP10K)", title="LONP1 by cell type") +
  th + theme(axis.text.x=element_text(angle=40, hjust=1), legend.position="none")
sv(p6c, "Fig6C_violin_celltype", 5.4, 3.9)

# 6D epithelial by stage: violin + per-sample pseudobulk points
epi <- d[d$cell_type=="Epithelial",]
pb <- aggregate(LONP1~sample+stage, epi, mean)
kw <- kruskal.test(epi$LONP1, epi$stage)
p6d <- ggplot(epi, aes(stage, LONP1, fill=stage)) +
  geom_violin(scale="width", linewidth=0.3, adjust=2) +
  geom_point(data=pb, aes(stage, LONP1), size=2.6, shape=21, fill="white", stroke=0.9,
             position=position_jitter(width=0.08, seed=1)) +
  scale_fill_manual(values=c("#4DBBD5","#91D1C2","#F39B7F","#E64B35")) +
  scale_x_discrete(labels=paste0(stageLab[levels(epi$stage)], "\n(n=", table(epi$stage)[levels(epi$stage)], ")")) +
  labs(x=NULL, y="LONP1, log1p(CP10K)",
       title="Epithelial cells across disease stages") +
  th + theme(legend.position="none") +
  annotate("text", x=2.5, y=max(epi$LONP1)*0.97,
           label=sprintf("Kruskal-Wallis p = %.1e (per cell)\ncircles = per-sample means", kw$p.value), size=3.2)
sv(p6d, "Fig6D_epithelial_stage", 5.0, 4.0)

# 6E cell-type composition per stage (关键：组织构成论据)
comp <- as.data.frame(prop.table(table(d$stage, d$cell_type), 1))
names(comp) <- c("stage","cell_type","frac")
p6e <- ggplot(comp, aes(stage, frac, fill=cell_type)) +
  geom_col(width=0.72, color="white", linewidth=0.2) +
  scale_fill_manual(values=ctCols, name=NULL) +
  scale_y_continuous(labels=percent) +
  scale_x_discrete(labels=stageLab[levels(comp$stage)]) +
  labs(x=NULL, y="Fraction of cells", title="Cellular composition by stage") +
  th + theme(axis.text.x=element_text(angle=25, hjust=1), legend.text=element_text(size=8.5))
sv(p6e, "Fig6E_composition", 5.6, 4.0)
write.csv(pb, file.path(rd,"19_GSE208653_epithelial_pseudobulk_LONP1.csv"), row.names=FALSE)
cat("FIG6_DONE\n")
