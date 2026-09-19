# 05_fig_pancancer.R —— 泛癌LONP1表达概览（TCGA肿瘤 vs GTEx+TCGA正常，Xena Toil）
# 口径：TcgaTargetGtex_rsem_gene_tpm (log2(TPM+0.001))；每癌种：肿瘤=TCGA Primary Tumor，
#      正常=同_primary_site的GTEx Normal Tissue + TCGA Solid Tissue Normal；正常n>=5才检验
suppressMessages({library(ggplot2); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE); base <- args[1]
rd <- file.path(base,"raw_data"); fig <- file.path(base,"figures","子图原文件"); dir.create(fig, showWarnings=FALSE, recursive=TRUE)
hub <- "https://toil.xenahubs.net/data/"
xq <- function(q){
  tb <- tempfile(); tf <- tempfile(); writeLines(q, tb)
  system(paste0("curl -s -X POST -H 'Content-Type: text/plain' --data-binary @", tb, " '", hub, "' -o ", tf))
  fromJSON(tf)
}
ds <- "TcgaTargetGTEX_phenotype.txt"
samples <- xq(sprintf('((fn [dataset limit] (map :value (query {:select [:value] :from [:dataset] :join [:field [:= :dataset.id :dataset_id] :code [:= :field.id :field_id]] :limit limit :where [:and [:= :dataset.name dataset] [:= :field.name "sampleID"]]}))) "%s" 30000)', ds))
sArr <- paste0('[', paste0('"', samples, '"', collapse=' '), ']')
fields <- c("_primary_site","_sample_type","_study","primary disease or tissue")
fArr <- paste0('[', paste0('"', fields, '"', collapse=' '), ']')
vals <- xq(sprintf('((fn [dataset samples probes] (fetch [{:table dataset :columns probes :samples samples}])) "%s" %s %s)', ds, sArr, fArr))
codes <- xq(sprintf('((fn [dataset fields] (query {:select [:P.name [#sql/call [:group_concat :value :order :ordering :separator #sql/call [:chr 9]] :code]] :from [[{:select [:field.id :field.name] :from [:field] :join [{:table [[[:name :varchar fields]] :T]} [:= :T.name :field.name]] :where [:= :dataset_id {:select [:id] :from [:dataset] :where [:= :name dataset]}]} :P]] :left-join [:code [:= :P.id :field_id]] :group-by [:P.id]})) "%s" %s)', ds, fArr))
dec <- function(idx, fname){
  cl <- strsplit(codes$code[codes$name==fname], "\t")[[1]]
  out <- rep(NA_character_, length(idx)); ok <- !is.na(idx); out[ok] <- cl[idx[ok]+1]; out
}
ph <- data.frame(sample=samples,
  site=dec(vals[1,], "_primary_site"), stype=dec(vals[2,], "_sample_type"),
  study=dec(vals[3,], "_study"), disease=dec(vals[4,], "primary disease or tissue"))
# LONP1表达（分批取）
expr <- rep(NA_real_, length(samples))
bs <- 5000
for(i in seq(1, length(samples), bs)){
  idx <- i:min(i+bs-1, length(samples))
  sA <- paste0('[', paste0('"', samples[idx], '"', collapse=' '), ']')
  r <- xq(sprintf('((fn [dataset samples genes] (let [probemap (:probemap (car (query {:select [:probemap] :from [:dataset] :where [:= :name dataset]}))) get-probes (fn [gene] (xena-query {:select ["name" "position"] :from [probemap] :where [:in :any "genes" [gene]]})) avg (fn [scores] (mean scores 0)) sfg (fn [gene] (let [probes (get-probes gene) pn (probes "name") scores (fetch [{:table dataset :samples samples :columns pn}])] {:gene gene :scores (if (car pn) (avg scores) [[]])}))] (map sfg genes))) "TcgaTargetGtex_rsem_gene_tpm" %s ["LONP1"])', sA))
  expr[idx] <- unlist(r$scores[[1]])
}
ph$LONP1 <- expr
# TCGA癌种缩写映射
ab <- c("Adrenocortical Cancer"="ACC","Bladder Urothelial Carcinoma"="BLCA","Breast Invasive Carcinoma"="BRCA",
"Cervical & Endocervical Cancer"="CESC","Cholangiocarcinoma"="CHOL",
"Colon Adenocarcinoma"="COAD","Diffuse Large B-Cell Lymphoma"="DLBC","Esophageal Carcinoma"="ESCA",
"Glioblastoma Multiforme"="GBM","Head & Neck Squamous Cell Carcinoma"="HNSC","Kidney Chromophobe"="KICH",
"Kidney Clear Cell Carcinoma"="KIRC","Kidney Papillary Cell Carcinoma"="KIRP","Acute Myeloid Leukemia"="LAML",
"Brain Lower Grade Glioma"="LGG","Liver Hepatocellular Carcinoma"="LIHC","Lung Adenocarcinoma"="LUAD",
"Lung Squamous Cell Carcinoma"="LUSC","Mesothelioma"="MESO","Ovarian Serous Cystadenocarcinoma"="OV",
"Pancreatic Adenocarcinoma"="PAAD","Pheochromocytoma & Paraganglioma"="PCPG","Prostate Adenocarcinoma"="PRAD",
"Rectum Adenocarcinoma"="READ","Sarcoma"="SARC","Skin Cutaneous Melanoma"="SKCM","Stomach Adenocarcinoma"="STAD",
"Testicular Germ Cell Tumor"="TGCT","Thyroid Carcinoma"="THCA","Thymoma"="THYM",
"Uterine Corpus Endometrioid Carcinoma"="UCEC","Uterine Carcinosarcoma"="UCS","Uveal Melanoma"="UVM")
norm_site <- function(s){ s <- tolower(s); s[s=="cervix uteri"] <- "cervix"; s }
ph$siteN <- norm_site(ph$site)
tcgaT <- ph[!is.na(ph$study) & ph$study=="TCGA" & ph$stype %in% c("Primary Tumor","Primary Solid Tumor") & ph$disease %in% names(ab) & !is.na(ph$LONP1),]
tcgaT$type <- ab[tcgaT$disease]
normals <- ph[!is.na(ph$study) & ((ph$study=="GTEX" & ph$stype=="Normal Tissue") | (ph$study=="TCGA" & ph$stype=="Solid Tissue Normal")) & !is.na(ph$LONP1),]
rows <- list(); stats <- list()
for(tp in sort(unique(tcgaT$type))){
  tt <- tcgaT[tcgaT$type==tp,]
  sites <- unique(tt$siteN)
  nn <- normals[normals$siteN %in% sites,]
  rows[[tp]] <- rbind(data.frame(type=tp, grp="Tumor", sample=tt$sample, LONP1=tt$LONP1),
                      if(nrow(nn)) data.frame(type=tp, grp="Normal", sample=nn$sample, LONP1=nn$LONP1))
  pv <- if(nrow(nn)>=5) suppressWarnings(wilcox.test(tt$LONP1, nn$LONP1)$p.value) else NA
  stats[[tp]] <- data.frame(type=tp, nT=nrow(tt), nN=nrow(nn),
    medT=median(tt$LONP1), medN=ifelse(nrow(nn), median(nn$LONP1), NA), p=pv)
}
dat <- do.call(rbind, rows); st <- do.call(rbind, stats)
st$padj <- p.adjust(st$p, "BH")
st$sig <- ifelse(is.na(st$padj), "", ifelse(st$padj<0.001,"***", ifelse(st$padj<0.01,"**", ifelse(st$padj<0.05,"*",""))))
write.csv(dat, file.path(rd,"15_pancancer_LONP1_per_sample.csv"), row.names=FALSE)
write.csv(st, file.path(rd,"16_pancancer_LONP1_stats.csv"), row.names=FALSE)
ymax <- max(dat$LONP1, na.rm=TRUE)
p <- ggplot(dat, aes(x=type, y=LONP1, fill=grp)) +
  geom_boxplot(width=0.68, outlier.size=0.25, outlier.alpha=0.25, linewidth=0.3,
               position=position_dodge2(preserve="single")) +
  scale_fill_manual(values=c(Tumor="#E64B35", Normal="#4DBBD5"), name=NULL) +
  geom_text(data=st, aes(x=type, y=ymax+0.25, label=sig), inherit.aes=FALSE, size=3.2) +
  labs(x=NULL, y=expression(bold("LONP1 expression, log"[2]*"(TPM+0.001)")),
       title="LONP1 across TCGA cancer types (tumors vs GTEx+TCGA normal)") +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=55, hjust=1, size=8.5, color="black",
        face=ifelse(sort(unique(dat$type))=="CESC","bold","plain")),
        axis.text.y=element_text(color="black"),
        axis.title=element_text(face="bold"), legend.position="top",
        plot.title=element_text(face="bold", hjust=0.5, size=12)) +
  coord_cartesian(ylim=c(2, ymax+0.5))  # 少数零表达样本(log2≈-10)超出显示范围，原始值见raw_data/15
ggsave(file.path(fig,"Fig1E_pancancer_LONP1.pdf"), p, width=11, height=4.6)
ggsave(file.path(fig,"Fig1E_pancancer_LONP1.png"), p, width=11, height=4.6, dpi=300)
cat("PANCAN_DONE types:", nrow(st), "\n")
print(st[st$type %in% c("CESC","OV","UCEC","BRCA"),])
