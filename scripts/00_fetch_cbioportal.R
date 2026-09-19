# 00_fetch_cbioportal.R
# 从cBioPortal API取TCGA-CESC (PanCancer Atlas)数据，落盘raw_data/
# 复现口径：study=cesc_tcga_pan_can_atlas_2018, profile=rna_seq_v2_mrna (RSEM),
#          sampleList=cesc_tcga_pan_can_atlas_2018_all, LONP1 entrez=9361
# 取数日期见文件时间戳；所有下游图表数值可由本脚本+01~05脚本完全复现
suppressMessages({library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
base <- args[1]  # 分析根目录
rd <- file.path(base, "raw_data")
api <- "https://www.cbioportal.org/api"
post <- function(path, body){
  tf <- tempfile(); tb <- tempfile()
  writeLines(toJSON(body, auto_unbox=TRUE), tb)
  system(paste0("curl -s -X POST -H 'Content-Type: application/json' --data @", tb,
                " '", api, path, "' -o ", tf))
  fromJSON(tf)
}
get <- function(path){ tf <- tempfile(); system(paste0("curl -s '", api, path, "' -o ", tf)); fromJSON(tf) }
prof <- "cesc_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"
sl <- "cesc_tcga_pan_can_atlas_2018_all"

# 1) 线粒体动态/质控基因panel + 免疫标志物 表达
mito <- c("DNM1L","MFN1","MFN2","OPA1","FIS1","MFF","MIEF1","MIEF2","PINK1","PRKN",
          "PPARGC1A","TFAM","NRF1","CLPP","HSPD1","HSPA9","AFG3L2","YME1L1","SPG7","LONP1")
imm <- c("CD8A","CD8B","CD4","IL7R","CD19","MS4A1","CD79A","NCAM1","KLRD1","NKG7",
         "NOS2","IRF5","PTGS2","CD163","MRC1","MS4A4A","CD86","CSF1R","CEACAM8","ITGAM",
         "CCR7","ITGAX","CD1C","NRP1","FOXP3","CCR8","IL2RA","PDCD1","CTLA4","LAG3",
         "HAVCR2","TIGIT","CD274","PDCD1LG2","GZMA","GZMB","PRF1","IFNG")
gsym <- unique(c(mito, imm))
g <- post("/genes/fetch?geneIdType=HUGO_GENE_SYMBOL", gsym)
gm <- setNames(g$hugoGeneSymbol, g$entrezGeneId)
d <- post(paste0("/molecular-profiles/", prof, "/molecular-data/fetch?projection=SUMMARY"),
          list(sampleListId=sl, entrezGeneIds=as.integer(names(gm))))
expr <- data.frame(sample=d$sampleId, gene=gm[as.character(d$entrezGeneId)], value=d$value)
w <- reshape(expr, idvar="sample", timevar="gene", direction="wide")
names(w) <- sub("^value\\.", "", names(w))
write.csv(w, file.path(rd, "08_CESC_gene_panel_expression_RSEM.csv"), row.names=FALSE)

# 2) CNA (log2CNA + GISTIC)
cna <- post("/molecular-profiles/cesc_tcga_pan_can_atlas_2018_log2CNA/molecular-data/fetch?projection=SUMMARY",
            list(sampleListId=sl, entrezGeneIds=list(9361L)))
gis <- post("/molecular-profiles/cesc_tcga_pan_can_atlas_2018_gistic/molecular-data/fetch?projection=SUMMARY",
            list(sampleListId=sl, entrezGeneIds=list(9361L)))
cn <- merge(data.frame(sample=cna$sampleId, log2CNA=cna$value),
            data.frame(sample=gis$sampleId, gistic=gis$value), all=TRUE)
write.csv(cn, file.path(rd, "07_CESC_LONP1_CNA.csv"), row.names=FALSE)

# 3) SUBTYPE（病人级临床）
cl <- get("/studies/cesc_tcga_pan_can_atlas_2018/clinical-data?clinicalDataType=PATIENT&projection=SUMMARY&pageSize=100000")
st <- cl[cl$clinicalAttributeId=="SUBTYPE", c("patientId","value")]
names(st) <- c("patient","subtype")
write.csv(st, file.path(rd, "10_CESC_subtype.csv"), row.names=FALSE)

# 4) 全基因组共表达（服务器端Spearman, |rho|>=0.3）
tf <- tempfile(); tb <- tempfile()
writeLines('{"entrezGeneId":9361,"sampleListId":"cesc_tcga_pan_can_atlas_2018_all"}', tb)
system(paste0("curl -s -X POST -H 'Content-Type: application/json' --data @", tb,
  " '", api, "/molecular-profiles/co-expressions/fetch?molecularProfileIdA=", prof,
  "&molecularProfileIdB=", prof, "&threshold=0.3' -o ", tf))
co <- fromJSON(tf)
co <- co[!is.na(co$spearmansCorrelation),]
g2 <- post("/genes/fetch?geneIdType=ENTREZ_GENE_ID", as.character(co$geneticEntityId))
gm2 <- setNames(g2$hugoGeneSymbol, g2$entrezGeneId)
out <- data.frame(gene=gm2[as.character(co$geneticEntityId)],
                  spearman_rho=co$spearmansCorrelation, p=co$pValue)
out <- out[order(-out$spearman_rho),]
write.csv(out, file.path(rd, "06_coexpression_LONP1_CESC.csv"), row.names=FALSE)
cat("FETCH_DONE  panel:", nrow(w), "cna:", nrow(cn), "subtype:", nrow(st), "coexp:", nrow(out), "\n")
