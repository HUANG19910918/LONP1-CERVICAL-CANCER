# 03_coloc.R —— coloc.abf（coloc 5.2.3）：暴露 eQTL（quant）× 宫颈癌 GWAS（cc）
# 先验：p1=1e-4, p2=1e-4, p12=1e-5（默认）；敏感性另报 p12=5e-6。强共定位证据：PP.H4 > 0.8；H3 = 两个独立信号。
# MAF：欧洲结局用 1000G EUR，BBJ 用 1000G EAS（缺失时用结局自身 af）；暴露 N 取文件 n 列；结局 s = 病例/(病例+对照)。
# 区段内暴露最小 P > 1e-4 时只记录、不解读（无 eQTL 信号）。
# 输出：<结果目录>/coloc_<标签>.csv
# 用法: Rscript 03_coloc.R <coloc 输入目录> <结果 csv>
suppressMessages({library(coloc); library(data.table)})
args <- commandArgs(trailingOnly = TRUE); ind <- args[1]; outf <- args[2]
cc <- list(FinnGen_CC = c(878, 222078), FinnGen_CIS = c(3701, 219058), FinnGen_SCC = c(222, 222078), FinnGen_ADC = c(142, 222078),
           UKB_C53 = c(273, 420258), `UKB_180.1` = c(299, 198725), `UKB_180.3` = c(2186, 198725), BBJ_CeC = c(967, 60614))
res <- list()
for (f in list.files(ind, pattern = "\\.tsv$", full.names = TRUE)) {
  p <- strsplit(sub("\\.tsv$", "", basename(f)), "__")[[1]]; gene <- p[1]; ex <- p[2]; on <- p[3]
  d <- fread(f); d <- d[is.finite(bx) & is.finite(sx) & is.finite(by) & is.finite(sy) & sx > 0 & sy > 0]
  if (ex == "eQTLGen_blood") next  # 显著 cis 文件只含显著 SNP，不能用于 coloc（待完整文件）
  d[, maf := if (on == "BBJ_CeC") fifelse(is.na(eas_af), af_out, eas_af) else fifelse(is.na(eur_af), af_out, eur_af)]
  d[, maf := pmin(maf, 1 - maf)]; d <- d[!is.na(maf) & maf > 0.005]; d <- d[!duplicated(d$key)]; maf <- d$maf
  if (nrow(d) < 50) next
  n1 <- cc[[on]][1]; n0 <- cc[[on]][2]
  D1 <- list(beta = d$bx, varbeta = d$sx^2, N = median(d$n), MAF = maf, type = "quant", snp = d$key, sdY = 1)
  D2 <- list(beta = d$by, varbeta = d$sy^2, N = n1 + n0, s = n1 / (n1 + n0), type = "cc", snp = d$key)
  r <- tryCatch(suppressMessages(coloc.abf(D1, D2)), error = function(e) NULL); if (is.null(r)) { message('coloc failed ', basename(f)); next }
  r2 <- suppressMessages(coloc.abf(D1, D2, p12 = 5e-6))
  top <- d[which.max(abs(d$bx / d$sx))]
  res[[length(res) + 1]] <- data.table(gene = gene, exposure = ex, outcome = on, nsnps = nrow(d),
                                        min_p_exposure = min(d$px), min_p_outcome = min(d$py),
                                        PP.H0 = r$summary[["PP.H0.abf"]], PP.H1 = r$summary[["PP.H1.abf"]], PP.H2 = r$summary[["PP.H2.abf"]],
                                        PP.H3 = r$summary[["PP.H3.abf"]], PP.H4 = r$summary[["PP.H4.abf"]],
                                        PP.H4_p12_5e6 = r2$summary[["PP.H4.abf"]], lead_eqtl = top$key,
                                        lead_eqtl_p_outcome = top$py)
}
out <- rbindlist(res); fwrite(out, outf)
print(out[, .(gene, exposure, outcome, nsnps, min_p_exposure = signif(min_p_exposure, 2), min_p_outcome = signif(min_p_outcome, 2),
              H3 = round(PP.H3, 3), H4 = round(PP.H4, 3))])
