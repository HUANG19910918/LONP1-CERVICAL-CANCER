# 07_MR_package_crosscheck.R —— 用 MendelianRandomization 包在同一套协调后工具变量上重算 IVW / MR-Egger / 加权中位数 / Cochran Q，
# 并与 02_harmonize_mr.py 的 Python 实现逐项比较。
# 模型对应：Python "IVW (fixed)" → mr_ivw(model="fixed")；"IVW (multiplicative RE)" → mr_ivw(model="random")；
#           单 SNP → Wald ratio（by/bx, SE=sy/|bx|，两边公式一致，直接重算）；
#           MR-Egger → mr_egger(model="random")（残差 SE 下限 1，与 Python 的 max(1, σ²) 一致）；
#           加权中位数 → mr_median(weighting="weighted", iterations=10000, seed=1)（bootstrap SE 两边随机误差不同，单独报告）。
# 输出：results/mr_package_crosscheck.csv、results/mr_package_crosscheck_summary.txt
# 用法: Rscript 07_MR_package_crosscheck.R <03_MR_coloc 目录>
suppressMessages({library(MendelianRandomization); library(data.table)})
base <- commandArgs(trailingOnly = TRUE)[1]; res <- file.path(base, "results")
M <- fread(file.path(res, "mr_all_genes_combined.csv"))
fs <- list.files(res, pattern = "^instruments_.*\\.csv$", full.names = TRUE); fs <- fs[file.size(fs) > 10]
I <- rbindlist(lapply(fs, fread), fill = TRUE)
out <- list()
for (k in seq_len(nrow(M))) {
  m <- M[k]; h <- I[gene == m$gene & exposure == m$exposure & outcome == m$outcome]
  if (nrow(h) != m$n_snp) { message("SNP count mismatch ", m$gene, " ", m$exposure, " ", m$outcome); next }
  o <- data.table(gene = m$gene, exposure = m$exposure, outcome = m$outcome, n_snp = m$n_snp, method = m$method,
                  py_b = log(m$OR), py_se = (log(m$OR_hi) - log(m$OR_lo)) / (2 * 1.959964), py_p = m$p)
  if (nrow(h) == 1) {
    o[, `:=`(r_b = h$by / h$bx, r_se = h$sy / abs(h$bx))]; o[, r_p := 2 * pnorm(-abs(r_b / r_se))]
  } else {
    inp <- mr_input(bx = h$bx, bxse = h$sx, by = h$by, byse = h$sy, snps = h$rsid)
    fit <- mr_ivw(inp, model = if (grepl("multiplicative", m$method)) "random" else "fixed")
    o[, `:=`(r_b = fit@Estimate, r_se = fit@StdError, r_p = fit@Pvalue, py_Q = m$Q, r_Q = fit@Heter.Stat[1], py_Q_p = m$Q_p, r_Q_p = fit@Heter.Stat[2])]
    if (nrow(h) >= 3) {
      eg <- mr_egger(inp); wm <- mr_median(inp, weighting = "weighted", iterations = 10000, seed = 1)
      o[, `:=`(py_egger_b = m$egger_b, r_egger_b = eg@Estimate, py_egger_se = m$egger_se, r_egger_se = eg@StdError.Est,
               py_egger_int = m$egger_int, r_egger_int = eg@Intercept, py_egger_int_p = m$egger_int_p, r_egger_int_p = eg@Pvalue.Int,
               py_wm_b = m$wm_b, r_wm_b = wm@Estimate, py_wm_se = m$wm_se, r_wm_se = wm@StdError, r_wm_p = wm@Pvalue)]
    }
  }
  out[[length(out) + 1]] <- o
}
R <- rbindlist(out, fill = TRUE); fwrite(R, file.path(res, "mr_package_crosscheck.csv"))
md <- function(a, b) { d <- abs(a - b); d <- d[is.finite(d)]; if (length(d)) sprintf("%.2e (n=%d)", max(d), length(d)) else "NA" }
rel <- function(a, b) { d <- abs(a - b) / pmax(abs(b), 1e-12); d <- d[is.finite(d)]; if (length(d)) sprintf("%.2e", max(d)) else "NA" }
L <- c(sprintf("MendelianRandomization %s; estimates compared: %d (Wald %d, IVW %d, Egger/WM %d)", packageVersion("MendelianRandomization"),
               nrow(R), sum(R$n_snp == 1), sum(R$n_snp >= 2), sum(R$n_snp >= 3)),
       paste("IVW/Wald beta  max |diff|:", md(R$py_b, R$r_b), " max rel:", rel(R$py_b, R$r_b)),
       paste("IVW/Wald SE    max |diff|:", md(R$py_se, R$r_se), " max rel:", rel(R$py_se, R$r_se)),
       paste("IVW/Wald P     max |diff|:", md(R$py_p, R$r_p)),
       paste("Cochran Q      max |diff|:", md(R$py_Q, R$r_Q), "; Q P max |diff|:", md(R$py_Q_p, R$r_Q_p)),
       paste("Egger slope    max |diff|:", md(R$py_egger_b, R$r_egger_b), "; SE:", md(R$py_egger_se, R$r_egger_se)),
       paste("Egger intercept max |diff|:", md(R$py_egger_int, R$r_egger_int), "; intercept P:", md(R$py_egger_int_p, R$r_egger_int_p)),
       paste("Weighted median beta max |diff|:", md(R$py_wm_b, R$r_wm_b), "; SE (bootstrap):", md(R$py_wm_se, R$r_wm_se)),
       paste("FDR<0.05 set / significance at P<0.05 identical:", all((R$py_p < 0.05) == (R$r_p < 0.05))))
writeLines(L, file.path(res, "mr_package_crosscheck_summary.txt")); cat(L, sep = "\n")
