# 06b_copykat.R —— CopyKAT 1.2.5 逐样本运行，交叉验证 inferCNV 恶性判定
# 参数：id.type="S"，ngene.chr=5，win.size=25，KS.cut=0.1，distance="euclidean"，genome="hg20"，
#       norm.cell.names = 同样本参考细胞（免疫+间质），set.seed(1)；每样本最多 1500 个上皮 + 500 个参考细胞（seed=1 抽样）。
# 注意：CopyKAT 1.2.5 在 "run with known normal" 路径下把全部细胞标为 diploid（试跑 GSE208653 SCC_5：3315 diploid / 226 not.defined，
#       但其 CNA 矩阵热图显示明显的 3p 缺失/3q 扩增克隆）。因此另按 CopyKAT 自身的层次聚类（r$hclustering）切成 2 簇，
#       含参考细胞比例更高的一簇判为 diploid，另一簇判为 aneuploid（与 CopyKAT 原文的判定逻辑一致），写入 pred_fixed 列。
# 输出：<样本>/copykat_prediction.txt（cell.names, copykat.pred, pred_fixed, is_ref）；运行后删除 CopyKAT 的大中间文件。
# 用法: Rscript 06b_copykat.R <06a 导出目录> [n.cores] [第几份/共几份，如 1/2]
suppressMessages({library(copykat); library(Matrix)})
args <- commandArgs(trailingOnly = TRUE); root <- args[1]; nc <- if (length(args) > 1) as.integer(args[2]) else 8
dirs <- sort(list.dirs(root, recursive = FALSE))
if (length(args) > 2) { pp <- as.integer(strsplit(args[3], "/")[[1]]); dirs <- dirs[seq_along(dirs) %% pp[2] == (pp[1] - 1)] }
for (d in dirs) {
  outf <- file.path(d, "copykat_prediction.txt")
  if (file.exists(outf) && "pred_fixed" %in% names(read.delim(outf, nrows = 1))) next
  m <- readMM(gzfile(file.path(d, "matrix.mtx.gz")))
  rownames(m) <- make.unique(read.csv(gzfile(file.path(d, "genes.tsv.gz")), header = FALSE)[[1]])
  colnames(m) <- read.csv(gzfile(file.path(d, "barcodes.tsv.gz")), header = FALSE)[[1]]
  ref <- read.csv(file.path(d, "ref_cells.txt"), header = FALSE)[[1]]
  set.seed(1)
  epi <- setdiff(colnames(m), ref); ref <- intersect(ref, colnames(m))
  if (length(epi) > 1500) epi <- sample(epi, 1500); if (length(ref) > 500) ref <- sample(ref, 500)
  m <- m[, c(epi, ref)]
  old <- setwd(d); set.seed(1)
  r <- tryCatch(copykat(rawmat = as.matrix(m), id.type = "S", ngene.chr = 5, win.size = 25, KS.cut = 0.1,
                        sam.name = "ck", distance = "euclidean", norm.cell.names = ref, output.seg = FALSE,
                        plot.genes = FALSE, genome = "hg20", n.cores = nc),
                error = function(e) { message("FAILED ", basename(d), ": ", conditionMessage(e)); NULL })
  setwd(old)
  if (!is.null(r)) {
    # 与 inferCNV 相同的单细胞规则，作用于 CopyKAT 的 CNA 矩阵（bin × cell，去掉前 3 列坐标）：
    #   ck_score = 各 bin CNA 值平方均值；ck_corr = 与上皮细胞中 ck_score 最高 5% 的平均谱的 Pearson 相关；
    #   pred_fixed：两项均 > 参考细胞 95 分位 → aneuploid，均 <= 95 分位 → diploid，其余 unresolved。
    # （k=2 层次聚类切分被弃用：它主要按上皮/非上皮表达差异分簇，正常样本 N_HPV_NEG_1 中 1366/1500 上皮被判 aneuploid。）
    X <- as.matrix(r$CNAmat[, -(1:3)]); cn <- colnames(X)
    key <- function(x) gsub("[^A-Za-z0-9]", "", x)
    isref <- key(cn) %in% key(ref)
    score <- colMeans(X^2)
    top <- !isref & score >= quantile(score[!isref], 0.95)
    prof <- rowMeans(X[, top, drop = FALSE])
    corr <- as.vector(cor(X, prof))
    ts <- quantile(score[isref], 0.95); tc <- quantile(corr[isref], 0.95)
    pf <- ifelse(score > ts & corr > tc, "aneuploid", ifelse(score <= ts & corr <= tc, "diploid", "unresolved"))
    pr <- r$prediction; mi <- match(key(pr$cell.names), key(cn))
    pr$is_ref <- key(pr$cell.names) %in% key(ref)
    pr$ck_score <- score[mi]; pr$ck_corr <- corr[mi]; pr$pred_fixed <- pf[mi]
    write.table(pr, outf, sep = "\t", quote = FALSE, row.names = FALSE)
    e <- !pr$is_ref
    cat(sprintf("done %s: epi=%d aneuploid=%d diploid=%d unresolved=%d\n", basename(d), sum(e),
                sum(pr$pred_fixed[e] == "aneuploid", na.rm = TRUE), sum(pr$pred_fixed[e] == "diploid", na.rm = TRUE),
                sum(pr$pred_fixed[e] == "unresolved", na.rm = TRUE)))
  }
  unlink(Sys.glob(file.path(d, "ck_copykat_*results*.txt")))
  rm(m, r); gc()
}
