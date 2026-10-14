# 05_combine_mr.py —— 合并 16 个基因的 MR 结果，按"每个暴露来源 × 结局"内跨基因做 BH-FDR 与 Bonferroni 校正
# 主分析定义：暴露 = eQTLGen 血液（最大样本），结局 = FinnGen 宫颈癌（发现），BBJ 作重复；其余为次要/敏感性。
# 用法: python -P 05_combine_mr.py <03_MR_coloc 目录>
import sys, os, glob
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
base = sys.argv[1]; res = os.path.join(base, "results")
A = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(res, "mr_*.csv"))) if "allgenes" not in f])
A = A[A.outcome != "-"].dropna(subset=["OR"]).reset_index(drop=True)
A["fdr"] = np.nan; A["bonf"] = np.nan
for (e, o), g in A.groupby(["exposure", "outcome"]):
    A.loc[g.index, "fdr"] = multipletests(g.p, method="fdr_bh")[1]; A.loc[g.index, "bonf"] = np.minimum(g.p * len(g), 1)
A.to_csv(os.path.join(res, "mr_all_genes_combined.csv"), index=False)
print("暴露×结局×基因 估计数:", len(A), "; 有工具变量的基因数（eQTLGen）:", A[A.exposure == "eQTLGen_blood"].gene.nunique())
print("全部估计中 名义 P<0.05:", int((A.p < 0.05).sum()), "; FDR<0.05:", int((A.fdr < 0.05).sum()))
pd.set_option("display.width", 220)
cols = ["gene", "exposure", "outcome", "n_snp", "OR", "OR_lo", "OR_hi", "p", "fdr"]
print("\n[主分析] eQTLGen 血液 → FinnGen 宫颈癌"); print(A[(A.exposure == "eQTLGen_blood") & (A.outcome == "FinnGen_CC")].sort_values("p")[cols].round(4).to_string(index=False))
print("\n[重复] eQTLGen 血液 → BBJ 宫颈癌"); print(A[(A.exposure == "eQTLGen_blood") & (A.outcome == "BBJ_CeC")].sort_values("p")[cols].round(4).to_string(index=False))
print("\n名义 P<0.05 的全部估计："); print(A[A.p < 0.05].sort_values("p")[cols + ["method", "egger_int_p", "Q_p"]].round(4).to_string(index=False))
