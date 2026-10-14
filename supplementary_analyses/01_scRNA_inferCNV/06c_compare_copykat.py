# 06c_compare_copykat.py —— inferCNV（主判定 q95）与 CopyKAT（同一 q95 规则作用于 CopyKAT CNA 矩阵）的逐样本一致性
# 只比较两法都给出明确判定的上皮细胞（inferCNV malignant/nonmalignant × CopyKAT aneuploid/diploid）；
# 报告：一致率、Cohen's kappa、各法的恶性比例；正常组织样本中两法的"恶性"比例（假阳性校准）。
# 另做 LONP1 敏感性：以"两法一致判恶性"的细胞为恶性上皮，重复 C1（恶性 vs 正常上皮参照）的 pseudobulk 比较。
# 用法: python -P 06c_compare_copykat.py <01_scRNA_inferCNV 目录> <cnv 目录> <copykat 目录> <qc 目录>
import sys, os, glob, re
import numpy as np, pandas as pd, scanpy as sc
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu, spearmanr
base, cnvd, ckd, qcd = sys.argv[1:5]; res = os.path.join(base, "results")
key = lambda s: re.sub(r"[^A-Za-z0-9]", "", s)
rows, units = [], []
NORMAL_REF = {"GSM6360680_N_HPV_NEG_1", "GSM6360681_N_HPV_NEG_2", "GSM6360682_N_1", "GSM6360683_N_2", "sample4_patient1", "sample5_patient2"}
for pf in sorted(glob.glob(os.path.join(ckd, "*", "copykat_prediction.txt"))):
    k = os.path.basename(os.path.dirname(pf))
    ck = pd.read_csv(pf, sep="\t"); ck = ck[~ck.is_ref]; ck["k"] = ck["cell.names"].map(key)
    d = pd.read_csv(os.path.join(cnvd, k + ".csv.gz")); d["k"] = d.cell.map(key)
    m = d.merge(ck[["k", "pred_fixed", "ck_score"]], on="k")
    m = m[m.cell_type == "Epithelial"]
    both = m[m.epi_class.isin(["malignant", "nonmalignant"]) & m.pred_fixed.isin(["aneuploid", "diploid"])]
    a = (both.epi_class == "malignant").values; b = (both.pred_fixed == "aneuploid").values
    po = (a == b).mean() if len(both) else np.nan
    pe = a.mean() * b.mean() + (1 - a.mean()) * (1 - b.mean()) if len(both) else np.nan
    kap = (po - pe) / (1 - pe) if len(both) and pe < 1 else np.nan
    rows.append({"dataset": d.dataset.iloc[0], "sample": d["sample"].iloc[0], "group": d.group.iloc[0], "epi_compared": len(m),
                 "both_definite": len(both), "agreement": round(po, 3), "kappa": round(kap, 3),
                 "pct_mal_infercnv": round(100 * (m.epi_class == "malignant").mean(), 1),
                 "pct_aneuploid_copykat": round(100 * (m.pred_fixed == "aneuploid").mean(), 1),
                 "pct_concordant_malignant": round(100 * ((m.epi_class == "malignant") & (m.pred_fixed == "aneuploid")).mean(), 1),
                 "spearman_cnvscore_vs_ckscore": round(spearmanr(m.cnv_score, m.ck_score, nan_policy="omit").statistic, 3) if m.ck_score.notna().sum() > 20 else np.nan})
    # LONP1 pseudobulk with concordant malignant cells
    grp = d.group.iloc[0]; sm = d["sample"].iloc[0]
    if grp == "Cancer":
        sel = m[(m.epi_class == "malignant") & (m.pred_fixed == "aneuploid")].cell; cat = "malignant_concordant"
    elif sm in NORMAL_REF:
        sel = m[(m.epi_class != "malignant") & (m.pred_fixed != "aneuploid")].cell; cat = "normal_concordant"
    else:
        continue
    ad = sc.read_h5ad(os.path.join(qcd, k + ".h5ad"))[sel.values]
    tot = float(ad.X.sum()); lon = float(ad[:, "LONP1"].X.sum())
    units.append({"dataset": d.dataset.iloc[0], "sample": sm, "category": cat, "n_cells": len(sel),
                  "log2CPM_LONP1": np.log2(lon / max(tot, 1) * 1e6 + 1)})
R = pd.DataFrame(rows); R.to_csv(os.path.join(res, "copykat_vs_infercnv.csv"), index=False)
print(R.to_string(index=False))
c = R[R.group == "Cancer"]; n = R[R.group.isin(["Normal", "AdjacentNormal"])]
print(f"\n癌样本 n={len(c)}: 一致率中位 {c.agreement.median():.3f}，kappa 中位 {c.kappa.median():.3f}；"
      f"恶性比例 inferCNV 中位 {c.pct_mal_infercnv.median():.1f}% vs CopyKAT {c.pct_aneuploid_copykat.median():.1f}%")
print(f"上皮细胞 inferCNV cnv_score 与 CopyKAT ck_score 的样本内 Spearman：癌样本中位 {c.spearman_cnvscore_vs_ckscore.median():.3f}（范围 {c.spearman_cnvscore_vs_ckscore.min():.3f}–{c.spearman_cnvscore_vs_ckscore.max():.3f}）")
print(f"正常/癌旁 n={len(n)}: inferCNV 恶性比例 {n.pct_mal_infercnv.min():.1f}–{n.pct_mal_infercnv.max():.1f}%；CopyKAT {n.pct_aneuploid_copykat.min():.1f}–{n.pct_aneuploid_copykat.max():.1f}%")
U = pd.DataFrame(units); U.to_csv(os.path.join(res, "pseudobulk_LONP1_concordant.csv"), index=False)
Ui = U[U.n_cells >= 100].copy(); Ui["is_mal"] = (Ui.category == "malignant_concordant").astype(int)
if Ui.is_mal.nunique() == 2:
    mw = mannwhitneyu(Ui[Ui.is_mal == 1].log2CPM_LONP1, Ui[Ui.is_mal == 0].log2CPM_LONP1)
    ols = smf.ols("log2CPM_LONP1 ~ is_mal + C(dataset)", Ui).fit()
    print(f"\n[敏感性 C1，两法一致] 恶性 n={Ui.is_mal.sum()} 中位 {Ui[Ui.is_mal == 1].log2CPM_LONP1.median():.3f} vs 正常 n={(1 - Ui.is_mal).sum()} "
          f"中位 {Ui[Ui.is_mal == 0].log2CPM_LONP1.median():.3f}; Mann–Whitney p={mw.pvalue:.3g}; OLS(+dataset) diff={ols.params['is_mal']:.3f} p={ols.pvalues['is_mal']:.3g}")
