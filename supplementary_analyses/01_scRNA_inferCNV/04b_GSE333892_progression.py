# 04b_GSE333892_progression.py —— GSE333892 队列内 正常→SIL→癌 的上皮 LONP1 比较 + 六个数据集合并的恶性 vs 正常上皮比较
# pseudobulk 定义同 04 号脚本（原始 counts 求和 → CPM → log2(CPM+1)；每单元 >= MIN_CELLS 细胞，默认 100）。
# GSE333892 类别：normal_epi = CON 样本中 inferCNV 非 malignant 的上皮；SIL_epi = SIL 样本全部上皮；
#                 malignant_epi = CC 样本 inferCNV malignant 上皮。CON_02 无上皮（已排除）。
# 注意：CON_01–05 的原始名为 P1/P2/P3/P4/P6.Normal，其中 P2/P3/P4 与 CC_09/10/11（P2/P3/P4.Tumor）为同患者癌旁，
#       故 CON 组并非全部健康供体；配对比较（P3、P4）单列。
# 统计：三组 Kruskal–Wallis；趋势 = 样本级 Spearman（组序 0/1/2 vs log2CPM）；两两 Mann–Whitney；
#       细胞级 MixedLM lonp1 ~ 组序（样本随机截距）。
# 合并 C1：六数据集正常上皮（GSE208653 4 + E-MTAB 2 + GSE333892 CON）vs 恶性上皮，OLS log2CPM ~ is_mal + dataset。
# 用法: python -P 04b_GSE333892_progression.py <01_scRNA_inferCNV 目录> <qc 目录> <cnv 目录> [MIN_CELLS]
import sys, os, glob, warnings
import numpy as np, pandas as pd, scanpy as sc
import statsmodels.formula.api as smf
from scipy.stats import kruskal, mannwhitneyu, spearmanr
warnings.filterwarnings("ignore")
base, qcd, cnvd = sys.argv[1:4]; MIN = int(sys.argv[4]) if len(sys.argv) > 4 else 100
res = os.path.join(base, "results"); log = []
def L(s=""): print(s); log.append(str(s))
units, cells = [], []
for f in sorted(glob.glob(os.path.join(cnvd, "GSE333892__*.csv.gz"))):
    k = os.path.basename(f)[:-7]; d = pd.read_csv(f).set_index("cell")
    a = sc.read_h5ad(os.path.join(qcd, k + ".h5ad"))
    tot = np.asarray(a.X.sum(1)).ravel(); lon = np.asarray(a[:, "LONP1"].X.todense()).ravel()
    c = pd.DataFrame({"tot": tot, "lon": lon}, index=a.obs_names).join(d, how="inner")
    c["lonp1"] = np.log1p(c.lon / c.tot * 1e4); e = c[c.cell_type == "Epithelial"]
    g = d.group.iloc[0]; sm = d["sample"].iloc[0]
    cat = {"Normal": ("normal_epi", e[e.epi_class != "malignant"]), "SIL": ("SIL_epi", e),
           "Cancer": ("malignant_epi", e[e.epi_class == "malignant"])}[g]
    x = cat[1]
    units.append({"sample": sm, "group": g, "category": cat[0], "n_cells": len(x),
                  "log2CPM_LONP1": np.log2(x.lon.sum() / max(x.tot.sum(), 1) * 1e6 + 1) if len(x) else np.nan,
                  "mean_lognorm": x.lonp1.mean() if len(x) else np.nan,
                  "SIL_pct_malignant_like": 100 * (e.epi_class == "malignant").mean() if g == "SIL" else np.nan})
    cells.append(x.assign(sample=sm, category=cat[0])[["sample", "category", "lonp1"]])
U = pd.DataFrame(units); U["included"] = U.n_cells >= MIN
U.to_csv(os.path.join(res, "GSE333892_pseudobulk_LONP1.csv"), index=False)
Ui = U[U.included].copy(); order = {"normal_epi": 0, "SIL_epi": 1, "malignant_epi": 2}; Ui["ord"] = Ui.category.map(order)
L(f"MIN_CELLS={MIN}; 纳入单元: " + ", ".join(f"{k}={int(v)}" for k, v in Ui.category.value_counts().items()) +
  f"; 排除（细胞不足）: {', '.join(U[~U.included]['sample'])}")
L(Ui.groupby("category").log2CPM_LONP1.describe().round(3).to_string())
kw = kruskal(*[Ui[Ui.category == c].log2CPM_LONP1 for c in order])
tr = spearmanr(Ui.ord, Ui.log2CPM_LONP1)
L(f"Kruskal–Wallis p={kw.pvalue:.3g}; 趋势 Spearman rho={tr.statistic:.3f} p={tr.pvalue:.3g}")
for a_, b_ in [("normal_epi", "SIL_epi"), ("SIL_epi", "malignant_epi"), ("normal_epi", "malignant_epi")]:
    p = mannwhitneyu(Ui[Ui.category == a_].log2CPM_LONP1, Ui[Ui.category == b_].log2CPM_LONP1).pvalue
    L(f"  {a_} vs {b_}: 中位 {Ui[Ui.category == a_].log2CPM_LONP1.median():.3f} vs {Ui[Ui.category == b_].log2CPM_LONP1.median():.3f}, MWU p={p:.3g}")
CL = pd.concat(cells); CL = CL[CL["sample"].isin(Ui["sample"])]; CL["ord"] = CL.category.map(order)
try:
    r = smf.mixedlm("lonp1 ~ ord", CL, groups=CL["sample"]).fit(method="powell")
    L(f"细胞级 MixedLM lonp1 ~ 组序: beta={r.params['ord']:.4f} SE={r.bse['ord']:.4f} p={r.pvalues['ord']:.3g}; cells={len(CL)}")
except Exception as ex:
    L(f"MixedLM 失败: {ex}")
L("\n[同患者配对] P3: CON_03 vs CC_10；P4: CON_04 vs CC_11（不受 MIN 限制，描述性）")
for n_, t_ in [("CON_03", "CC_10"), ("CON_04", "CC_11")]:
    a_ = U[U["sample"] == n_].iloc[0]; b_ = U[U["sample"] == t_].iloc[0]
    L(f"  {n_} 正常上皮 n={a_.n_cells} log2CPM={a_.log2CPM_LONP1:.3f}  vs  {t_} 恶性上皮 n={b_.n_cells} log2CPM={b_.log2CPM_LONP1:.3f}")
# 合并 C1
P = pd.read_csv(os.path.join(res, "pseudobulk_LONP1_units.csv"))
P = P[P.included & P.category.isin(["malignant_epi", "normal_epi"]) & (P.dataset != "GSE333892")][["dataset", "sample", "category", "log2CPM_LONP1"]]
Q = Ui[Ui.category.isin(["malignant_epi", "normal_epi"])].assign(dataset="GSE333892")[["dataset", "sample", "category", "log2CPM_LONP1"]]
A = pd.concat([P, Q]); A["is_mal"] = (A.category == "malignant_epi").astype(int)
A.to_csv(os.path.join(res, "pseudobulk_LONP1_pooled_6datasets.csv"), index=False)
ols = smf.ols("log2CPM_LONP1 ~ is_mal + C(dataset)", A).fit()
mw = mannwhitneyu(A[A.is_mal == 1].log2CPM_LONP1, A[A.is_mal == 0].log2CPM_LONP1)
L(f"\n[合并 C1，6 数据集] 恶性 n={A.is_mal.sum()} 中位 {A[A.is_mal == 1].log2CPM_LONP1.median():.3f} vs 正常 n={(1 - A.is_mal).sum()} 中位 "
  f"{A[A.is_mal == 0].log2CPM_LONP1.median():.3f}; OLS(+dataset) diff={ols.params['is_mal']:.3f} (SE {ols.bse['is_mal']:.3f}, 95% CI "
  f"{ols.conf_int().loc['is_mal', 0]:.3f}–{ols.conf_int().loc['is_mal', 1]:.3f}) p={ols.pvalues['is_mal']:.3g}; Mann–Whitney p={mw.pvalue:.3g}")
for ds, g in A.groupby("dataset"):
    if g.is_mal.nunique() == 2:
        L(f"  {ds}: 恶性 n={g.is_mal.sum()} 中位 {g[g.is_mal == 1].log2CPM_LONP1.median():.3f} vs 正常 n={(1 - g.is_mal).sum()} 中位 {g[g.is_mal == 0].log2CPM_LONP1.median():.3f}"
          f"; MWU p={mannwhitneyu(g[g.is_mal == 1].log2CPM_LONP1, g[g.is_mal == 0].log2CPM_LONP1).pvalue:.3g}")
open(os.path.join(res, "04b_GSE333892_progression_log.txt"), "w", encoding="utf-8").write("\n".join(log) + "\n")
