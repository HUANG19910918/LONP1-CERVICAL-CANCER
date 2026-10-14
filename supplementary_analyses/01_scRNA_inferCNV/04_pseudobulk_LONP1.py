# 04_pseudobulk_LONP1.py —— LONP1 在恶性/非恶性/正常上皮中的 pseudobulk 与混合模型比较（补充分析 1，阶段 5）
# 单元：样本 × 细胞类别；pseudobulk = 该单元全部细胞原始 counts 求和 → CPM（以该样本全部基因为分母）→ log2(CPM+1)。
# 纳入阈值：每个 pseudobulk 单元 >= MIN_CELLS 个细胞（默认 100，记录于输出）。
# 类别：
#   malignant_epi   = 癌样本中 inferCNV 判为 malignant 的上皮
#   nonmal_epi_tumor= 癌样本中 inferCNV 判为 nonmalignant 的上皮
#   normal_epi      = 正常上皮参照：GSE208653 4 个正常 + E-MTAB 癌旁 sample4/5（patient3 癌旁 sample6 无上皮，排除）；
#                     inferCNV 判为 malignant 的细胞不计入
#   HSIL_epi        = GSE208653 HSIL 上皮（描述性）
# 比较：
#   C1 恶性上皮 vs 正常上皮：样本级 OLS log2CPM ~ group + dataset（固定效应）；另做 Mann–Whitney；
#      细胞级 MixedLM log1p(CP10K) ~ group，随机截距 = 样本。
#   C2 同样本内 恶性 vs 非恶性上皮：配对 Wilcoxon（pseudobulk）+ 细胞级 MixedLM（随机截距 = 样本）。
#   C3 恶性细胞中 19p13.3 增益 vs 非增益：同样本配对 Wilcoxon（pseudobulk）+ 细胞级 MixedLM；另报 LONP1 与 seg_19p13_3 的样本内 Spearman。
#   C4 E-MTAB 配对：patient1/2 癌 恶性上皮 vs 同患者癌旁正常上皮（描述性，n=2）。
# 用法: python -P 04_pseudobulk_LONP1.py <01_scRNA_inferCNV 目录> <qc h5ad 目录> <cnv 目录> [MIN_CELLS]
import sys, os, glob, warnings
import numpy as np, pandas as pd, scanpy as sc
import statsmodels.formula.api as smf
from scipy.stats import wilcoxon, mannwhitneyu, spearmanr
warnings.filterwarnings("ignore")
base, qcd, cnvd = sys.argv[1:4]; MIN_CELLS = int(sys.argv[4]) if len(sys.argv) > 4 else 100
res = os.path.join(base, "results")
NORMAL_REF = {("GSE208653", s) for s in ["GSM6360680_N_HPV_NEG_1", "GSM6360681_N_HPV_NEG_2", "GSM6360682_N_1", "GSM6360683_N_2"]} | \
             {("E-MTAB-11948", "sample4_patient1"), ("E-MTAB-11948", "sample5_patient2")}
log = []
def L(s=""): print(s); log.append(str(s))

units, cells = [], []
for f in sorted(glob.glob(os.path.join(cnvd, "*.csv.gz"))):
    key = os.path.basename(f)[:-7]
    d = pd.read_csv(f).set_index("cell")
    a = sc.read_h5ad(os.path.join(qcd, key + ".h5ad"))
    ds, sm, grp = d.dataset.iloc[0], d["sample"].iloc[0], d.group.iloc[0]
    if "LONP1" not in a.var_names: continue
    tot = np.asarray(a.X.sum(1)).ravel(); lon = np.asarray(a[:, "LONP1"].X.todense()).ravel()
    cell_df = pd.DataFrame({"tot": tot, "lon": lon}, index=a.obs_names).join(d, how="inner")
    cell_df["lonp1_lognorm"] = np.log1p(cell_df.lon / cell_df.tot * 1e4)
    e = cell_df[cell_df.cell_type == "Epithelial"]
    cats = {}
    if grp == "Cancer":
        cats["malignant_epi"] = e[e.epi_class == "malignant"]
        cats["nonmal_epi_tumor"] = e[e.epi_class == "nonmalignant"]
        m = e[e.epi_class == "malignant"]
        cats["mal_gain19p13_3"] = m[m.gain19p13_3]; cats["mal_nogain19p13_3"] = m[~m.gain19p13_3]
    elif (ds, sm) in NORMAL_REF:
        cats["normal_epi"] = e[e.epi_class != "malignant"]
    elif grp == "HSIL":
        cats["HSIL_epi"] = e
    for c, x in cats.items():
        units.append({"dataset": ds, "sample": sm, "group": grp, "hpv": d.hpv.iloc[0], "category": c, "n_cells": len(x),
                      "log2CPM_LONP1": np.log2(x.lon.sum() / max(x.tot.sum(), 1) * 1e6 + 1) if len(x) else np.nan,
                      "mean_lognorm": x.lonp1_lognorm.mean() if len(x) else np.nan,
                      "pct_detected": 100 * (x.lon > 0).mean() if len(x) else np.nan,
                      "median_seg19p13_3": x.seg_19p13_3.median() if len(x) else np.nan})
        if len(x): cells.append(x.assign(category=c)[["dataset", "sample", "category", "lonp1_lognorm", "seg_19p13_3", "cnv_score"]])
U = pd.DataFrame(units); CELLS = pd.concat(cells)
U["included"] = U.n_cells >= MIN_CELLS
U.to_csv(os.path.join(res, "pseudobulk_LONP1_units.csv"), index=False)
L(f"MIN_CELLS per pseudobulk unit = {MIN_CELLS}")
L(U.groupby("category").agg(units=("sample", "size"), included=("included", "sum"), cells=("n_cells", "sum")).to_string())
Ui = U[U.included]

def mixed(df, formula):
    # 依次尝试多种优化器；只接受收敛且 SE 有限的结果，返回 (beta, SE, p, 方法)
    for m in ("lbfgs", "powell", "nm", "bfgs"):
        try:
            r = smf.mixedlm(formula, df, groups=df["sample"]).fit(method=m, reml=True)
            k = [q for q in r.params.index if q not in ("Intercept", "Group Var")][0]
            if r.converged and np.isfinite(r.bse[k]) and r.bse[k] < 10:
                return r.params[k], r.bse[k], r.pvalues[k], m
        except Exception:
            pass
    return np.nan, np.nan, np.nan, "未收敛"

# C1
a1 = Ui[Ui.category.isin(["malignant_epi", "normal_epi"])].copy()
a1["is_mal"] = (a1.category == "malignant_epi").astype(int)
L("\n[C1] 恶性上皮（癌样本）vs 正常上皮参照 —— pseudobulk log2CPM")
L(a1.groupby("category").log2CPM_LONP1.describe().round(3).to_string())
mw = mannwhitneyu(a1[a1.is_mal == 1].log2CPM_LONP1, a1[a1.is_mal == 0].log2CPM_LONP1)
ols = smf.ols("log2CPM_LONP1 ~ is_mal + C(dataset)", a1).fit()
L(f"  Mann–Whitney p={mw.pvalue:.3g}; OLS(+dataset) diff={ols.params['is_mal']:.3f} (SE {ols.bse['is_mal']:.3f}) p={ols.pvalues['is_mal']:.3g}; "
  f"n_mal={int(a1.is_mal.sum())}, n_norm={int((1 - a1.is_mal).sum())}")
for ds in ["GSE208653", "E-MTAB-11948"]:
    s = a1[a1.dataset == ds]
    if s.is_mal.nunique() == 2:
        L(f"  {ds} 内部: mal mean={s[s.is_mal == 1].log2CPM_LONP1.mean():.3f} (n={s.is_mal.sum()}) vs normal {s[s.is_mal == 0].log2CPM_LONP1.mean():.3f} (n={(1 - s.is_mal).sum()})")
c1 = CELLS[CELLS["sample"].isin(a1["sample"]) & CELLS.category.isin(["malignant_epi", "normal_epi"])].assign(is_mal=lambda x: (x.category == "malignant_epi").astype(int))
b, se, p, mm = mixed(c1, "lonp1_lognorm ~ is_mal")
L(f"  细胞级 MixedLM（样本随机截距，{mm}）: beta={b:.3f} SE={se:.3f} p={p:.3g}; cells={len(c1)}")
mw2 = mannwhitneyu(a1[a1.is_mal == 1].mean_lognorm, a1[a1.is_mal == 0].mean_lognorm)
L(f"  样本级 细胞均值 log1p(CP10K): 恶性 中位 {a1[a1.is_mal == 1].mean_lognorm.median():.3f} vs 正常 {a1[a1.is_mal == 0].mean_lognorm.median():.3f}; Mann–Whitney p={mw2.pvalue:.3g}")

# C2
L("\n[C2] 同样本内 恶性 vs 非恶性上皮（癌样本）")
w = Ui[Ui.category.isin(["malignant_epi", "nonmal_epi_tumor"])].pivot_table(index="sample", columns="category", values="log2CPM_LONP1").dropna()
if len(w) >= 3:
    st = wilcoxon(w.malignant_epi, w.nonmal_epi_tumor)
    L(f"  配对样本 n={len(w)}; 中位差(恶性−非恶性)={np.median(w.malignant_epi - w.nonmal_epi_tumor):.3f}; 恶性更高的样本 {int((w.malignant_epi > w.nonmal_epi_tumor).sum())}/{len(w)}; Wilcoxon p={st.pvalue:.3g}")
else:
    L(f"  可配对样本不足（n={len(w)}）")
c2 = CELLS[CELLS["sample"].isin(w.index) & CELLS.category.isin(["malignant_epi", "nonmal_epi_tumor"])].assign(is_mal=lambda x: (x.category == "malignant_epi").astype(int))
b, se, p, mm = mixed(c2, "lonp1_lognorm ~ is_mal"); L(f"  细胞级 MixedLM（{mm}）: beta={b:.3f} SE={se:.3f} p={p:.3g}; cells={len(c2)}")
w2m = Ui[Ui.category.isin(["malignant_epi", "nonmal_epi_tumor"])].pivot_table(index="sample", columns="category", values="mean_lognorm").dropna()
if len(w2m) >= 3:
    st = wilcoxon(w2m.malignant_epi, w2m.nonmal_epi_tumor)
    L(f"  配对 细胞均值 log1p(CP10K): 中位差={np.median(w2m.malignant_epi - w2m.nonmal_epi_tumor):.3f}; 恶性更高 {int((w2m.malignant_epi > w2m.nonmal_epi_tumor).sum())}/{len(w2m)}; Wilcoxon p={st.pvalue:.3g}")

# C3
L("\n[C3] 恶性细胞中 19p13.3 增益 vs 非增益")
w3 = Ui[Ui.category.isin(["mal_gain19p13_3", "mal_nogain19p13_3"])].pivot_table(index="sample", columns="category", values="log2CPM_LONP1").dropna()
if len(w3) >= 3:
    st = wilcoxon(w3.mal_gain19p13_3, w3.mal_nogain19p13_3)
    L(f"  配对样本 n={len(w3)}; 中位差(增益−非增益)={np.median(w3.mal_gain19p13_3 - w3.mal_nogain19p13_3):.3f}; 增益更高 {int((w3.mal_gain19p13_3 > w3.mal_nogain19p13_3).sum())}/{len(w3)}; Wilcoxon p={st.pvalue:.3g}")
else:
    L(f"  可配对样本不足（n={len(w3)}）")
c3 = CELLS[CELLS.category.isin(["mal_gain19p13_3", "mal_nogain19p13_3"])].assign(gain=lambda x: (x.category == "mal_gain19p13_3").astype(int))
b, se, p, mm = mixed(c3, "lonp1_lognorm ~ gain"); L(f"  细胞级 MixedLM（{mm}）: beta={b:.3f} SE={se:.3f} p={p:.3g}; cells={len(c3)}")
MIN_C3 = 30
w3b = U[U.category.isin(["mal_gain19p13_3", "mal_nogain19p13_3"]) & (U.n_cells >= MIN_C3)].pivot_table(index="sample", columns="category", values="mean_lognorm").dropna()
if len(w3b) >= 3:
    st = wilcoxon(w3b.mal_gain19p13_3, w3b.mal_nogain19p13_3)
    L(f"  敏感性（每单元 >= {MIN_C3} 细胞，细胞均值）: 配对 n={len(w3b)}; 中位差={np.median(w3b.mal_gain19p13_3 - w3b.mal_nogain19p13_3):.3f}; 增益更高 {int((w3b.mal_gain19p13_3 > w3b.mal_nogain19p13_3).sum())}/{len(w3b)}; Wilcoxon p={st.pvalue:.3g}")
else:
    L(f"  敏感性（>= {MIN_C3} 细胞）配对样本 n={len(w3b)}，不足 3")
mal = CELLS[CELLS.category == "malignant_epi"]
rs = mal.groupby("sample").apply(lambda x: spearmanr(x.lonp1_lognorm, x.seg_19p13_3).statistic if len(x) >= MIN_CELLS else np.nan).dropna()
L(f"  样本内 LONP1 ~ seg_19p13_3 Spearman: 中位 {rs.median():.3f}，范围 {rs.min():.3f}–{rs.max():.3f}，>0 的样本 {int((rs > 0).sum())}/{len(rs)}")

# C4
L("\n[C4] E-MTAB-11948 同患者 癌恶性上皮 vs 癌旁正常上皮")
em = U[U.dataset == "E-MTAB-11948"]
for pt in ["patient1", "patient2", "patient3"]:
    x = em[em["sample"].str.endswith(pt)]
    L("  " + pt + ": " + "; ".join(f"{r.category} n={r.n_cells} log2CPM={r.log2CPM_LONP1:.3f}" for r in x.itertuples()))
L("\n[HSIL]"); L(U[U.category == "HSIL_epi"][["sample", "n_cells", "log2CPM_LONP1"]].round(3).to_string(index=False))
open(os.path.join(res, "04_pseudobulk_LONP1_log.txt"), "w", encoding="utf-8").write("\n".join(log) + "\n")
