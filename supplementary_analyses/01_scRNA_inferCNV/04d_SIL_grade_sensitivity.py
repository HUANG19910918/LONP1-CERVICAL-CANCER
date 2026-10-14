# 04d_SIL_grade_sensitivity.py —— GSE333892 SIL 组按病理级别的敏感性分析（倒 U 型是否由高级别 SIL 驱动）
# GEO 未提供 SIL 分级；只能由 orig.ident 推断：
#   明确高级别（HSIL/CIN3）：SIL_11 "CIN.III.1"、SIL_12/13/14 "hsil.*"、SIL_19 "HSIL.2"
#   明确低级别：SIL_05 "ascus..."（ASCUS）
#   未知：其余 13 个（"CIN1/3/4/6/7.<姓名缩写>" 中的数字更像病例编号——存在 "CIN4"，非合法分级；其余仅姓名缩写）
# 另用客观替代指标：SIL 上皮中 inferCNV 恶性样克隆比例 >=20%（CNV-high SIL） vs <20%。
# 比较：各亚组 vs 正常、vs 癌（Mann–Whitney），以及去除高级别/CNV-high 后倒 U 型是否仍成立。
# 用法: python -P 04d_SIL_grade_sensitivity.py <01_scRNA_inferCNV 目录>
import sys, os
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, kruskal, spearmanr
base = sys.argv[1]; res = os.path.join(base, "results")
U = pd.read_csv(os.path.join(res, "GSE333892_pseudobulk_LONP1.csv")); U = U[U.included]
HIGH = {"SIL_11", "SIL_12", "SIL_13", "SIL_14", "SIL_19"}; LOW = {"SIL_05"}
S = U[U.category == "SIL_epi"].copy()
S["label_grade"] = np.where(S["sample"].isin(HIGH), "HSIL_labelled", np.where(S["sample"].isin(LOW), "ASCUS", "unknown"))
S["cnv_class"] = np.where(S.SIL_pct_malignant_like >= 20, "CNV_high", "CNV_low")
N = U[U.category == "normal_epi"].log2CPM_LONP1; T = U[U.category == "malignant_epi"].log2CPM_LONP1
out = []
def L(s): print(s); out.append(s)
L(S[["sample", "label_grade", "cnv_class", "SIL_pct_malignant_like", "n_cells", "log2CPM_LONP1"]].sort_values("log2CPM_LONP1").round(2).to_string(index=False))
L(f"\n正常 n={len(N)} 中位 {N.median():.3f}；癌恶性 n={len(T)} 中位 {T.median():.3f}")
for col, levels in [("label_grade", ["HSIL_labelled", "unknown"]), ("cnv_class", ["CNV_high", "CNV_low"])]:
    for lv in levels:
        v = S[S[col] == lv].log2CPM_LONP1
        if len(v) < 2: continue
        L(f"{col}={lv}: n={len(v)} 中位 {v.median():.3f}; vs 正常 MWU p={mannwhitneyu(v, N).pvalue:.3g}; vs 癌 MWU p={mannwhitneyu(v, T).pvalue:.3g}")
    a, b = [S[S[col] == lv].log2CPM_LONP1 for lv in levels]
    L(f"  {levels[0]} vs {levels[1]}: MWU p={mannwhitneyu(a, b).pvalue:.3g}")
low_only = S[(S.label_grade != "HSIL_labelled") & (S.cnv_class == "CNV_low")].log2CPM_LONP1
L(f"\n去除明确 HSIL 且去除 CNV-high 后的 SIL: n={len(low_only)} 中位 {low_only.median():.3f}; vs 正常 p={mannwhitneyu(low_only, N).pvalue:.3g}; "
  f"vs 癌 p={mannwhitneyu(low_only, T).pvalue:.3g}; 三组 Kruskal p={kruskal(N, low_only, T).pvalue:.3g}")
r = spearmanr(S.SIL_pct_malignant_like, S.log2CPM_LONP1)
L(f"SIL 内 LONP1 与 CNV 克隆比例的 Spearman rho={r.statistic:.3f} p={r.pvalue:.3g}（n={len(S)}）")
S.to_csv(os.path.join(res, "GSE333892_SIL_grade_sensitivity.csv"), index=False)
open(os.path.join(res, "04d_SIL_grade_sensitivity_log.txt"), "w", encoding="utf-8").write("\n".join(out) + "\n")
