# Schoenfeld residual test of the proportional hazards assumption
# for the four Cox models reported in the manuscript (TCGA-CESC).
import os, sys
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

rd = sys.argv[1]
if os.path.isdir(os.path.join(rd, "raw_data")):   # 允许传入分析根目录或 raw_data 目录
    rd = os.path.join(rd, "raw_data")
d2 = pd.read_csv(os.path.join(rd, "02_TCGA_CESC_clinical_LONP1_survival.csv"))
L = []
def P(s):
    print(s, flush=True); L.append(s)

P("Proportional hazards assumption: Schoenfeld residual tests (lifelines proportional_hazard_test)")
P("Data: 02_TCGA_CESC_clinical_LONP1_survival.csv; models reproduce 09_figs_mito_surv_immune_plotnine.py")
P("")

def report(cf, df, tag):
    for tr in ("rank", "km"):
        res = proportional_hazard_test(cf, df, time_transform=tr)
        s = res.summary.reset_index()
        cols = [c for c in s.columns if c in ("index", "test_statistic", "p")]
        P(f"[{tag}] time_transform={tr}")
        for _, r in s.iterrows():
            P(f"    {r['index']:<20s} chi2={r['test_statistic']:.3f}  p={r['p']:.4f}")
        # global: combine per-covariate chi2 (df = n covariates)
        from scipy.stats import chi2 as chi2d
        gstat = s.test_statistic.sum(); gdf = len(s)
        P(f"    GLOBAL               chi2={gstat:.3f}  df={gdf}  p={chi2d.sf(gstat, gdf):.4f}")
    P("")

# --- univariate models (OS / DSS / PFS), continuous LONP1 ---
for tk, ek, lab in [("OS_months","OS_event","Univariate OS"),
                    ("DSS_months","DSS_event","Univariate DSS"),
                    ("PFS_months","PFS_event","Univariate PFS")]:
    x = d2.dropna(subset=[tk, ek]); x = x[x[tk] > 0].copy()
    df = x[[tk, ek, "LONP1_log2RSEM"]]
    cf = CoxPHFitter().fit(df, tk, ek)
    P(f"=== {lab}: n={len(df)}, events={int(df[ek].sum())} ===")
    report(cf, df, lab)

# --- multivariable OS model ---
x = d2.dropna(subset=["OS_months","OS_event"]); x = x[x.OS_months > 0].copy()
x["Tg"] = x.path_T.astype(str).str.extract(r"^(T[1-4])")[0]
x = x.dropna(subset=["Tg","age"])
x = x[x.path_N.isin(["N0","N1"])].copy()
x["T2"] = (x.Tg == "T2").astype(int)
x["T34"] = x.Tg.isin(["T3","T4"]).astype(int)
x["N1"] = (x.path_N == "N1").astype(int)
df = x[["OS_months","OS_event","LONP1_log2RSEM","age","T2","T34","N1"]]
cf = CoxPHFitter().fit(df, "OS_months", "OS_event")
P(f"=== Multivariable OS: n={len(df)}, events={int(df.OS_event.sum())} ===")
report(cf, df, "Multivariable OS")

out = os.path.join(rd, "47_PH_assumption_schoenfeld.txt")
open(out, "w", encoding="utf-8").write("\n".join(L) + "\n")
print("WROTE", out)
