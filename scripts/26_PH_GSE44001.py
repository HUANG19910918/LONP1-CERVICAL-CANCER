import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lifelines import CoxPHFitter, CoxTimeVaryingFitter
from lifelines.statistics import proportional_hazard_test
from lifelines.utils import to_episodic_format
from scipy.stats import chi2 as chi2d
rd = sys.argv[1]
if os.path.isdir(os.path.join(rd, "raw_data")):   # 允许传入分析根目录或 raw_data 目录
    rd = os.path.join(rd, "raw_data")
L=[]
def P(s): print(s,flush=True); L.append(str(s))
d = pd.read_csv(os.path.join(rd,"44_GSE44001_LONP1_DFS.csv"))
x = d.dropna(subset=["DFS_months","DFS_event","LONP1_log2"]); x = x[x.DFS_months>0].copy()
x["stage_adv"] = (~x.stage.isin(["IA1","IA2","IB1"])).astype(int)
P("GSE44001 external cohort: proportional hazards check")
P(f"n={len(x)}, events={int(x.DFS_event.sum())}")
P("")
# univariate
df = x[["DFS_months","DFS_event","LONP1_log2"]].reset_index(drop=True)
df["lon_c"] = df.LONP1_log2 - df.LONP1_log2.mean()
cf = CoxPHFitter().fit(df[["DFS_months","DFS_event","lon_c"]], "DFS_months","DFS_event")
s0 = cf.summary.loc["lon_c"]
P(f"[univariate DFS] HR={np.exp(s0['coef']):.2f} ({np.exp(s0['coef lower 95%']):.2f}-{np.exp(s0['coef upper 95%']):.2f}), p={s0['p']:.3f}")
for tr in ("rank","km"):
    r = proportional_hazard_test(cf, df[["DFS_months","DFS_event","lon_c"]], time_transform=tr).summary
    P(f"   Schoenfeld ({tr}): chi2={r.test_statistic.iloc[0]:.3f}, p={r.p.iloc[0]:.4f}")
# multivariable
xm = x.dropna(subset=["diameter_cm"]).reset_index(drop=True)
dm = xm[["DFS_months","DFS_event","LONP1_log2","stage_adv","diameter_cm"]].copy()
dm["LONP1_log2"] = dm.LONP1_log2 - dm.LONP1_log2.mean()
cfm = CoxPHFitter().fit(dm, "DFS_months","DFS_event")
P("")
P(f"[multivariable DFS] n={len(dm)}, events={int(dm.DFS_event.sum())}")
for tr in ("rank","km"):
    r = proportional_hazard_test(cfm, dm, time_transform=tr).summary.reset_index()
    for _, row in r.iterrows():
        P(f"   Schoenfeld ({tr}) {row['index']:<14} chi2={row['test_statistic']:.3f}  p={row['p']:.4f}")
    P(f"   Schoenfeld ({tr}) GLOBAL         chi2={r.test_statistic.sum():.3f} df={len(r)} p={chi2d.sf(r.test_statistic.sum(), len(r)):.4f}")
# time-varying for LONP1
P("")
dd = df[["DFS_months","DFS_event","lon_c"]].copy(); dd["id"]=np.arange(len(dd))
long = to_episodic_format(dd, duration_col="DFS_months", event_col="DFS_event", id_col="id")
long = long[long["stop"]>0].copy()
C = float(np.log(df.loc[df.DFS_event==1,"DFS_months"]).mean())
long["gt"] = np.log(long["stop"]) - C
long["lon_x_gt"] = long["lon_c"]*long["gt"]
ctv = CoxTimeVaryingFitter().fit(long[["id","start","stop","DFS_event","lon_c","lon_x_gt"]],
        id_col="id", event_col="DFS_event", start_col="start", stop_col="stop", show_progress=False)
s = ctv.summary; V = ctv.variance_matrix_
b0,b1 = s.loc["lon_c","coef"], s.loc["lon_x_gt","coef"]
P(f"[time-varying DFS]  c={C:.3f} (t={np.exp(C):.1f} mo)")
P(f"   b1 (LONP1 x log t) = {b1:+.4f} (se {s.loc['lon_x_gt','se(coef)']:.4f}, p={s.loc['lon_x_gt','p']:.4f})")
P(f"   {'month':>6} {'HR(t)':>8} {'95% CI':>18} {'n_at_risk':>10}")
for t in [6,12,24,36,60,120]:
    g = np.log(t)-C; est = b0+b1*g
    var = V.iloc[0,0] + (g**2)*V.iloc[1,1] + 2*g*V.iloc[0,1]
    se = np.sqrt(max(var,0))
    P(f"   {t:>6} {np.exp(est):>8.2f} {'%.2f-%.2f'%(np.exp(est-1.96*se),np.exp(est+1.96*se)):>18} {int((df.DFS_months>=t).sum()):>10}")
open(os.path.join(rd,"50_PH_GSE44001.txt"),"w",encoding="utf-8").write("\n".join(L)+"\n")
print("WROTE 50_PH_GSE44001.txt")
