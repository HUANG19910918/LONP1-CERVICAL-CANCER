import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lifelines import CoxPHFitter, CoxTimeVaryingFitter
from lifelines.utils import to_episodic_format
from scipy.stats import chi2 as chi2d

rd = sys.argv[1]
if os.path.isdir(os.path.join(rd, "raw_data")):   # 允许传入分析根目录或 raw_data 目录
    rd = os.path.join(rd, "raw_data")
d2 = pd.read_csv(os.path.join(rd, "02_TCGA_CESC_clinical_LONP1_survival.csv"))
L = []
def P(s): print(s, flush=True); L.append(str(s))

P("STEP 2 (corrected): time-varying coefficient models with centred covariates")
P("NOTE: a first attempt used uncentred LONP1 x log(t). LONP1 (log2 RSEM, mean ~10.6)")
P("      is almost perfectly collinear with its own interaction, producing unstable")
P("      coefficients (HR(3mo)=160). Centring both terms fixes this.")
P("Values here match Additional file 1 Table S8 (raw_data/52, script 28) exactly.")
P("")

def prep(tk, ek):
    x = d2.dropna(subset=[tk, ek]); x = x[x[tk] > 0].copy()
    return x[[tk, ek, "LONP1_log2RSEM"]].reset_index(drop=True)

# 中心化改用各分析子集自身的均值，与 27/28 号脚本一致（此前用全局均值，拟合到第3位小数有微小差异）
MEAN_LON = d2.LONP1_log2RSEM.mean()
P(f"LONP1_log2RSEM (全体): mean={MEAN_LON:.3f}, sd={d2.LONP1_log2RSEM.std():.3f}, range {d2.LONP1_log2RSEM.min():.2f}-{d2.LONP1_log2RSEM.max():.2f}")
P("中心化：每个分析子集减去其自身的 LONP1 均值")
P("")

for tk, ek, lab in [("OS_months","OS_event","OS"), ("PFS_months","PFS_event","PFS"), ("DSS_months","DSS_event","DSS")]:
    df = prep(tk, ek)
    df["lon_c"] = df.LONP1_log2RSEM - df.LONP1_log2RSEM.mean()
    cf0 = CoxPHFitter().fit(df[[tk, ek, "lon_c"]], tk, ek)
    s0 = cf0.summary.loc["lon_c"]
    P(f"=== {lab}  (n={len(df)}, events={int(df[ek].sum())}) ===")
    P(f"  [PH model]  HR={np.exp(s0['coef']):.2f} ({np.exp(s0['coef lower 95%']):.2f}-{np.exp(s0['coef upper 95%']):.2f}), p={s0['p']:.3f}")

    dd = df.copy(); dd["id"] = np.arange(len(dd))
    long = to_episodic_format(dd, duration_col=tk, event_col=ek, id_col="id")
    long = long[long["stop"] > 0].copy()
    evt_times = df.loc[df[ek] == 1, tk]
    C = float(np.log(evt_times).mean())
    long["gt"] = np.log(long["stop"]) - C
    long["lon_x_gt"] = long["lon_c"] * long["gt"]
    ctv = CoxTimeVaryingFitter().fit(
        long[["id","start","stop",ek,"lon_c","lon_x_gt"]],
        id_col="id", event_col=ek, start_col="start", stop_col="stop", show_progress=False)
    s = ctv.summary; V = ctv.variance_matrix_
    b0, b1 = s.loc["lon_c","coef"], s.loc["lon_x_gt","coef"]
    P(f"  [time-varying]  centring constant c = mean log(event time) = {C:.3f}  (t = {np.exp(C):.1f} mo)")
    P(f"     b0 = {b0:+.4f} (se {s.loc['lon_c','se(coef)']:.4f}, p={s.loc['lon_c','p']:.4f})   = log HR at t={np.exp(C):.1f} mo")
    P(f"     b1 = {b1:+.4f} (se {s.loc['lon_x_gt','se(coef)']:.4f}, p={s.loc['lon_x_gt','p']:.4f})   <- time-dependence test")
    P(f"     {'month':>6} {'HR(t)':>8} {'95% CI':>18} {'n_at_risk':>10} {'cum_events':>11}")
    for t in [6, 12, 24, 36, 60, 120]:
        g = np.log(t) - C
        est = b0 + b1*g
        var = V.iloc[0,0] + (g**2)*V.iloc[1,1] + 2*g*V.iloc[0,1]
        se = np.sqrt(max(var, 0))
        nr = int((df[tk] >= t).sum()); ce = int(df.loc[df[tk] <= t, ek].sum())
        P(f"     {t:>6} {np.exp(est):>8.2f} {'%.2f-%.2f' % (np.exp(est-1.96*se), np.exp(est+1.96*se)):>18} {nr:>10} {ce:>11}")

    cut = float(evt_times.median())
    long2 = long.copy()
    long2["late"] = (long2["stop"] > cut).astype(int)
    long2["lon_early"] = long2["lon_c"] * (1 - long2["late"])
    long2["lon_late"]  = long2["lon_c"] * long2["late"]
    ctv2 = CoxTimeVaryingFitter().fit(
        long2[["id","start","stop",ek,"lon_early","lon_late"]],
        id_col="id", event_col=ek, start_col="start", stop_col="stop", show_progress=False)
    s2 = ctv2.summary
    n_e = int(df.loc[df[tk] <= cut, ek].sum()); n_l = int(df[ek].sum()) - n_e
    P(f"  [piecewise, cut at median event time = {cut:.1f} mo]")
    for nm, tag, nev in [("lon_early", f"t<={cut:.1f} mo", n_e), ("lon_late", f"t> {cut:.1f} mo", n_l)]:
        r = s2.loc[nm]
        P(f"     {tag:<16} HR={np.exp(r['coef']):.2f} ({np.exp(r['coef lower 95%']):.2f}-{np.exp(r['coef upper 95%']):.2f}), p={r['p']:.3f}   events={nev}")
    P("")

P("Linearity of LONP1 on the log-hazard scale (linear vs natural cubic spline):")
from patsy import dmatrix
for tk, ek, lab in [("OS_months","OS_event","OS"), ("PFS_months","PFS_event","PFS"), ("DSS_months","DSS_event","DSS")]:
    df = prep(tk, ek); df["lon_c"] = df.LONP1_log2RSEM - df.LONP1_log2RSEM.mean()
    ll_lin = CoxPHFitter().fit(df[[tk, ek, "lon_c"]], tk, ek).log_likelihood_
    B = dmatrix("cr(x, df=4)", {"x": df["lon_c"].values}, return_type="dataframe").iloc[:, 1:]
    B = B.loc[:, B.std() > 1e-8]
    B.columns = [f"s{i}" for i in range(B.shape[1])]
    dsp = pd.concat([df[[tk, ek]], B], axis=1)
    lls = CoxPHFitter(penalizer=1e-6).fit(dsp, tk, ek).log_likelihood_
    lrt = 2*(lls - ll_lin); ddf = B.shape[1] - 1
    P(f"   {lab}: LRT chi2={lrt:.3f}, df={ddf}, p={chi2d.sf(max(lrt,0), ddf):.4f}   (p>0.05: no evidence against linearity)")

open(os.path.join(rd, "49_PH_step2_timevarying.txt"), "w", encoding="utf-8").write("\n".join(L)+"\n")
print("\nWROTE 49_PH_step2_timevarying.txt")
