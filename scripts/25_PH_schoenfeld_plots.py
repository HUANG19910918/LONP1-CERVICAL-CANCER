# Scaled Schoenfeld residual plots for the univariate Cox models (LONP1, TCGA-CESC).
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from statsmodels.nonparametric.smoothers_lowess import lowess

rd, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)
d2 = pd.read_csv(os.path.join(rd, "02_TCGA_CESC_clinical_LONP1_survival.csv"))
for tk, ek, lab in [("OS_months","OS_event","OS"), ("PFS_months","PFS_event","PFS"), ("DSS_months","DSS_event","DSS")]:
    x = d2.dropna(subset=[tk, ek]); x = x[x[tk] > 0].copy()
    df = x[[tk, ek, "LONP1_log2RSEM"]].reset_index(drop=True)
    cf = CoxPHFitter().fit(df, tk, ek)
    sch = cf.compute_residuals(df, "scaled_schoenfeld")
    beta = cf.params_["LONP1_log2RSEM"]
    t = df.loc[sch.index, tk].values
    y = sch["LONP1_log2RSEM"].values + beta
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, tr in zip(axes, ["rank", "km"]):
        xt = pd.Series(t).rank().values if tr == "rank" else t
        ax.scatter(xt, y, s=18, alpha=0.6, color="#3C5488")
        lo = lowess(y, xt, frac=0.8); ax.plot(lo[:, 0], lo[:, 1], color="#E64B35", lw=2)
        ax.axhline(beta, ls="--", color="grey", lw=1); ax.axhline(0, ls=":", color="black", lw=0.8)
        pv = proportional_hazard_test(cf, df, time_transform=tr).summary["p"].iloc[0]
        ax.set_title(f"{lab}: scaled Schoenfeld, {tr} transform (p={pv:.4f})", fontsize=10)
        ax.set_xlabel("rank(time)" if tr == "rank" else "time (months)")
        ax.set_ylabel("beta(t) for LONP1 (log HR)")
    fig.tight_layout(); fig.savefig(os.path.join(outdir, f"PH_schoenfeld_{lab}.png"), dpi=130); plt.close(fig)
    print("saved", lab)
