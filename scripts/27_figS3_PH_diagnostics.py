# 27_figS3_PH_diagnostics.py —— Figure S3: 比例风险假设诊断与时变系数模型
# 复现口径:
#   A: 单因素Cox(LONP1连续)缩放Schoenfeld残差 beta(t) vs 随访时间, 按终点分面, lowess趋势, 标注Schoenfeld检验p(rank变换)
#   B: 时变系数模型 log HR(t)=b0+b1*(log t - c) 的 HR(t) 及95%CI(delta法), 含外部队列GSE44001作对照
# 输入: raw_data/02_TCGA_CESC_clinical_LONP1_survival.csv, raw_data/44_GSE44001_LONP1_DFS.csv
# 输出: figures/【1】图片源文件/(1) 子图原文件/FigS3A_schoenfeld.png|pdf, FigS3B_HR_over_time.png|pdf
#       raw_data/51_figS3_source_data.csv
# 用法: python3 27_figS3_PH_diagnostics.py <分析根目录>
import sys, os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lifelines import CoxPHFitter, CoxTimeVaryingFitter
from lifelines.statistics import proportional_hazard_test
from lifelines.utils import to_episodic_format
from statsmodels.nonparametric.smoothers_lowess import lowess
from plotnine import *
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import setup_font, TH, NPG

base = sys.argv[1]
rd = os.path.join(base, "raw_data")
fig = os.path.join(base, "figures", "【1】图片源文件", "(1) 子图原文件")
os.makedirs(fig, exist_ok=True)
FAM = setup_font(base); T = TH(FAM)
def sv(p, name, w, h):
    p.save(os.path.join(fig, name + ".pdf"), width=w, height=h, verbose=False)
    p.save(os.path.join(fig, name + ".png"), width=w, height=h, dpi=300, verbose=False)

d2 = pd.read_csv(os.path.join(rd, "02_TCGA_CESC_clinical_LONP1_survival.csv"))
g44 = pd.read_csv(os.path.join(rd, "44_GSE44001_LONP1_DFS.csv"))

ENDPOINTS = [("OS_months","OS_event","TCGA-CESC: overall survival"),
             ("DSS_months","DSS_event","TCGA-CESC: disease-specific survival"),
             ("PFS_months","PFS_event","TCGA-CESC: progression-free survival")]

def prep(df, tk, ek, xcol):
    x = df.dropna(subset=[tk, ek, xcol]); x = x[x[tk] > 0].copy()
    out = x[[tk, ek, xcol]].rename(columns={xcol: "lon"}).reset_index(drop=True)
    out["lon"] = out["lon"] - out["lon"].mean()
    return out

# ---------------- Panel A: scaled Schoenfeld residuals ----------------
rowsA, annA = [], []
for tk, ek, lab in ENDPOINTS:
    df = prep(d2, tk, ek, "LONP1_log2RSEM")
    cf = CoxPHFitter().fit(df, tk, ek)
    beta = cf.params_["lon"]
    sch = cf.compute_residuals(df, "scaled_schoenfeld")
    t = df.loc[sch.index, tk].values
    y = sch["lon"].values + beta
    rowsA.append(pd.DataFrame({"t": t, "beta": y, "endpoint": lab}))
    lo = lowess(y, np.log(t), frac=0.8)
    annA.append(pd.DataFrame({"t": np.exp(lo[:, 0]), "beta": lo[:, 1], "endpoint": lab}))
    p_rank = proportional_hazard_test(cf, df, time_transform="rank").summary["p"].iloc[0]
    p_km = proportional_hazard_test(cf, df, time_transform="km").summary["p"].iloc[0]
    print(f"{lab}: beta={beta:.4f}, Schoenfeld p(rank)={p_rank:.4f}, p(km)={p_km:.4f}")
    annA[-1].attrs["p"] = p_rank

A = pd.concat(rowsA); Atr = pd.concat(annA)
labs_order = [e[2] for e in ENDPOINTS]
A["endpoint"] = pd.Categorical(A.endpoint, labs_order)
Atr["endpoint"] = pd.Categorical(Atr.endpoint, labs_order)
ptxt = []
for tk, ek, lab in ENDPOINTS:
    df = prep(d2, tk, ek, "LONP1_log2RSEM")
    cf = CoxPHFitter().fit(df, tk, ek)
    pr = proportional_hazard_test(cf, df, time_transform="rank").summary["p"].iloc[0]
    pk = proportional_hazard_test(cf, df, time_transform="km").summary["p"].iloc[0]
    ptxt.append({"endpoint": lab, "t": 1.3, "beta": 5.6,
                 "lab": f"Schoenfeld p = {pr:.3f} (rank), {pk:.3f} (KM)"})
PT = pd.DataFrame(ptxt); PT["endpoint"] = pd.Categorical(PT.endpoint, labs_order)

pA = (ggplot(A, aes("t", "beta"))
      + geom_hline(yintercept=0, linetype="dashed", color="#808080", size=0.4)
      + geom_point(size=1.5, alpha=0.5, color=NPG["navy"], shape="o", stroke=0)
      + geom_line(data=Atr, mapping=aes("t", "beta"), color=NPG["red"], size=1.0)
      + geom_text(data=PT, mapping=aes("t", "beta", label="lab"), ha="left", size=7.5, family=FAM)
      + scale_x_log10(breaks=[1, 3, 6, 12, 24, 60, 120], labels=["1","3","6","12","24","60","120"])
      + facet_wrap("~endpoint", nrow=1)
      + labs(x="Time (months, log scale)", y="Scaled Schoenfeld residual\n(time-specific log HR for LONP1)")
      + T + theme(figure_size=(10.5, 3.4), strip_text=element_text(size=9, family=FAM),
                  panel_spacing=0.05))
sv(pA, "FigS3A_schoenfeld", 10.5, 3.4)

# ---------------- Panel B: HR(t) from time-varying coefficient model ----------------
def hr_curve(df, tk, ek, lab, tmin, tmax):
    dd = df.copy(); dd["id"] = np.arange(len(dd))
    long = to_episodic_format(dd, duration_col=tk, event_col=ek, id_col="id")
    long = long[long["stop"] > 0].copy()
    C = float(np.log(df.loc[df[ek] == 1, tk]).mean())
    long["gt"] = np.log(long["stop"]) - C
    long["lon_x_gt"] = long["lon"] * long["gt"]
    ctv = CoxTimeVaryingFitter().fit(long[["id","start","stop",ek,"lon","lon_x_gt"]],
            id_col="id", event_col=ek, start_col="start", stop_col="stop", show_progress=False)
    s = ctv.summary; V = ctv.variance_matrix_
    b0, b1 = s.loc["lon","coef"], s.loc["lon_x_gt","coef"]
    ts = np.exp(np.linspace(np.log(tmin), np.log(tmax), 120))
    g = np.log(ts) - C
    est = b0 + b1*g
    var = V.iloc[0,0] + (g**2)*V.iloc[1,1] + 2*g*V.iloc[0,1]
    se = np.sqrt(np.maximum(var, 0))
    return pd.DataFrame({"t": ts, "HR": np.exp(est), "lo": np.exp(est-1.96*se),
                         "hi": np.exp(est+1.96*se), "endpoint": lab}), s.loc["lon_x_gt","p"]

rowsB, pB_txt = [], []
for tk, ek, lab in ENDPOINTS:
    df = prep(d2, tk, ek, "LONP1_log2RSEM")
    cur, p1 = hr_curve(df, tk, ek, lab, 6, 60)
    rowsB.append(cur); pB_txt.append({"endpoint": lab, "lab": f"time-dependence p = {p1:.3f}"})
dfg = prep(g44, "DFS_months", "DFS_event", "LONP1_log2")
labg = "GSE44001: disease-free survival"
cur, p1 = hr_curve(dfg, "DFS_months", "DFS_event", labg, 6, 60)
rowsB.append(cur); pB_txt.append({"endpoint": labg, "lab": f"time-dependence p = {p1:.3f}"})

order_b = labs_order + [labg]
B = pd.concat(rowsB); B["endpoint"] = pd.Categorical(B.endpoint, order_b)
PB = pd.DataFrame(pB_txt); PB["t"] = 6.3; PB["HR"] = 3.6
PB["endpoint"] = pd.Categorical(PB.endpoint, order_b)

pB = (ggplot(B, aes("t", "HR"))
      + geom_hline(yintercept=1, linetype="dashed", color="#808080", size=0.4)
      + geom_ribbon(aes(ymin="lo", ymax="hi"), fill=NPG["blue"], alpha=0.25)
      + geom_line(color=NPG["navy"], size=0.9)
      + geom_text(data=PB, mapping=aes("t", "HR", label="lab"), ha="left", size=7.5, family=FAM)
      + scale_x_log10(breaks=[6, 12, 24, 36, 60], labels=["6","12","24","36","60"])
      + scale_y_log10(breaks=[0.1, 0.25, 0.5, 1, 2, 4], labels=["0.1","0.25","0.5","1","2","4"])
      + facet_wrap("~endpoint", nrow=1)
      + labs(x="Time (months, log scale)", y="Hazard ratio for LONP1\n(per log2 unit, 95% CI)")
      + T + theme(figure_size=(12.5, 3.4), strip_text=element_text(size=9, family=FAM),
                  panel_spacing=0.05))
sv(pB, "FigS3B_HR_over_time", 12.5, 3.4)

B.round(4).to_csv(os.path.join(rd, "51_figS3_source_data.csv"), index=False)
print("Figure S3 panels written to:", fig)
