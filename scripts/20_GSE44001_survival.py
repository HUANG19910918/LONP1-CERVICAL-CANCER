# 20_GSE44001_survival.py —— 外部队列 GSE44001（早期宫颈癌 n=300，Illumina DASL GPL14951，无病生存DFS）验证 LONP1 预后
# 复现口径：
#   数据: GSE44001_series_matrix.txt.gz（GEO FTP；作者已 log2 + quantile 归一化）；LONP1 探针 ILMN_1766125（GPL14951 上唯一 LONP1 探针）
#   临床: 样本 characteristics 中 Stage / largest diameter / DFS(months) / status_of_dfs
#   统计: 单因素 Cox（连续 log2 表达，lifelines）；中位数分组 KM + log-rank；多因素 Cox（LONP1 + 分期 IB1参照 + 肿瘤最大径）
# 输出: raw_data/44_GSE44001_LONP1_DFS.csv, 45_GSE44001_survival_stats.txt；figures/FigS2_GSE44001_KM.png/pdf
# 用法: python3 20_GSE44001_survival.py <分析根目录> <series_matrix.gz>
import sys, os, gzip, re
import numpy as np, pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from plotnine import *
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import setup_font, TH, NPG
base, sm = sys.argv[1], sys.argv[2]
rd = os.path.join(base, "raw_data"); fig = os.path.join(base, "figures", "子图原文件"); os.makedirs(fig, exist_ok=True)
FAM = setup_font(base); T = TH(FAM)
PROBE = "ILMN_1766125"
ids = None; chars = []; expr = None
with gzip.open(sm, "rt") as f:
    for line in f:
        if line.startswith('"ID_REF"'): ids = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"): chars.append([x.strip('"') for x in line.rstrip("\n").split("\t")[1:]])
        elif line.startswith(f'"{PROBE}"'): expr = [float(x) for x in line.rstrip("\n").split("\t")[1:]]
d = pd.DataFrame({"gsm": ids, "LONP1_log2": expr})
for row in chars:
    key = row[0].split(":")[0].strip()
    d[key] = [x.split(":", 1)[1].strip() if ":" in x else np.nan for x in row]
d = d.rename(columns={"Stage": "stage", "largest diameter": "diameter_cm", "disease_free_survival_(dfs)_(months)": "DFS_months", "status_of_dfs": "DFS_event"})
d["DFS_months"] = pd.to_numeric(d.DFS_months, errors="coerce"); d["DFS_event"] = pd.to_numeric(d.DFS_event, errors="coerce"); d["diameter_cm"] = pd.to_numeric(d.diameter_cm, errors="coerce")
d.to_csv(os.path.join(rd, "44_GSE44001_LONP1_DFS.csv"), index=False)
x = d.dropna(subset=["DFS_months", "DFS_event", "LONP1_log2"]); x = x[x.DFS_months > 0].copy()
L = []
def P(s): print(s); L.append(s)
P(f"GSE44001: n={len(x)} with DFS, events={int(x.DFS_event.sum())}; LONP1 log2 median={x.LONP1_log2.median():.2f} range={x.LONP1_log2.min():.2f}-{x.LONP1_log2.max():.2f}; stages: {x.stage.value_counts().to_dict()}")
cf = CoxPHFitter().fit(x[["DFS_months", "DFS_event", "LONP1_log2"]], "DFS_months", "DFS_event"); s = cf.summary.iloc[0]
P(f"Univariate Cox DFS (per log2 unit): HR={s['exp(coef)']:.2f} (95% CI {s['exp(coef) lower 95%']:.2f}-{s['exp(coef) upper 95%']:.2f}), p={s['p']:.3f}")
# per SD
sd = x.LONP1_log2.std(); P(f"  per-SD HR={np.exp(cf.params_.iloc[0]*sd):.2f} (SD={sd:.2f})")
med = x.LONP1_log2.median(); x["grp"] = np.where(x.LONP1_log2 > med, "High", "Low")
lr = logrank_test(x.loc[x.grp == "High", "DFS_months"], x.loc[x.grp == "Low", "DFS_months"], x.loc[x.grp == "High", "DFS_event"], x.loc[x.grp == "Low", "DFS_event"])
P(f"Median split (cutoff {med:.2f}; High n={(x.grp=='High').sum()}, Low n={(x.grp=='Low').sum()}): log-rank p={lr.p_value:.3f}")
# multivariable: stage groups (IA/IB1 vs IB2/II+), diameter
x["stage_adv"] = (~x.stage.isin(["IA1", "IA2", "IB1"])).astype(int)
xm = x.dropna(subset=["diameter_cm"])
cf2 = CoxPHFitter().fit(xm[["DFS_months", "DFS_event", "LONP1_log2", "stage_adv", "diameter_cm"]], "DFS_months", "DFS_event")
for cov, r in cf2.summary.iterrows():
    P(f"Multivariable (n={len(xm)}, events={int(xm.DFS_event.sum())}) {cov}: HR={r['exp(coef)']:.2f} ({r['exp(coef) lower 95%']:.2f}-{r['exp(coef) upper 95%']:.2f}) p={r['p']:.3f}")
open(os.path.join(rd, "45_GSE44001_survival_stats.txt"), "w").write("\n".join(L) + "\n")
# KM figure
rows, cens = [], []
for g in ["Low", "High"]:
    kmf = KaplanMeierFitter().fit(x.loc[x.grp == g, "DFS_months"], x.loc[x.grp == g, "DFS_event"])
    sf = kmf.survival_function_.reset_index(); sf.columns = ["t", "s"]; sf["g"] = g; rows.append(sf)
    cc = x.loc[(x.grp == g) & (x.DFS_event == 0)]; cens.append(pd.DataFrame({"t": cc.DFS_months.values, "s": kmf.survival_function_at_times(cc.DFS_months.values).values, "g": g}))
dd = pd.concat(rows); ce = pd.concat(cens)
dd["g"] = pd.Categorical(dd.g, ["Low", "High"]); ce["g"] = pd.Categorical(ce.g, ["Low", "High"])
ann = f"Log-rank p = {lr.p_value:.3f}\nCox (continuous) HR = {s['exp(coef)']:.2f}\n(95% CI {s['exp(coef) lower 95%']:.2f}-{s['exp(coef) upper 95%']:.2f}), p = {s['p']:.3f}"
p = (ggplot(dd, aes("t", "s", color="g")) + geom_step(size=0.9) + geom_point(data=ce, shape="+", size=2.0, show_legend=False)
     + scale_color_manual(values={"Low": NPG["blue"], "High": NPG["red"]}, labels=[f"Low (n={(x.grp=='Low').sum()})", f"High (n={(x.grp=='High').sum()})"])
     + scale_y_continuous(limits=(0, 1), labels=lambda l: [f"{v*100:.0f}%" for v in l])
     + labs(x="Time (months)", y="Disease-free survival probability", color="", title=f"GSE44001 (early-stage cervical cancer, n={len(x)})\nDFS by LONP1 (median split)")
     + annotate("text", x=dd.t.max() * 0.55, y=0.16, label=ann, size=8, ha="left", family=FAM)
     + T + theme(legend_position=(0.22, 0.22), legend_background=element_blank(), plot_title=element_text(size=10.5, weight="bold", ha="center", family=FAM)))
p.save(os.path.join(fig, "FigS2_GSE44001_KM.png"), width=4.3, height=4.2, dpi=300, verbose=False)
p.save(os.path.join(fig, "FigS2_GSE44001_KM.pdf"), width=4.3, height=4.2, verbose=False)
print("DONE20")
