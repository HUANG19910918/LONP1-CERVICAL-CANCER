import os, sys
import numpy as np, pandas as pd
rd = sys.argv[1]
d2 = pd.read_csv(os.path.join(rd,"02_TCGA_CESC_clinical_LONP1_survival.csv"))
L=[]
def P(s): print(s,flush=True); L.append(str(s))

P("STEP 1: data integrity / coding / sample consistency")
P(f"rows={len(d2)}  cols={list(d2.columns)}")
P("")
for ek in ["OS_event","DSS_event","PFS_event"]:
    P(f"{ek}: unique={sorted(d2[ek].dropna().unique().tolist())}  n_missing={d2[ek].isna().sum()}  n_events={int(d2[ek].fillna(0).sum())}")
P("")
for tk in ["OS_months","DSS_months","PFS_months"]:
    s = d2[tk]
    P(f"{tk}: n_missing={s.isna().sum()}  min={s.min():.3f}  median={s.median():.2f}  max={s.max():.2f}  n_zero={(s==0).sum()}  n_neg={(s<0).sum()}")
P("")
# analysis sets
def aset(tk,ek):
    x = d2.dropna(subset=[tk,ek]); x = x[x[tk]>0]
    return set(x.patient_id)
sOS, sDSS, sPFS = aset("OS_months","OS_event"), aset("DSS_months","DSS_event"), aset("PFS_months","PFS_event")
x = d2.dropna(subset=["OS_months","OS_event"]); x = x[x.OS_months>0].copy()
x["Tg"] = x.path_T.astype(str).str.extract(r"^(T[1-4])")[0]
x = x.dropna(subset=["Tg","age"]); x = x[x.path_N.isin(["N0","N1"])]
sMV = set(x.patient_id)
P(f"analysis sets: OS={len(sOS)}  DSS={len(sDSS)}  PFS={len(sPFS)}  multivariableOS={len(sMV)}")
P(f"MV subset of OS? {sMV.issubset(sOS)};  OS\\MV dropped={len(sOS-sMV)}")
P(f"OS vs PFS identical? {sOS==sPFS};  OS\\PFS={len(sOS-sPFS)}  PFS\\OS={len(sPFS-sOS)}")
P(f"OS vs DSS: OS\\DSS={len(sOS-sDSS)}  DSS\\OS={len(sDSS-sOS)}")
P("")
# why patients dropped from MV model
d = d2.dropna(subset=["OS_months","OS_event"]); d = d[d.OS_months>0].copy()
d["Tg"] = d.path_T.astype(str).str.extract(r"^(T[1-4])")[0]
P("reasons for exclusion from multivariable model (within the OS set, n=%d):" % len(d))
P(f"  missing/unparsable T stage: {d.Tg.isna().sum()}   (raw path_T values: {sorted(set(d.path_T.astype(str)))[:12]} ...)")
P(f"  missing age:                {d.age.isna().sum()}")
P(f"  path_N not in N0/N1:        {(~d.path_N.isin(['N0','N1'])).sum()}   (values: {d.path_N.astype(str).value_counts().to_dict()})")
P("")
# risk set over time
P("risk set over follow-up (OS and PFS):")
P(f"{'month':>7} {'at_risk_OS':>11} {'cumEv_OS':>9} {'at_risk_PFS':>12} {'cumEv_PFS':>10}")
xo = d2.dropna(subset=["OS_months","OS_event"]); xo = xo[xo.OS_months>0]
xp = d2.dropna(subset=["PFS_months","PFS_event"]); xp = xp[xp.PFS_months>0]
for t in [0,6,12,24,36,48,60,84,120,180,240]:
    nro=(xo.OS_months>=t).sum(); ceo=int(xo.loc[xo.OS_months<=t,"OS_event"].sum())
    nrp=(xp.PFS_months>=t).sum(); cep=int(xp.loc[xp.PFS_months<=t,"PFS_event"].sum())
    P(f"{t:>7} {nro:>11} {ceo:>9} {nrp:>12} {cep:>10}")
P("")
P(f"OS: median follow-up (all) = {xo.OS_months.median():.1f} mo; median event time = {xo.loc[xo.OS_event==1,'OS_months'].median():.1f} mo")
P(f"PFS: median follow-up (all) = {xp.PFS_months.median():.1f} mo; median event time = {xp.loc[xp.PFS_event==1,'PFS_months'].median():.1f} mo")
open(os.path.join(rd,"48_PH_step1_data_checks.txt"),"w",encoding="utf-8").write("\n".join(L)+"\n")
print("WROTE 48_PH_step1_data_checks.txt")
