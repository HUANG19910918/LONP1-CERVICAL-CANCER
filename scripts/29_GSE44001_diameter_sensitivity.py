# 29_GSE44001_diameter_sensitivity.py —— GSE44001 多因素Cox：剔除直径记录为0的病例后的敏感性分析
# 背景: 44_GSE44001_LONP1_DFS.csv 中 14 例 largest diameter 记录为 0 cm (IB1 9例, IA2 5例),
#       按 FIGO 定义 IB1 直径应 >0, 疑为缺失被记为 0; 本脚本评估其对"每增加1 cm"HR 的影响。
# 输出: raw_data/53_GSE44001_diameter_sensitivity.txt
# 用法: python3 29_GSE44001_diameter_sensitivity.py <分析根目录>
import sys, os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lifelines import CoxPHFitter

rd = os.path.join(sys.argv[1], "raw_data")
g = pd.read_csv(os.path.join(rd, "44_GSE44001_LONP1_DFS.csv"))
L = []
def P(s): print(s, flush=True); L.append(str(s))

x = g.dropna(subset=["DFS_months","DFS_event","LONP1_log2","diameter_cm"]); x = x[x.DFS_months > 0].copy()
x["stage_adv"] = (~x.stage.isin(["IA1","IA2","IB1"])).astype(int)

P("GSE44001 multivariable Cox: sensitivity to tumour diameters recorded as 0 cm")
z = x[x.diameter_cm == 0]
P(f"diameter = 0 cm: n={len(z)} ({z.stage.value_counts().to_dict()}), events={int(z.DFS_event.sum())}")
P("")
LAB = {"LONP1_log2":"LONP1 (per log2 unit)", "stage_adv":"FIGO IB2-IIA vs IA-IB1", "diameter_cm":"Largest diameter (per cm)"}

def run(df, tag):
    cf = CoxPHFitter().fit(df[["DFS_months","DFS_event","LONP1_log2","stage_adv","diameter_cm"]],
                           "DFS_months", "DFS_event")
    P(f"=== {tag}: n={len(df)}, events={int(df.DFS_event.sum())} ===")
    for cov in ["LONP1_log2","stage_adv","diameter_cm"]:
        r = cf.summary.loc[cov]
        P(f"   {LAB[cov]:<28} HR={r['exp(coef)']:.3f} ({r['exp(coef) lower 95%']:.3f}-{r['exp(coef) upper 95%']:.3f})  p={r['p']:.4f}")
    P("")
    return cf

run(x, "All patients (as published)")
run(x[x.diameter_cm > 0].copy(), "Excluding diameter = 0 cm")

# also: treat 0 as missing and compare complete-case LONP1 estimate
P("Conclusion aid: change in the diameter HR and in the LONP1 HR between the two models above.")
open(os.path.join(rd, "53_GSE44001_diameter_sensitivity.txt"), "w", encoding="utf-8").write("\n".join(L)+"\n")
print("WROTE 53_GSE44001_diameter_sensitivity.txt")
