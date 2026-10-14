# 06d_fig_copykat.py —— CopyKAT 交叉验证图（PDF + 600 dpi PNG；Arial）
# 左：每样本恶性比例 inferCNV vs CopyKAT（同一 q95 规则），按组着色；右：癌样本两法一致率与 Cohen's kappa
# 用法: python -P 06d_fig_copykat.py <01_scRNA_inferCNV 目录>
import sys, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
base = sys.argv[1]; res = os.path.join(base, "results"); fig = os.path.join(base, "figures", "publication")
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "font.size": 7, "axes.linewidth": 0.7, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.labelsize": 8, "legend.fontsize": 6.5, "legend.frameon": False})
R = pd.read_csv(os.path.join(res, "copykat_vs_infercnv.csv"))
G = {"Cancer": "#DC0000", "Normal": "#4DBBD5", "AdjacentNormal": "#3C5488", "HSIL": "#F39B7F"}
f, axs = plt.subplots(1, 2, figsize=(5.4, 2.6), gridspec_kw={"width_ratios": [1, 1.4]})
ax = axs[0]
for g, x in R.groupby("group"):
    ax.scatter(x.pct_mal_infercnv, x.pct_aneuploid_copykat, s=16, c=G[g], edgecolor="white", lw=0.4,
               label={"AdjacentNormal": "Adjacent normal"}.get(g, g))
ax.plot([0, 100], [0, 100], color="grey", lw=0.6, ls="--"); ax.set_xlim(-3, 103); ax.set_ylim(-3, 103)
ax.set_xlabel("Malignant epithelium, inferCNV (%)"); ax.set_ylabel("Aneuploid epithelium, CopyKAT (%)")
ax.legend(loc="upper left", handletextpad=0.1)
ax = axs[1]
c = R[R.group == "Cancer"].sort_values("kappa", ascending=False, na_position="last").reset_index(drop=True)
lab = [s.split("_", 1)[1] if s.startswith("GSM") else s.replace("_patient", " P") for s in c["sample"]]
x = np.arange(len(c))
ax.bar(x - 0.2, c.agreement, 0.38, color="#8491B4", label="Agreement")
ax.bar(x + 0.2, c.kappa.fillna(0), 0.38, color="#00A087", label="Cohen's κ")
ax.set_xticks(x); ax.set_xticklabels(lab, rotation=90, fontsize=5.5); ax.axhline(0, color="black", lw=0.5)
ax.set_ylabel("inferCNV–CopyKAT concordance\n(cells called by both methods)"); ax.set_ylim(-0.3, 1.05)
ax.set_title(f"Tumors: median agreement {c.agreement.median():.2f}, κ {c.kappa.median():.2f}", fontsize=6.5)
ax.legend(loc="lower left", ncol=2)
f.tight_layout()
for e in ("pdf", "png"): f.savefig(os.path.join(fig, f"Fig_S1E_copykat_crossvalidation.{e}"), dpi=600)
print("saved Fig_S1E_copykat_crossvalidation")
