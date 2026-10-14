# 04e_fig_SIL_grade.py —— GSE333892 SIL 分级敏感性图（PDF + 600 dpi PNG；Arial）
# 用法: python -P 04e_fig_SIL_grade.py <01_scRNA_inferCNV 目录>
import sys, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
base = sys.argv[1]; res = os.path.join(base, "results"); fig = os.path.join(base, "figures", "publication")
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "font.size": 7, "axes.linewidth": 0.7, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.labelsize": 8, "legend.fontsize": 6, "legend.frameon": False})
U = pd.read_csv(os.path.join(res, "GSE333892_pseudobulk_LONP1.csv")); U = U[U.included]
S = pd.read_csv(os.path.join(res, "GSE333892_SIL_grade_sensitivity.csv"))
groups = [("Normal", U[U.category == "normal_epi"].log2CPM_LONP1, "#4DBBD5"),
          ("SIL\nCNV-low,\nnot HSIL-labelled", S[(S.label_grade != "HSIL_labelled") & (S.cnv_class == "CNV_low")].log2CPM_LONP1, "#F39B7F"),
          ("SIL\nHSIL-labelled\nor CNV-high", S[(S.label_grade == "HSIL_labelled") | (S.cnv_class == "CNV_high")].log2CPM_LONP1, "#E64B35"),
          ("Cancer\n(malignant)", U[U.category == "malignant_epi"].log2CPM_LONP1, "#DC0000")]
rng = np.random.default_rng(3)
f, ax = plt.subplots(figsize=(3.2, 2.6))
for i, (lab, v, c) in enumerate(groups):
    ax.boxplot(v, positions=[i], widths=0.5, showfliers=False, medianprops=dict(color="black", lw=1), boxprops=dict(lw=0.6),
               whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    ax.scatter(i + rng.uniform(-0.13, 0.13, len(v)), v, s=9, c=c, edgecolor="white", lw=0.3, zorder=3)
    ax.text(i, 2.25, f"n={len(v)}", ha="center", fontsize=6)
ax.set_xticks(range(4)); ax.set_xticklabels([g[0] for g in groups], fontsize=6); ax.set_ylim(2.1, 5.6)
ax.set_ylabel("LONP1 pseudobulk log2(CPM+1)"); ax.set_title("GSE333892: SIL by inferred grade / CNV status", fontsize=7)
f.tight_layout()
for e in ("pdf", "png"): f.savefig(os.path.join(fig, f"Fig_S1G_SIL_grade_sensitivity.{e}"), dpi=600)
print("saved")
