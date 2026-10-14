# 04c_fig_progression.py —— GSE333892 正常→SIL→癌上皮 LONP1 + 六数据集合并恶性 vs 正常上皮（PDF + 600 dpi PNG；Arial）
# 用法: python -P 04c_fig_progression.py <01_scRNA_inferCNV 目录>
import sys, os, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
base = sys.argv[1]; res = os.path.join(base, "results"); fig = os.path.join(base, "figures", "publication")
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "font.size": 7, "axes.linewidth": 0.7, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.labelsize": 8, "legend.fontsize": 6, "legend.frameon": False})
NPG = {"GSE208653": "#E64B35", "E-MTAB-11948": "#3C5488", "GSE197461": "#00A087", "GSE236738": "#F39B7F", "GSE333892": "#8491B4"}
U = pd.read_csv(os.path.join(res, "GSE333892_pseudobulk_LONP1.csv")); U = U[U.included]
A = pd.read_csv(os.path.join(res, "pseudobulk_LONP1_pooled_6datasets.csv"))
L = open(os.path.join(res, "04b_GSE333892_progression_log.txt"), encoding="utf-8").read()
rng = np.random.default_rng(2)
f, axs = plt.subplots(1, 2, figsize=(5.0, 2.6), gridspec_kw={"width_ratios": [1.1, 1]})
ax = axs[0]; cats = [("normal_epi", "Normal", "#4DBBD5"), ("SIL_epi", "SIL", "#F39B7F"), ("malignant_epi", "Cancer\n(malignant)", "#DC0000")]
for i, (c, lab, col) in enumerate(cats):
    v = U[U.category == c].log2CPM_LONP1
    ax.boxplot(v, positions=[i], widths=0.5, showfliers=False, medianprops=dict(color="black", lw=1), boxprops=dict(lw=0.6),
               whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    ax.scatter(i + rng.uniform(-0.14, 0.14, len(v)), v, s=9, c=col, edgecolor="white", lw=0.3, zorder=3)
    ax.text(i, ax.get_ylim()[0] if False else 2.2, f"n={len(v)}", ha="center", fontsize=6)
ax.set_xticks(range(3)); ax.set_xticklabels([l for _, l, _ in cats]); ax.set_ylabel("LONP1 pseudobulk log2(CPM+1)")
ax.set_ylim(2.1, 5.9)
kw = re.search(r"Kruskal–Wallis p=([\d.e-]+)", L).group(1)
p1 = re.search(r"normal_epi vs SIL_epi: .*MWU p=([\d.e-]+)", L).group(1); p2 = re.search(r"SIL_epi vs malignant_epi: .*MWU p=([\d.e-]+)", L).group(1)
p3 = re.search(r"normal_epi vs malignant_epi: .*MWU p=([\d.e-]+)", L).group(1)
for (x0, x1, y, p) in [(0, 1, 5.25, p1), (1, 2, 5.45, p2), (0, 2, 5.7, p3)]:
    ax.plot([x0, x0, x1, x1], [y - 0.05, y, y, y - 0.05], color="black", lw=0.6)
    ax.text((x0 + x1) / 2, y + 0.02, f"P={float(p):.1g}" if float(p) < 0.01 else f"P={float(p):.2f}", ha="center", va="bottom", fontsize=6)
ax.set_title(f"GSE333892 epithelium (Kruskal–Wallis P={float(kw):.1g})", fontsize=7)
ax = axs[1]
dss = [d for d in ["GSE208653", "E-MTAB-11948", "GSE197461", "GSE236738", "GSE333892"] if d in set(A.dataset)]
for j, (flag, lab) in enumerate([(0, "Normal\nepithelium"), (1, "Malignant\nepithelium")]):
    g = A[A.is_mal == flag]
    ax.scatter(j + rng.uniform(-0.16, 0.16, len(g)), g.log2CPM_LONP1, s=10, c=[NPG[d] for d in g.dataset], edgecolor="white", lw=0.3, zorder=3)
    ax.hlines(g.log2CPM_LONP1.median(), j - 0.28, j + 0.28, color="black", lw=1.2, zorder=4)
ax.set_xticks([0, 1]); ax.set_xticklabels(["Normal\nepithelium", "Malignant\nepithelium"]); ax.set_xlim(-0.6, 1.6)
m = re.search(r"\[合并 C1，6 数据集\].*diff=([-\d.]+) \(SE [\d.]+, 95% CI ([-\d.]+)–([-\d.]+)\) p=([\d.e-]+)", L)
ax.set_title(f"Pooled, 5 cohorts (OLS + cohort)\nΔ={float(m.group(1)):.2f} (95% CI {float(m.group(2)):.2f}–{float(m.group(3)):.2f}), P={float(m.group(4)):.2f}", fontsize=6.5)
for d in dss: ax.scatter([], [], s=10, c=NPG[d], label=d)
ax.legend(loc="lower right", handletextpad=0.1, fontsize=5.5); ax.set_ylabel("LONP1 pseudobulk log2(CPM+1)")
f.tight_layout()
for e in ("pdf", "png"): f.savefig(os.path.join(fig, f"Fig_S1F_LONP1_progression_pooled.{e}"), dpi=600)
print("saved")
