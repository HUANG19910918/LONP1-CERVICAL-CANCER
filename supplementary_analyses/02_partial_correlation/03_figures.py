# 03_figures.py —— 分析 2 出版级图（PDF 矢量 + 600 dpi PNG；Arial；配色沿用 v8 figstyle.py 的 NPG 调色板）
# (a) 19p13.3 富集 OR + 95% CI（Fisher 条件 MLE，scipy odds_ratio）：原始 / 偏相关 A / 偏相关 B × 两种截断（rho>=0.4、top164）
# (b) 原始 rho vs 偏相关 A rho 散点，19p13.3 基因单独着色，标注 v8 线粒体模块基因
# (c) 三层 fgsea 线粒体/蛋白稳态相关通路 NES 点阵图（点大小 = -log10 padj）
# 用法: python -P 03_figures.py <02_partial_correlation 目录>
import sys, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats.contingency import odds_ratio
from scipy.stats import fisher_exact

base = sys.argv[1]; res = os.path.join(base, "results"); fig = os.path.join(base, "figures"); os.makedirs(fig, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 8, "axes.linewidth": 0.8,
                     "axes.spines.top": False, "axes.spines.right": False, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
                     "axes.labelsize": 9, "axes.titlesize": 9, "legend.fontsize": 7, "legend.frameon": False})
NPG = {"red": "#E64B35", "blue": "#4DBBD5", "green": "#00A087", "navy": "#3C5488", "salmon": "#F39B7F", "purple": "#8491B4",
       "grey": "#B0B0B0"}
LAYERS = [("rho_raw", "Raw Spearman"), ("rho_partial_A", "Partial A\n(LONP1 CN + gene CN)"),
          ("rho_partial_B", "Partial B\n(A + 19p arm CN + purity)")]
def save(f, name):
    f.savefig(os.path.join(fig, name + ".pdf")); f.savefig(os.path.join(fig, name + ".png"), dpi=600); plt.close(f)
    print("saved", name)

t = pd.read_csv(os.path.join(res, "LONP1_partial_spearman_genomewide.csv"))
ann = t[t.cytoband.notna()].copy()
rows = []
for col, lab in LAYERS:
    v = ann[col]
    for cut, sel in [("rho ≥ 0.4", v >= 0.4), ("Top 164 genes", v.rank(ascending=False, method="first") <= 164)]:
        a = int((sel & ann.on_19p13_3).sum()); b = int((sel & ~ann.on_19p13_3).sum())
        c = int((~sel & ann.on_19p13_3).sum()); d = int((~sel & ~ann.on_19p13_3).sum())
        r = odds_ratio([[a, b], [c, d]]); ci = r.confidence_interval(0.95)
        p = fisher_exact([[a, b], [c, d]]).pvalue
        rows.append({"layer": col, "label": lab, "cut": cut, "a": a, "n_sel": a + b, "OR_cond": r.statistic,
                     "OR_sample": (a * d) / (b * c), "lo": ci.low, "hi": ci.high, "p": p})
E = pd.DataFrame(rows); E.to_csv(os.path.join(res, "fig_a_19p13.3_OR_CI.csv"), index=False)
print(E.round(3).to_string(index=False))

# (a)
f, ax = plt.subplots(figsize=(3.6, 2.6))
cols = {"rho ≥ 0.4": NPG["red"], "Top 164 genes": NPG["navy"]}
for k, (cut, g) in enumerate(E.groupby("cut", sort=False)):
    y = np.arange(len(LAYERS)) + (k - 0.5) * 0.28
    ax.errorbar(g.OR_sample, y, xerr=[g.OR_sample - g.lo, g.hi - g.OR_sample], fmt="o", ms=4, color=cols[cut],
                ecolor=cols[cut], elinewidth=0.9, capsize=2, label=cut)
    for yi, (o, n, a) in zip(y, g[["OR_sample", "n_sel", "a"]].values):
        ax.text(g.hi.max() * 1.25, yi, f"{o:.1f} ({int(a)}/{int(n)})", va="center", fontsize=6.5, color=cols[cut])
ax.set_xscale("log"); ax.axvline(1, color="grey", lw=0.6, ls="--")
ax.set_yticks(range(len(LAYERS))); ax.set_yticklabels([l for _, l in LAYERS]); ax.invert_yaxis()
ax.set_xlabel("19p13.3 enrichment odds ratio (95% CI)"); ax.set_xlim(0.8, E.hi.max() * 6)
ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, handletextpad=0.3)
f.tight_layout(); save(f, "Fig_P2a_19p13.3_enrichment_OR")

# (b)
lab = ["CLPP", "TIMM13", "TRAP1", "PHB", "PHB2", "HSPD1", "HSPE1", "HSPA9", "ATF4", "ATF5", "SHMT2", "THOP1"]
f, ax = plt.subplots(figsize=(3.2, 3.0))
o = t[~t.on_19p13_3.astype(bool)]; i = t[t.on_19p13_3.astype(bool)]
ax.scatter(o.rho_raw, o.rho_partial_A, s=1.5, c=NPG["grey"], alpha=0.35, lw=0, rasterized=True, label=f"Other genes (n={len(o):,})")
ax.scatter(i.rho_raw, i.rho_partial_A, s=5, c=NPG["red"], alpha=0.85, lw=0, label=f"19p13.3 genes (n={len(i)})")
lim = [-0.6, 0.8]; ax.plot(lim, lim, color="black", lw=0.6, ls="--"); ax.set_xlim(lim); ax.set_ylim(lim)
from scipy.stats import pearsonr
OFF = {"THOP1": (-34, 8), "CLPP": (6, -2), "TIMM13": (6, -8), "TRAP1": (6, -6), "PHB": (-30, 6), "PHB2": (-32, 4),
       "SHMT2": (-34, 4), "HSPD1": (8, -2), "ATF4": (-28, -2), "HSPE1": (6, -9), "HSPA9": (6, -8), "ATF5": (-26, 4)}
for g in lab:
    x = t[t.gene == g]
    if len(x):
        xr, yr = x.rho_raw.iloc[0], x.rho_partial_A.iloc[0]
        ax.scatter(xr, yr, s=12, facecolor="none", edgecolor=NPG["navy"], lw=0.8)
        dx, dy = OFF.get(g, (4, 3))
        ax.annotate("PHB1" if g == "PHB" else g, (xr, yr), xytext=(dx, dy), textcoords="offset points", fontsize=6,
                    color=NPG["navy"], arrowprops=dict(arrowstyle="-", lw=0.4, color=NPG["navy"], shrinkA=0, shrinkB=2))
ax.set_xlabel("LONP1 co-expression, raw Spearman ρ"); ax.set_ylabel("Partial Spearman ρ\n(adjusted for LONP1 and gene copy number)")
ax.legend(loc="lower right", markerscale=3, handletextpad=0.2)
f.tight_layout(); save(f, "Fig_P2b_raw_vs_partial_rho")

# (c)
fg = pd.read_csv(os.path.join(res, "fgsea_GOBP_Reactome_by_layer.csv"))
P = [("REACTOME_TRANSLATION", "Translation (R)"), ("GOBP_MITOCHONDRIAL_TRANSLATION", "Mitochondrial translation (GO)"),
     ("GOBP_MITOCHONDRIAL_GENE_EXPRESSION", "Mitochondrial gene expression (GO)"),
     ("GOBP_OXIDATIVE_PHOSPHORYLATION", "Oxidative phosphorylation (GO)"),
     ("REACTOME_RESPIRATORY_ELECTRON_TRANSPORT", "Respiratory electron transport (R)"),
     ("GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY", "Respiratory chain complex assembly (GO)"),
     ("REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT", "Mitochondrial protein import (R)"),
     ("REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION", "Mitochondrial protein degradation (R)"),
     ("GOBP_PROTEIN_FOLDING", "Protein folding (GO)"),
     ("REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY", "GCN2 response to amino acid deficiency (R)")]
sub = fg[fg.pathway.isin([p for p, _ in P])]
sub.to_csv(os.path.join(res, "fig_c_fgsea_selected.csv"), index=False)
f = plt.figure(figsize=(5.2, 3.4)); ax = f.add_axes([0.47, 0.10, 0.22, 0.86]); cax = f.add_axes([0.73, 0.55, 0.02, 0.3])
norm = plt.Normalize(2.0, 4.0); cmap = plt.cm.Reds
for yi, (p, nm) in enumerate(P):
    for xi, (col, _) in enumerate(LAYERS):
        r = sub[(sub.pathway == p) & (sub.stat == col)]
        if not len(r): continue
        nes, padj = r.NES.iloc[0], r.padj.iloc[0]
        s = min(-np.log10(max(padj, 1e-100)), 60)
        ax.scatter(xi, yi, s=8 + s * 2.2, c=[cmap(norm(nes))], edgecolor="black", lw=0.4)
ax.set_xticks(range(3)); ax.set_xticklabels(["Raw", "Partial\nA", "Partial\nB"]); ax.set_xlim(-0.6, 2.6)
ax.set_yticks(range(len(P))); ax.set_yticklabels([n for _, n in P]); ax.invert_yaxis()
ax.tick_params(axis="y", length=0)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); cb = f.colorbar(sm, cax=cax); cb.set_label("NES", fontsize=7)
cb.ax.tick_params(labelsize=6)
for s in (10, 30, 60):
    ax.scatter([], [], s=8 + s * 2.2, c="white", edgecolor="black", lw=0.4, label=f"{s}" if s < 60 else "≥60")
ax.legend(title="−log10 padj", loc="upper left", bbox_to_anchor=(1.28, 0.42), labelspacing=1.1, title_fontsize=7)
save(f, "Fig_P2c_fgsea_mito_pathways")
