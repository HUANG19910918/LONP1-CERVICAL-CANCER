# 06_figures_allgenes.py —— 16 个基因 cis-eQTL MR 汇总图与 coloc 汇总热图（PDF + 600 dpi PNG；Arial）
# D：森林图——eQTLGen 血液 → FinnGen 宫颈癌（发现）与 BBJ 宫颈癌（重复），每基因一行；标注工具变量数与 FDR
# E：coloc 热图——每基因 × 结局取 4 个 GTEx 组织中最大的 PP.H4
# 用法: python -P 06_figures_allgenes.py <03_MR_coloc 目录>
import sys, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
base = sys.argv[1]; res = os.path.join(base, "results"); fig = os.path.join(base, "figures")
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "font.size": 7, "axes.linewidth": 0.7, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.labelsize": 8, "legend.fontsize": 6, "legend.frameon": False})
A = pd.read_csv(os.path.join(res, "mr_all_genes_combined.csv"))
UPR = {"TRAP1", "PHB2", "HSPD1", "HSPE1", "ATF4"}
def save(f, n):
    for e in ("pdf", "png"): f.savefig(os.path.join(fig, f"{n}.{e}"), dpi=600)
    plt.close(f); print("saved", n)
# D
e = A[A.exposure == "eQTLGen_blood"]
genes = ["LONP1"] + sorted([g for g in e.gene.unique() if g != "LONP1"], key=lambda g: e[(e.gene == g) & (e.outcome == "FinnGen_CC")].p.min())
missing = [g for g in pd.read_csv(os.path.join(res, "genes_coordinates.csv")).gene if g not in set(e.gene)]
f, axs = plt.subplots(1, 2, figsize=(5.6, 0.22 * len(genes) + 0.9), sharey=True)
for ax, (o, ttl) in zip(axs, [("FinnGen_CC", "FinnGen cervical cancer (878 cases)"), ("BBJ_CeC", "BBJ cervical cancer (967 cases)")]):
    for i, g in enumerate(genes):
        r = e[(e.gene == g) & (e.outcome == o)]
        if not len(r): continue
        r = r.iloc[0]; col = "#DC0000" if r.fdr < 0.05 else ("#3C5488" if g == "LONP1" else "#7E7E7E")
        ax.errorbar(r.OR, -i, xerr=[[r.OR - r.OR_lo], [min(r.OR_hi, 20) - r.OR]], fmt="s" if g == "LONP1" else "o", ms=3, color=col,
                    elinewidth=0.8, capsize=1.5)
        ax.text(21, -i, f"{int(r.n_snp)} | {r.OR:.2f} | FDR {r.fdr:.2f}", va="center", fontsize=5)
    ax.axvline(1, color="grey", lw=0.6, ls="--"); ax.set_xscale("log"); ax.set_xlim(0.1, 20)
    ax.set_xticks([0.1, 0.5, 1, 2, 5, 20]); ax.set_xticklabels(["0.1", "0.5", "1", "2", "5", "20"]); ax.set_title(ttl, fontsize=6.5)
    ax.set_xlabel("OR per 1-SD expression (eQTLGen blood)")
axs[0].set_yticks([-i for i in range(len(genes))])
axs[0].set_yticklabels([("* " if g in UPR else "") + g for g in genes]); axs[0].tick_params(axis="y", length=0)
axs[0].get_yticklabels()[0].set_fontweight("bold")
axs[1].text(21, 1.0, "SNPs | OR | FDR", fontsize=5, va="center")
foot = "* prespecified UPRmt module; red = FDR<0.05 across genes. " + (f"No eQTLGen instrument: {', '.join(missing)}." if missing else "")
f.text(0.01, 0.005, foot, fontsize=5)
f.tight_layout(rect=(0, 0.03, 0.93, 1)); save(f, "Fig_MR_D_allgenes_forest")
# E
C = pd.concat([pd.read_csv(os.path.join(res, "coloc_all_genes_GTEx.csv")), pd.read_csv(os.path.join(res, "coloc_all_genes_eQTLGen.csv"))])
oc = ["FinnGen_CC", "FinnGen_SCC", "FinnGen_ADC", "FinnGen_CIS", "UKB_C53", "UKB_180.3", "BBJ_CeC"]
P = C[C.outcome.isin(oc)].groupby(["gene", "outcome"])["PP.H4"].max().unstack()[oc]
P = P.loc[["LONP1"] + sorted([g for g in P.index if g != "LONP1"])]
f, ax = plt.subplots(figsize=(3.6, 0.17 * len(P) + 1.0))
im = ax.imshow(P.values, aspect="auto", cmap="Reds", vmin=0, vmax=1)
for i in range(len(P)):
    for j in range(len(oc)):
        v = P.values[i, j]
        if np.isfinite(v): ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=4.8, color="black" if v < 0.6 else "white")
ax.set_xticks(range(len(oc))); ax.set_xticklabels(["FG\nCC", "FG\nSCC", "FG\nADC", "FG\nCIS", "UKB\nC53", "UKB\nCIN", "BBJ\nCC"], fontsize=6)
ax.set_yticks(range(len(P))); ax.set_yticklabels(P.index, fontsize=6); ax.get_yticklabels()[0].set_fontweight("bold")
cb = f.colorbar(im, ax=ax, fraction=0.06, pad=0.03); cb.set_label("max PP.H4 (eQTLGen + 4 GTEx tissues)", fontsize=6); cb.ax.tick_params(labelsize=5)
ax.set_title("coloc: 16 genes × cervical outcomes", fontsize=6.5)
f.tight_layout(); save(f, "Fig_MR_E_allgenes_coloc")
