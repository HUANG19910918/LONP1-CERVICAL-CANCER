# 04_figures_LONP1.py —— LONP1 cis-eQTL MR + coloc 出版级图（PDF + 600 dpi PNG；Arial；NPG 配色）
# A：MR 森林图（暴露组织 × 结局；Wald ratio / IVW 的 OR 与 95% CI）
# B：LONP1 ±1 Mb 区段图：GTEx 血液 LONP1 eQTL、FinnGen 宫颈癌、BBJ 宫颈癌的 −log10 P（GRCh37 坐标）；工具变量 SNP 标出
# C：coloc PP.H4 热图：LONP1（4 个组织）与窗口内 24 个有 cis 信号的邻近基因（GTEx 血液）× 结局
# 用法: python -P 04_figures_LONP1.py <03_MR_coloc 目录>
import sys, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
base = sys.argv[1]; res = os.path.join(base, "results"); fig = os.path.join(base, "figures"); os.makedirs(fig, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "font.size": 7, "axes.linewidth": 0.7, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.labelsize": 8, "legend.fontsize": 6, "legend.frameon": False})
NPG = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148"]
def save(f, n):
    for e in ("pdf", "png"): f.savefig(os.path.join(fig, f"{n}.{e}"), dpi=600)
    plt.close(f); print("saved", n)
OUTL = {"FinnGen_CC": "FinnGen cervical cancer (878)", "FinnGen_SCC": "FinnGen SCC (222)", "FinnGen_ADC": "FinnGen ADC (142)",
        "FinnGen_CIS": "FinnGen carcinoma in situ (3,701)", "UKB_C53": "UKB C53 (273)", "UKB_180.1": "UKB phecode 180.1 (299)",
        "UKB_180.3": "UKB CIN 180.3 (2,186)", "BBJ_CeC": "BBJ cervical cancer (967, EAS)"}
EXPL = {"eQTLGen_blood": "eQTLGen blood", "GTEx_blood": "GTEx blood", "GTEx_uterus": "GTEx uterus", "GTEx_vagina": "GTEx vagina",
        "GTEx_fibroblast": "GTEx fibroblast"}
M = pd.read_csv(os.path.join(res, "mr_LONP1.csv")); M = M[M.outcome != "-"].dropna(subset=["OR"])
# A
outs = [o for o in OUTL if o in set(M.outcome)]; exps = [e for e in EXPL if e in set(M.exposure)]
f, ax = plt.subplots(figsize=(4.2, 0.18 * len(outs) * len(exps) + 0.8))
y = 0; yt, yl = [], []
for o in outs:
    for j, e in enumerate(exps):
        r = M[(M.outcome == o) & (M.exposure == e)]
        if not len(r): continue
        r = r.iloc[0]
        ax.errorbar(r.OR, -y, xerr=[[r.OR - r.OR_lo], [r.OR_hi - r.OR]], fmt="o", ms=3, color=NPG[j], elinewidth=0.8, capsize=1.5)
        ax.text(4.3, -y, f"{r.OR:.2f} ({r.OR_lo:.2f}–{r.OR_hi:.2f})  P={r.p:.2f}", va="center", fontsize=5.5)
        y += 1
    yt.append(-(y - (len(exps) + 1) / 2)); yl.append(OUTL[o]); y += 0.8
ax.axvline(1, color="grey", lw=0.6, ls="--"); ax.set_xscale("log"); ax.set_xlim(0.2, 4)
ax.set_xticks([0.25, 0.5, 1, 2, 4]); ax.set_xticklabels(["0.25", "0.5", "1", "2", "4"])
ax.set_yticks(yt); ax.set_yticklabels(yl); ax.tick_params(axis="y", length=0)
ax.set_xlabel("OR per 1-SD increase in genetically predicted LONP1 expression")
for j, e in enumerate(exps): ax.errorbar([], [], fmt="o", ms=3, color=NPG[j], label=EXPL[e])
ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5.5, handletextpad=0.1)
f.tight_layout(); save(f, "Fig_MR_A_LONP1_forest")
# B
I = pd.read_csv(os.path.join(res, "instruments_LONP1.csv"))
tracks = [("GTEx_blood", "FinnGen_CC", "x", "GTEx blood LONP1 eQTL"), ("GTEx_blood", "FinnGen_CC", "y", OUTL["FinnGen_CC"]),
          ("GTEx_blood", "BBJ_CeC", "y", OUTL["BBJ_CeC"])]
f, axs = plt.subplots(len(tracks), 1, figsize=(4.2, 3.6), sharex=True)
lon37 = (5691845, 5720583)
for ax, (e, o, w, lab) in zip(axs, tracks):
    d = pd.read_csv(os.path.join(base, "coloc_input", f"LONP1__{e}__{o}.tsv"), sep="\t")
    d["pos"] = d.key.str.split(":").str[1].astype(int) / 1e6
    p = d.px if w == "x" else d.py
    ax.scatter(d.pos, -np.log10(p), s=2, c="#3C5488" if w == "x" else "#E64B35", lw=0, rasterized=True)
    ins = I[(I.exposure == e)].key.unique()
    h = d[d.key.isin(ins)]
    ax.scatter(h.pos, -np.log10(h.px if w == "x" else h.py), s=18, marker="D", facecolor="none", edgecolor="black", lw=0.7)
    ax.axvspan(lon37[0] / 1e6, lon37[1] / 1e6, color="grey", alpha=0.2, lw=0)
    ax.axhline(-np.log10(5e-8), color="grey", lw=0.5, ls="--")
    ax.set_ylabel("−log10 P"); ax.set_title(lab, fontsize=6.5, loc="left"); ax.set_ylim(0, max(8.5, float(-np.log10(p.min())) * 1.1))
axs[-1].set_xlabel("Chromosome 19 position (Mb, GRCh37)")
axs[0].text(lon37[0] / 1e6, axs[0].get_ylim()[1] * 0.92, " LONP1", fontsize=6)
f.tight_layout(); save(f, "Fig_MR_B_LONP1_regional")
# C
C1 = pd.read_csv(os.path.join(res, "coloc_LONP1.csv")); C2 = pd.read_csv(os.path.join(res, "coloc_LONP1_neighbours.csv"))
nb = pd.read_csv(os.path.join(res, "neighbour_genes_LONP1_window.csv")).set_index("ensg")
C2["label"] = [nb.loc[g, "symbol"] if isinstance(nb.loc[g, "symbol"], str) else g for g in C2.gene]; C2["order"] = [nb.loc[g, "start38"] for g in C2.gene]
C1["label"] = "LONP1 · " + C1.exposure.map(EXPL).str.replace("GTEx ", ""); C1["order"] = -1
oc = ["FinnGen_CC", "FinnGen_CIS", "UKB_180.3", "BBJ_CeC"]
A = pd.concat([C1[C1.outcome.isin(oc)], C2[C2.outcome.isin(oc)]])
P = A.pivot_table(index=["order", "label"], columns="outcome", values="PP.H4").reset_index().sort_values(["order", "label"])
P = P[P.label != "LONP1"]  # 邻近基因中的 LONP1 行与 C1 的 GTEx 血液重复
f, ax = plt.subplots(figsize=(3.0, 0.14 * len(P) + 1.0))
im = ax.imshow(P[oc].values, aspect="auto", cmap="Reds", vmin=0, vmax=1)
for i in range(len(P)):
    for j in range(len(oc)):
        v = P[oc].values[i, j]
        if np.isfinite(v): ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=4.5, color="black" if v < 0.6 else "white")
ax.set_xticks(range(len(oc))); ax.set_xticklabels(["FinnGen\nCC", "FinnGen\nCIS", "UKB\nCIN", "BBJ\nCC"], fontsize=6)
ax.set_yticks(range(len(P))); ax.set_yticklabels(P.label, fontsize=5.5)
for i, l in enumerate(P.label):
    if l.startswith("LONP1"): ax.get_yticklabels()[i].set_fontweight("bold")
cb = f.colorbar(im, ax=ax, fraction=0.06, pad=0.03); cb.set_label("PP.H4 (shared causal variant)", fontsize=6); cb.ax.tick_params(labelsize=5)
ax.set_title("coloc: LONP1 and 19p13.3 neighbours (±1 Mb)", fontsize=6.5)
f.tight_layout(); save(f, "Fig_MR_C_coloc_heatmap")
