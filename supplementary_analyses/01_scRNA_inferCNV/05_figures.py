# 05_figures.py —— 分析 1（inferCNV + pseudobulk LONP1）出版级图（PDF 矢量 + 600 dpi PNG；Arial；NPG 配色）
# A 各样本上皮细胞 inferCNV 分类构成（主判定 q95；黑点 = q99 敏感性判定的恶性比例）
# B LONP1 pseudobulk log2CPM：正常上皮参照 / HSIL / 恶性上皮（点 = 样本，颜色 = 数据集；E-MTAB 同患者连线）
# C 同一肿瘤内 恶性 vs 非恶性上皮（左：pseudobulk log2CPM；右：细胞均值 log1p(CP10K)）
# D 恶性细胞 19p13.3 增益 vs 非增益（每单元 >=30 细胞的样本，细胞均值）+ 样本内 LONP1 与 19p13.3 CNV 的 Spearman
# 用法: python -P 05_figures.py <01_scRNA_inferCNV 目录> <cnv 目录> <qc 目录>
import sys, os, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

base, cnvd, qcd = sys.argv[1:4]
res = os.path.join(base, "results"); fig = os.path.join(base, "figures", "publication"); os.makedirs(fig, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "font.size": 7, "axes.linewidth": 0.7, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.labelsize": 8, "legend.fontsize": 6.5, "legend.frameon": False,
                     "xtick.major.width": 0.7, "ytick.major.width": 0.7})
NPG = {"GSE208653": "#E64B35", "E-MTAB-11948": "#3C5488", "GSE197461": "#00A087", "GSE236738": "#F39B7F", "GSE333892": "#8491B4"}
CLS = {"malignant": "#DC0000", "unresolved": "#D9D9D9", "nonmalignant": "#4DBBD5"}
def save(f, n):
    f.savefig(os.path.join(fig, n + ".pdf")); f.savefig(os.path.join(fig, n + ".png"), dpi=600); plt.close(f); print("saved", n)
def short(s):
    return s.split("_", 1)[1] if s.startswith("GSM") else s.replace("_patient", " P")

S = pd.read_csv(os.path.join(res, "infercnv_summary.csv"))
s95 = S[S.rule == "q95"].copy(); s99 = S[S.rule == "q99"].set_index(["dataset", "sample"])
order_g = {"Normal": 0, "AdjacentNormal": 1, "HSIL": 2, "SIL": 2, "Cancer": 3}
s95["og"] = s95.group.map(order_g); s95 = s95.sort_values(["og", "dataset", "sample"]).reset_index(drop=True)

# A
f, ax = plt.subplots(figsize=(7.0, 2.5))
x = np.arange(len(s95)); bottom = np.zeros(len(s95))
for c in ["malignant", "unresolved", "nonmalignant"]:
    v = 100 * s95[c] / s95.n_epi; ax.bar(x, v, bottom=bottom, color=CLS[c], width=0.8, lw=0, label=c.capitalize()); bottom += v
q99 = [s99.loc[(d, s), "pct_malignant"] for d, s in zip(s95.dataset, s95["sample"])]
ax.scatter(x, q99, s=6, c="black", zorder=3, label="Malignant (q99 rule)")
ax.set_xticks(x); ax.set_xticklabels([short(s) for s in s95["sample"]], rotation=90, fontsize=5.5)
for t, d in zip(ax.get_xticklabels(), s95.dataset): t.set_color(NPG[d])
ax.set_ylabel("Epithelial cells (%)"); ax.set_ylim(0, 100); ax.set_xlim(-0.6, len(s95) - 0.4)
gpos = s95.groupby("group", sort=False).apply(lambda g: (g.index.min(), g.index.max()))
for g, (a, b) in gpos.items():
    ax.annotate({"Normal": "Normal", "AdjacentNormal": "Adjacent", "HSIL": "HSIL", "Cancer": "Cancer"}.get(g, g), ((a + b) / 2, 103),
                ha="center", fontsize=7, annotation_clip=False)
    if a > 0: ax.axvline(a - 0.5, color="black", lw=0.5)
ax.legend(ncol=4, loc="upper left", bbox_to_anchor=(0, 1.32))
for i, d in enumerate(NPG):
    if d in set(s95.dataset): ax.text(1.005, 0.9 - i * 0.1, d, color=NPG[d], transform=ax.transAxes, fontsize=6)
f.tight_layout(); save(f, "Fig_S1A_infercnv_classification")

U = pd.read_csv(os.path.join(res, "pseudobulk_LONP1_units.csv")); Ui = U[U.included]
# B
f, ax = plt.subplots(figsize=(2.6, 2.6))
cats = [("normal_epi", "Normal\nepithelium"), ("HSIL_epi", "HSIL\nepithelium"), ("malignant_epi", "Malignant\nepithelium")]
rng = np.random.default_rng(1)
for i, (c, _) in enumerate(cats):
    d = Ui[Ui.category == c]
    jx = i + rng.uniform(-0.15, 0.15, len(d))
    ax.scatter(jx, d.log2CPM_LONP1, s=14, c=[NPG[k] for k in d.dataset], edgecolor="white", lw=0.4, zorder=3)
    ax.hlines(d.log2CPM_LONP1.median(), i - 0.28, i + 0.28, color="black", lw=1.2, zorder=4)
em = Ui[(Ui.dataset == "E-MTAB-11948") & Ui.category.isin(["normal_epi", "malignant_epi"])].copy()
em["pt"] = em["sample"].str.extract(r"(patient\d)")
for pt, g in em.groupby("pt"):
    if g.category.nunique() == 2:
        y0 = g[g.category == "normal_epi"].log2CPM_LONP1.iloc[0]; y1 = g[g.category == "malignant_epi"].log2CPM_LONP1.iloc[0]
        ax.plot([0, 2], [y0, y1], color=NPG["E-MTAB-11948"], lw=0.5, ls=":", zorder=2)
ax.set_xticks(range(3)); ax.set_xticklabels([l for _, l in cats]); ax.set_xlim(-0.5, 2.5)
ax.set_ylabel("LONP1 pseudobulk log2(CPM+1)")
L = open(os.path.join(res, "04_pseudobulk_LONP1_log.txt"), encoding="utf-8").read()
import re
m = re.search(r"OLS\(\+dataset\) diff=([-\d.]+) \(SE ([\d.]+)\) p=([\d.e-]+)", L)
ax.text(0.02, 0.98, f"Malignant vs normal: Δ={float(m.group(1)):.2f}, P={float(m.group(3)):.2f}\n(OLS adjusted for dataset)", transform=ax.transAxes,
        va="top", fontsize=6)
for i, d in enumerate([k for k in NPG if k in set(Ui.dataset)]):
    ax.scatter([], [], s=14, c=NPG[d], label=d)
ax.legend(loc="lower right", handletextpad=0.1, fontsize=5.5)
f.tight_layout(); save(f, "Fig_S1B_LONP1_malignant_vs_normal")

# C
f, axs = plt.subplots(1, 2, figsize=(3.6, 2.4))
for ax, val, lab in [(axs[0], "log2CPM_LONP1", "Pseudobulk log2(CPM+1)"), (axs[1], "mean_lognorm", "Mean per-cell log1p(CP10K)")]:
    w = Ui[Ui.category.isin(["malignant_epi", "nonmal_epi_tumor"])].pivot_table(index=["dataset", "sample"], columns="category", values=val).dropna()
    for (d, s), r in w.iterrows():
        ax.plot([0, 1], [r.nonmal_epi_tumor, r.malignant_epi], color=NPG[d], lw=0.8, marker="o", ms=3)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Non-\nmalignant", "Malignant"]); ax.set_xlim(-0.3, 1.3); ax.set_ylabel(lab)
    from scipy.stats import wilcoxon
    p = wilcoxon(w.malignant_epi, w.nonmal_epi_tumor).pvalue
    ax.set_title(f"n={len(w)} tumors, Wilcoxon P={p:.2f}", fontsize=6.5)
f.tight_layout(); save(f, "Fig_S1C_LONP1_within_tumour")

# D
cells = []
for f_ in sorted(glob.glob(os.path.join(cnvd, "*.csv.gz"))):
    d = pd.read_csv(f_)
    if d.group.iloc[0] != "Cancer": continue
    m = d[(d.cell_type == "Epithelial") & (d.epi_class == "malignant")]
    if len(m) < 100: continue
    import scanpy as sc
    a = sc.read_h5ad(os.path.join(qcd, os.path.basename(f_)[:-7] + ".h5ad"))
    a = a[m.cell.values]
    tot = np.asarray(a.X.sum(1)).ravel(); lon = np.asarray(a[:, "LONP1"].X.todense()).ravel()
    m = m.assign(lonp1=np.log1p(lon / tot * 1e4))
    cells.append(m[["dataset", "sample", "lonp1", "seg_19p13_3", "gain19p13_3"]])
C = pd.concat(cells); C.to_csv(os.path.join(res, "fig_D_malignant_cells_LONP1_19p.csv.gz"), index=False)
f, axs = plt.subplots(1, 2, figsize=(4.2, 2.6), gridspec_kw={"width_ratios": [1, 1]})
from scipy.stats import wilcoxon
ax = axs[0]; pairs = []
for (d, s_), g in C.groupby(["dataset", "sample"]):
    a0 = g[~g.gain19p13_3].lonp1; a1 = g[g.gain19p13_3].lonp1
    if len(a1) >= 30 and len(a0) >= 30: pairs.append((d, s_, a0.mean(), a1.mean()))
P_ = pd.DataFrame(pairs, columns=["dataset", "sample", "nogain", "gain"]); P_.to_csv(os.path.join(res, "fig_D_paired_tumours.csv"), index=False)
for r in P_.itertuples():
    ax.plot([0, 1], [r.nogain, r.gain], color=NPG[r.dataset], lw=0.7, marker="o", ms=2.5, alpha=0.9)
wp = wilcoxon(P_.gain, P_.nogain).pvalue
ax.set_xticks([0, 1]); ax.set_xticklabels(["No 19p13.3\ngain", "19p13.3\ngain"]); ax.set_xlim(-0.3, 1.3)
ax.set_ylabel("LONP1 mean log1p(CP10K)\nmalignant cells per tumor")
ax.set_title(f"{int((P_.gain > P_.nogain).sum())}/{len(P_)} tumors higher with gain\npaired Wilcoxon P={wp:.3f}", fontsize=6.5)
lg = L[L.index("[C3]"):L.index("[C4]")]; mb = re.search(r"beta=([\d.-]+) SE=([\d.]+) p=([\d.e-]+)", lg)
ax.text(0.02, 0.02, f"Cell-level mixed model β={float(mb.group(1)):.2f}, P={float(mb.group(3)):.0e}", transform=ax.transAxes, fontsize=5.5)
for dd in [k for k in NPG if k in set(P_.dataset)]: ax.plot([], [], color=NPG[dd], lw=0.8, label=dd)
ax.legend(loc="upper left", fontsize=5, handlelength=1.2)
ax = axs[1]
rs = C.groupby(["dataset", "sample"])[["lonp1", "seg_19p13_3"]].apply(lambda g: spearmanr(g.lonp1, g.seg_19p13_3).statistic).dropna().reset_index(name="rho")
ax.scatter(rng.uniform(-0.12, 0.12, len(rs)), rs.rho, s=12, c=[NPG[d] for d in rs.dataset], edgecolor="white", lw=0.4)
ax.axhline(0, color="grey", lw=0.6, ls="--"); ax.hlines(rs.rho.median(), -0.25, 0.25, color="black", lw=1.2)
ax.set_xticks([]); ax.set_xlim(-0.5, 0.5); ax.set_ylabel("Within-tumor Spearman ρ\nLONP1 vs 19p13.3 CNV")
ax.set_title(f"{int((rs.rho > 0).sum())}/{len(rs)} tumors ρ>0\nWilcoxon vs 0 P={wilcoxon(rs.rho).pvalue:.3f}", fontsize=6.5)
f.tight_layout(); save(f, "Fig_S1D_LONP1_19p13.3_gain")
