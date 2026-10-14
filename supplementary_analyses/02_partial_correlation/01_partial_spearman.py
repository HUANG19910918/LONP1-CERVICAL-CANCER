# 01_partial_spearman.py —— LONP1 共表达的拷贝数偏相关（补充分析 2）
# 口径沿用 v8 15_coexpression_positional.py / 31_chr19_positional_stats.py：
#   表达 = cBioPortal cesc_tcga_pan_can_atlas_2018 data_mrna_seq_v2_rsem.txt（RSEM），Spearman；
#   基因位置/细胞带 = v8 raw_data/30_coexpression_LONP1_CESC_genomewide.csv（GRCh38 + UCSC hg38 cytoBand）；
#   19p13.3 富集 = Fisher 精确检验（双侧，同 31 号脚本），v8 rho>=0.4 时 OR=82.2。
# 新增：偏 Spearman = 全部变量取秩 → 对协变量秩做 OLS 取残差 → 残差 Pearson，df = n-2-k。
#   模型 A（主）：协变量 = LONP1 拷贝数 + g 自身拷贝数（data_log2_cna.txt，连续 log2 ratio）
#   模型 A_thr：同上，拷贝数用 GISTIC2 thresholded（data_cna.txt，-2..2）
#   模型 B（敏感性）：A + 19p 臂水平拷贝数（19p 上全部基因 log2 CN 的样本均值）+ ABSOLUTE 纯度
#   模型 B_arm：A + data_armlevel_cna.txt 的 19p 臂调用（Loss/Unchanged/Gain → -1/0/1）+ 纯度
# 用法: python -I 01_partial_spearman.py <02_partial_correlation 目录> <v8 干实验分析目录>
import sys, os
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, rankdata, t as tdist

out_dir, v8 = sys.argv[1], sys.argv[2]
dd = os.path.join(out_dir, "data"); res = os.path.join(out_dir, "results"); os.makedirs(res, exist_ok=True)
log = []
def L(s=""):
    print(s); log.append(str(s))

def read_cbio(path):
    d = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    d = d.dropna(subset=["Entrez_Gene_Id"])
    d["Entrez_Gene_Id"] = d["Entrez_Gene_Id"].astype(int)
    d = d.drop_duplicates("Entrez_Gene_Id").set_index("Entrez_Gene_Id")
    sym = d["Hugo_Symbol"]
    m = d.drop(columns=[c for c in ("Hugo_Symbol", "Cytoband") if c in d.columns]).apply(pd.to_numeric, errors="coerce")
    return m, sym

expr, esym = read_cbio(os.path.join(dd, "data_mrna_seq_v2_rsem.txt"))
cnl, _ = read_cbio(os.path.join(dd, "data_log2_cna.txt"))
cnt, _ = read_cbio(os.path.join(dd, "data_cna.txt"))
L(f"expression: {expr.shape[0]} genes x {expr.shape[1]} samples; log2 CNA: {cnl.shape}; GISTIC thr: {cnt.shape}")

# ---- 位置注释（v8 CSV 30） ----
pos = pd.read_csv(os.path.join(v8, "raw_data", "30_coexpression_LONP1_CESC_genomewide.csv"))
pos = pos.drop_duplicates("entrez").set_index("entrez")
LON = 9361

# ---- 1. 复现 v8 原始 Spearman（全部有表达的样本） ----
E = expr.dropna(axis=0, how="any")
E = E[(E.std(axis=1) > 0)]
samples_all = list(E.columns)
def spearman_vec(x, M):
    rx = rankdata(x); rx = (rx - rx.mean()) / rx.std()
    R = np.apply_along_axis(rankdata, 1, M)
    R = (R - R.mean(1, keepdims=True)) / R.std(1, keepdims=True)
    return (R * rx).mean(1)
def pval(r, n, k=0):
    df = n - 2 - k
    tt = r * np.sqrt(df / np.clip(1 - r ** 2, 1e-12, None))
    return 2 * tdist.sf(np.abs(tt), df)
rho_raw = pd.Series(spearman_vec(E.loc[LON].values, E.values), index=E.index)
n_raw = E.shape[1]
cmp = pd.DataFrame({"rho_new": rho_raw}).join(pos[["gene", "rho"]], how="inner")
cmp = cmp[cmp.index != LON]
L(f"[复现] 样本 n={n_raw}；与 v8 共有基因 {len(cmp)}；rho 相关 r={np.corrcoef(cmp.rho_new, cmp.rho)[0,1]:.5f}；"
  f"最大绝对差 {np.abs(cmp.rho_new - cmp.rho).max():.4f}；CLPP v8={pos.loc[8192,'rho']:.3f} 新={rho_raw.get(8192, np.nan):.3f}")

# ---- 2. 偏相关共同样本 ----
pur = pd.read_csv(os.path.join(v8, "raw_data", "_cache", "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"), sep="\t")
pur = pur[pur["call status"] == "called"]
pur["sid"] = pur["array"].str[:15]
pur = pur.drop_duplicates("sid").set_index("sid")["purity"]
arm = pd.read_csv(os.path.join(dd, "data_armlevel_cna.txt"), sep="\t")
arm19 = arm[arm["NAME"].astype(str).str.strip() == "19p"].iloc[0]
arm19 = arm19[[c for c in arm19.index if str(c).startswith("TCGA")]].map(
    lambda v: {"Loss": -1, "Unchanged": 0, "Gain": 1}.get(str(v).strip(), np.nan)).astype(float)
L(f"19p 臂调用分布: {arm19.value_counts(dropna=False).to_dict()}")

S = [s for s in samples_all if s in cnl.columns and s in cnt.columns]
L(f"表达 ∩ 拷贝数 样本 n={len(S)}")
genes = [g for g in E.index if g in cnl.index and g in cnt.index and g != LON]
L(f"表达 ∩ 拷贝数 基因 {len(genes)}（无拷贝数而未纳入偏相关: {E.shape[0] - 1 - len(genes)}）")

on19p = pos.cytoband.astype(str).str.startswith("19p")
p19_genes = [g for g in pos.index[on19p] if g in cnl.index]
arm19_cont = cnl.loc[p19_genes, S].mean(0)

def ranks(v): return rankdata(v)
def partial_all(S_, covfun, k_label):
    Ev = E.loc[:, S_]
    y = ranks(Ev.loc[LON].values)
    out = {}
    base = covfun(S_)  # n x k0 shared covariates (ranks)
    Rg = np.apply_along_axis(ranks, 1, Ev.loc[genes].values)
    for i, g in enumerate(genes):
        cg = cnl.loc[g, S_].values if k_label != "thr" else cnt.loc[g, S_].values
        if np.isnan(cg).any():
            out[g] = np.nan; continue
        X = np.column_stack([np.ones(len(S_)), base, ranks(cg)])
        beta_y, *_ = np.linalg.lstsq(X, y, rcond=None); ry = y - X @ beta_y
        beta_x, *_ = np.linalg.lstsq(X, Rg[i], rcond=None); rx = Rg[i] - X @ beta_x
        den = np.sqrt((ry ** 2).sum() * (rx ** 2).sum())
        out[g] = (ry * rx).sum() / den if den > 0 else np.nan
    return pd.Series(out), base.shape[1] + 1

def cov_A(S_): return ranks(cnl.loc[LON, S_].values)[:, None]
def cov_Athr(S_): return ranks(cnt.loc[LON, S_].values)[:, None]
def cov_B(S_): return np.column_stack([ranks(cnl.loc[LON, S_].values), ranks(arm19_cont[S_].values), ranks(pur[S_].values)])
def cov_Barm(S_): return np.column_stack([ranks(cnl.loc[LON, S_].values), ranks(arm19[S_].values), ranks(pur[S_].values)])

S_pur = [s for s in S if s in pur.index and s in arm19.index and not np.isnan(arm19[s])]
L(f"有 ABSOLUTE 纯度 + 19p 臂调用的样本 n={len(S_pur)}")

tab = pd.DataFrame(index=genes)
tab["gene"] = [esym.get(g) for g in genes]
tab["cytoband"] = [pos.cytoband.get(g) if g in pos.index else None for g in genes]
tab["on_19p13_3"] = tab.cytoband.eq("19p13.3")
tab["rho_raw_allS"] = rho_raw.reindex(genes).values
tab["rho_raw"] = spearman_vec(E.loc[LON, S].values, E.loc[genes, S].values)
for nm, f, S_, lab in [("A", cov_A, S, "log2"), ("A_thr", cov_Athr, S, "thr"),
                       ("B", cov_B, S_pur, "log2"), ("B_arm", cov_Barm, S_pur, "log2")]:
    r, k = partial_all(S_, f, lab)
    tab[f"rho_partial_{nm}"] = r.reindex(genes).values
    tab[f"p_partial_{nm}"] = pval(tab[f"rho_partial_{nm}"].values, len(S_), k)
    L(f"模型 {nm}: n={len(S_)}, 协变量数 k={k}, 有效基因 {tab[f'rho_partial_{nm}'].notna().sum()}")
# LONP1 与自身拷贝数
from scipy.stats import spearmanr
L(f"LONP1 表达 ~ LONP1 log2 CN: Spearman rho={spearmanr(E.loc[LON,S], cnl.loc[LON,S]).statistic:.3f} (n={len(S)})")
L(f"19p 臂连续值 vs 臂调用 Spearman={spearmanr(arm19_cont[S_pur], arm19[S_pur]).statistic:.3f}")

# ---- 3. 19p13.3 富集 ----
ann = tab[tab.cytoband.notna()].copy()
L(f"有位置注释基因 N={len(ann)}，其中 19p13.3 = {ann.on_19p13_3.sum()}")
def fisher_at(sel):
    a = int((sel & ann.on_19p13_3).sum()); b = int((sel & ~ann.on_19p13_3).sum())
    c = int((~sel & ann.on_19p13_3).sum()); d = int((~sel & ~ann.on_19p13_3).sum())
    o, p = fisher_exact([[a, b], [c, d]])
    return a, a + b, o, p
rows = []
for col in ["rho_raw_allS", "rho_raw", "rho_partial_A", "rho_partial_A_thr", "rho_partial_B", "rho_partial_B_arm"]:
    v = ann[col]
    for thr in (0.3, 0.4, 0.5):
        a, n, o, p = fisher_at(v >= thr)
        rows.append({"layer": col, "cut": f"rho>={thr}", "n_sel": n, "n_19p13.3": a, "OR": o, "p": p})
    for N in (50, 100, 164):
        sel = v.rank(ascending=False, method="first") <= N
        a, n, o, p = fisher_at(sel)
        rows.append({"layer": col, "cut": f"top{N}", "n_sel": n, "n_19p13.3": a, "OR": o, "p": p})
enr = pd.DataFrame(rows)
enr.to_csv(os.path.join(res, "19p13.3_enrichment_by_layer.csv"), index=False)
L("\n[19p13.3 富集，Fisher 双侧]")
L(enr.to_string(index=False, float_format=lambda x: f"{x:.3g}"))

# ---- 4. 三层 top 基因 ----
def top(col, N=30, excl=False):
    x = tab if not excl else tab[~tab.on_19p13_3]
    return x.sort_values(col, ascending=False).head(N)[["gene", "cytoband", col]]
L("\n[top30 原始]"); L(top("rho_raw").to_string(index=False))
L("\n[top30 剔除 19p13.3]"); L(top("rho_raw", excl=True).to_string(index=False))
L("\n[top30 偏相关 A]"); L(top("rho_partial_A").to_string(index=False))
L("\n[top30 偏相关 B（+19p 臂 + 纯度）]"); L(top("rho_partial_B").to_string(index=False))
key = ["CLPP", "TIMM13", "HSPD1", "HSPE1", "HSPA9", "TRAP1", "PHB1", "PHB2", "ATF4", "ATF5", "DDIT3", "SHMT2", "MTHFD2", "PMPCB"]
L("\n[v8 线粒体模块基因]")
L(tab[tab.gene.isin(key)][["gene", "cytoband", "rho_raw", "rho_partial_A", "rho_partial_A_thr", "rho_partial_B", "rho_partial_B_arm"]]
  .round(3).to_string(index=False))
tab.index.name = "entrez"
tab.sort_values("rho_partial_A", ascending=False).to_csv(os.path.join(res, "LONP1_partial_spearman_genomewide.csv"))
open(os.path.join(res, "01_partial_spearman_log.txt"), "w", encoding="utf-8").write("\n".join(log) + "\n")
print("DONE")
