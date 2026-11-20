# 01_control_genes.py —— 分析 1：同区域对照基因 / 其他染色体线粒体基因 / 随机表达匹配基因
# 对每个"种子基因"走与 LONP1 完全相同的流程：
#   原始 Spearman；模型 A（种子基因 CN + 伙伴基因 CN）；模型 B（A + 19p 臂连续 CN + ABSOLUTE 纯度）
# 口径沿用 02_partial_correlation/01_partial_spearman.py（全部取秩 → 对协变量秩 OLS 残差 → 残差 Pearson）。
# 富集读数：基因集 AUC（集合内基因的 rho 高于集合外基因的概率；0.5 = 无富集），
#   分别在 全部基因 / 剔除 19p13.3 / 剔除整条 19 号染色体 的背景下计算（后者即分析 3 的口径）。
# 基因集：若提供 GMT（由 00_dump_genesets.R 从 msigdbr 导出）则用 MSigDB；否则用按基因符号前缀定义的集合。
# 用法: python -I 01_control_genes.py <02_partial_correlation 目录> <干实验分析目录> <输出目录> [genesets.gmt]
import sys, os, re, glob
import numpy as np, pandas as pd
from scipy.stats import rankdata, fisher_exact

here = os.path.dirname(os.path.abspath(__file__))
if len(sys.argv) > 3:
    import glob
here = os.path.dirname(os.path.abspath(__file__))
if len(sys.argv) > 3:
    base, v8, out = sys.argv[1], sys.argv[2], sys.argv[3]
    gmt = sys.argv[4] if len(sys.argv) > 4 else None
else:  # no arguments: locate everything relative to this script (…/supplementary-analyses/06_control_genes/)
    out = here
    base = os.path.join(os.path.dirname(here), "02_partial_correlation")
    root = os.path.dirname(os.path.dirname(here))
    hits = glob.glob(os.path.join(root, "*", "raw_data", "30_coexpression_LONP1_CESC_genomewide.csv"))
    assert hits, "cannot find raw_data/30_coexpression_LONP1_CESC_genomewide.csv"
    v8 = os.path.dirname(os.path.dirname(hits[0]))
    gmt = os.path.join(here, "genesets_selected.gmt")
os.makedirs(out, exist_ok=True)
dd = os.path.join(base, "data")
log = []
def L(s=""):
    print(s, flush=True); log.append(str(s))

def read_cbio(path):
    d = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    d = d.dropna(subset=["Entrez_Gene_Id"]); d["Entrez_Gene_Id"] = d["Entrez_Gene_Id"].astype(int)
    d = d.drop_duplicates("Entrez_Gene_Id").set_index("Entrez_Gene_Id")
    sym = d["Hugo_Symbol"]
    m = d.drop(columns=[c for c in ("Hugo_Symbol", "Cytoband") if c in d.columns]).apply(pd.to_numeric, errors="coerce")
    return m, sym

expr, esym = read_cbio(os.path.join(dd, "data_mrna_seq_v2_rsem.txt"))
cnl, _ = read_cbio(os.path.join(dd, "data_log2_cna.txt"))
pos = pd.read_csv(os.path.join(v8, "raw_data", "30_coexpression_LONP1_CESC_genomewide.csv")).drop_duplicates("entrez").set_index("entrez")
pur = pd.read_csv(os.path.join(v8, "raw_data", "_cache", "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"), sep="\t")
pur = pur[pur["call status"] == "called"]; pur["sid"] = pur["array"].str[:15]
pur = pur.drop_duplicates("sid").set_index("sid")["purity"]
LON = 9361

E = expr.dropna(axis=0, how="any"); E = E[E.std(axis=1) > 0]
S = [s for s in E.columns if s in cnl.columns]
C = cnl.loc[:, S]
genes = [g for g in E.index if g in C.index and not C.loc[g].isna().any()]
E = E.loc[genes, S]; C = C.loc[genes]
sym = esym.reindex(genes).astype(str)
cyto = pos.cytoband.reindex(genes).astype(str)
chrom = pos.chr.reindex(genes).astype(str)
cyto.loc[LON] = "19p13.3"; chrom.loc[LON] = "19"   # LONP1 itself is absent from the partner table
on19p133 = (cyto == "19p13.3").values
on19 = (chrom == "19").values | cyto.str.startswith("19").values
# model B sample set = 02_partial_correlation (n=240): ABSOLUTE purity AND a non-missing 19p arm-level call
# (the call itself is not a covariate; the filter keeps models B and B_arm on the same samples)
arm = pd.read_csv(os.path.join(dd, "data_armlevel_cna.txt"), sep="\t")
a19 = arm[arm["NAME"].astype(str).str.strip() == "19p"].iloc[0]
a19 = a19[[c for c in a19.index if str(c).startswith("TCGA")]].map(
    lambda v: {"Loss": -1, "Unchanged": 0, "Gain": 1}.get(str(v).strip(), np.nan)).astype(float)
S_pur = [s for s in S if s in pur.index and not np.isnan(pur[s]) and s in a19.index and not np.isnan(a19[s])]
ip = np.array([S.index(s) for s in S_pur])
p19 = [g for g in pos.index[pos.cytoband.astype(str).str.startswith("19p")] if g in cnl.index]
arm19 = cnl.loc[p19, S].mean(axis=0).values
L(f"genes with expression+CN: {len(genes)}; samples n={len(S)} (model A), n={len(S_pur)} (model B); 19p13.3 genes: {on19p133.sum()}; chr19 genes: {on19.sum()}")

def rk(M):
    return np.apply_along_axis(rankdata, -1, M).astype(np.float64)
def zs(M):
    M = M - M.mean(-1, keepdims=True); sd = np.sqrt((M ** 2).sum(-1, keepdims=True)); sd[sd == 0] = np.nan
    return M / sd

Ev, Cv = E.values, C.values
RE_A, RC_A = rk(Ev), rk(Cv)
RE_B, RC_B = rk(Ev[:, ip]), rk(Cv[:, ip])
ZE_A = zs(RE_A)
logE = np.log2(Ev + 1); mean_expr = logE.mean(1)
cn_expr_rho = np.nansum(zs(RE_A) * zs(RC_A), axis=1)
gi = {g: i for i, g in enumerate(genes)}

def resid(M, Z):
    Q, _ = np.linalg.qr(Z)
    return M - (M @ Q) @ Q.T

def partial(i, RE, RC, extra):
    """partial Spearman of seed i with every gene, adjusting for [seed CN, extra..., partner CN]."""
    n = RE.shape[1]
    Z = np.column_stack([np.ones(n), RC[i]] + [rankdata(e) for e in extra])
    x = resid(RE[i][None, :], Z)[0]
    Y = resid(RE, Z); K = resid(RC, Z)
    kn = np.sqrt((K ** 2).sum(1)); ok = kn > 1e-6 * np.sqrt(n) * n   # partner CN not collinear with the fixed covariates
    xn = x / np.sqrt((x ** 2).sum())
    Yn = Y / np.sqrt((Y ** 2).sum(1, keepdims=True))
    Kn = np.zeros_like(K); Kn[ok] = K[ok] / kn[ok, None]
    a = Yn @ xn; b = Kn @ xn; c = (Yn * Kn).sum(1)
    r = np.where(ok, (a - b * c) / np.sqrt(np.clip((1 - b ** 2) * (1 - c ** 2), 1e-12, None)), a)
    r[i] = np.nan
    return r

# ---------- gene sets ----------
up = sym.str.upper().values
sets = {}
if gmt and os.path.exists(gmt):
    for line in open(gmt, encoding="utf8"):
        f = line.rstrip("\n").split("\t")
        sets[f[0]] = set(f[2:])
    L(f"gene sets from GMT: {len(sets)}")
else:
    rx = {"mitoribosome_MRPL_MRPS": r"^MRP[LS]\d+$",
          "OXPHOS_subunits": r"^(NDUF[ABCSV]\d+[A-Z]?\d*|SDH[ABCD]|UQCR(C[12]|FS1|B|H|Q|10|11)|CYC1|COX(4I1|5A|5B|6A1|6B1|6C|7A2|7B|7C|8A)|ATP5\w+)$",
          "mito_import_TIMM_TOMM": r"^(TIMM|TOMM)\d+[A-Z]?\d*$",
          "cytosolic_ribosome_RPL_RPS": r"^RP[LS]\d+[AXY]?\d*$"}
    for k, pat in rx.items():
        sets[k] = set(sym[sym.str.upper().str.match(pat)].values)
    L("gene sets (symbol-prefix definitions): " + ", ".join(f"{k} n={len(v)}" for k, v in sets.items()))
mito_any = set().union(*[v for k, v in sets.items() if "cytosolic" not in k and "CYTOPLASMIC" not in k and "RIBOSOME_BIOGENESIS" not in k])
if "GOCC_MITOCHONDRION" in sets:
    mito_any = set(sets["GOCC_MITOCHONDRION"])
MANUAL_19P133_MITO = {"LONP1", "CLPP", "NDUFA11", "NDUFS7", "UQCR11", "ATP5F1D", "ATP5D", "TIMM13", "POLRMT", "MRPL54", "THOP1",
                      "MICOS13", "C19orf70", "GPX4", "SLC25A23", "SLC25A41", "PET100", "TIMM44", "MRPL4", "NDUFA7", "ECSIT",
                      "GADD45GIP1", "TMEM160", "ATP8B3", "DOHH", "PTCD", "OXLD1", "FASTK", "MRPS", "AURKAIP1", "MARCHF2", "TRMT1",
                      "ACSBG2", "ALKBH7", "NDUFB7", "TUFM", "GTPBP3", "PPAN", "QTRT1", "DNM2", "TIMM29", "FDX2", "FDX1L", "RDH13", "PCSK4"}
is_mito = np.array([(s in mito_any) or (s in MANUAL_19P133_MITO) for s in sym.values])
set_idx = {k: np.array([s in v for s in sym.values]) for k, v in sets.items()}
set_idx = {k: v for k, v in set_idx.items() if 10 <= v.sum() <= 600}

def auc(r, member, universe):
    m = universe & ~np.isnan(r)
    rr = rankdata(r[m]); mem = member[m]
    n1 = mem.sum(); n0 = (~mem).sum()
    if n1 < 5: return np.nan
    return (rr[mem].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

# ---------- seeds ----------
rng = np.random.default_rng(1)
iL = gi[LON]
seeds = {}
for i in np.where(on19p133)[0]:
    seeds[i] = "19p13.3"
PANEL = ["AFG3L2", "YME1L1", "SPG7", "PMPCB", "HSPD1", "HSPE1", "HSPA9", "TRAP1", "DNAJA3", "PHB", "PHB1", "PHB2", "TFAM", "CLPX",
         "LRPPRC", "GRPEL1", "TUFM", "MRPL12", "MRPS12", "NDUFS1", "SDHA", "UQCRC1", "COX4I1", "ATP5A1", "ATP5F1A", "TIMM23",
         "TOMM20", "CS", "SSBP1", "POLG", "TIMM44", "HSCB", "IMMT", "CHCHD4"]
for s in PANEL:
    hit = np.where((sym.values == s) & ~on19)[0]
    if len(hit): seeds.setdefault(hit[0], "mito_other_chr")
cand = np.where(~on19 & ~is_mito & (np.abs(mean_expr - mean_expr[iL]) <= 1.0) & ~np.isin(np.arange(len(genes)), list(seeds)))[0]
for i in rng.choice(cand, size=min(300, len(cand)), replace=False):
    seeds[i] = "random_expr_matched_non_mito"
L(f"seeds: " + str(pd.Series(seeds).value_counts().to_dict()) + f"; LONP1 mean log2(RSEM+1)={mean_expr[iL]:.2f}, CN-expr rho={cn_expr_rho[iL]:.3f}")

U_all = np.ones(len(genes), bool); U_x133 = ~on19p133; U_x19 = ~on19
rows = []
check = None
for n_done, (i, kind) in enumerate(seeds.items()):
    layers = {"raw": ZE_A @ ZE_A[i], "A": partial(i, RE_A, RC_A, []), "B": partial(i, RE_B, RC_B, [arm19[ip], pur[S_pur].values])}
    layers["raw"][i] = np.nan
    if i == iL: check = layers
    row = dict(entrez=genes[i], gene=sym.values[i], seed_class=kind, cytoband=cyto.values[i], is_mito=bool(is_mito[i]),
               mean_log2_expr=mean_expr[i], cn_expr_rho=cn_expr_rho[i])
    for lay, r in layers.items():
        sel = r >= 0.4
        row[f"n_rho04_{lay}"] = int(np.nansum(sel))
        a = int((sel & on19p133).sum()); b = int((sel & ~on19p133).sum())
        c = int((~sel & on19p133 & ~np.isnan(r)).sum()); d = int((~sel & ~on19p133 & ~np.isnan(r)).sum())
        row[f"n19p133_rho04_{lay}"] = a
        row[f"OR19p133_{lay}"] = fisher_exact([[a, b], [c, d]])[0] if (a + b) > 0 else np.nan
        for k, mem in set_idx.items():
            mm = mem.copy(); mm[i] = False
            row[f"AUC_{k}|{lay}|all"] = auc(r, mm, U_all)
            row[f"AUC_{k}|{lay}|excl19p13.3"] = auc(r, mm, U_x133)
            row[f"AUC_{k}|{lay}|exclChr19"] = auc(r, mm, U_x19)
    rows.append(row)
    if n_done % 50 == 0: print("seed", n_done, "/", len(seeds), flush=True)
T = pd.DataFrame(rows)
T.to_csv(os.path.join(out, "control_genes_seed_table.csv"), index=False)

# ---------- validation against the published LONP1 table ----------
old = pd.read_csv(os.path.join(base, "results", "LONP1_partial_spearman_genomewide.csv")).set_index("entrez")
for lay, col in [("raw", "rho_raw"), ("A", "rho_partial_A"), ("B", "rho_partial_B")]:
    new = pd.Series(check[lay], index=genes).reindex(old.index)
    m = new.notna() & old[col].notna()
    L(f"[validation] LONP1 {lay}: n={m.sum()} r={np.corrcoef(new[m], old.loc[m, col])[0, 1]:.5f} max|diff|={np.abs(new[m] - old.loc[m, col]).max():.4f}")

# ---------- summary ----------
lon = T[T.entrez == LON].iloc[0]
T19 = T[T.seed_class == "19p13.3"]
matched = T19[(~T19.is_mito) & (np.abs(T19.mean_log2_expr - lon.mean_log2_expr) <= 1.0) & (np.abs(T19.cn_expr_rho - lon.cn_expr_rho) <= 0.10)]
L(f"\n19p13.3 seeds: {len(T19)}; non-mito matched controls (expr +/-1 log2, CN-expr rho +/-0.10): {len(matched)}: {', '.join(matched.gene)}")
def pct(v, ref):
    ref = ref.dropna(); return 100.0 * (ref < v).mean() if len(ref) else np.nan
groups = {"19p13.3 all": T19[T19.entrez != LON], "19p13.3 non-mito": T19[(~T19.is_mito)], "19p13.3 matched non-mito": matched,
          "mito other chr": T[T.seed_class == "mito_other_chr"], "random matched non-mito": T[T.seed_class == "random_expr_matched_non_mito"]}
summ = []
for col in [c for c in T.columns if c.startswith("AUC_")] + [c for c in T.columns if c.startswith("OR19p133_")]:
    r = dict(readout=col, LONP1=lon[col])
    for gname, G in groups.items():
        r[f"{gname} median"] = G[col].median(); r[f"{gname} IQR"] = f"{G[col].quantile(.25):.3f}-{G[col].quantile(.75):.3f}"
        r[f"LONP1 pct in {gname}"] = pct(lon[col], G[col])
    summ.append(r)
SM = pd.DataFrame(summ); SM.to_csv(os.path.join(out, "control_genes_summary.csv"), index=False)
pd.set_option("display.width", 300); pd.set_option("display.max_columns", 50); pd.set_option("display.max_rows", 500)
L(SM[["readout", "LONP1"] + [c for c in SM.columns if c.endswith("median") or c.startswith("LONP1 pct")]].round(3).to_string(index=False))
open(os.path.join(out, "01_control_genes_log.txt"), "w", encoding="utf-8").write("\n".join(log) + "\n")
print("DONE")
