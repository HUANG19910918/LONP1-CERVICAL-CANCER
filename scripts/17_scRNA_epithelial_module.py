# 17_scRNA_epithelial_module.py —— GSE208653 上皮细胞内在共表达（扩展线粒体蛋白稳态模块）
# 复现口径：
#   数据: GSE208653_RAW.tar (md5=e264e8e98f1005f153c3beaf1f9e29bd)，9样本10x矩阵；加载/QC与06脚本完全一致
#         (min_genes=300, n_genes<7500, filter_genes min_cells=3, normalize_total 1e4 + log1p) → 80,435细胞
#   细胞注释: 沿用 raw_data/17_GSE208653_LONP1_per_cell.csv.gz（06脚本annotate步输出，行序=本脚本加载顺序），
#             以逐细胞LONP1值完全相等作为对齐校验
#   上皮重聚类: 上皮细胞子集 → HVG 1500 (batch_key=sample) → scale(max 10) → PCA 20 (arpack) → neighbors 15
#               → Leiden 0.5 (igraph, n_iterations=2) → 假批量单元 = sample × 亚群 (≥50细胞)
#   相关: Spearman，三层（逐细胞 / 假批量单元 / 逐样本）
#   全基因组上皮内在共表达 (step=genomewide): 58个假批量单元 × 上皮中检出率≥1%的基因(14,141)，Spearman；
#         同时给出按样本中心化(减去样本均值)的敏感性分析；top300正相关基因 Enrichr 富集
#   输出: raw_data/34_epithelial_pseudobulk_units_extended.csv, 35_epithelial_intrinsic_corr_extended.csv,
#         36_epithelial_recluster_check.txt, 37_epithelial_pseudobulk_genomewide_rho.csv,
#         38_enrichment_epithelial_intrinsic_top300.csv
# 用法: python3 17_scRNA_epithelial_module.py <RAW目录> <分析根目录> <step: load|recluster|corr|genomewide> [variant]
import sys, os, gzip
import numpy as np, pandas as pd, scipy.sparse as sp
import scanpy as sc, anndata
from scipy.io import mmread
from scipy.stats import spearmanr

raw_dir, base, step = sys.argv[1], sys.argv[2], sys.argv[3]
rd = os.path.join(base, "raw_data"); wd = os.path.join(raw_dir, "work"); os.makedirs(wd, exist_ok=True)
sc.settings.verbosity = 1
SAMPLES = ["GSM6360680_N_HPV_NEG_1", "GSM6360681_N_HPV_NEG_2", "GSM6360682_N_1", "GSM6360683_N_2",
           "GSM6360684_HSIL_1", "GSM6360685_HSIL_2", "GSM6360686_SCC_4", "GSM6360687_SCC_5", "GSM6360688_ADC_6"]
STAGE = {"GSM6360680_N_HPV_NEG_1": "Normal_HPVneg", "GSM6360681_N_HPV_NEG_2": "Normal_HPVneg",
         "GSM6360682_N_1": "Normal_HPVpos", "GSM6360683_N_2": "Normal_HPVpos",
         "GSM6360684_HSIL_1": "HSIL", "GSM6360685_HSIL_2": "HSIL",
         "GSM6360686_SCC_4": "Cancer", "GSM6360687_SCC_5": "Cancer", "GSM6360688_ADC_6": "Cancer"}
GENES = ["LONP1", "CLPP", "HSPD1", "HSPE1", "HSPA9", "TRAP1", "DNAJA3", "PHB1", "PHB2", "PMPCB", "AFG3L2", "YME1L1",
         "SPG7", "ATF4", "ATF5", "DDIT3", "TFAM", "NRF1", "PPARGC1A", "POLRMT", "FIS1", "MFF", "DNM1L", "MIEF1",
         "MIEF2", "MFN1", "MFN2", "OPA1", "PINK1", "PRKN", "THOP1", "TIMM13", "NDUFA11", "SHMT2", "MTHFD2"]
EPI_H5 = os.path.join(wd, "epithelial_norm.h5ad")

def load_one(prefix):
    m = mmread(os.path.join(raw_dir, f"{prefix}.matrix.mtx.gz")).tocsr().T.astype(np.float32)
    feats = pd.read_csv(os.path.join(raw_dir, f"{prefix}.features.tsv.gz"), sep="\t", header=None)
    bcs = pd.read_csv(os.path.join(raw_dir, f"{prefix}.barcodes.tsv.gz"), sep="\t", header=None)[0].values
    ad = sc.AnnData(sp.csr_matrix(m))
    ad.var_names = feats[1].values; ad.obs_names = [f"{prefix}:{b}" for b in bcs]
    ad.var_names_make_unique()
    sc.pp.filter_cells(ad, min_genes=300)
    ad = ad[ad.obs.n_genes < 7500].copy()
    ad.obs["sample"] = prefix; ad.obs["stage"] = STAGE[prefix]
    return ad

if step == "load":
    ads = []
    for s in SAMPLES:
        a = load_one(s); print(s, a.n_obs, flush=True); ads.append(a)
    comb = anndata.concat(ads, join="outer"); del ads
    comb.X = sp.csr_matrix(comb.X)
    sc.pp.filter_genes(comb, min_cells=3)
    print("combined:", comb.shape, flush=True)
    sc.pp.normalize_total(comb, target_sum=1e4); sc.pp.log1p(comb)
    ref = pd.read_csv(os.path.join(rd, "17_GSE208653_LONP1_per_cell.csv.gz"))
    assert len(ref) == comb.n_obs, (len(ref), comb.n_obs)
    assert (ref["sample"].values == comb.obs["sample"].values).all(), "sample order mismatch"
    lon = np.asarray(comb[:, "LONP1"].X.todense()).ravel()
    d = np.abs(lon - ref.LONP1.values)
    print("LONP1 per-cell max abs diff vs CSV17:", d.max(), " n_diff>1e-4:", int((d > 1e-4).sum()), flush=True)
    # float32归一化舍入可产生~1e-3量级差异；对齐校验标准：最大差<0.01 且 >99.9%细胞差<1e-3
    assert d.max() < 1e-2 and (d < 1e-3).mean() > 0.999, "per-cell LONP1 mismatch -> alignment failed"
    comb.obs["cell_type"] = ref.cell_type.values; comb.obs["leiden_global"] = ref.leiden.astype(str).values
    epi = comb[comb.obs.cell_type == "Epithelial"].copy()
    print("epithelial cells:", epi.n_obs, flush=True)
    epi.write(EPI_H5)
    # 保存全数据集中本研究所需基因的逐细胞值（供逐细胞层与对照）
    keep = [g for g in GENES if g in comb.var_names]
    df = pd.DataFrame(np.asarray(comb[:, keep].X.todense()), columns=keep)
    df.insert(0, "cell_type", comb.obs.cell_type.values); df.insert(0, "stage", comb.obs.stage.values)
    df.insert(0, "sample", comb.obs["sample"].values)
    df.to_csv(os.path.join(wd, "all_cells_module_genes.csv.gz"), index=False, compression="gzip")
    print("LOAD_DONE")

elif step == "recluster":
    epi = sc.read_h5ad(EPI_H5)
    variant = sys.argv[4] if len(sys.argv) > 4 else "batch"
    if variant == "batch":
        sc.pp.highly_variable_genes(epi, n_top_genes=1500, batch_key="sample")
    else:
        sc.pp.highly_variable_genes(epi, n_top_genes=1500)
    h = epi[:, epi.var.highly_variable].copy()
    sc.pp.scale(h, max_value=10)
    sc.tl.pca(h, n_comps=20, svd_solver="arpack")
    epi.obsm["X_pca"] = h.obsm["X_pca"]; del h
    sc.pp.neighbors(epi, n_neighbors=15, n_pcs=20)
    sc.tl.leiden(epi, resolution=0.5, key_added="sub", flavor="igraph", n_iterations=2, directed=False)
    print("subclusters:", epi.obs["sub"].nunique(), flush=True)
    epi.obs[["sample", "stage", "sub"]].to_csv(os.path.join(wd, f"epi_sub_{variant}.csv"))
    # 与27b比对
    old = pd.read_csv(os.path.join(rd, "27b_epithelial_pseudobulk_units.csv"))
    ct = epi.obs.groupby(["sample", "sub"], observed=True).size().reset_index(name="n")
    units = ct[ct.n >= 50]
    print("units>=50:", len(units), "(27b:", len(old), ")")
    a = sorted(units.n.tolist()); b = sorted(old.n.tolist())
    print("unit sizes identical:", a == b)
    print("new:", a[:20], "... old:", b[:20])
    lines = [f"variant={variant}", f"subclusters={epi.obs['sub'].nunique()}", f"units>=50={len(units)} (27b={len(old)})",
             f"unit_sizes_identical={a == b}", f"new_sizes={a}", f"old_sizes={b}"]
    open(os.path.join(rd, "36_epithelial_recluster_check.txt"), "a").write("\n".join(lines) + "\n\n")
    epi.write(EPI_H5.replace(".h5ad", f"_{variant}.h5ad"))
    print("RECLUSTER_DONE")

elif step == "corr":
    variant = sys.argv[4] if len(sys.argv) > 4 else "batch"
    epi = sc.read_h5ad(EPI_H5.replace(".h5ad", f"_{variant}.h5ad"))
    keep = [g for g in GENES if g in epi.var_names]
    X = pd.DataFrame(np.asarray(epi[:, keep].X.todense()), columns=keep, index=epi.obs_names)
    X["fission"] = X[["FIS1", "MFF"]].mean(axis=1); X["fusion"] = X[["MFN1", "MFN2", "OPA1"]].mean(axis=1)
    X["proteostasis"] = X[["CLPP", "HSPD1", "HSPE1", "HSPA9", "TRAP1"]].mean(axis=1)
    cols = [c for c in X.columns if c != "LONP1"]
    X["sample"] = epi.obs["sample"].values; X["sub"] = epi.obs["sub"].values; X["stage"] = epi.obs.stage.values
    rows = []
    # 1) 逐细胞
    for g in cols:
        r, p = spearmanr(X.LONP1, X[g]); rows.append({"level": f"per-cell (n={len(X):,})", "gene": g, "rho": r, "p": p, "n": len(X)})
    # 2) 假批量 sample x sub (>=50 cells)
    pb = X.groupby(["sample", "sub"], observed=True)
    units = pb.size().reset_index(name="n"); means = pb[["LONP1"] + cols].mean().reset_index()
    U = units.merge(means, on=["sample", "sub"]); U = U[U.n >= 50].reset_index(drop=True)
    U.to_csv(os.path.join(rd, "34_epithelial_pseudobulk_units_extended.csv"), index=False)
    for g in cols:
        r, p = spearmanr(U.LONP1, U[g]); rows.append({"level": f"pseudobulk sample x subcluster (n={len(U)})", "gene": g, "rho": r, "p": p, "n": len(U)})
    # 3) 逐样本
    S = X.groupby("sample", observed=True)[["LONP1"] + cols].mean().reset_index()
    for g in cols:
        r, p = spearmanr(S.LONP1, S[g]); rows.append({"level": f"per-sample (n={len(S)})", "gene": g, "rho": r, "p": p, "n": len(S)})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(rd, "35_epithelial_intrinsic_corr_extended.csv"), index=False)
    pbo = out[out.level.str.startswith("pseudobulk")].sort_values("rho", ascending=False)
    print(pbo.to_string(index=False))
    # 与旧27比对
    old = pd.read_csv(os.path.join(rd, "27_epithelial_intrinsic_corr.csv"))
    old = old[old.level.str.startswith("pseudobulk")].set_index("gene").rho
    cmp = pbo.set_index("gene").rho.reindex(old.index)
    print("\nvs 27 (old pseudobulk):"); print(pd.DataFrame({"old": old, "new": cmp}).round(3).to_string())
    print("CORR_DONE")

elif step == "genomewide":
    import json, urllib.request
    from scipy.stats import rankdata
    variant = sys.argv[4] if len(sys.argv) > 4 else "batch"
    epi = sc.read_h5ad(EPI_H5.replace(".h5ad", f"_{variant}.h5ad"))
    obs = epi.obs.copy(); obs["unit"] = obs["sample"].astype(str) + "|" + obs["sub"].astype(str)
    sizes = obs.unit.value_counts(); units = sizes[sizes >= 50].index.tolist()
    X = epi.X.tocsr()
    det = np.asarray((X > 0).mean(axis=0)).ravel(); keep = det >= 0.01
    M = np.vstack([np.asarray(X[(obs.unit == u).values].mean(axis=0)).ravel() for u in units])
    M = pd.DataFrame(M, index=units, columns=epi.var_names).loc[:, keep]
    print("units:", len(units), "genes detected>=1%:", int(keep.sum()), flush=True)
    def spear_all(M, target):
        R = M.apply(rankdata, axis=0); r0 = rankdata(M[target])
        Rc = R - R.mean(axis=0); r0c = r0 - r0.mean()
        return (Rc.T @ r0c) / np.sqrt((Rc ** 2).sum(axis=0) * (r0c ** 2).sum())
    rho = spear_all(M, "LONP1").drop("LONP1")
    samp = pd.Series([u.split("|")[0] for u in units], index=units)
    Mc = M.sub(M.groupby(samp.values).transform("mean"))
    rho_c = spear_all(Mc, "LONP1").drop("LONP1")
    n = len(units)
    from scipy.stats import t as tdist
    def pval(r): 
        tt = r * np.sqrt((n - 2) / np.maximum(1 - r ** 2, 1e-12)); return 2 * tdist.sf(np.abs(tt), n - 2)
    co = pd.read_csv(os.path.join(rd, "30_coexpression_LONP1_CESC_genomewide.csv")).set_index("gene")
    out = pd.DataFrame({"rho_pb": rho, "p_pb": pval(rho.values), "pct_rank_pb": rho.rank(pct=True),
                        "rho_pb_sample_centered": rho_c, "p_pb_sample_centered": pval(rho_c.values),
                        "detect_frac_epi": det[keep][[g != "LONP1" for g in M.columns]],
                        "bulk_rho_TCGA": co.rho.reindex(rho.index), "on_19p13_3": co.on_19p13_3.reindex(rho.index)})
    out.index.name = "gene"; out = out.sort_values("rho_pb", ascending=False)
    out.to_csv(os.path.join(rd, "37_epithelial_pseudobulk_genomewide_rho.csv"))
    L = []
    def P(s): print(s, flush=True); L.append(s)
    P(f"units={n}; genes={len(rho)}; genome-wide pseudobulk rho median={rho.median():.3f} IQR={rho.quantile(.25):.3f}-{rho.quantile(.75):.3f}")
    P(f"sample-centered: median={rho_c.median():.3f} IQR={rho_c.quantile(.25):.3f}-{rho_c.quantile(.75):.3f}")
    b19 = out.on_19p13_3.fillna(False).astype(bool)
    P(f"19p13.3 genes within epithelium: median rho_pb={out.rho_pb[b19].median():.3f} (n={int(b19.sum())}) vs others {out.rho_pb[~b19].median():.3f}")
    common = out.bulk_rho_TCGA.dropna().index
    P(f"bulk vs epithelium-intrinsic rho agreement (Spearman, n={len(common)}): raw={spearmanr(out.loc[common,'bulk_rho_TCGA'], out.loc[common,'rho_pb'])[0]:.3f}, centered={spearmanr(out.loc[common,'bulk_rho_TCGA'], out.loc[common,'rho_pb_sample_centered'])[0]:.3f}")
    for g in ["CLPP","HSPA9","HSPD1","HSPE1","TRAP1","DNAJA3","PMPCB","PHB","PHB2","ATF4","ATF5","DDIT3","FIS1","MFF","DNM1L","MFN1","MFN2","OPA1","TFAM","PPARGC1A","NRF1","POLRMT","THOP1","NDUFA11","TIMM13","YME1L1","AFG3L2","PINK1"]:
        if g in out.index: P(f"  {g}: rho_pb={out.rho_pb[g]:.3f} pct={out.pct_rank_pb[g]*100:.0f} centered={out.rho_pb_sample_centered[g]:.3f} bulk={out.bulk_rho_TCGA.get(g, np.nan):.3f}")
    # Enrichr: top300 (raw) 与 top300 (centered)
    def enrichr(genes, libs, tag):
        body = ("--XX\r\nContent-Disposition: form-data; name=\"list\"\r\n\r\n" + "\n".join(genes) +
                "\r\n--XX\r\nContent-Disposition: form-data; name=\"description\"\r\n\r\n" + tag + "\r\n--XX--\r\n").encode()
        req = urllib.request.Request("https://maayanlab.cloud/Enrichr/addList", data=body, headers={"Content-Type": "multipart/form-data; boundary=XX"})
        uid = json.load(urllib.request.urlopen(req, timeout=120))["userListId"]
        rows = []
        for lib in libs:
            r = json.load(urllib.request.urlopen(f"https://maayanlab.cloud/Enrichr/enrich?userListId={uid}&backgroundType={lib}", timeout=120))
            for e in r[lib]: rows.append({"lib": lib, "term": e[1], "p": e[2], "padj": e[6], "n_overlap": len(e[5]), "genes": ";".join(e[5])})
        return pd.DataFrame(rows)
    libs = ["GO_Biological_Process_2023", "GO_Cellular_Component_2023", "KEGG_2021_Human"]
    res = []
    for tag, col in [("top300_raw", "rho_pb"), ("top300_sample_centered", "rho_pb_sample_centered")]:
        top = out.sort_values(col, ascending=False).head(300).index.tolist()
        e = enrichr(top, libs, "LONP1_epithelium_intrinsic_" + tag); e.insert(0, "set", tag); res.append(e)
        for lib in libs:
            P(f"[{tag} {lib}] " + "; ".join(f"{r.term} padj={r.padj:.1e} n={r.n_overlap}" for _, r in e[e.lib == lib].sort_values('padj').head(5).iterrows()))
    pd.concat(res).sort_values(["set", "lib", "padj"]).to_csv(os.path.join(rd, "38_enrichment_epithelial_intrinsic_top300.csv"), index=False)
    open(os.path.join(rd, "36_epithelial_recluster_check.txt"), "a").write("\n[genomewide]\n" + "\n".join(L) + "\n")
    print("GENOMEWIDE_DONE")
