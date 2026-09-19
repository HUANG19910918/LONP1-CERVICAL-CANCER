# 06_scRNA_GSE208653.py —— GSE208653宫颈癌进展队列单细胞分析：LONP1细胞类型定位
# 复现口径：
#   数据：GSE208653_RAW.tar (md5=e264e8e98f1005f153c3beaf1f9e29bd)，9样本10x矩阵
#   （作者已用CellRanger v5.0.1+GRCh38定量，>200基因/细胞、线粒体UMI<20%过滤，且矩阵中已移除线粒体基因）
#   本脚本QC：min_genes=300, min_cells=3, n_genes_by_counts<7500(疑似双细胞上限)
#   归一化: total-count 1e4 + log1p；HVG 2000 (seurat_v3, batch_key=sample)
#   降维聚类: PCA 30PC → harmony? 否(内存限制，改用标准PCA+BBKNN? 否) → neighbors(15) → leiden(res=0.5) → UMAP
#   注释marker: 上皮EPCAM/KRT17/KRT5/KRT8; T/NK CD3D/CD2/NKG7; B MS4A1/CD79A; 浆细胞IGHG1/MZB1;
#             髓系LYZ/CD68/AIF1; 肥大细胞TPSAB1/CPA3; 内皮PECAM1/VWF; 成纤维COL1A1/DCN; 平滑肌ACTA2/MYH11
#   用法: python3 06_scRNA_GSE208653.py <workdir> <step>   step∈{load,cluster,umap,annotate}
import sys, os, glob, gzip
import numpy as np
import scanpy as sc
import pandas as pd
import scipy.sparse as sp

wd = sys.argv[1]; step = sys.argv[2]
os.chdir(wd)
sc.settings.verbosity = 1
H5 = "combined.h5ad"

SAMPLES = {
 "GSM6360680_N_HPV_NEG_1": ("NO_HPV", "Normal_HPVneg"),
 "GSM6360681_N_HPV_NEG_2": ("NO_HPV", "Normal_HPVneg"),
 "GSM6360682_N_1":         ("N_HPV",  "Normal_HPVpos"),
 "GSM6360683_N_2":         ("N_HPV",  "Normal_HPVpos"),
 "GSM6360684_HSIL_1":      ("HSIL",   "HSIL"),
 "GSM6360685_HSIL_2":      ("HSIL",   "HSIL"),
 "GSM6360686_SCC_4":       ("CA",     "Cancer"),
 "GSM6360687_SCC_5":       ("CA",     "Cancer"),
 "GSM6360688_ADC_6":       ("CA",     "Cancer"),
}

def load_one(prefix):
    from scipy.io import mmread
    m = mmread(f"{prefix}.matrix.mtx.gz").tocsr().T.astype(np.float32)  # cells x genes
    feats = pd.read_csv(f"{prefix}.features.tsv.gz", sep="\t", header=None)
    bcs = pd.read_csv(f"{prefix}.barcodes.tsv.gz", sep="\t", header=None)[0].values
    ad = sc.AnnData(m)
    ad.var_names = feats[1].values  # gene symbols
    ad.obs_names = [f"{prefix}:{b}" for b in bcs]
    ad.var_names_make_unique()
    sc.pp.filter_cells(ad, min_genes=300)
    ad = ad[ad.obs.n_genes < 7500].copy()
    ad.obs["sample"] = prefix
    ad.obs["stage"] = SAMPLES[prefix][1]
    return ad

if step.startswith("s"):  # s0..s8 逐样本
    i = int(step[1:]); k = list(SAMPLES)[i]
    out = f"s_{i}.h5ad"
    if os.path.exists(out): print("SKIP", out)
    else:
        ad = load_one(k); ad.write(out)
        print("SAMPLE_DONE", k, "cells:", ad.n_obs)

elif step == "merge":
    import anndata
    ads = [sc.read_h5ad(f"s_{i}.h5ad") for i in range(9)]
    comb = anndata.concat(ads, join="outer")
    comb.X = sp.csr_matrix(comb.X)
    sc.pp.filter_genes(comb, min_cells=3)
    comb.write(H5)
    print("LOAD_DONE cells:", comb.n_obs, "genes:", comb.n_vars)

elif step == "norm":
    # 原始counts保留于s_*.h5ad（磁盘），combined转为归一化log值
    ad = sc.read_h5ad(H5)
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, n_top_genes=2000, batch_key="sample")
    ad.write(H5)
    print("NORM_DONE")

elif step == "pca":
    ad = sc.read_h5ad(H5)
    adh = ad[:, ad.var.highly_variable].copy()
    sc.pp.scale(adh, max_value=10)
    sc.tl.pca(adh, n_comps=30, svd_solver="arpack")
    np.save("pca.npy", adh.obsm["X_pca"])
    print("PCA_DONE")

elif step == "cluster":
    ad = sc.read_h5ad(H5)
    ad.obsm["X_pca"] = np.load("pca.npy")
    sc.pp.neighbors(ad, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(ad, resolution=0.5, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
    ad.write(H5)
    print("CLUSTER_DONE nclusters:", ad.obs.leiden.nunique())

elif step == "umap":
    ad = sc.read_h5ad(H5)
    sc.tl.umap(ad)
    ad.write(H5)
    print("UMAP_DONE")

elif step == "annotate":
    ad = sc.read_h5ad(H5)
    markers = {
      "Epithelial": ["EPCAM","KRT17","KRT5","KRT8","KRT18"],
      "T_NK": ["CD3D","CD3E","CD2","NKG7","TRAC"],
      "B": ["MS4A1","CD79A","CD19"],
      "Plasma": ["MZB1","IGHG1","JCHAIN"],
      "Myeloid": ["LYZ","CD68","AIF1","CD14"],
      "Mast": ["TPSAB1","TPSB2","CPA3"],
      "Endothelial": ["PECAM1","VWF","CLDN5"],
      "Fibroblast": ["COL1A1","COL1A2","DCN","LUM"],
      "SmoothMuscle": ["ACTA2","MYH11","TAGLN"],
    }
    raw = ad  # X已是log1p(CP10K)，全基因保留
    # score each cluster by mean marker expression → assign argmax
    cl = ad.obs["leiden"]
    rows = []
    assign = {}
    for c in sorted(cl.unique(), key=int):
        idx = (cl==c).values
        best, bestv = None, -1
        sc_row = {"cluster": c, "n": int(idx.sum())}
        for ct, gs in markers.items():
            gs2 = [g for g in gs if g in raw.var_names]
            v = float(np.asarray(raw[idx, gs2].X.mean()))
            sc_row[ct] = round(v,3)
            if v > bestv: best, bestv = ct, v
        assign[c] = best
        rows.append(sc_row)
    ad.obs["cell_type"] = cl.map(assign).astype(str)
    ad.write(H5)
    pd.DataFrame(rows).to_csv("cluster_marker_scores.csv", index=False)
    # LONP1 stats
    lon = np.asarray(raw[:, "LONP1"].X.todense()).ravel()
    df = pd.DataFrame({"cell_type": ad.obs.cell_type.values, "stage": ad.obs.stage.values,
                       "sample": ad.obs["sample"].values, "LONP1": lon,
                       "umap1": ad.obsm["X_umap"][:,0], "umap2": ad.obsm["X_umap"][:,1],
                       "leiden": ad.obs.leiden.values})
    df.to_csv("LONP1_per_cell.csv.gz", index=False, compression="gzip")
    print("ANNOT_DONE")
    print(pd.crosstab(ad.obs.cell_type, ad.obs.stage))
else:
    print("unknown step")
