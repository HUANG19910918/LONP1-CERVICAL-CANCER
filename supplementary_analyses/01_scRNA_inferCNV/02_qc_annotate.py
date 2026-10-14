# 02_qc_annotate.py —— 各数据集逐样本统一 QC + 细胞类型注释（补充分析 1，阶段 3a）
# QC 口径沿用 v8 06_scRNA_GSE208653.py：min_genes=300、n_genes<7500、min_cells=3；
#   新增：线粒体 UMI 比例 <20%（有 MT- 基因的数据集；GSE208653/GSE197461 矩阵中作者已移除线粒体基因）；
#         scanpy 内置 scrublet 逐样本去双细胞（threshold=0.25）（predicted_doublet=True 剔除）。
# 注释：逐样本 normalize_total(1e4)+log1p → HVG 2000 → PCA 30 → neighbors 15 → leiden(res=1.0)；
#   簇级 marker 均值取最大者（marker 集同 v8 annotate 步骤）。
# 输出：D:\LONP1_supp_data\work\qc\<dataset>__<sample>.h5ad（X=原始 counts，obs 含 QC 指标与 cell_type）；
#       results/qc_summary.csv、results/celltype_by_sample.csv
# 用法: python -P 02_qc_annotate.py <01_scRNA_inferCNV 目录> <输出 h5ad 目录>
import sys, os, glob
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.io import mmread
import scipy.sparse as sp

base, outh = sys.argv[1], sys.argv[2]
raw = os.path.join(base, "raw"); res = os.path.join(base, "results")
os.makedirs(outh, exist_ok=True); os.makedirs(res, exist_ok=True)
sc.settings.verbosity = 0

MARKERS = {
    "Epithelial": ["EPCAM", "KRT17", "KRT5", "KRT8", "KRT18"],
    "T_NK": ["CD3D", "CD3E", "CD2", "NKG7", "TRAC"],
    "B": ["MS4A1", "CD79A", "CD19"],
    "Plasma": ["MZB1", "IGHG1", "JCHAIN"],
    "Myeloid": ["LYZ", "CD68", "AIF1", "CD14"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "SmoothMuscle": ["ACTA2", "MYH11", "TAGLN"],
}

# 样本表：dataset, sample, group, hpv, mtx, barcodes, features, symbol 列
S = []
G208 = {"N_HPV_NEG_1": ("Normal", "HPV-"), "N_HPV_NEG_2": ("Normal", "HPV-"), "N_1": ("Normal", "HPV+"), "N_2": ("Normal", "HPV+"),
        "HSIL_1": ("HSIL", "HPV+"), "HSIL_2": ("HSIL", "HPV+"), "SCC_4": ("Cancer", "HPV+"), "SCC_5": ("Cancer", "HPV+"), "ADC_6": ("Cancer", "HPV+")}
for m in sorted(glob.glob(os.path.join(raw, "GSE208653", "extracted", "GSM*.matrix.mtx.gz"))):
    p = m.replace("matrix.mtx.gz", ""); nm = os.path.basename(p).rstrip(".").split("_", 1)[1]
    S.append(("GSE208653", os.path.basename(p).rstrip("."), *G208[nm], m, p + "barcodes.tsv.gz", p + "features.tsv.gz", 1))
G197 = {"SCC_1": "HPV58", "SCC_2": "HPV16", "SCC_3": "HPV16", "ADC_1": "HPV16", "ADC_2": "HPV18", "ADC_3": "HPV45", "ADC_4": "HPV-", "ADC_5": "HPV-"}
for m in sorted(glob.glob(os.path.join(raw, "GSE197461", "extracted", "GSM*_matrix.mtx.gz"))):
    p = m.replace("matrix.mtx.gz", ""); nm = os.path.basename(p).rstrip("_").split("_", 1)[1]
    S.append(("GSE197461", os.path.basename(p).rstrip("_"), "Cancer", G197[nm], m, p + "barcodes.tsv.gz", p + "features.tsv.gz", 1))
e = os.path.join(raw, "E-MTAB-11948")
for i in range(1, 7):
    feat = os.path.join(e, "sample1features.tsv.gz" if i == 5 else f"sample{i}features.tsv.gz")  # sample5 features 误传 → 用 sample1（解压内容一致）
    S.append(("E-MTAB-11948", f"sample{i}_patient{(i - 1) % 3 + 1}", "Cancer" if i <= 3 else "AdjacentNormal", "NA",
              os.path.join(e, f"sample{i}matrix.mtx.gz"), os.path.join(e, f"sample{i}barcodes.tsv.gz"), feat, 1))
for m in sorted(glob.glob(os.path.join(raw, "GSE236738", "extracted", "GSM757477[79]_*.matrix.mtx.gz")) +
               glob.glob(os.path.join(raw, "GSE236738", "extracted", "GSM7574781_*.matrix.mtx.gz"))):  # 只用治疗前 TN
    p = m.replace("matrix.mtx.gz", "")
    if all(os.path.exists(p + x) for x in ("barcodes.tsv.gz", "features.tsv.gz")):
        S.append(("GSE236738", os.path.basename(p).rstrip("."), "Cancer", "NA", m, p + "barcodes.tsv.gz", p + "features.tsv.gz", 1))
for m in sorted(glob.glob(os.path.join(raw, "GSE173231", "extracted", "GSM*_matrix.mtx.gz"))):
    p = m.replace("matrix.mtx.gz", "")
    S.append(("GSE173231", os.path.basename(p).rstrip("_"), "NormalBiopsy", "NA", m, p + "barcodes.tsv.gz", p + "genes.tsv.gz", 1))
# GSE333892：整体 rds 的 counts 层为 log1p(CP10K)，由 00b_GSE333892_recover_counts.R 用 expm1(data)×nCount_RNA/1e4 还原为整数 counts，
# 按 isolate 逐样本导出（作者已做 min.features=200、线粒体过滤、DoubletDetection 去双细胞）。
G333 = {"CC": "Cancer", "CON": "Normal", "SIL": "SIL"}
d333 = os.environ.get("GSE333892_DIR", r"D:\LONP1_supp_data\GSE333892\per_sample")
for p in sorted(glob.glob(os.path.join(d333, "*", "matrix.mtx.gz"))):
    dd = os.path.dirname(p); iso = os.path.basename(dd)
    if os.path.exists(os.path.join(dd, "metadata.csv")):
        S.append(("GSE333892", iso, G333[iso.split("_")[0]], "NA", p, os.path.join(dd, "barcodes.tsv.gz"), os.path.join(dd, "features.tsv.gz"), 0))

only = set(sys.argv[3].split(",")) if len(sys.argv) > 3 else None
qc_rows, ct_rows = [], []
for ds, sm, grp, hpv, m, b, f, col in S:
    if only and ds not in only: continue
    out = os.path.join(outh, f"{ds}__{sm}.h5ad")
    if os.path.exists(out):
        a = sc.read_h5ad(out)
    else:
        X = mmread(m).tocsr().T.astype(np.float32)
        genes = pd.read_csv(f, sep="\t", header=None)[col].astype(str).values
        bcs = pd.read_csv(b, sep="\t", header=None)[0].astype(str).values
        a = ad.AnnData(sp.csr_matrix(X)); a.var_names = genes; a.var_names_make_unique()
        a.obs_names = [f"{ds}:{sm}:{x}" for x in bcs]
        n0 = a.n_obs
        a.var["mt"] = a.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
        has_mt = bool(a.var.mt.sum() > 0)
        keep = (a.obs.n_genes_by_counts >= 300) & (a.obs.n_genes_by_counts < 7500)
        if has_mt: keep &= a.obs.pct_counts_mt < 20
        a = a[keep.values].copy()
        n1 = a.n_obs
        try:
            sc.pp.scrublet(a, threshold=0.25, random_state=0)
            dbl = int(a.obs.predicted_doublet.sum()); a = a[~a.obs.predicted_doublet.values].copy()
        except Exception as ex:
            print("scrublet failed", sm, ex); dbl = -1
        sc.pp.filter_genes(a, min_cells=3)
        a.obs["dataset"], a.obs["sample"], a.obs["group"], a.obs["hpv"] = ds, sm, grp, hpv
        # 注释
        w = a.copy(); sc.pp.normalize_total(w, target_sum=1e4); sc.pp.log1p(w)
        sc.pp.highly_variable_genes(w, n_top_genes=2000); h = w[:, w.var.highly_variable].copy()
        sc.pp.scale(h, max_value=10); sc.tl.pca(h, n_comps=30, random_state=0)
        sc.pp.neighbors(h, n_neighbors=15, random_state=0)
        sc.tl.leiden(h, resolution=1.0, key_added="leiden", flavor="igraph", n_iterations=2, directed=False, random_state=0)
        cl = h.obs.leiden.values; a.obs["leiden"] = cl
        assign = {}
        for c in np.unique(cl):
            idx = cl == c; best, bv = None, -1
            for ct, gs in MARKERS.items():
                gs2 = [g for g in gs if g in w.var_names]
                v = float(np.asarray(w[idx, gs2].X.mean())) if gs2 else -1
                if v > bv: best, bv = ct, v
            assign[c] = best
        a.obs["cell_type"] = pd.Categorical([assign[c] for c in cl])
        a.uns["qc"] = {"n_raw": n0, "n_after_filter": n1, "n_doublet_removed": dbl, "has_mt": has_mt}
        a.write(out, compression="gzip")
    q = a.uns["qc"]
    qc_rows.append({"dataset": ds, "sample": sm, "group": grp, "hpv": hpv, "cells_raw": q["n_raw"], "after_filter": q["n_after_filter"],
                    "doublets_removed": q["n_doublet_removed"], "cells_final": a.n_obs, "genes": a.n_vars, "mt_filter": q["has_mt"],
                    "median_genes": float(a.obs.n_genes_by_counts.median()), "LONP1_detected": "LONP1" in a.var_names})
    vc = a.obs.cell_type.value_counts()
    ct_rows.append({"dataset": ds, "sample": sm, "group": grp, **vc.to_dict()})
    print(f"{ds:13s} {sm:28s} {grp:14s} raw={q['n_raw']:6d} final={a.n_obs:6d} epi={vc.get('Epithelial', 0):6d}", flush=True)
pd.DataFrame(qc_rows).to_csv(os.path.join(res, "qc_summary.csv"), index=False)
pd.DataFrame(ct_rows).fillna(0).to_csv(os.path.join(res, "celltype_by_sample.csv"), index=False)
print("DONE")
