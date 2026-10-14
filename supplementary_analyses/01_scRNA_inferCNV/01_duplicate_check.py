# 01_duplicate_check.py —— 复旦妇产来源数据集间的样本查重
# 方法: (1) 每样本 barcode 集合的 Jaccard 重叠；(2) 伪 bulk（全部细胞 counts 求和 → CPM → log1p）
#       在共同基因上的 Pearson 相关，取高变的前 2000 个基因。r>0.98 且 barcode 大量重叠视为同一样本。
# 用法: python -P 01_duplicate_check.py <01_scRNA_inferCNV 目录>
import sys, os, glob, gzip
import numpy as np, pandas as pd
from scipy.io import mmread

base = sys.argv[1]; raw = os.path.join(base, "raw"); out = os.path.join(base, "results"); os.makedirs(out, exist_ok=True)

def load(mtx, bc, feat, col=1):
    m = mmread(mtx).tocsr()  # genes x cells
    genes = pd.read_csv(feat, sep="\t", header=None)[col].values
    bcs = pd.read_csv(bc, sep="\t", header=None)[0].astype(str).str.replace(r"-\d+$", "", regex=True).values
    s = np.asarray(m.sum(1)).ravel()
    pb = pd.Series(s, index=genes).groupby(level=0).sum()
    return pb, set(bcs), m.shape[1]

samples = {}
for d, pat in [("GSE197461", "GSM*_matrix.mtx.gz"), ("GSE208653", "GSM*.matrix.mtx.gz")]:
    for mtx in sorted(glob.glob(os.path.join(raw, d, "extracted", pat))):
        pre = mtx.replace("matrix.mtx.gz", "")
        samples[f"{d}:{os.path.basename(pre).rstrip('._')}"] = (mtx, pre + "barcodes.tsv.gz", pre + "features.tsv.gz")
g = os.path.join(raw, "GSE168652", "extracted")
for t in ("GSM5155196_Tumor/Tumor", "GSM5155197_Normal/Normal"):
    p = os.path.join(g, t)
    samples[f"GSE168652:{t.split('/')[0]}"] = (p + "/matrix.mtx.gz", p + "/barcodes.tsv.gz", p + "/features.tsv.gz")
e = os.path.join(raw, "E-MTAB-11948")
for i in range(1, 7):
    bc = os.path.join(e, f"sample{i}barcodes.tsv.gz")
    if not os.path.exists(bc): print("skip E-MTAB sample", i, "(barcodes missing)"); continue
    samples[f"E-MTAB-11948:sample{i}"] = (os.path.join(e, f"sample{i}matrix.mtx.gz"), bc, os.path.join(e, "sample1features.tsv.gz"))

pbs, bcs, ncell = {}, {}, {}
for k, (m, b, f) in samples.items():
    pbs[k], bcs[k], ncell[k] = load(m, b, f); print("loaded", k, ncell[k])
P = pd.DataFrame(pbs).fillna(0)
P = P.loc[(P > 0).all(1)]
P = np.log1p(P / P.sum(0) * 1e6)
hv = P.var(1).sort_values(ascending=False).index[:2000]
R = np.corrcoef(P.loc[hv].T.values); R = pd.DataFrame(R, index=P.columns, columns=P.columns)
keys = list(samples)
J = pd.DataFrame([[len(bcs[a] & bcs[b]) / len(bcs[a] | bcs[b]) for b in keys] for a in keys], index=keys, columns=keys)
R.round(4).to_csv(os.path.join(out, "dupcheck_pseudobulk_pearson.csv")); J.round(4).to_csv(os.path.join(out, "dupcheck_barcode_jaccard.csv"))
rows = []
for i, a in enumerate(keys):
    for b in keys[i + 1:]:
        if a.split(":")[0] == b.split(":")[0]: continue
        rows.append({"a": a, "b": b, "pearson_hv2000": R.loc[a, b], "barcode_jaccard": J.loc[a, b],
                     "shared_barcodes": len(bcs[a] & bcs[b]), "cells_a": ncell[a], "cells_b": ncell[b]})
res = pd.DataFrame(rows).sort_values("pearson_hv2000", ascending=False)
res.to_csv(os.path.join(out, "dupcheck_cross_dataset_pairs.csv"), index=False)
print(f"common genes used: {P.shape[0]}; HV genes: {len(hv)}")
print(res.head(20).to_string(index=False, float_format=lambda x: f"{x:.4f}"))
print("\nwithin-dataset max off-diagonal pearson (reference):")
for d in sorted(set(k.split(":")[0] for k in keys)):
    ks = [k for k in keys if k.startswith(d)]
    if len(ks) > 1:
        sub = R.loc[ks, ks].values; np.fill_diagonal(sub, np.nan); print(f"  {d}: {np.nanmax(sub):.4f}")
