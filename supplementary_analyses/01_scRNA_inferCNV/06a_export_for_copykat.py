# 06a_export_for_copykat.py —— 为 CopyKAT 交叉验证导出逐样本 counts（10x mtx）
# 每样本取全部上皮细胞（>3000 时随机抽 3000，seed=0）+ 参考细胞（免疫+间质，>800 时随机抽 800）；
# 只导出已有 inferCNV 结果的样本（癌、正常、癌旁、HSIL/SIL），参考细胞名单另存 ref_cells.txt。
# 用法: python -P 06a_export_for_copykat.py <qc 目录> <cnv 目录> <输出目录> [数据集逗号列表]
import sys, os, glob, gzip
import numpy as np, pandas as pd, scanpy as sc
from scipy.io import mmwrite
qcd, cnvd, outd = sys.argv[1:4]; only = set(sys.argv[4].split(",")) if len(sys.argv) > 4 else None
rng = np.random.default_rng(0)
for f in sorted(glob.glob(os.path.join(cnvd, "*.csv.gz"))):
    key = os.path.basename(f)[:-7]
    if only and key.split("__")[0] not in only: continue
    od = os.path.join(outd, key)
    if os.path.exists(os.path.join(od, "ref_cells.txt")): continue
    d = pd.read_csv(f, usecols=["cell", "cell_type", "is_ref"])
    epi = d.cell[d.cell_type == "Epithelial"].values; ref = d.cell[d.is_ref].values
    if len(epi) > 3000: epi = rng.choice(epi, 3000, replace=False)
    if len(ref) > 800: ref = rng.choice(ref, 800, replace=False)
    a = sc.read_h5ad(os.path.join(qcd, key + ".h5ad"))
    a = a[np.concatenate([epi, ref])]
    os.makedirs(od, exist_ok=True)
    mmwrite(os.path.join(od, "matrix.mtx"), a.X.T.tocsr().astype(np.int32))
    with open(os.path.join(od, "matrix.mtx"), "rb") as i, gzip.open(os.path.join(od, "matrix.mtx.gz"), "wb") as o: o.write(i.read())
    os.remove(os.path.join(od, "matrix.mtx"))
    pd.Series(a.var_names).to_csv(os.path.join(od, "genes.tsv.gz"), index=False, header=False)
    pd.Series(a.obs_names).to_csv(os.path.join(od, "barcodes.tsv.gz"), index=False, header=False)
    pd.Series(ref).to_csv(os.path.join(od, "ref_cells.txt"), index=False, header=False)
    print("exported", key, len(epi), len(ref), flush=True)
