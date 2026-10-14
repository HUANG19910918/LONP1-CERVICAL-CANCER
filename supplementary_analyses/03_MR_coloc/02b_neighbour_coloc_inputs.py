# 02b_neighbour_coloc_inputs.py —— LONP1 ±1 Mb 内其他基因（GTEx v8 血液）的 coloc 输入，作为阴性/对照
# 纳入：窗口内 GTEx 血液有任一 SNP P<1e-4 的基因（有 cis-eQTL 信号），与各结局在 LONP1 窗口内配对（同 02 号脚本的协调规则）。
# 用法: python -P 02b_neighbour_coloc_inputs.py <03_MR_coloc 目录> <regions 目录> <ref 目录>
import sys, os, importlib.util
import numpy as np, pandas as pd
base, regd, refd = sys.argv[1:4]
spec = importlib.util.spec_from_file_location("hm", os.path.join(base, "02_harmonize_mr.py"))
src = open(os.path.join(base, "02_harmonize_mr.py"), encoding="utf-8").read()
src = src[:src.index("EQ = None\nif eqtlgen_f")]  # 只取函数定义部分
sys.argv = [sys.argv[0], base, regd, refd, "-", "LONP1"]
g = {}; exec(compile(src, "02_harmonize_mr.py", "exec"), g)
gdir = os.path.join(regd, "LONP1"); kg = g["kg_table"](gdir)
a = pd.read_csv(os.path.join(gdir, "GTEx_blood_allgenes.tsv"), sep="\t", header=0, low_memory=False)
if "molecular_trait_id" not in a.columns: a = pd.read_csv(os.path.join(gdir, "GTEx_blood_allgenes.tsv"), sep="\t", header=None, names=g["GTEX_COLS"], low_memory=False)
a = a[a.type == "SNP"]
sig = a.groupby("molecular_trait_id").pvalue.min(); genes = sig[sig < 1e-4].index.tolist()
print(f"LONP1 窗口内 GTEx 血液基因 {a.molecular_trait_id.nunique()} 个，有 P<1e-4 cis 信号的 {len(genes)} 个")
a["key"] = [g["lift38"](c, x) for c, x in zip(a.chromosome, a.position)]
outs = {k: g["load_outcome"](k, gdir) for k in ["FinnGen_CC", "FinnGen_CIS", "UKB_180.3", "BBJ_CeC"]}
cid = os.path.join(base, "coloc_input_neighbours"); os.makedirs(cid, exist_ok=True)
meta = []
for ensg in genes:
    e = a[a.molecular_trait_id == ensg]
    ex = pd.DataFrame({"key": e.key, "rsid": e.rsid, "ea": e.alt, "oa": e.ref, "beta": e.beta, "se": e.se, "p": e.pvalue, "n": 670}).dropna(subset=["key"])
    exa = g["align"](ex, kg, "EUR")
    for on, o in outs.items():
        oa = g["align"](o, kg, g["OUT"][on][4])
        m = exa[["key", "rsid", "beta", "se", "p", "n", "eur_af", "eas_af"]].merge(oa[["key", "beta", "se", "p", "af_out"]], on="key", suffixes=("_x", "_y"))
        m.rename(columns={"beta_x": "bx", "se_x": "sx", "p_x": "px", "beta_y": "by", "se_y": "sy", "p_y": "py"}).to_csv(
            os.path.join(cid, f"{ensg}__GTEx_blood__{on}.tsv"), sep="\t", index=False)
    meta.append({"ensg": ensg, "min_p": float(sig[ensg]), "n_snp": len(exa)})
pd.DataFrame(meta).to_csv(os.path.join(base, "results", "neighbour_genes_LONP1_window.csv"), index=False)
print("written", len(meta))
