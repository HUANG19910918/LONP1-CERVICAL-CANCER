# 02c_eqtlgen_full_coloc_inputs.py —— 从 eQTLGen 完整 cis 文件（全部 SNP–基因对）提取 16 个基因，生成 coloc 输入
# Zscore→beta/se 同 02 号脚本（f 取 1000G EUR ALT 频率，N=NrSamples）；协调规则复用 02 号脚本的 align/load_outcome。
# 输出：coloc_input_eqtlgen/<gene>__eQTLGen_full__<outcome>.tsv；同时缓存提取结果 D:\LONP1_supp_data\mr\eqtlgen\eqtlgen_full_16genes.tsv.gz
# 用法: python -P 02c_eqtlgen_full_coloc_inputs.py <03_MR_coloc 目录> <regions 目录> <ref 目录> <eQTLGen 完整文件>
import sys, os
import numpy as np, pandas as pd
base, regd, refd, full = sys.argv[1:5]
src = open(os.path.join(base, "02_harmonize_mr.py"), encoding="utf-8").read()
src = src[:src.index("EQ = None\nif eqtlgen_f")]
sys.argv = [sys.argv[0], base, regd, refd, "-", "LONP1"]
g = {}; exec(compile(src, "02_harmonize_mr.py", "exec"), g)
coords = g["coords"]
cache = os.path.join(os.path.dirname(full), "eqtlgen_full_16genes.tsv.gz")
if os.path.exists(cache):
    E = pd.read_csv(cache, sep="\t")
else:
    want_sym = set(coords.alias_GRCh37) | set(coords.index); want_ens = set(coords.ensg_GRCh38)
    parts = []
    for ch in pd.read_csv(full, sep="\t", chunksize=2_000_000, usecols=["Pvalue", "SNP", "SNPChr", "SNPPos", "AssessedAllele", "OtherAllele",
                                                                     "Zscore", "Gene", "GeneSymbol", "NrSamples"]):
        parts.append(ch[ch.Gene.isin(want_ens) | ch.GeneSymbol.isin(want_sym)])
    E = pd.concat(parts); E.to_csv(cache, sep="\t", index=False)
print("eQTLGen 完整文件提取行数:", len(E), "; 基因:", E.GeneSymbol.value_counts().to_dict())
cid = os.path.join(base, "coloc_input_eqtlgen"); os.makedirs(cid, exist_ok=True)
for gene in coords.index:
    gdir = os.path.join(regd, gene)
    e = E[(E.Gene == coords.loc[gene, "ensg_GRCh38"]) | (E.GeneSymbol == coords.loc[gene, "alias_GRCh37"])].copy()
    if not len(e): print("no eQTLGen rows", gene); continue
    kg = g["kg_table"](gdir)
    e["key"] = e.SNPChr.astype(str) + ":" + e.SNPPos.astype(str)
    e = e.join(kg[["eur_af", "kalt", "kref"]], on="key", how="inner")
    f = np.where(e.AssessedAllele == e.kalt, e.eur_af, 1 - e.eur_af); n = e.NrSamples.values
    den = np.sqrt(2 * f * (1 - f) * (n + e.Zscore ** 2))
    ex = pd.DataFrame({"key": e.key.values, "rsid": e.SNP.values, "ea": e.AssessedAllele.values, "oa": e.OtherAllele.values,
                       "beta": e.Zscore.values / den, "se": 1 / den, "p": e.Pvalue.values, "n": n}).replace([np.inf, -np.inf], np.nan).dropna()
    exa = g["align"](ex, kg, "EUR")
    for on in g["OUT"]:
        o = g["load_outcome"](on, gdir)
        if o is None or not len(o): continue
        oa = g["align"](o, kg, g["OUT"][on][4])
        m = exa[["key", "rsid", "beta", "se", "p", "n", "eur_af", "eas_af"]].merge(oa[["key", "beta", "se", "p", "af_out"]], on="key", suffixes=("_x", "_y"))
        m.rename(columns={"beta_x": "bx", "se_x": "sx", "p_x": "px", "beta_y": "by", "se_y": "sy", "p_y": "py"}).to_csv(
            os.path.join(cid, f"{gene}__eQTLGen_full__{on}.tsv"), sep="\t", index=False)
    print(gene, "eQTLGen SNPs aligned:", len(exa))
