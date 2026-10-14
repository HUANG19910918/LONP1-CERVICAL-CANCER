# 03_infercnv.py —— 逐样本 inferCNV（infercnvpy），划分恶性/非恶性上皮，推断 19p13.3 区段 CNV（补充分析 1，阶段 4）
# 参考细胞：同一样本内的免疫 + 间质细胞（T_NK、B、Plasma、Myeloid、Mast、Endothelial、Fibroblast、SmoothMuscle）。
# 基因坐标：v8 raw_data/30_coexpression_LONP1_CESC_genomewide.csv（Ensembl GRCh38 + UCSC hg38 细胞带）。
# infercnvpy.tl.infercnv: window_size=250, step=10, exclude_chromosomes=X/Y, calculate_gene_values=True；输入 log1p(CP10K)。
# 每个细胞的指标：
#   cnv_score = 全部窗口 CNV 值的均方；cnv_corr = 与本样本 CNV score 最高 5% 上皮细胞平均谱的 Pearson 相关；
#   seg_19p13_3 = 19p13.3 上基因的 gene-level CNV 均值；arm_19p = 19p 全部基因均值；lonp1_cnv = LONP1 基因 CNV 值。
# 恶性判定（Tirosh/Puram 式双阈值，阈值来自本样本参考细胞）：
#   上皮细胞 cnv_score > 参考细胞 99 分位 且 cnv_corr > 参考细胞 99 分位 → malignant；
#   两项都 <= 95 分位 → nonmalignant；其余 → unresolved（不进入比较）。
#   正常组织样本（GSE208653 正常、E-MTAB 癌旁）同样计算，作为"正常上皮中被误判为恶性"的比例校准。
# 19p13.3 增益：seg_19p13_3 > 参考细胞 99 分位 → gain19p13.3。
# 输出：D:\LONP1_supp_data\work\cnv\<样本>.csv.gz（每细胞指标）、results/infercnv_summary.csv、figures/infercnv_<样本>.png
# 用法: python -P 03_infercnv.py <01_scRNA_inferCNV 目录> <qc h5ad 目录> <cnv 输出目录> <v8 干实验分析目录>
import sys, os, glob, warnings
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc
import infercnvpy as cnv
warnings.filterwarnings("ignore")
def main():
    base, qcd, outd, v8 = sys.argv[1:5]
    res = os.path.join(base, "results"); figd = os.path.join(base, "figures"); os.makedirs(outd, exist_ok=True); os.makedirs(figd, exist_ok=True)
    only = set(sys.argv[5].split(",")) if len(sys.argv) > 5 else None
    REF = ["T_NK", "B", "Plasma", "Myeloid", "Mast", "Endothelial", "Fibroblast", "SmoothMuscle"]

    pos = pd.read_csv(os.path.join(v8, "raw_data", "30_coexpression_LONP1_CESC_genomewide.csv"))
    pos = pos.dropna(subset=["chr", "start", "end"]).drop_duplicates("gene").set_index("gene")
    # v8 共表达表不含 LONP1 自身 → 手工补 GRCh38 坐标（Ensembl ENSG00000196365，chr19:5,691,834-5,720,463，19p13.3）
    pos.loc["LONP1", ["chr", "start", "end", "cytoband"]] = ["19", 5691834, 5720463, "19p13.3"]
    pos = pos[pos.chr.astype(str).isin([str(i) for i in range(1, 23)])]
    pos["chromosome"] = "chr" + pos.chr.astype(str)
    b19p133 = set(pos.index[pos.cytoband.eq("19p13.3")]); b19p = set(pos.index[pos.cytoband.astype(str).str.startswith("19p")])

    rows = []
    for f in sorted(glob.glob(os.path.join(qcd, "*.h5ad"))):
        key = os.path.basename(f)[:-5]
        if only and key.split("__")[0] not in only: continue
        out = os.path.join(outd, key + ".csv.gz")
        a = sc.read_h5ad(f)
        grp = a.obs.group.iloc[0]
        n_epi = int((a.obs.cell_type == "Epithelial").sum()); n_ref = int(a.obs.cell_type.isin(REF).sum())
        if grp == "NormalBiopsy" or n_epi < 20 or n_ref < 50:
            print(f"SKIP {key}: group={grp} epi={n_epi} ref={n_ref}"); continue
        if os.path.exists(out):
            d = pd.read_csv(out)
        else:
            a.var["chromosome"] = pos.chromosome.reindex(a.var_names).values
            a.var["start"] = pos.start.reindex(a.var_names).values
            a.var["end"] = pos.end.reindex(a.var_names).values
            a = a[:, a.var.chromosome.notna()].copy()
            a.layers["counts"] = a.X.copy()
            sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
            a.obs["cnv_ref"] = np.where(a.obs.cell_type.isin(REF), "ref", a.obs.cell_type.astype(str))
            cnv.tl.infercnv(a, reference_key="cnv_ref", reference_cat=["ref"], window_size=250, step=10,
                            exclude_chromosomes=None, calculate_gene_values=True, n_jobs=8, chunksize=1000)
            X = a.obsm["X_cnv"]; X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
            score = (X ** 2).mean(1)
            epi = (a.obs.cell_type == "Epithelial").values; isref = (a.obs.cnv_ref == "ref").values
            top = epi & (score >= np.quantile(score[epi], 0.95))
            prof = X[top].mean(0)
            Xc = X - X.mean(1, keepdims=True); pc = prof - prof.mean()
            corr = (Xc @ pc) / (np.linalg.norm(Xc, axis=1) * np.linalg.norm(pc) + 1e-12)
            G = a.layers["gene_values_cnv"]; G = G.toarray() if hasattr(G, "toarray") else np.asarray(G)
            gi = {g: i for i, g in enumerate(a.var_names)}
            def segmean(gs):
                idx = [gi[g] for g in gs if g in gi]
                return np.nanmean(G[:, idx], axis=1) if idx else np.full(a.n_obs, np.nan)
            d = pd.DataFrame({"cell": a.obs_names, "dataset": a.obs.dataset.values, "sample": a.obs["sample"].values,
                              "group": a.obs.group.values, "hpv": a.obs.hpv.values, "cell_type": a.obs.cell_type.astype(str).values,
                              "is_ref": isref, "cnv_score": score, "cnv_corr": corr,
                              "seg_19p13_3": segmean(b19p133), "arm_19p": segmean(b19p),
                              "lonp1_cnv": G[:, gi["LONP1"]] if "LONP1" in gi else np.nan,
                              "n_genes_19p13_3": sum(g in gi for g in b19p133)})
            r = d[d.is_ref]
            t99s, t95s = r.cnv_score.quantile(0.99), r.cnv_score.quantile(0.95)
            t99c, t95c = r.cnv_corr.quantile(0.99), r.cnv_corr.quantile(0.95)
            e = d.cell_type.eq("Epithelial")
            d["epi_class"] = np.where(~e, "non_epithelial",
                             np.where((d.cnv_score > t99s) & (d.cnv_corr > t99c), "malignant",
                             np.where((d.cnv_score <= t95s) & (d.cnv_corr <= t95c), "nonmalignant", "unresolved")))
            d["gain19p13_3"] = d.seg_19p13_3 > r.seg_19p13_3.quantile(0.99)
            d["loss19p13_3"] = d.seg_19p13_3 < r.seg_19p13_3.quantile(0.01)
            d.attrs = {}
            d.to_csv(out, index=False)
            try:
                cnv.pl.chromosome_heatmap(a[np.isin(a.obs.cnv_ref, ["ref", "Epithelial"])], groupby="cnv_ref", show=False,
                                          save=False, figsize=(12, 5))
                import matplotlib.pyplot as plt
                plt.savefig(os.path.join(figd, f"infercnv_{key}.png"), dpi=110, bbox_inches="tight"); plt.close("all")
            except Exception as ex:
                print("plot failed", key, ex)
        e = d[d.cell_type == "Epithelial"]
        vc = e.epi_class.value_counts()
        rows.append({"dataset": d.dataset.iloc[0], "sample": d["sample"].iloc[0], "group": d.group.iloc[0], "hpv": d.hpv.iloc[0],
                     "n_ref": int(d.is_ref.sum()), "n_epi": len(e), "malignant": int(vc.get("malignant", 0)),
                     "nonmalignant": int(vc.get("nonmalignant", 0)), "unresolved": int(vc.get("unresolved", 0)),
                     "pct_malignant": round(100 * vc.get("malignant", 0) / max(len(e), 1), 1),
                     "pct_gain19p13_3_in_malignant": round(100 * e[e.epi_class == "malignant"].gain19p13_3.mean(), 1) if vc.get("malignant", 0) else np.nan,
                     "median_seg19p13_3_malignant": e[e.epi_class == "malignant"].seg_19p13_3.median() if vc.get("malignant", 0) else np.nan,
                     "median_seg19p13_3_ref": d[d.is_ref].seg_19p13_3.median()})
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(res, "infercnv_summary.csv"), index=False)
    print("DONE")


if __name__ == "__main__":
    main()
