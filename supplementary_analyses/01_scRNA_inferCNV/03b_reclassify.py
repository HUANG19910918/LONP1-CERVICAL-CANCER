# 03b_reclassify.py —— inferCNV 恶性判定阈值的主分析/敏感性分析（不重跑 inferCNV，只改写判定列）
# 背景：03 号脚本原判定为双阈值 >参考细胞 99 分位，癌样本中"未定"细胞过多（如 E-MTAB sample3 7,112/10,853）。
# 阈值比较（results/infercnv_threshold_comparison.csv）：正常组织上皮的"恶性"比例在 q95 下为 1.5–9.2%，q99 下为 0.3–2.8%。
# 主判定（epi_class）：cnv_score 与 cnv_corr 均 > 本样本参考细胞 95 分位 → malignant；均 <= 95 分位 → nonmalignant；其余 unresolved。
# 敏感性（epi_class_q99）：同上，阈值为 99 分位（即 03 号脚本原判定）。
# 用法: python -P 03b_reclassify.py <cnv 目录> <results 目录>
import sys, os, glob
import numpy as np, pandas as pd
cnvd, res = sys.argv[1:3]
def classify(d, q):
    r = d[d.is_ref]; ts, tc = r.cnv_score.quantile(q), r.cnv_corr.quantile(q)
    e = d.cell_type.eq("Epithelial")
    return np.where(~e, "non_epithelial", np.where((d.cnv_score > ts) & (d.cnv_corr > tc), "malignant",
                    np.where((d.cnv_score <= ts) & (d.cnv_corr <= tc), "nonmalignant", "unresolved")))
rows, cmp = [], []
for f in sorted(glob.glob(os.path.join(cnvd, "*.csv.gz"))):
    d = pd.read_csv(f)
    if "epi_class_q99" not in d.columns:
        d["epi_class_q99"] = classify(d, 0.99)
    d["epi_class"] = classify(d, 0.95)
    d.to_csv(f, index=False)
    e = d[d.cell_type == "Epithelial"]
    for col in ("epi_class", "epi_class_q99"):
        vc = e[col].value_counts(); m = e[e[col] == "malignant"]
        rows.append({"dataset": d.dataset.iloc[0], "sample": d["sample"].iloc[0], "group": d.group.iloc[0], "rule": "q95" if col == "epi_class" else "q99",
                     "n_epi": len(e), "malignant": int(vc.get("malignant", 0)), "nonmalignant": int(vc.get("nonmalignant", 0)),
                     "unresolved": int(vc.get("unresolved", 0)), "pct_malignant": round(100 * vc.get("malignant", 0) / max(len(e), 1), 1),
                     "pct_gain19p13_3_in_malignant": round(100 * m.gain19p13_3.mean(), 1) if len(m) else np.nan})
S = pd.DataFrame(rows); S.to_csv(os.path.join(res, "infercnv_summary.csv"), index=False)
p = S.pivot_table(index=["dataset", "sample", "group", "n_epi"], columns="rule", values="pct_malignant").reset_index()
p.to_csv(os.path.join(res, "infercnv_threshold_comparison.csv"), index=False)
print(S[S.rule == "q95"].drop(columns="rule").to_string(index=False))
print("\n正常组织上皮被判恶性的比例（q95 / q99）:")
print(p[p.group.isin(["Normal", "AdjacentNormal"])].to_string(index=False))
