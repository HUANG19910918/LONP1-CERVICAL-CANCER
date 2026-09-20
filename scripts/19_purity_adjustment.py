# 19_purity_adjustment.py —— TCGA-CESC 肿瘤纯度（ABSOLUTE）与 LONP1：组织构成效应的纯度层面检验
# 复现口径：
#   纯度: TCGA PanCanAtlas ABSOLUTE 结果 TCGA_mastercalls.abs_tables_JSedit.fixed.txt
#         (https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5, GDC PanCanAtlas publications page)，取 call status=="called"
#   表达: raw_data/08 (RSEM→log2(RSEM+1))；CNA: raw_data/07；上皮分数: raw_data/28 (TCGA样本, Xena)；免疫签名: raw_data/08 成员基因
#   统计: Spearman；纯度校正 = 对两变量分别用 OLS(原值) 去除纯度后取残差再算 Spearman
#         (注意: 这不是先取秩的标准偏 Spearman, 正文与图注按实际算法表述)；OLS 含纯度协变量
# 输出: raw_data/42_purity_merged.csv, 43_purity_stats.txt；figures/FigS1_purity.png/pdf
# 用法: python3 19_purity_adjustment.py <分析根目录> <ABSOLUTE文件>
import sys, os
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
import statsmodels.api as sm
from plotnine import *
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import setup_font, TH, NPG
base, absf = sys.argv[1], sys.argv[2]
rd = os.path.join(base, "raw_data"); fig = os.path.join(base, "figures", "子图原文件"); os.makedirs(fig, exist_ok=True)
FAM = setup_font(base); T = TH(FAM)

ab = pd.read_csv(absf, sep="\t"); ab = ab[ab["call status"] == "called"][["array", "purity", "ploidy"]].rename(columns={"array": "sample"})
d8 = pd.read_csv(os.path.join(rd, "08_CESC_gene_panel_expression_RSEM.csv"))
d7 = pd.read_csv(os.path.join(rd, "07_CESC_LONP1_CNA.csv"))
d28 = pd.read_csv(os.path.join(rd, "28_TCGA_GTEx_epithelial_markers.csv")); d28 = d28[d28.source == "TCGA"][["sample_id", "epi_score"]].rename(columns={"sample_id": "sample"})
z = lambda v: (v - v.mean()) / v.std()
lg = lambda g: np.log2(d8[g] + 1)
m = pd.DataFrame({"sample": d8["sample"], "LONP1_log2": lg("LONP1"),
                  "M1_sig": np.mean([z(lg(g)) for g in ["NOS2", "IRF5", "PTGS2"]], axis=0),
                  "DC_sig": np.mean([z(lg(g)) for g in ["ITGAX", "CD1C", "NRP1"]], axis=0),
                  "CD4_sig": np.mean([z(lg(g)) for g in ["CD4", "IL7R"]], axis=0)})
m = m.merge(ab, on="sample", how="inner").merge(d7, on="sample", how="left").merge(d28, on="sample", how="left")
m.to_csv(os.path.join(rd, "42_purity_merged.csv"), index=False)
L = []
def P(s): print(s); L.append(s)
P(f"CESC samples with ABSOLUTE purity and LONP1 expression: n={len(m)}; purity median={m.purity.median():.2f} IQR={m.purity.quantile(.25):.2f}-{m.purity.quantile(.75):.2f}")
r, p = spearmanr(m.purity, m.LONP1_log2); P(f"LONP1 ~ purity: Spearman rho={r:.3f} p={p:.2e}")
mm = m.dropna(subset=["epi_score"]); r, p = spearmanr(mm.purity, mm.epi_score); P(f"epithelial marker score ~ purity: rho={r:.3f} p={p:.2e} (n={len(mm)})")
# CNA-expression with purity adjustment
mc = m.dropna(subset=["log2CNA"])
r0, p0 = spearmanr(mc.log2CNA, mc.LONP1_log2)
def resid(y, X):
    X = sm.add_constant(X); return y - sm.OLS(y, X).fit().predict(X)
rr, pr = spearmanr(resid(mc.LONP1_log2, mc[["purity"]]), resid(mc.log2CNA, mc[["purity"]]))
P(f"LONP1 ~ log2CNA: raw rho={r0:.3f} p={p0:.2e}; purity-adjusted partial rho={rr:.3f} p={pr:.2e} (n={len(mc)})")
ols = sm.OLS(mc.LONP1_log2, sm.add_constant(mc[["log2CNA", "purity"]])).fit()
P(f"OLS LONP1 ~ log2CNA + purity: beta_CNA={ols.params['log2CNA']:.3f} (p={ols.pvalues['log2CNA']:.2e}), beta_purity={ols.params['purity']:.3f} (p={ols.pvalues['purity']:.2e}), R2={ols.rsquared:.3f}")
# immune signatures with purity adjustment
for sig in ["M1_sig", "DC_sig", "CD4_sig"]:
    r0, p0 = spearmanr(m.LONP1_log2, m[sig]); rs, ps = spearmanr(m.purity, m[sig])
    rr, pr = spearmanr(resid(m.LONP1_log2, m[["purity"]]), resid(m[sig], m[["purity"]]))
    P(f"{sig}: LONP1 raw rho={r0:.3f} (p={p0:.2e}); sig~purity rho={rs:.3f}; purity-adjusted partial rho={rr:.3f} (p={pr:.2e})")
# purity tertiles
m["ptert"] = pd.qcut(m.purity, 3, labels=["Low", "Mid", "High"])
P("LONP1 by purity tertile (mean log2): " + ", ".join(f"{k}={v:.2f}" for k, v in m.groupby("ptert", observed=True).LONP1_log2.mean().items()))
open(os.path.join(rd, "43_purity_stats.txt"), "w").write("\n".join(L) + "\n")

# Figure S1: A purity vs LONP1; B purity vs epithelial score; C CNA vs LONP1 residuals (purity-adjusted)
r1, p1 = spearmanr(m.purity, m.LONP1_log2)
pA = (ggplot(m, aes("purity", "LONP1_log2")) + geom_point(size=1.4, alpha=0.55, color=NPG["navy"])
      + geom_smooth(method="lm", color=NPG["red"], fill=NPG["salmon"], size=0.8)
      + labs(x="Tumor purity (ABSOLUTE)", y="LONP1 mRNA, log2(RSEM+1)", title=f"LONP1 vs tumor purity (n={len(m)})\nSpearman rho = {r1:.2f}, p = {p1:.1e}")
      + T + theme(plot_title=element_text(size=10.5, weight="bold", ha="center", family=FAM)))
pA.save(os.path.join(fig, "FigS1A_purity_LONP1.png"), width=3.8, height=3.8, dpi=300, verbose=False)
pA.save(os.path.join(fig, "FigS1A_purity_LONP1.pdf"), width=3.8, height=3.8, verbose=False)
r2, p2 = spearmanr(mm.purity, mm.epi_score)
pB = (ggplot(mm, aes("purity", "epi_score")) + geom_point(size=1.4, alpha=0.55, color=NPG["navy"])
      + geom_smooth(method="lm", color=NPG["red"], fill=NPG["salmon"], size=0.8)
      + labs(x="Tumor purity (ABSOLUTE)", y="Epithelial marker score (z)", title=f"Epithelial score vs purity (n={len(mm)})\nSpearman rho = {r2:.2f}, p = {p2:.1e}")
      + T + theme(plot_title=element_text(size=10.5, weight="bold", ha="center", family=FAM)))
pB.save(os.path.join(fig, "FigS1B_purity_episcore.png"), width=3.8, height=3.8, dpi=300, verbose=False)
pB.save(os.path.join(fig, "FigS1B_purity_episcore.pdf"), width=3.8, height=3.8, verbose=False)
mc = mc.assign(rx=resid(mc.log2CNA, mc[["purity"]]), ry=resid(mc.LONP1_log2, mc[["purity"]]))
pC = (ggplot(mc, aes("rx", "ry")) + geom_point(size=1.4, alpha=0.55, color=NPG["navy"])
      + geom_smooth(method="lm", color=NPG["red"], fill=NPG["salmon"], size=0.8)
      + labs(x="LONP1 log2 CNA ratio (purity-adj. residual)", y="LONP1 mRNA (purity-adj. residual)",
             title=f"Copy number vs expression, purity-adjusted\nSpearman rho of purity residuals = {spearmanr(mc.rx, mc.ry)[0]:.2f}, p = {spearmanr(mc.rx, mc.ry)[1]:.1e}")
      + T + theme(plot_title=element_text(size=10.5, weight="bold", ha="center", family=FAM)))
pC.save(os.path.join(fig, "FigS1C_CNA_purity_adjusted.png"), width=3.8, height=3.8, dpi=300, verbose=False)
pC.save(os.path.join(fig, "FigS1C_CNA_purity_adjusted.pdf"), width=3.8, height=3.8, verbose=False)
print("DONE19")
