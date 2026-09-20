# 30_regen_figS1C.py —— 重新生成 Figure S1C，修正标题中的统计量名称
# 原标题误称 "partial Spearman rho"；实际算法为对原值回归去除纯度后取残差再算 Spearman，
# 并非先取秩的标准偏 Spearman，故改为 "Spearman rho of purity residuals"。
# 输入: raw_data/42_purity_merged.csv
# 输出: figures/【1】图片源文件/(1) 子图原文件/FigS1C_CNA_purity_adjusted.png|pdf
# 用法: python3 30_regen_figS1C.py <分析根目录>
import sys, os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr
from plotnine import *
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import setup_font, TH, NPG

base = sys.argv[1]
rd = os.path.join(base, "raw_data")
fig = os.path.join(base, "figures", "【1】图片源文件", "(1) 子图原文件")
os.makedirs(fig, exist_ok=True)
FAM = setup_font(base); T = TH(FAM)

m = pd.read_csv(os.path.join(rd, "42_purity_merged.csv"))
mc = m.dropna(subset=["log2CNA", "LONP1_log2", "purity"]).copy()
def resid(y, X):
    X = sm.add_constant(X); return y - sm.OLS(y, X).fit().predict(X)
mc = mc.assign(rx=resid(mc.log2CNA, mc[["purity"]]), ry=resid(mc.LONP1_log2, mc[["purity"]]))
r, p = spearmanr(mc.rx, mc.ry)
print(f"n={len(mc)}  Spearman rho of purity residuals = {r:.4f}, p = {p:.3e}")

pC = (ggplot(mc, aes("rx", "ry")) + geom_point(size=1.4, alpha=0.55, color=NPG["navy"])
      + geom_smooth(method="lm", color=NPG["red"], fill=NPG["salmon"], size=0.8)
      + labs(x="LONP1 log2 CNA ratio (purity-adj. residual)", y="LONP1 mRNA (purity-adj. residual)",
             title=f"Copy number vs expression, purity-adjusted\nSpearman rho of purity residuals = {r:.2f}, p = {p:.1e}")
      + T + theme(plot_title=element_text(size=10.5, weight="bold", ha="center", family=FAM)))
pC.save(os.path.join(fig, "FigS1C_CNA_purity_adjusted.png"), width=3.8, height=3.8, dpi=300, verbose=False)
pC.save(os.path.join(fig, "FigS1C_CNA_purity_adjusted.pdf"), width=3.8, height=3.8, verbose=False)
print("regenerated FigS1C in:", fig)
