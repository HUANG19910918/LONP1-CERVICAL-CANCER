# 12_upgrade_panels.py —— 图型升级：气泡富集图(3B/3C/3D)、marker气泡矩阵(6F)、
#                         边际分布散点(2A)、raincloud(1A)
# 用法: python3 12...py <目录> <step>  step∈{enrichdot, markerdot, marginal, raincloud}
import sys, os
import pandas as pd, numpy as np
from plotnine import *
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import setup_font, TH, NPG

base = sys.argv[1]; step = sys.argv[2] if len(sys.argv)>2 else "all"
rd = os.path.join(base,"raw_data"); fig = os.path.join(base,"figures", "子图原文件"); os.makedirs(fig, exist_ok=True)
FAM = setup_font(base); T = TH(FAM)
def sv(p,name,w,h):
    p.save(os.path.join(fig,name+".png"), width=w, height=h, dpi=300, verbose=False)
    p.save(os.path.join(fig,name+".pdf"), width=w, height=h, verbose=False)
print("font:", FAM)

if step in ("all","enrichdot"):
    # clusterProfiler风格dotplot：x=-log10 padj，点大小=基因数，颜色=padj
    def dot(e, title, name, w=5.8, h=3.4):
        e = e.copy(); e["term"] = e.term.str.replace(r" \(GO:\d+\)","",regex=True)
        e = e.sort_values("padj")
        e["term"] = pd.Categorical(e.term, e.term.tolist()[::-1])
        e["nl"] = -np.log10(e.padj)
        p = (ggplot(e, aes("nl","term",size="n_overlap",color="padj"))
             + geom_point()
             + scale_color_gradient(low="#B2182B", high="#FDDBC7", name="adj. p",
                                    labels=lambda l:[f"{v:.0e}" for v in l])
             + scale_size_continuous(range=(2.2,6.5), name="genes")
             + labs(x="-log10 adjusted p", y="", title=title)
             + T + theme(axis_text_y=element_text(size=8.5, family=FAM),
                         legend_key_size=9, legend_text=element_text(size=7.5),
                         legend_title=element_text(size=8.5)))
        sv(p, name, w, h)
    ea = pd.read_csv(os.path.join(rd,"11_enrichment_enrichr_pos.csv"))
    cc = ea[ea.lib=="GO_Cellular_Component_2023"].head(4)
    bp = ea[ea.lib=="GO_Biological_Process_2023"].iloc[[0,5,7,9]]
    dot(pd.concat([cc,bp]), "GO enrichment (LONP1-positive genes)", "Fig3B_GO_enrichment")
    dot(ea[ea.lib=="KEGG_2021_Human"].head(8), "KEGG enrichment (LONP1-positive genes)", "Fig3C_KEGG_enrichment")
    en = pd.read_csv(os.path.join(rd,"21_enrichment_enrichr_neg.csv"))
    en2 = pd.concat([en[en.lib=="KEGG_2021_Human"].head(5), en[en.lib=="GO_Biological_Process_2023"].head(3)])
    dot(en2, "Enrichment of LONP1-negative genes", "Fig3D_neg_enrichment")
    print("ENRICHDOT done")

if step in ("all","markerdot"):
    d = pd.read_csv(os.path.join(rd,"26_marker_dotplot_data.csv"))
    ct_order = ["Epithelial","T_NK","B","Plasma","Myeloid","Mast","Endothelial","Fibroblast","SmoothMuscle"]
    gene_order = ["EPCAM","KRT5","KRT8","KRT17","CD3D","CD2","NKG7","MS4A1","CD79A","MZB1","JCHAIN",
                  "LYZ","CD68","AIF1","TPSAB1","CPA3","PECAM1","VWF","COL1A1","DCN","ACTA2","TAGLN"]
    d = d[d.gene.isin(gene_order)].copy()
    d["cell_type"] = pd.Categorical(d.cell_type, ct_order[::-1])
    d["gene"] = pd.Categorical(d.gene, gene_order)
    # 每基因按最大值归一化颜色，展示特异性
    d["scaled"] = d.groupby("gene")["mean_expr"].transform(lambda v: v/v.max())
    p = (ggplot(d, aes("gene","cell_type",size="pct_pos",color="scaled"))
         + geom_point()
         + scale_color_gradient(low="#F0F0F0", high="#B2182B", name="scaled\nmean expr")
         + scale_size_continuous(range=(0.5,5.2), name="fraction\nexpressing",
                                 labels=lambda l:[f"{v*100:.0f}%" for v in l])
         + labs(x="", y="", title="Canonical marker expression across annotated cell types")
         + T + theme(axis_text_x=element_text(rotation=60, ha="right", size=8.5, family=FAM),
                     axis_text_y=element_text(size=9, family=FAM),
                     legend_key_size=9, legend_text=element_text(size=7.5),
                     legend_title=element_text(size=8)))
    sv(p, "Fig6F_marker_dotplot", 7.6, 3.6)
    print("MARKERDOT done")

if step in ("all","marginal"):
    # Fig2A升级：散点+边际直方图（matplotlib gridspec）
    import matplotlib
    import matplotlib.pyplot as plt
    from scipy.stats import spearmanr
    d7 = pd.read_csv(os.path.join(rd,"07_CESC_LONP1_CNA.csv"))
    d8 = pd.read_csv(os.path.join(rd,"08_CESC_gene_panel_expression_RSEM.csv"))
    d8["expr"] = np.log2(d8.LONP1+1)
    m = d7.merge(d8[["sample","expr"]]).dropna(subset=["log2CNA","expr"])
    rho, pv = spearmanr(m.log2CNA, m.expr)
    f = plt.figure(figsize=(4.4,4.4), dpi=300)
    gs = f.add_gridspec(2,2, width_ratios=(5,1.1), height_ratios=(1.1,5), hspace=0.06, wspace=0.06)
    axS = f.add_subplot(gs[1,0]); axT = f.add_subplot(gs[0,0], sharex=axS); axR = f.add_subplot(gs[1,1], sharey=axS)
    axS.scatter(m.log2CNA, m.expr, s=13, alpha=0.55, color=NPG["navy"], edgecolors="none")
    b, a = np.polyfit(m.log2CNA, m.expr, 1)
    xs = np.linspace(m.log2CNA.min(), m.log2CNA.max(), 50)
    axS.plot(xs, a+b*xs, color=NPG["red"], lw=1.6)
    axS.set_xlabel("LONP1 copy number, log2(CNA ratio)", fontsize=10, fontweight="bold", family=FAM)
    axS.set_ylabel("LONP1 mRNA, log2(RSEM+1)", fontsize=10, fontweight="bold", family=FAM)
    axS.text(0.03, 0.95, f"Spearman rho = {rho:.2f}\np = {pv:.2e}", transform=axS.transAxes,
             fontsize=8.5, va="top", family=FAM)
    axT.hist(m.log2CNA, bins=32, color=NPG["blue"], alpha=0.75)
    axR.hist(m.expr, bins=32, orientation="horizontal", color=NPG["salmon"], alpha=0.85)
    for a2 in (axT, axR):
        a2.axis("off")
    for spine in ("top","right"): axS.spines[spine].set_visible(False)
    axS.tick_params(labelsize=9)
    for lab in axS.get_xticklabels()+axS.get_yticklabels(): lab.set_family(FAM)
    axT.set_title(f"TCGA-CESC (n={len(m)})", fontsize=11, fontweight="bold", family=FAM)
    for ext in ("png","pdf"):
        f.savefig(os.path.join(fig, f"Fig2A_CNA_expression_scatter.{ext}"), bbox_inches="tight", facecolor="white", dpi=300)
    print("MARGINAL done")

if step in ("all","raincloud"):
    d = pd.read_csv(os.path.join(rd,"01_TCGA_GTEx_LONP1_expression_tumor_vs_normal.csv"))
    from scipy.stats import mannwhitneyu
    w = mannwhitneyu(d.loc[d.group=="Tumor","LONP1_log2TPM"], d.loc[d.group=="Normal","LONP1_log2TPM"])
    d["g"] = pd.Categorical(d.group, ["Normal","Tumor"])
    rngj = np.random.default_rng(1)
    d["xj"] = d.g.cat.codes + 1 + 0.16 + rngj.uniform(-0.055, 0.055, len(d))
    n_n, n_t = (d.group=="Normal").sum(), (d.group=="Tumor").sum()
    p = (ggplot(d, aes("g","LONP1_log2TPM",fill="g"))
         + geom_violin(style="left", size=0.35, alpha=0.65, position=position_nudge(x=-0.12), width=0.9, adjust=1.2)
         + geom_boxplot(width=0.14, outlier_shape="", alpha=0.9, size=0.5)
         + geom_point(aes(x="xj"), size=0.9, alpha=0.4, color="#555555")
         + scale_fill_manual(values=[NPG["blue"],NPG["red"]])
         + scale_x_discrete(labels=[f"Normal\n(n={n_n})", f"Tumor\n(n={n_t})"])
         + labs(x="", y="LONP1 expression, log2(TPM+0.001)", title="TCGA-CESC + GTEx")
         + annotate("text", x=1.5, y=d.LONP1_log2TPM.max()+0.25, label=f"p = {w.pvalue:.2e}", size=10, family=FAM)
         + T + theme(legend_position="none"))
    sv(p, "Fig1A_TCGA_GTEx_expression", 3.4, 3.9)
    print("RAINCLOUD done")
