# 07b_fig_scRNA_plotnine.py —— Figure 6（ggplot2语法，plotnine渲染）
# 说明：沙盒无法安装R时的等价替代；07_fig_scRNA.R为R标准版，两者图形语法一致
# 输入: raw_data/17_GSE208653_LONP1_per_cell.csv.gz
import sys, os
import pandas as pd, numpy as np
from scipy.stats import kruskal
from plotnine import *
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import setup_font, TH as THfn

base = sys.argv[1]
FAM = setup_font(base)
rd = os.path.join(base, "raw_data"); fig = os.path.join(base, "figures", "子图原文件"); os.makedirs(fig, exist_ok=True)
d = pd.read_csv(os.path.join(rd, "17_GSE208653_LONP1_per_cell.csv.gz"))
stage_order = ["Normal_HPVneg","Normal_HPVpos","HSIL","Cancer"]
stage_lab = {"Normal_HPVneg":"Normal HPV-","Normal_HPVpos":"Normal HPV+","HSIL":"HSIL","Cancer":"Cancer"}
d["stage"] = pd.Categorical(d["stage"], categories=stage_order, ordered=True)
ctCols = {"Epithelial":"#E64B35","T_NK":"#4DBBD5","Myeloid":"#00A087","Plasma":"#3C5488",
  "B":"#F39B7F","Fibroblast":"#8491B4","Endothelial":"#91D1C2","Mast":"#DC0000","SmoothMuscle":"#7E6148"}
TH = THfn(FAM)

def sv(p, name, w, h):
    p.save(os.path.join(fig, name+".png"), width=w, height=h, dpi=300, verbose=False)
    p.save(os.path.join(fig, name+".pdf"), width=w, height=h, verbose=False)

step = sys.argv[2] if len(sys.argv)>2 else "all"

if step in ("all","a"):
    p = (ggplot(d, aes("umap1","umap2",color="cell_type"))
         + geom_point(size=0.35, alpha=0.8, stroke=0)
         + scale_color_manual(values=ctCols)
         + guides(color=guide_legend(override_aes={"size":3,"alpha":1}))
         + labs(x="UMAP1", y="UMAP2", title="GSE208653: 80,435 cells, 9 samples", color="")
         + TH + theme(legend_text=element_text(size=9)))
    sv(p, "Fig6A_UMAP_celltype", 6.0, 4.4)
    print("6A done")

if step in ("all","b"):
    d2 = d.sort_values("LONP1")
    p = (ggplot(d2, aes("umap1","umap2",color="LONP1"))
         + geom_point(size=0.35, alpha=0.8, stroke=0)
         + scale_color_gradient(low="#e0e0e0", high="#B2182B", name="LONP1")
         + labs(x="UMAP1", y="UMAP2", title="LONP1 expression") + TH)
    sv(p, "Fig6B_UMAP_LONP1", 5.2, 4.4)
    print("6B done")

if step in ("all","c"):
    ord_ct = d.groupby("cell_type").LONP1.mean().sort_values(ascending=False).index.tolist()
    d["ctF"] = pd.Categorical(d["cell_type"], categories=ord_ct, ordered=True)
    p = (ggplot(d, aes("ctF","LONP1",fill="ctF"))
         + geom_violin(style="full", size=0.3, adjust=2, scale="width")
         + stat_summary(fun_y=np.mean, geom="point", size=1.6, color="black")
         + scale_fill_manual(values=ctCols)
         + labs(x="", y="LONP1, log1p(CP10K)", title="LONP1 by cell type")
         + TH + theme(axis_text_x=element_text(rotation=40, ha="right"), legend_position="none"))
    sv(p, "Fig6C_violin_celltype", 5.4, 3.9)
    print("6C done")

if step in ("all","d"):
    epi = d[d.cell_type=="Epithelial"].copy()
    pb = epi.groupby(["sample","stage"], observed=True).LONP1.mean().reset_index()
    kw = kruskal(*[epi[epi.stage==s].LONP1 for s in stage_order])
    ncounts = epi.stage.value_counts().reindex(stage_order)
    xlabels = [f"{stage_lab[s]}\n(n={ncounts[s]})" for s in stage_order]
    p = (ggplot(epi, aes("stage","LONP1",fill="stage"))
         + geom_violin(style="full", size=0.3, adjust=2, scale="width")
         + geom_jitter(data=pb, mapping=aes("stage","LONP1"), width=0.08, size=2.6,
                       fill="white", color="black", stroke=0.6, random_state=1)
         + scale_fill_manual(values=["#4DBBD5","#91D1C2","#F39B7F","#E64B35"])
         + scale_x_discrete(labels=xlabels)
         + labs(x="", y="LONP1, log1p(CP10K)", title="Epithelial cells across disease stages")
         + annotate("text", x=2.5, y=float(epi.LONP1.max())*0.97,
                    label=f"Kruskal-Wallis p = {kw.pvalue:.1e} (per cell)\ncircles = per-sample means", size=9)
         + TH + theme(legend_position="none"))
    sv(p, "Fig6D_epithelial_stage", 5.0, 4.0)
    print("6D done")

if step in ("all","e"):
    comp = (d.groupby(["stage","cell_type"], observed=True).size().reset_index(name="n"))
    comp["frac"] = comp.groupby("stage", observed=True)["n"].transform(lambda x: x/x.sum())
    comp["stageL"] = comp["stage"].map(stage_lab)
    comp["stageL"] = pd.Categorical(comp["stageL"], categories=[stage_lab[s] for s in stage_order], ordered=True)
    p = (ggplot(comp, aes("stageL","frac",fill="cell_type"))
         + geom_col(width=0.72, color="white", size=0.2)
         + scale_fill_manual(values=ctCols, name="")
         + scale_y_continuous(labels=lambda l: [f"{v*100:.0f}%" for v in l])
         + labs(x="", y="Fraction of cells", title="Cellular composition by stage")
         + TH + theme(axis_text_x=element_text(rotation=25, ha="right"), legend_text=element_text(size=8)))
    sv(p, "Fig6E_composition", 5.6, 4.0)
    print("6E done")
