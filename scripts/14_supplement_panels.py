# 14_supplement_panels.py —— Fig3F(bulk vs 上皮内在rho对比) + Fig6G(组成效应一致性检验)
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

if step in ("all","f3f"):
    bulk = pd.read_csv(os.path.join(rd,"12_mito_panel_spearman.csv"))[["gene","rho"]].rename(columns={"rho":"Bulk (TCGA-CESC, n=294)"})
    epi = pd.read_csv(os.path.join(rd,"27_epithelial_intrinsic_corr.csv"))
    epi = epi[epi.level.str.startswith("pseudobulk")][["gene","rho"]].rename(columns={"rho":"Epithelium-intrinsic (scRNA pseudobulk, n=63)"})
    genes = ["CLPP","HSPA9","HSPD1","FIS1","MFF","DNM1L","MFN1","MFN2","OPA1"]
    m = bulk.merge(epi, on="gene"); m = m[m.gene.isin(genes)]
    m["gene"] = pd.Categorical(m.gene, genes[::-1])
    L = m.melt(id_vars="gene", var_name="setting", value_name="rho")
    p = (ggplot(L, aes("gene","rho"))
         + geom_line(aes(group="gene"), color="#BBBBBB", size=0.7)
         + geom_point(aes(color="setting"), size=3)
         + geom_hline(yintercept=0, linetype="dashed", color="#808080", size=0.4)
         + scale_color_manual(values={"Bulk (TCGA-CESC, n=294)":NPG["navy"],
                                      "Epithelium-intrinsic (scRNA pseudobulk, n=63)":NPG["red"]})
         + coord_flip()
         + labs(x="", y="Spearman correlation with LONP1", color="",
                title="Bulk vs epithelium-intrinsic correlations")
         + T + theme(legend_position="bottom", legend_text=element_text(size=8, family=FAM),
                     plot_title=element_text(size=11, weight="bold", ha="center", family=FAM)))
    sv(p, "Fig3F_bulk_vs_intrinsic", 5.0, 4.3)
    print("3F done")

if step in ("all","f6g"):
    d = pd.read_csv(os.path.join(rd,"28_TCGA_GTEx_epithelial_markers.csv"))
    d["g"] = pd.Categorical(d.group, ["Normal","Tumor"])
    ann = ("OLS: tumor effect 0.47 (p=0.005) ->\n0.33 (p=0.15) after adjusting for\n"
           "epithelial score (collinearity rho=0.33)")
    p = (ggplot(d, aes("epi_score","LONP1_log2TPM",color="g"))
         + geom_point(size=1.6, alpha=0.6)
         + geom_smooth(method="lm", size=0.8, se=True, alpha=0.15)
         + scale_color_manual(values={"Normal":NPG["blue"],"Tumor":NPG["red"]}, name="")
         + labs(x="Epithelial marker score (mean z of EPCAM/KRT5/KRT8/KRT17)",
                y="LONP1, log2(TPM+0.001)",
                title="Composition consistency check (TCGA+GTEx)")
         + annotate("text", x=d.epi_score.min()+0.2, y=7.3, label=ann, size=7.6, ha="left", family=FAM)
         + T + theme(legend_position=(0.85,0.22), legend_background=element_blank(),
                     plot_title=element_text(size=11, weight="bold", ha="center", family=FAM)))
    sv(p, "Fig6G_composition_check", 5.2, 4.0)
    print("6G done")
