# 09_figs_mito_surv_immune_plotnine.py —— Fig3A-E, Fig4A-E, Fig5A-D（统一字体+新panel）
# 数据: raw_data/02,08,11,12,13,14,21   用法: python3 09...py <目录> <step>
import sys, os
import pandas as pd, numpy as np
from scipy.stats import spearmanr, mannwhitneyu
from plotnine import *
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import setup_font, TH, NPG

base = sys.argv[1]; step = sys.argv[2] if len(sys.argv)>2 else "all"
rd = os.path.join(base,"raw_data"); fig = os.path.join(base,"figures", "子图原文件"); os.makedirs(fig, exist_ok=True)
FAM = setup_font(base); T = TH(FAM)
def sv(p,name,w,h):
    p.save(os.path.join(fig,name+".png"), width=w, height=h, dpi=300, verbose=False)
    p.save(os.path.join(fig,name+".pdf"), width=w, height=h, verbose=False)
def pfmt(p): return f"p = {p:.2e}" if p<0.001 else f"p = {p:.3f}"
print("font:", FAM)

if step in ("all","f3a"):
    r = pd.read_csv(os.path.join(rd,"12_mito_panel_spearman.csv")).sort_values("rho")
    r["gene"] = pd.Categorical(r.gene, r.gene.tolist())
    r["xlab"] = np.where(r.rho>0, r.rho+0.045, r.rho-0.045)
    p = (ggplot(r, aes("rho","gene",color="class"))
         + geom_segment(aes(x=0, xend="rho", yend="gene"), size=0.7)
         + geom_point(size=2.6) + geom_vline(xintercept=0, linetype="dashed", color="#808080", size=0.4)
         + geom_text(aes(x="xlab", label="sig"), size=9, color="black", family=FAM)
         + scale_color_manual(values={"Protease/chaperone":NPG["red"],"Fission":NPG["salmon"],
                                      "Fusion":NPG["blue"],"Mitophagy/biogenesis":NPG["purple"]})
         + labs(x="Spearman correlation with LONP1", y="", title="Mitochondrial dynamics & quality control", color="")
         + T + theme(legend_position="right"))
    sv(p,"Fig3A_mito_panel_lollipop",6.2,4.6); print("3A")

if step in ("all","f3bc"):
    ea = pd.read_csv(os.path.join(rd,"11_enrichment_enrichr_pos.csv"))
    def bar(e, title, color, name, w=5.6, h=3.6):
        # 横向条形图：用 x=类别 + coord_flip，规避plotnine横向geom_col方向bug
        e = e.copy(); e["term"] = e.term.str.replace(r" \(GO:\d+\)","",regex=True)
        e["term"] = pd.Categorical(e.term, e.term.tolist()[::-1])
        e["nl"] = -np.log10(e.padj)
        p = (ggplot(e, aes("term","nl")) + geom_col(fill=color, alpha=0.85, width=0.7)
             + geom_text(aes(label="n_overlap"), va="bottom", nudge_y=0.12, size=8, family=FAM)
             + labs(x="", y="-log10 adjusted p", title=title)
             + scale_y_continuous(expand=(0,0,0.14,0))
             + coord_flip()
             + T + theme(axis_text_y=element_text(size=8.5, family=FAM)))
        sv(p, name, w, h)
    cc = ea[ea.lib=="GO_Cellular_Component_2023"].head(4)
    bp = ea[ea.lib=="GO_Biological_Process_2023"].iloc[[0,5,7,9]]
    bar(pd.concat([cc,bp]), "GO enrichment (LONP1-positive genes)", NPG["red"], "Fig3B_GO_enrichment")
    kk = ea[ea.lib=="KEGG_2021_Human"].head(8)
    bar(kk, "KEGG enrichment (LONP1-positive genes)", NPG["navy"], "Fig3C_KEGG_enrichment")
    print("3BC")

if step in ("all","f3de"):
    en = pd.read_csv(os.path.join(rd,"21_enrichment_enrichr_neg.csv"))
    en2 = pd.concat([en[en.lib=="KEGG_2021_Human"].head(5), en[en.lib=="GO_Biological_Process_2023"].head(3)])
    en2 = en2.copy(); en2["term"] = en2.term.str.replace(r" \(GO:\d+\)","",regex=True)
    en2 = en2.sort_values("padj")
    en2["term"] = pd.Categorical(en2.term, en2.term.tolist()[::-1])
    en2["nl"] = -np.log10(en2.padj)
    p = (ggplot(en2, aes("term","nl")) + geom_col(fill=NPG["blue"], alpha=0.85, width=0.7)
         + geom_text(aes(label="n_overlap"), va="bottom", nudge_y=0.04, size=8, family=FAM)
         + geom_hline(yintercept=-np.log10(0.05), linetype="dashed", color="#808080", size=0.4)
         + labs(x="", y="-log10 adjusted p", title="Enrichment of LONP1-negative genes")
         + scale_y_continuous(expand=(0,0,0.14,0))
         + coord_flip()
         + T + theme(axis_text_y=element_text(size=8.5, family=FAM)))
    sv(p,"Fig3D_neg_enrichment",5.6,3.4)
    d8 = pd.read_csv(os.path.join(rd,"08_CESC_gene_panel_expression_RSEM.csv"))
    z = lambda v: (v-v.mean())/v.std()
    lg = np.log2(d8.LONP1+1)
    fis = np.mean([z(np.log2(d8[g]+1)) for g in ["FIS1","MFF"]], axis=0)
    fus = np.mean([z(np.log2(d8[g]+1)) for g in ["MFN1","MFN2","OPA1"]], axis=0)
    bal = fis - fus
    rho, pv = spearmanr(lg, bal)
    df = pd.DataFrame({"x":lg, "y":bal})
    pd.DataFrame({"sample":d8["sample"],"LONP1_log2":lg,"fission_z":fis,"fusion_z":fus,"balance":bal}).to_csv(
        os.path.join(rd,"23_fission_fusion_balance.csv"), index=False)
    p = (ggplot(df, aes("x","y")) + geom_point(size=1.4, alpha=0.55, color=NPG["navy"])
         + geom_smooth(method="lm", color=NPG["red"], fill=NPG["salmon"], size=0.8)
         + labs(x="LONP1 mRNA, log2(RSEM+1)", y="Fission-fusion balance score (z)",
                title=f"Fission-fusion balance vs LONP1\n(Spearman rho = {rho:.2f}, p = {pv:.1e})")
         + T + theme(plot_title=element_text(size=10.5, weight="bold", ha="center", family=FAM)))
    sv(p,"Fig3E_fission_fusion_balance",3.8,3.8)
    print("3DE", f"rho={rho:.3f} p={pv:.2e}")

if step in ("all","f4"):
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
    d2 = pd.read_csv(os.path.join(rd,"02_TCGA_CESC_clinical_LONP1_survival.csv"))
    def km_panel(tk, ek, label, name):
        x = d2.dropna(subset=[tk,ek]); x = x[x[tk]>0].copy()
        med = x.LONP1_log2RSEM.median()
        x["grp"] = np.where(x.LONP1_log2RSEM>med, "High", "Low")
        lr = logrank_test(x.loc[x.grp=="High",tk], x.loc[x.grp=="Low",tk],
                          x.loc[x.grp=="High",ek], x.loc[x.grp=="Low",ek])
        cf = CoxPHFitter().fit(x[[tk,ek,"LONP1_log2RSEM"]], tk, ek)
        s = cf.summary.iloc[0]
        rows, cens = [], []
        for g in ["Low","High"]:
            kmf = KaplanMeierFitter().fit(x.loc[x.grp==g,tk], x.loc[x.grp==g,ek])
            sf = kmf.survival_function_.reset_index()
            sf.columns = ["t","s"]; sf["g"] = g; rows.append(sf)
            cc = x.loc[(x.grp==g)&(x[ek]==0)]
            ct = cc[tk].values
            si = kmf.survival_function_at_times(ct).values
            cens.append(pd.DataFrame({"t":ct,"s":si,"g":g}))
        dd = pd.concat(rows); ce = pd.concat(cens)
        dd["g"] = pd.Categorical(dd.g, ["Low","High"]); ce["g"] = pd.Categorical(ce.g, ["Low","High"])
        nl, nh = (x.grp=="Low").sum(), (x.grp=="High").sum()
        ann = (f"Log-rank p = {lr.p_value:.3f}\nCox (continuous) HR = {s['exp(coef)']:.2f}\n"
               f"(95% CI {s['exp(coef) lower 95%']:.2f}-{s['exp(coef) upper 95%']:.2f}), p = {s['p']:.3f}")
        p = (ggplot(dd, aes("t","s",color="g")) + geom_step(size=0.9)
             + geom_point(data=ce, shape="+", size=2.0, show_legend=False)
             + scale_color_manual(values={"Low":NPG["blue"],"High":NPG["red"]},
                                  labels=[f"Low (n={nl})", f"High (n={nh})"])
             + labs(color="")
             + scale_y_continuous(limits=(0,1), labels=lambda l:[f"{v*100:.0f}%" for v in l])
             + labs(x="Time (months)", y=f"{label} probability", title=f"{label} by LONP1 (median split)")
             + annotate("text", x=dd.t.max()*0.55, y=0.16, label=ann, size=8, ha="left", family=FAM)
             + T + theme(legend_position=(0.78,0.88), legend_background=element_blank()))
        sv(p, name, 4.3, 4.0)
        return dict(ep=label, n=len(x), ev=int(x[ek].sum()), hr=s['exp(coef)'],
                    lo=s['exp(coef) lower 95%'], hi=s['exp(coef) upper 95%'], p=s['p'])
    uni = [km_panel("OS_months","OS_event","Overall survival","Fig4A_KM_OS"),
           km_panel("DSS_months","DSS_event","Disease-specific survival","Fig4B_KM_DSS"),
           km_panel("PFS_months","PFS_event","Progression-free survival","Fig4C_KM_PFS")]
    fu = pd.DataFrame(uni)
    fu["lab"] = fu.apply(lambda r: f"{r.ep.split()[0] if r.ep!='Overall survival' else 'OS'} (n={r.n}, ev={r.ev})", axis=1)
    fu.loc[fu.ep=="Disease-specific survival","lab"] = f"DSS (n={fu.iloc[1].n}, ev={fu.iloc[1].ev})"
    fu.loc[fu.ep=="Progression-free survival","lab"] = f"PFS (n={fu.iloc[2].n}, ev={fu.iloc[2].ev})"
    fu["lab"] = pd.Categorical(fu.lab, fu.lab.tolist()[::-1])
    fu["txt"] = fu.apply(lambda r: f"HR {r.hr:.2f} ({r.lo:.2f}-{r.hi:.2f}), p={r.p:.2f}", axis=1)
    p = (ggplot(fu, aes("hr","lab")) + geom_vline(xintercept=1, linetype="dashed", color="#808080")
         + geom_errorbarh(aes(xmin="lo",xmax="hi"), height=0.15, size=0.7, color=NPG["navy"])
         + geom_point(size=3, shape="s", color=NPG["red"])
         + geom_text(aes(label="txt"), nudge_y=0.28, size=8, family=FAM)
         + scale_x_log10(limits=(0.4,2.2))
         + labs(x="Hazard ratio (per log2 unit LONP1)", y="", title="Univariate Cox regression")
         + T)
    sv(p,"Fig4D_forest_univariate",5.4,3.4)
    # Fig4E 多因素Cox（OS ~ LONP1 + age + T期 + N期）
    x = d2.dropna(subset=["OS_months","OS_event"]); x = x[x.OS_months>0].copy()
    x["Tg"] = x.path_T.astype(str).str.extract(r"^(T[1-4])")[0]
    x = x.dropna(subset=["Tg","age"])
    x = x[x.path_N.isin(["N0","N1"])].copy()
    x["T2"] = (x.Tg=="T2").astype(int); x["T34"] = x.Tg.isin(["T3","T4"]).astype(int)
    x["N1"] = (x.path_N=="N1").astype(int)
    cf = CoxPHFitter().fit(x[["OS_months","OS_event","LONP1_log2RSEM","age","T2","T34","N1"]], "OS_months","OS_event")
    sm = cf.summary.reset_index()
    lab_map = {"LONP1_log2RSEM":"LONP1 (per log2)","age":"Age (per year)","T2":"T2 vs T1","T34":"T3/T4 vs T1","N1":"N1 vs N0"}
    sm["lab"] = sm["covariate"].map(lab_map)
    sm = sm.iloc[::-1]
    sm["lab"] = pd.Categorical(sm.lab, sm.lab.tolist())
    sm["HR"] = np.exp(sm["coef"]); sm["HRlo"] = np.exp(sm["coef lower 95%"]); sm["HRhi"] = np.exp(sm["coef upper 95%"])
    sm["txt"] = sm.apply(lambda r: f"HR {r.HR:.2f} ({r.HRlo:.2f}-{r.HRhi:.2f}), p={r['p']:.3f}", axis=1)
    sm.to_csv(os.path.join(rd,"24_multivariable_cox_OS.csv"), index=False)
    ann = f"n={len(x)}, events={int(x.OS_event.sum())}"
    p = (ggplot(sm, aes("HR","lab")) + geom_vline(xintercept=1, linetype="dashed", color="#808080")
         + geom_errorbarh(aes(xmin="HRlo", xmax="HRhi"), height=0.15, size=0.7, color=NPG["navy"])
         + geom_point(size=3, shape="s", color=NPG["red"])
         + geom_text(aes(label="txt"), nudge_y=0.3, size=8, family=FAM)
         + scale_x_log10()
         + labs(x="Hazard ratio (multivariable, OS)", y="", title=f"Multivariable Cox ({ann})")
         + T)
    sv(p,"Fig4E_forest_multivariable",5.6,3.8)
    print("F4 done. multivariable LONP1:", sm[sm.covariate=="LONP1_log2RSEM"][["coef","p"]].round(4).to_dict("records"))

if step in ("all","f5"):
    d14 = pd.read_csv(os.path.join(rd,"14_immune_signature_spearman.csv")).sort_values("rho")
    d14["sig"] = pd.Categorical(d14.sig, d14.sig.tolist())
    d14["col"] = np.where(d14.padj>=0.05, "ns", np.where(d14.rho<0, "neg", "pos"))
    d14["labx"] = np.where(d14.rho<0, 0.012, -0.012)
    stars = np.where(d14.padj<0.001,"***",np.where(d14.padj<0.01,"**",np.where(d14.padj<0.05,"*","")))
    d14["sigL"] = d14.sig.astype(str) + np.where(stars!="", " "+pd.Series(stars, index=d14.index), "")
    d14["sigL"] = pd.Categorical(d14.sigL, d14.sigL.tolist())
    p = (ggplot(d14, aes("sigL","rho",fill="col")) + geom_col(width=0.68, alpha=0.9)
         + geom_hline(yintercept=0, color="#666666", size=0.4)
         + scale_fill_manual(values={"neg":NPG["blue"],"pos":NPG["red"],"ns":"#cccccc"})
         + labs(x="", y="Spearman correlation with LONP1",
                title="Immune marker signatures (TCGA-CESC)\n(* FDR<0.05, ** <0.01, *** <0.001)")
         + coord_flip()
         + T + theme(legend_position="none", plot_title=element_text(size=10, weight="bold", ha="center", family=FAM)))
    sv(p,"Fig5A_immune_bar",5.4,4.4)
    d8 = pd.read_csv(os.path.join(rd,"08_CESC_gene_panel_expression_RSEM.csv"))
    z = lambda v:(v-v.mean())/v.std(); lg = np.log2(d8.LONP1+1)
    def sig_scatter(genes, label, name):
        s = np.mean([z(np.log2(d8[g]+1)) for g in genes], axis=0)
        rho, pv = spearmanr(lg, s)
        df = pd.DataFrame({"x":lg,"y":s})
        p = (ggplot(df, aes("x","y")) + geom_point(size=1.4, alpha=0.55, color=NPG["navy"])
             + geom_smooth(method="lm", color=NPG["red"], fill=NPG["salmon"], size=0.8)
             + labs(x="LONP1 mRNA, log2(RSEM+1)", y=f"{label} signature (z-score)",
                    title=f"LONP1 vs {label}\n(rho = {rho:.2f}, p = {pv:.1e})") + T
             + theme(plot_title=element_text(size=10.5, weight="bold", ha="center", family=FAM)))
        sv(p, name, 3.8, 3.8)
    sig_scatter(["NOS2","IRF5","PTGS2"], "M1 macrophage", "Fig5B_M1_scatter")
    sig_scatter(["ITGAX","CD1C","NRP1"], "Dendritic cell", "Fig5C_DC_scatter")
    cps = ["PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT"]
    rows = []
    for g in cps:
        rho, pv = spearmanr(lg, np.log2(d8[g]+1)); rows.append({"gene":g,"rho":rho,"p":pv})
    cp = pd.DataFrame(rows); cp["padj"] = cp.p*len(cp)  # Bonferroni上界，保守
    cp = cp.sort_values("rho"); cp["gene"] = pd.Categorical(cp.gene, cp.gene.tolist())
    cp.round(4).to_csv(os.path.join(rd,"25_checkpoint_gene_corr.csv"), index=False)
    cp["col"] = np.where(np.minimum(cp.padj,1)>=0.05, "ns", np.where(cp.rho<0, "neg", "pos"))
    p = (ggplot(cp, aes("gene","rho",fill="col")) + geom_col(width=0.65, alpha=0.9)
         + geom_hline(yintercept=0, color="#666666", size=0.4)
         + scale_fill_manual(values={"neg":NPG["blue"],"pos":NPG["red"],"ns":"#cccccc"})
         + labs(x="", y="Spearman correlation with LONP1", title="Immune checkpoint genes")
         + coord_flip()
         + T + theme(legend_position="none"))
    sv(p,"Fig5D_checkpoint_genes",4.2,3.6)
    print("F5 done")
