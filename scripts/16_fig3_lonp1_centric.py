# 16_fig3_lonp1_centric.py —— Figure 3（LONP1中心版）：共表达全景 / 19p13.3位置效应 / 线粒体模块 /
#                              bulk vs 上皮内在 / 去19p13.3富集 / 负相关富集
# 数据: raw_data/30 (全基因组共表达+位置), 31 (去19p13.3富集), 21 (负相关富集), 37/38 (上皮内在全基因组假批量rho与富集, 脚本17)
# 用法: python3 16_fig3_lonp1_centric.py <分析根目录> [step]
import sys, os
import pandas as pd, numpy as np
from plotnine import *
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import setup_font, TH, NPG

base = sys.argv[1]; step = sys.argv[2] if len(sys.argv) > 2 else "all"
rd = os.path.join(base, "raw_data"); fig = os.path.join(base, "figures", "子图原文件"); os.makedirs(fig, exist_ok=True)
FAM = setup_font(base); T = TH(FAM)
def sv(p, name, w, h):
    p.save(os.path.join(fig, name + ".png"), width=w, height=h, dpi=300, verbose=False)
    p.save(os.path.join(fig, name + ".pdf"), width=w, height=h, verbose=False)
print("font:", FAM)
co = pd.read_csv(os.path.join(rd, "30_coexpression_LONP1_CESC_genomewide.csv"))
co["pos"] = np.where(co.on_19p13_3, "19p13.3 (LONP1 locus)", "Other loci")
POSCOL = {"19p13.3 (LONP1 locus)": NPG["red"], "Other loci": "#B0B0B0"}

# ---------- 3A 全基因组共表达秩图 ----------
if step in ("all", "a"):
    d = co.dropna(subset=["rho"]).sort_values("rho", ascending=False).reset_index(drop=True)
    d["rank"] = np.arange(1, len(d) + 1)
    npos = int((d.rho >= 0.3).sum()); nneg = int((d.rho <= -0.3).sum())
    top = d.head(6); bot = d.tail(5).iloc[::-1]
    top_lab = "Top positive partners:\n" + "\n".join(f"{g} ({r:.2f})" + ("" if b else "") for g, r, b in zip(top.gene, top.rho, top.on_19p13_3))
    bot_lab = "Top negative partners:\n" + "\n".join(f"{g} ({r:.2f})" for g, r in zip(bot.gene, bot.rho))
    d_other = d[d.pos == "Other loci"]; d_19 = d[d.pos != "Other loci"]
    p = (ggplot()
         + geom_point(d_other, aes("rank", "rho"), color="#B8B8B8", size=0.5, alpha=0.5)
         + geom_point(d_19, aes("rank", "rho"), color=NPG["red"], size=1.0, alpha=0.85)
         + geom_hline(yintercept=[0.3, -0.3], linetype="dashed", color="#808080", size=0.4)
         + geom_hline(yintercept=0, color="#666666", size=0.3)
         + annotate("text", x=len(d) * 0.17, y=0.62, label=top_lab, size=7, family=FAM, ha="left", va="top", color=NPG["red"])
         + annotate("text", x=len(d) * 0.17, y=0.20, label="(all six on 19p13.3)", size=6.5, family=FAM, ha="left", color=NPG["red"])
         + annotate("text", x=len(d) * 0.98, y=-0.10, label=bot_lab, size=7, family=FAM, ha="right", va="top", color=NPG["navy"])
         + annotate("text", x=len(d) * 0.60, y=0.345, label=f"rho >= 0.3: n = {npos}", size=7.5, family=FAM)
         + annotate("text", x=len(d) * 0.40, y=-0.345, label=f"rho <= -0.3: n = {nneg}", size=7.5, family=FAM)
         + annotate("text", x=len(d) * 0.98, y=0.55, ha="right", size=7.5, family=FAM, color=NPG["red"],
                    label="red = genes on 19p13.3 (n = 188)")
         + scale_x_continuous(labels=lambda l: [f"{int(v/1000)}k" for v in l])
         + labs(x=f"Genes ranked by correlation (n = {len(d):,})", y="Spearman rho with LONP1",
                title="Genome-wide LONP1 co-expression (TCGA-CESC, n = 294)")
         + T + theme(plot_title=element_text(size=11, weight="bold", ha="center", family=FAM)))
    sv(p, "Fig3A_coexp_rank", 5.2, 4.0); print("3A")

# ---------- 3B chr19 位置效应 ----------
if step in ("all", "b"):
    c19 = co[co.chr.astype(str) == "19"].dropna(subset=["rho", "start"]).copy()
    c19["Mb"] = c19.start / 1e6
    lon = c19[c19.gene == "LONP1"]
    lon_mb = float(lon.Mb.iloc[0]) if len(lon) else 5.69
    # 19p13.3 边界 (hg38 cytoBand): 0-6.9 Mb
    band = pd.DataFrame({"xmin": [0], "xmax": [6.9], "ymin": [-0.6], "ymax": [0.8]})
    # 滑动中位数（1 Mb窗，步长0.25 Mb）
    xs = np.arange(0, c19.Mb.max(), 0.25); med = []
    for x in xs:
        w = c19[(c19.Mb >= x - 0.5) & (c19.Mb < x + 0.5)]
        med.append(w.rho.median() if len(w) >= 5 else np.nan)
    sm = pd.DataFrame({"Mb": xs, "med": med}).dropna()
    m19 = c19[c19.on_19p13_3].rho.median(); mo = co[(~co.on_19p13_3) & co.chr.notna()].rho.median()
    lab = c19[c19.gene.isin(["CLPP"])].copy()
    p = (ggplot()
         + geom_rect(band, aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax"), fill="#FDECEA", alpha=0.9)
         + annotate("text", x=3.45, y=0.76, label="19p13.3", size=8, family=FAM, color=NPG["red"], fontweight="bold")
         + geom_hline(yintercept=0, color="#666666", size=0.3)
         + geom_point(c19, aes("Mb", "rho"), color="#8C8C8C", size=0.9, alpha=0.55)
         + geom_line(sm, aes("Mb", "med"), color=NPG["navy"], size=0.9)
         + geom_vline(xintercept=lon_mb, linetype="dashed", color=NPG["red"], size=0.5)
         + geom_point(lab, aes("Mb", "rho"), color=NPG["red"], size=1.8)
         + geom_text(lab, aes("Mb", "rho", label="gene"), size=7, family=FAM, nudge_x=1.2, ha="left", color=NPG["red"])
         + annotate("text", x=lon_mb + 0.6, y=-0.5, label="LONP1", size=8, family=FAM, color=NPG["red"], ha="left", fontweight="bold")
         + annotate("text", x=57, y=0.70, ha="right", size=7.5, family=FAM,
                    label=f"median rho: 19p13.3 = {m19:.2f}; all other loci = {mo:.2f}\nline = 1-Mb sliding median")
         + coord_cartesian(ylim=(-0.6, 0.8), xlim=(0, 59))
         + labs(x="Position on chromosome 19 (Mb, GRCh38)", y="Spearman rho with LONP1",
                title="Positional co-expression around the LONP1 locus")
         + T + theme(plot_title=element_text(size=11, weight="bold", ha="center", family=FAM)))
    sv(p, "Fig3B_chr19_positional", 5.2, 4.0); print("3B")

# ---------- 3C 线粒体蛋白稳态/动态模块 ----------
MOD = [
    ("Matrix protease", ["CLPP", "AFG3L2", "YME1L1", "SPG7", "PMPCB"]),
    ("Chaperone", ["HSPD1", "HSPE1", "HSPA9", "TRAP1", "DNAJA3"]),
    ("UPRmt / ISR TF", ["ATF4", "ATF5", "DDIT3"]),
    ("Prohibitin scaffold", ["PHB1", "PHB2"]),
    ("Biogenesis / mtDNA", ["TFAM", "NRF1", "PPARGC1A", "POLRMT"]),
    ("Fission", ["DNM1L", "FIS1", "MFF", "MIEF1", "MIEF2"]),
    ("Fusion", ["MFN1", "MFN2", "OPA1"]),
    ("Mitophagy", ["PINK1", "PRKN"]),
]
if step in ("all", "c"):
    rows = []
    for cls, genes in MOD:
        for g in genes:
            x = co[co.gene == g]
            if len(x): rows.append({"gene": g, "cls": cls, "rho": x.rho.iloc[0], "p": x.p.iloc[0], "b19": bool(x.on_19p13_3.iloc[0])})
    m = pd.DataFrame(rows)
    m["sig"] = np.where(m.p < 0.001, "***", np.where(m.p < 0.01, "**", np.where(m.p < 0.05, "*", "")))
    m["lab"] = m.gene + np.where(m.b19, " (19p13.3)", "")
    m = m.sort_values("rho")
    m["lab"] = pd.Categorical(m.lab, m.lab.tolist())
    m["xlab"] = np.where(m.rho > 0, m.rho + 0.04, m.rho - 0.04)
    m.drop(columns=["xlab"]).round(5).to_csv(os.path.join(rd, "33_mito_module_spearman.csv"), index=False)
    cols = {"Matrix protease": NPG["red"], "Chaperone": NPG["salmon"], "UPRmt / ISR TF": NPG["darkred"],
            "Prohibitin scaffold": NPG["brown"], "Biogenesis / mtDNA": NPG["purple"],
            "Fission": NPG["green"], "Fusion": NPG["blue"], "Mitophagy": NPG["teal"]}
    p = (ggplot(m, aes("rho", "lab", color="cls"))
         + geom_segment(aes(x=0, xend="rho", yend="lab"), size=0.7)
         + geom_point(size=2.6)
         + geom_vline(xintercept=0, linetype="dashed", color="#808080", size=0.4)
         + geom_text(aes(x="xlab", label="sig"), size=8, color="black", family=FAM)
         + scale_color_manual(values=cols, breaks=list(cols.keys()))
         + scale_x_continuous(limits=(-0.32, 0.82))
         + labs(x="Spearman correlation with LONP1 (TCGA-CESC, n = 294)", y="",
                title="Mitochondrial proteostasis, biogenesis and dynamics genes", color="")
         + T + theme(legend_position="right", legend_key_size=8, legend_text=element_text(size=8, family=FAM),
                     axis_text_y=element_text(size=8.5, family=FAM),
                     plot_title=element_text(size=11, weight="bold", ha="center", family=FAM)))
    sv(p, "Fig3C_mito_module_lollipop", 6.4, 5.6); print("3C")

# ---------- 3D bulk vs 上皮内在（扩展模块；上皮内在=raw_data/37全基因组假批量rho，n=58单元） ----------
if step in ("all", "d"):
    ep = pd.read_csv(os.path.join(rd, "37_epithelial_pseudobulk_genomewide_rho.csv")).set_index("gene")
    q1, med, q3 = ep.rho_pb.quantile(0.25), ep.rho_pb.median(), ep.rho_pb.quantile(0.75)
    bq1, bq3 = co.rho.quantile(0.25), co.rho.quantile(0.75)
    groups = [("Proteostasis", ["CLPP", "HSPA9", "TRAP1", "DNAJA3", "PMPCB", "PHB1", "PHB2", "HSPE1", "HSPD1", "YME1L1", "AFG3L2"]),
              ("ISR TF", ["ATF4", "ATF5", "DDIT3"]),
              ("Biogenesis", ["POLRMT", "TFAM", "PPARGC1A", "NRF1"]),
              ("Dynamics", ["FIS1", "OPA1", "DNM1L", "MFF", "MFN1", "MFN2"])]
    alias = {"PHB1": "PHB"}  # scRNA注释中PHB1记为PHB
    rows = []
    for fac, genes in groups:
        for g in genes:
            ge = alias.get(g, g)
            if g in set(co.gene) and ge in ep.index:
                rows.append({"gene": g, "facet": fac, "Bulk tumor (TCGA-CESC, n = 294)": float(co[co.gene == g].rho.iloc[0]),
                             "Epithelium-intrinsic (scRNA pseudobulk, n = 58)": float(ep.rho_pb[ge]),
                             "pct": float(ep.pct_rank_pb[ge])})
    m = pd.DataFrame(rows)
    m["facet"] = pd.Categorical(m.facet, [g[0] for g in groups])
    order = [g for _, gs in groups for g in gs if g in set(m.gene)]
    m["gene"] = pd.Categorical(m.gene, order[::-1])
    m.round(4).to_csv(os.path.join(rd, "39_fig3D_bulk_vs_epithelium_module.csv"), index=False)
    L = m.melt(id_vars=["gene", "facet", "pct"], var_name="setting", value_name="rho")
    band = pd.DataFrame({"xmin": [q1], "xmax": [q3], "ymin": [-np.inf], "ymax": [np.inf]})
    p = (ggplot(L, aes("gene", "rho"))
         + geom_rect(band, aes(xmin="ymin", xmax="ymax", ymin="xmin", ymax="xmax"), fill="#F1E5DC", alpha=0.7, inherit_aes=False)
         + geom_hline(yintercept=med, linetype="dotted", color="#B09C85", size=0.6)
         + geom_hline(yintercept=0, linetype="dashed", color="#808080", size=0.4)
         + geom_line(aes(group="gene"), color="#BBBBBB", size=0.7)
         + geom_point(aes(color="setting"), size=2.7)
         + scale_color_manual(values={"Bulk tumor (TCGA-CESC, n = 294)": NPG["navy"],
                                      "Epithelium-intrinsic (scRNA pseudobulk, n = 58)": NPG["red"]})
         + facet_grid("facet ~ .", scales="free_y", space="free_y")
         + coord_flip()
         + labs(x="", y="Spearman correlation with LONP1", color="",
                title="Bulk vs epithelium-intrinsic correlations\n(band = epithelium genome-wide IQR; dotted = median)")
         + T + theme(legend_position="bottom", legend_direction="vertical",
                     legend_text=element_text(size=8, family=FAM),
                     axis_text_y=element_text(size=8.5, family=FAM),
                     strip_text=element_text(size=8.5, weight="bold", family=FAM),
                     strip_background=element_rect(fill="#EFEFEF", color="none"),
                     plot_title=element_text(size=10, weight="bold", ha="center", family=FAM)))
    sv(p, "Fig3D_bulk_vs_intrinsic", 4.8, 6.6); print("3D", f"epi IQR {q1:.2f}-{q3:.2f} med {med:.2f}; bulk IQR {bq1:.2f}-{bq3:.2f}")

# ---------- 3E 去19p13.3后富集 ----------
if step in ("all", "e"):
    e = pd.read_csv(os.path.join(rd, "31_enrichment_enrichr_pos_excl19p13.3.csv"))
    keep = {"GO_Cellular_Component_2023": ["Mitochondrial Inner Membrane (GO:0005743)", "Mitochondrial Membrane (GO:0031966)",
                                          "Mitochondrial Matrix (GO:0005759)", "Mitochondrial Ribosome (GO:0005761)"],
            "GO_Biological_Process_2023": ["Translation (GO:0006412)", "Cytoplasmic Translation (GO:0002181)",
                                          "Mitochondrial Translation (GO:0032543)", "Mitochondrial Gene Expression (GO:0140053)",
                                          "Aerobic Respiration (GO:0009060)", "Cellular Respiration (GO:0045333)"],
            "KEGG_2021_Human": ["Ribosome", "Oxidative phosphorylation"]}
    rows = []
    for lib, terms in keep.items():
        for t in terms:
            x = e[(e.lib == lib) & (e.term == t)]
            if len(x): rows.append({"lib": {"GO_Cellular_Component_2023": "GO CC", "GO_Biological_Process_2023": "GO BP", "KEGG_2021_Human": "KEGG"}[lib],
                                    "term": t.split(" (GO:")[0], "padj": x.padj.iloc[0], "n": x.n_overlap.iloc[0]})
    d = pd.DataFrame(rows); d["nl"] = -np.log10(d.padj)
    d["term"] = d.term.replace({"Oxidative phosphorylation": "Oxidative phosphorylation (n.s.)"})
    d = d.sort_values(["lib", "padj"], ascending=[True, False])
    d["term"] = pd.Categorical(d.term, d.term.tolist())
    d["lib"] = pd.Categorical(d.lib, ["GO CC", "GO BP", "KEGG"])
    n_in = 102
    p = (ggplot(d, aes("nl", "term", size="n", color="padj"))
         + geom_vline(xintercept=-np.log10(0.05), linetype="dashed", color="#808080", size=0.4)
         + geom_point()
         + scale_color_gradient(low="#B2182B", high="#FDDBC7", name="adj. p", labels=lambda l: [f"{v:.0e}" for v in l])
         + scale_size_continuous(range=(2.0, 6.5), name="genes")
         + facet_grid("lib ~ .", scales="free_y", space="free_y")
         + labs(x="-log10 adjusted p", y="",
                title=f"Enrichment of LONP1-positive genes\n(rho >= 0.4, 19p13.3 genes excluded; n = {n_in})")
         + T + theme(axis_text_y=element_text(size=8.5, family=FAM), legend_key_size=9,
                     legend_text=element_text(size=7.5), legend_title=element_text(size=8.5),
                     strip_text=element_text(size=8.5, weight="bold", family=FAM),
                     strip_background=element_rect(fill="#EFEFEF", color="none"),
                     plot_title=element_text(size=10.5, weight="bold", ha="center", family=FAM)))
    sv(p, "Fig3E_enrichment_excl19p13", 5.8, 4.2); print("3E")

# ---------- 3F 上皮内在共表达伙伴富集（top300, raw_data/38） ----------
if step in ("all", "f"):
    e = pd.read_csv(os.path.join(rd, "38_enrichment_epithelial_intrinsic_top300.csv"))
    e = e[e.set == "top300_raw"]
    keep = {"GO_Cellular_Component_2023": ["Mitochondrial Membrane (GO:0031966)", "Mitochondrial Inner Membrane (GO:0005743)",
                                          "Organelle Inner Membrane (GO:0019866)", "Endoplasmic Reticulum Membrane (GO:0005789)"],
            "GO_Biological_Process_2023": ["Mitochondrial Respiratory Chain Complex Assembly (GO:0033108)",
                                          "Mitochondrial Cytochrome C Oxidase Assembly (GO:0033617)",
                                          "Respiratory Chain Complex IV Assembly (GO:0008535)", "Protein Transport (GO:0015031)"]}
    rows = []
    for lib, terms in keep.items():
        for t in terms:
            x = e[(e.lib == lib) & (e.term == t)]
            if len(x): rows.append({"lib": {"GO_Cellular_Component_2023": "GO CC", "GO_Biological_Process_2023": "GO BP"}[lib],
                                    "term": t.split(" (GO:")[0].replace("Mitochondrial Respiratory Chain Complex Assembly", "Mito. respiratory chain complex assembly")
                                    .replace("Mitochondrial Cytochrome C Oxidase Assembly", "Mito. cytochrome c oxidase assembly"),
                                    "padj": x.padj.iloc[0], "n": x.n_overlap.iloc[0]})
    d = pd.DataFrame(rows); d["nl"] = -np.log10(d.padj)
    d = d.sort_values(["lib", "padj"], ascending=[True, False]); d["term"] = pd.Categorical(d.term, d.term.tolist())
    d["lib"] = pd.Categorical(d.lib, ["GO CC", "GO BP"])
    p = (ggplot(d, aes("nl", "term", size="n", color="padj"))
         + geom_vline(xintercept=-np.log10(0.05), linetype="dashed", color="#808080", size=0.4)
         + geom_point()
         + scale_color_gradient(low="#B2182B", high="#FDDBC7", name="adj. p", labels=lambda l: [f"{v:.0e}" for v in l])
         + scale_size_continuous(range=(2.0, 6.5), name="genes")
         + facet_grid("lib ~ .", scales="free_y", space="free_y")
         + labs(x="-log10 adjusted p", y="",
                title="Enrichment of epithelium-intrinsic LONP1 partners\n(top 300 by pseudobulk rho; KEGG: none significant)")
         + T + theme(axis_text_y=element_text(size=8.5, family=FAM), legend_key_size=9,
                     legend_text=element_text(size=7.5), legend_title=element_text(size=8.5),
                     strip_text=element_text(size=8.5, weight="bold", family=FAM),
                     strip_background=element_rect(fill="#EFEFEF", color="none"),
                     plot_title=element_text(size=10, weight="bold", ha="center", family=FAM)))
    sv(p, "Fig3F_epithelium_enrichment", 5.8, 3.6); print("3F")

# ---------- 3G 负相关富集（沿用21号数据） ----------
if step in ("all", "g"):
    en = pd.read_csv(os.path.join(rd, "21_enrichment_enrichr_neg.csv"))
    en2 = pd.concat([en[en.lib == "KEGG_2021_Human"].head(5), en[en.lib == "GO_Biological_Process_2023"].head(3)]).copy()
    en2["term"] = en2.term.str.replace(r" \(GO:\d+\)", "", regex=True)
    en2 = en2.sort_values("padj")
    en2["term"] = pd.Categorical(en2.term, en2.term.tolist()[::-1]); en2["nl"] = -np.log10(en2.padj)
    p = (ggplot(en2, aes("nl", "term", size="n_overlap", color="padj"))
         + geom_vline(xintercept=-np.log10(0.05), linetype="dashed", color="#808080", size=0.4)
         + geom_point()
         + scale_color_gradient(low="#2166AC", high="#D1E5F0", name="adj. p", labels=lambda l: [f"{v:.0e}" for v in l])
         + scale_size_continuous(range=(2.0, 6.5), name="genes")
         + labs(x="-log10 adjusted p", y="", title="Enrichment of LONP1-negative genes (bulk)\n(rho <= -0.35; n = 96)")
         + T + theme(axis_text_y=element_text(size=8.5, family=FAM), legend_key_size=9,
                     legend_text=element_text(size=7.5), legend_title=element_text(size=8.5),
                     plot_title=element_text(size=10.5, weight="bold", ha="center", family=FAM)))
    sv(p, "Fig3G_neg_enrichment", 5.0, 3.6); print("3G")
print("FIG3_DONE")
