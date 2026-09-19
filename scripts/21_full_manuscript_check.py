# 21_full_manuscript_check.py —— 稿件全文数值 vs raw_data 程序化复核（v6→v7 口径）
# 覆盖：摘要/Results/图注/Additional files 中所有数值。用法: python3 21_full_manuscript_check.py <分析根目录>
# 输出: raw_data/46_manuscript_number_check.txt
import sys, os, warnings, gzip
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind, spearmanr, kruskal, fisher_exact
import statsmodels.api as sm
warnings.filterwarnings("ignore")
base = sys.argv[1]; rd = os.path.join(base, "raw_data")
R = lambda f: pd.read_csv(os.path.join(rd, f))
chk = []
def C(name, val, exp, tol=0.006):
    ok = abs(float(val) - float(exp)) <= tol; chk.append((name, round(float(val), 5), exp, ok))
def CL(name, val, exp, tol=0.05): C(name, np.log10(float(val)), np.log10(float(exp)), tol)

# ---- R1 TCGA+GTEx (Results para 1, Fig 1B) ----
d1 = R("01_TCGA_GTEx_LONP1_expression_tumor_vs_normal.csv")
t, n = d1[d1.group == "Tumor"].LONP1_log2TPM, d1[d1.group == "Normal"].LONP1_log2TPM
C("n tumor", len(t), 304, 0); C("n normal", len(n), 13, 0)
C("mean T", t.mean(), 5.60); C("sd T", t.std(), 0.59); C("mean N", n.mean(), 5.13); C("sd N", n.std(), 0.33)
C("log2FC", t.mean() - n.mean(), 0.47); C("fold", 2 ** (t.mean() - n.mean()), 1.4, 0.05)
w = ttest_ind(t, n, equal_var=False); C("Welch t", w.statistic, 4.78, 0.01); CL("Welch p", w.pvalue, 2.3e-4); CL("MWU p", mannwhitneyu(t, n).pvalue, 7.2e-4)
# ---- GEO cohorts ----
g4 = R("04_GSE9750_LONP1_209017_s_at.csv"); ca = g4[g4.group == "cancer"].value_MAS5 if "cancer" in set(g4.group) else None
grpcol = g4.group.unique().tolist()
tum = g4[g4.group.str.lower().str.startswith("c") & ~g4.group.str.lower().str.contains("line")].value_MAS5
nor = g4[g4.group.str.lower().str.startswith("n")].value_MAS5
C("GSE9750 n cancer", len(tum), 31, 0); C("GSE9750 n normal", len(nor), 24, 0); C("GSE9750 fold", tum.mean() / nor.mean(), 1.9, 0.05); CL("GSE9750 p", mannwhitneyu(tum, nor).pvalue, 2.3e-5)
g5 = R("05_GSE7803_LONP1_209017_s_at.csv"); s5, n5, h5 = g5[g5.group == "scc"].value_log2, g5[g5.group == "normal"].value_log2, g5[g5.group == "hsil"].value_log2
C("GSE7803 n", len(s5) + len(n5) + len(h5), 38, 0); C("GSE7803 p", mannwhitneyu(s5, n5).pvalue, 0.052, 0.001)
g3 = R("03_GSE63514_LONP1_209017_s_at.csv"); c3, n3 = g3[g3.group == "Cancer"].value_log2norm, g3[g3.group == "Normal"].value_log2norm
C("GSE63514 n cancer", len(c3), 28, 0); C("GSE63514 n normal", len(n3), 24, 0); C("GSE63514 p", mannwhitneyu(c3, n3).pvalue, 0.98, 0.005)
C("GSE63514 KW p", kruskal(*[g.value_log2norm for _, g in g3.groupby("group")]).pvalue, 0.35, 0.005)
C("GSE63514 CIN1 n", (g3.group == "CIN1").sum(), 14, 0); C("CIN2 n", (g3.group == "CIN2").sum(), 22, 0); C("CIN3 n", (g3.group == "CIN3").sum(), 40, 0)
# ---- pan-cancer ----
p16 = R("16_pancancer_LONP1_stats.csv"); C("pancancer types", len(p16), 32, 0); CL("CESC padj", p16[p16.type == "CESC"].padj.iloc[0], 7.4e-4)
C("OV direction (tumor<normal)", float(p16[p16.type == "OV"].medT.iloc[0] < p16[p16.type == "OV"].medN.iloc[0]), 1, 0)
# ---- single cell ----
sc17 = pd.read_csv(os.path.join(rd, "17_GSE208653_LONP1_per_cell.csv.gz"))
C("cells", len(sc17), 80435, 0); ct = sc17.cell_type.value_counts()
C("epithelial n", ct["Epithelial"], 20363, 0); C("T_NK n", ct["T_NK"], 32027, 0); C("Myeloid n", ct["Myeloid"], 13136, 0); C("n cell types", len(ct), 9, 0)
m = sc17.groupby("cell_type").LONP1.mean()
C("epi mean", m["Epithelial"], 0.143, 0.0006); C("epi pct pos", (sc17[sc17.cell_type == "Epithelial"].LONP1 > 0).mean() * 100, 24.9, 0.06)
C("T_NK mean", m["T_NK"], 0.061, 0.0006); C("Myeloid mean", m["Myeloid"], 0.065, 0.0006)
C("min non-epi mean", m.drop("Epithelial").min(), 0.050, 0.0006); C("max non-epi mean", m.drop("Epithelial").max(), 0.127, 0.0006)
epi = sc17[sc17.cell_type == "Epithelial"]; ms = epi.groupby("stage").LONP1.mean()
C("epi cancer", ms["Cancer"], 0.113, 0.0006); C("epi HPVneg", ms["Normal_HPVneg"], 0.132, 0.0006); C("epi HPVpos", ms["Normal_HPVpos"], 0.142, 0.0006); C("epi HSIL", ms["HSIL"], 0.196, 0.0006)
normal_all = epi[epi.stage.str.startswith("Normal")].LONP1; C("epi normal pooled (abstract v3 0.136)", normal_all.mean(), 0.136, 0.001)
CL("epi cancer vs normal MWU", mannwhitneyu(epi[epi.stage == "Cancer"].LONP1, normal_all).pvalue, 1.0e-16, 0.3)
CL("epi KW stages", kruskal(*[g.LONP1 for _, g in epi.groupby("stage")]).pvalue, 3.3e-23, 0.3)
pb = R("19_GSE208653_epithelial_pseudobulk_LONP1.csv"); cn = pb[pb.stage == "Cancer"].mean_LONP1; nn = pb[pb.stage != "Cancer"].mean_LONP1
C("pseudobulk cancer lowest 3", float(cn.max() < nn.min()), 1, 0); C("pseudobulk MWU p (cancer vs 4 normal)", mannwhitneyu(cn, pb[pb.stage.str.startswith("Normal")].mean_LONP1).pvalue, 0.057, 0.003)
# composition regression
c29 = R("29_composition_regression.csv"); d28 = R("28_TCGA_GTEx_epithelial_markers.csv")
C("epi score normal", d28[d28.group == "Normal"].epi_score.mean(), -1.99, 0.01); C("epi score tumor", d28[d28.group == "Tumor"].epi_score.mean(), 0.09, 0.01)
C("beta unadj", c29[(c29.model == "unadjusted") & (c29.term == "tumor")].beta.iloc[0], 0.47); C("p unadj", c29[(c29.model == "unadjusted") & (c29.term == "tumor")].p.iloc[0], 0.005, 0.0005)
C("beta adj", c29[(c29.model == "adjusted") & (c29.term == "tumor")].beta.iloc[0], 0.33); C("p adj", c29[(c29.model == "adjusted") & (c29.term == "tumor")].p.iloc[0], 0.15, 0.005)
C("collinearity rho", spearmanr(d28.tumor, d28.epi_score)[0], 0.33)
# ---- CNA / GISTIC / subtype / T N ----
d7 = R("07_CESC_LONP1_CNA.csv"); d8 = R("08_CESC_gene_panel_expression_RSEM.csv"); d8["e"] = np.log2(d8.LONP1 + 1)
mc = d7.merge(d8[["sample", "e"]]).dropna(subset=["log2CNA"]); C("CNA n", len(mc), 290, 0)
r, p = spearmanr(mc.log2CNA, mc.e); C("CNA rho", r, 0.41); CL("CNA p", p, 3.1e-13)
mg = d7.merge(d8[["sample", "e"]]).dropna(subset=["gistic"]); gm = mg.groupby("gistic").e.agg(["mean", "count"])
C("GISTIC -1 mean", gm.loc[-1, "mean"], 10.68); C("GISTIC -1 n", gm.loc[-1, "count"], 94, 0); C("GISTIC 0 mean", gm.loc[0, "mean"], 11.02); C("GISTIC 0 n", gm.loc[0, "count"], 159, 0)
C("GISTIC 1 mean", gm.loc[1, "mean"], 11.45); C("GISTIC 2 mean", gm.loc[2, "mean"], 11.42); C("GISTIC 1+2 n", gm.loc[1, "count"] + gm.loc[2, "count"], 31, 0)
CL("GISTIC KW p", kruskal(*[g.e for _, g in mg.groupby("gistic")]).pvalue, 1.1e-8)
d10 = R("10_CESC_subtype.csv"); d8["patient"] = d8["sample"].str[:12]; ms8 = d8.merge(d10, on="patient")
ad = ms8[ms8.subtype.str.contains("Adeno", case=False)].e; sq = ms8[ms8.subtype.str.contains("Squam", case=False)].e
C("adeno n", len(ad), 43, 0); C("squam n", len(sq), 229, 0); C("adeno mean", ad.mean(), 11.12); C("adeno sd", ad.std(), 0.53); C("squam mean", sq.mean(), 10.91); C("squam sd", sq.std(), 0.59)
C("subtype p", ttest_ind(ad, sq, equal_var=False).pvalue, 0.022, 0.001)
d2 = R("02_TCGA_CESC_clinical_LONP1_survival.csv")
d2["Tg"] = d2.path_T.astype(str).str.extract(r"^(T[1-4])")[0]; d2["Tg"] = d2.Tg.replace({"T3": "T3/T4", "T4": "T3/T4"})
xT = d2.dropna(subset=["Tg"]); C("T stage KW p", kruskal(*[g.LONP1_log2RSEM for _, g in xT.groupby("Tg")]).pvalue, 0.49, 0.005)
xN = d2[d2.path_N.isin(["N0", "N1"])]; C("N MWU p", mannwhitneyu(xN[xN.path_N == "N0"].LONP1_log2RSEM, xN[xN.path_N == "N1"].LONP1_log2RSEM).pvalue, 0.999, 0.001)
# ---- survival ----
c13 = R("13_cox_forest_data.csv").set_index("ep")
for ep, hr, lo, hi, p in [("OS", 0.84, 0.54, 1.31, 0.45), ("DSS", 0.99, 0.59, 1.64, 0.96), ("PFS", 0.85, 0.55, 1.32, 0.47)]:
    C(f"{ep} HR", c13.loc[ep, "hr"], hr); C(f"{ep} lo", c13.loc[ep, "lo"], lo); C(f"{ep} hi", c13.loc[ep, "hi"], hi); C(f"{ep} p", c13.loc[ep, "p"], p, 0.006)
C("OS n", c13.loc["OS", "n"], 281, 0); C("OS ev", c13.loc["OS", "ev"], 68, 0); C("DSS n", c13.loc["DSS", "n"], 277, 0); C("DSS ev", c13.loc["DSS", "ev"], 51, 0)
from lifelines.statistics import logrank_test
x = d2.dropna(subset=["OS_months", "OS_event"]); x = x[x.OS_months > 0]; med = x.LONP1_log2RSEM.median(); hi = x.LONP1_log2RSEM > med
C("OS logrank p", logrank_test(x[hi].OS_months, x[~hi].OS_months, x[hi].OS_event, x[~hi].OS_event).p_value, 0.43, 0.006)
c24 = R("24_multivariable_cox_OS.csv").set_index("covariate")
C("MV LONP1 HR", c24.loc["LONP1_log2RSEM", "HR"], 0.98); C("MV LONP1 lo", c24.loc["LONP1_log2RSEM", "HRlo"], 0.54); C("MV LONP1 hi", c24.loc["LONP1_log2RSEM", "HRhi"], 1.81); C("MV LONP1 p", c24.loc["LONP1_log2RSEM", "p"], 0.96, 0.006)
C("MV N1 HR", c24.loc["N1", "HR"], 2.91); C("MV N1 lo", c24.loc["N1", "HRlo"], 1.40); C("MV N1 hi", c24.loc["N1", "HRhi"], 6.04); C("MV N1 p", c24.loc["N1", "p"], 0.004, 0.0006)
# ---- immune ----
d14 = R("14_immune_signature_spearman.csv").set_index("sig")
C("M1 rho", d14.loc["Macrophage M1", "rho"], -0.22); C("M1 FDR", d14.loc["Macrophage M1", "padj"], 0.0012, 0.0001); C("DC rho", d14.loc["Dendritic cell", "rho"], -0.22); C("DC FDR", d14.loc["Dendritic cell", "padj"], 0.0012, 0.0001)
C("CD4 rho", d14.loc["CD4/Th", "rho"], -0.16); C("CD4 FDR", d14.loc["CD4/Th", "padj"], 0.028, 0.001)
C("n sig significant", int((d14.padj < 0.05).sum()), 3, 0)
d25 = R("25_checkpoint_gene_corr.csv").set_index("gene"); C("CD274 rho", d25.loc["CD274", "rho"], -0.17); C("CD274 padj", d25.loc["CD274", "padj"], 0.023, 0.001); C("max |rho| checkpoint", d25.rho.abs().max(), 0.17, 0.005)
C("n checkpoint significant", int((d25.padj < 0.05).sum()), 1, 0)
# ---- co-expression (bulk) ----
co = R("30_coexpression_LONP1_CESC_genomewide.csv"); ann = co[co.chr.notna()]
C("genes", len(co), 17378, 0); C("|rho|>=0.3", (co.rho.abs() >= 0.3).sum(), 970, 0); C("pos", (co.rho >= 0.3).sum(), 628, 0); C("neg", (co.rho <= -0.3).sum(), 342, 0)
C("top20 19p", ann.head(20).on_19p13_3.sum(), 18, 0); C("top50 19p", ann.head(50).on_19p13_3.sum(), 35, 0)
for g, e in {"CLPP": 0.72, "THOP1": 0.68, "HDGFL2": 0.67, "WDR18": 0.67, "SGTA": 0.64, "NDUFA11": 0.61, "TIMM13": 0.60, "TRAP1": 0.41, "HSPD1": 0.33, "HSPE1": 0.30, "HSPA9": 0.23, "PHB1": 0.40, "PHB2": 0.39, "PMPCB": 0.22, "ATF4": 0.32, "DDIT3": 0.22, "ATF5": 0.14, "FIS1": 0.32, "MFF": 0.15, "MFN1": -0.23, "OPA1": -0.22, "MFN2": -0.13, "DNM1L": -0.02}.items():
    C("bulk rho " + g, co[co.gene == g].rho.iloc[0], e)
C("ATF5 p", co[co.gene == "ATF5"].p.iloc[0], 0.018, 0.0006); CL("TRAP1 p", co[co.gene == "TRAP1"].p.iloc[0], 2.5e-13)
C("19p13.3 n", ann.on_19p13_3.sum(), 188, 0); C("annotated n", len(ann), 17324, 0); C("bg frac %", ann.on_19p13_3.mean() * 100, 1.1, 0.05)
sel = ann[ann.rho >= 0.4]; rest = ann[ann.rho < 0.4]; orr, pv = fisher_exact([[sel.on_19p13_3.sum(), (~sel.on_19p13_3).sum()], [rest.on_19p13_3.sum(), (~rest.on_19p13_3).sum()]], alternative="greater")
C("rho>=0.4 n", len(sel), 164, 0); C("19p in sel", sel.on_19p13_3.sum(), 62, 0); C("OR", orr, 82, 0.5); CL("Fisher p", pv, 1.1e-81, 0.3)
C("19p frac>=0.3 %", (ann[ann.on_19p13_3].rho >= 0.3).mean() * 100, 44, 0.5); C("others frac>=0.3 %", (ann[~ann.on_19p13_3].rho >= 0.3).mean() * 100, 3.2, 0.05)
c19 = ann[ann.chr.astype(str) == "19"]; C("median 19p13.3", ann[ann.on_19p13_3].rho.median(), 0.24); C("median 19p", c19[c19.start < 24.2e6].rho.median(), 0.15); C("median 19q", c19[c19.start >= 26.2e6].rho.median(), -0.01); C("median rest genome", ann[ann.chr.astype(str) != "19"].rho.median(), -0.03)
C("excl input n", ((co.rho >= 0.4) & (~co.on_19p13_3)).sum(), 102, 0); C("neg input n", (co.rho <= -0.35).sum(), 96, 0)
e31 = R("31_enrichment_enrichr_pos_excl19p13.3.csv"); e11 = R("11_enrichment_enrichr_pos.csv"); e21 = R("21_enrichment_enrichr_neg.csv")
def padj(e, lib, t): return e[(e.lib == lib) & (e.term.str.startswith(t))].padj.iloc[0]
for lib, t, x in [("GO_Cellular_Component_2023", "Mitochondrial Membrane", 2.3e-12), ("GO_Cellular_Component_2023", "Mitochondrial Inner Membrane", 2.3e-12), ("GO_Cellular_Component_2023", "Mitochondrial Matrix", 5.9e-3), ("GO_Cellular_Component_2023", "Mitochondrial Ribosome", 2.6e-3), ("GO_Biological_Process_2023", "Translation (", 2.3e-12), ("GO_Biological_Process_2023", "Cytoplasmic Translation", 3.6e-6), ("GO_Biological_Process_2023", "Mitochondrial Translation", 8.4e-5), ("GO_Biological_Process_2023", "Mitochondrial Gene Expression", 1.0e-4), ("GO_Biological_Process_2023", "Aerobic Respiration", 8.7e-4), ("GO_Biological_Process_2023", "Cellular Respiration", 4.5e-3), ("KEGG_2021_Human", "Ribosome", 3.3e-8), ("KEGG_2021_Human", "Oxidative phosphorylation", 0.076)]:
    CL("excl " + t[:22], padj(e31, lib, t), x)
CL("all OXPHOS padj", padj(e11, "KEGG_2021_Human", "Oxidative phosphorylation"), 1.0e-3); C("all OXPHOS n", e11[(e11.lib == "KEGG_2021_Human") & (e11.term == "Oxidative phosphorylation")].n_overlap.iloc[0], 8, 0)
C("OXPHOS 4 named genes on 19p13.3 & rho>=0.4", sum(bool(co[(co.gene == g)].on_19p13_3.iloc[0]) and co[co.gene == g].rho.iloc[0] >= 0.4 for g in ["NDUFA11", "NDUFS7", "UQCR11", "ATP5F1D"]), 4, 0)
CL("neg Ras padj", padj(e21, "KEGG_2021_Human", "Ras signaling"), 3.0e-4)
# ---- epithelium-intrinsic ----
g37 = R("37_epithelial_pseudobulk_genomewide_rho.csv").set_index("gene"); b19 = g37.on_19p13_3.astype(object).fillna(False).astype(bool)
C("epi genes", len(g37), 14140, 0); C("epi median", g37.rho_pb.median(), 0.28); C("epi q1", g37.rho_pb.quantile(.25), 0.10); C("epi q3", g37.rho_pb.quantile(.75), 0.43)
C("epi 19p median", g37.rho_pb[b19].median(), 0.26); C("epi other median", g37.rho_pb[~b19].median(), 0.28)
for g, (r_, pc) in {"TRAP1": (0.49, 85), "DNAJA3": (0.47, 82), "PMPCB": (0.46, 81), "CLPP": (0.45, 79), "HSPA9": (0.45, 78), "PHB": (0.41, 71), "YME1L1": (0.41, 73), "HSPD1": (0.24, 44), "HSPE1": (0.31, 55), "PHB2": (0.22, 40), "ATF5": (0.47, 81), "ATF4": (0.16, 32), "DDIT3": (0.00, 16), "POLRMT": (0.62, 97), "PPARGC1A": (0.45, None), "TFAM": (0.44, None), "FIS1": (0.51, None), "OPA1": (0.47, None), "DNM1L": (0.31, None), "MFF": (0.28, None), "MFN1": (0.20, None), "MFN2": (0.05, None), "THOP1": (None, 43), "TIMM13": (None, 39)}.items():
    if r_ is not None: C("epi rho " + g, g37.rho_pb[g], r_)
    if pc is not None: C("epi pct " + g, g37.pct_rank_pb[g] * 100, pc, 0.5)
e38 = R("38_enrichment_epithelial_intrinsic_top300.csv"); raw = e38[e38.set == "top300_raw"]; cen = e38[e38.set == "top300_sample_centered"]
CL("epi complex assembly", padj(raw, "GO_Biological_Process_2023", "Mitochondrial Respiratory Chain Complex Assembly"), 9.6e-8); CL("epi mito membrane", padj(raw, "GO_Cellular_Component_2023", "Mitochondrial Membrane"), 6.5e-5); CL("epi inner membrane", padj(raw, "GO_Cellular_Component_2023", "Mitochondrial Inner Membrane"), 7.6e-3)
C("epi KEGG none", float(raw[raw.lib == "KEGG_2021_Human"].padj.min() > 0.05), 1, 0); CL("cen inner membrane", padj(cen, "GO_Cellular_Component_2023", "Mitochondrial Inner Membrane"), 1.7e-4); CL("cen complex assembly", padj(cen, "GO_Biological_Process_2023", "Mitochondrial Respiratory Chain Complex Assembly"), 1.4e-5)
common = g37.bulk_rho_TCGA.dropna().index; C("bulk-epi agreement", spearmanr(g37.loc[common, "bulk_rho_TCGA"], g37.loc[common, "rho_pb"])[0], 0.26); C("agreement n", len(common), 13385, 0)
C("units", R("34_epithelial_pseudobulk_units_extended.csv").shape[0], 58, 0); C("bulk q1", co.rho.quantile(.25), -0.13); C("bulk q3", co.rho.quantile(.75), 0.08)
# ---- purity ----
d42 = R("42_purity_merged.csv"); C("purity n", len(d42), 291, 0); r, p = spearmanr(d42.purity, d42.LONP1_log2); C("purity rho", r, 0.12); C("purity p", p, 0.045, 0.001)
C("purity median", d42.purity.median(), 0.67); C("purity q1", d42.purity.quantile(.25), 0.54); C("purity q3", d42.purity.quantile(.75), 0.79)
mm = d42.dropna(subset=["epi_score"]); r, p = spearmanr(mm.purity, mm.epi_score); C("epi~purity rho", r, 0.13); C("epi~purity p", p, 0.025, 0.001)
def resid(y, X): X = sm.add_constant(X); return y - sm.OLS(y, X).fit().predict(X)
mc2 = d42.dropna(subset=["log2CNA"]); r, p = spearmanr(resid(mc2.LONP1_log2, mc2[["purity"]]), resid(mc2.log2CNA, mc2[["purity"]])); C("CNA partial rho", r, 0.45); CL("CNA partial p", p, 1.2e-15)
ols = sm.OLS(mc2.LONP1_log2, sm.add_constant(mc2[["log2CNA", "purity"]])).fit(); C("OLS beta CNA", ols.params["log2CNA"], 0.97); CL("OLS p CNA", ols.pvalues["log2CNA"], 1.1e-19)
for sig, (r_, p_) in {"M1_sig": (-0.22, 1.6e-4), "DC_sig": (-0.19, 1.0e-3), "CD4_sig": (-0.13, 0.026)}.items():
    r, p = spearmanr(resid(d42.LONP1_log2, d42[["purity"]]), resid(d42[sig], d42[["purity"]])); C(f"{sig} partial rho", r, r_); (CL if p_ < 0.01 else C)(f"{sig} partial p", p, p_, 0.05 if p_ < 0.01 else 0.001)
C("DC~purity", spearmanr(d42.purity, d42.DC_sig)[0], -0.54); C("CD4~purity", spearmanr(d42.purity, d42.CD4_sig)[0], -0.58)
# ---- GSE44001 ----
from lifelines import CoxPHFitter
d44 = R("44_GSE44001_LONP1_DFS.csv"); x = d44.dropna(subset=["DFS_months", "DFS_event"]); x = x[x.DFS_months > 0]
C("GSE44001 n", len(x), 300, 0); C("GSE44001 events", x.DFS_event.sum(), 38, 0)
cf = CoxPHFitter().fit(x[["DFS_months", "DFS_event", "LONP1_log2"]], "DFS_months", "DFS_event"); s = cf.summary.iloc[0]
C("GSE44001 HR", s["exp(coef)"], 1.07); C("GSE44001 lo", s["exp(coef) lower 95%"], 0.72); C("GSE44001 hi", s["exp(coef) upper 95%"], 1.60); C("GSE44001 p", s["p"], 0.73, 0.006)
med = x.LONP1_log2.median(); hi = x.LONP1_log2 > med; C("GSE44001 logrank p", logrank_test(x[hi].DFS_months, x[~hi].DFS_months, x[hi].DFS_event, x[~hi].DFS_event).p_value, 0.83, 0.006)
x = x.assign(stage_adv=(~x.stage.isin(["IA1", "IA2", "IB1"])).astype(int)); cf2 = CoxPHFitter().fit(x[["DFS_months", "DFS_event", "LONP1_log2", "stage_adv", "diameter_cm"]], "DFS_months", "DFS_event"); s2 = cf2.summary
C("GSE44001 diam HR", s2.loc["diameter_cm", "exp(coef)"], 1.29); C("diam lo", s2.loc["diameter_cm", "exp(coef) lower 95%"], 1.10); C("diam hi", s2.loc["diameter_cm", "exp(coef) upper 95%"], 1.52); C("diam p", s2.loc["diameter_cm", "p"], 0.002, 0.0006)
C("MV LONP1 HR (S2 legend)", s2.loc["LONP1_log2", "exp(coef)"], 1.03); C("MV stage HR", s2.loc["stage_adv", "exp(coef)"], 1.96); C("MV stage p", s2.loc["stage_adv", "p"], 0.063, 0.001)

bad = [c for c in chk if not c[3]]
lines = [f"Full manuscript number check ({len(chk)} items): {len(bad)} mismatches"] + [f"  FAIL {c[0]}: computed={c[1]} manuscript={c[2]}" for c in bad] + ["", "All checked items:"] + [f"  {'OK ' if c[3] else 'BAD'} {c[0]}: {c[1]} vs {c[2]}" for c in chk]
open(os.path.join(rd, "46_manuscript_number_check.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines[:1 + len(bad)]))
