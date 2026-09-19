# patch_v3.py —— 由 build_manuscript.js (v2) 生成 build_manuscript_v3.js（LONP1中心框架改写）
# 只替换指定段落（按行首前缀定位），其余内容原样保留；改写依据见 实验记录.md R15
import sys, re, os
src = sys.argv[1]; dst = sys.argv[2]
lines = open(src, encoding="utf-8").read().split("\n")

def js(s):  # 转成 JS 双引号字符串内容
    return s.replace("\\", "\\\\").replace('"', '\\"')

def replace_line(prefix, new_text, kind="P"):
    """prefix: 行内唯一前缀（去掉缩进后）; new_text: 新段落纯文本"""
    hits = [i for i, l in enumerate(lines) if l.strip().startswith(prefix)]
    assert len(hits) == 1, (prefix, hits)
    i = hits[0]
    indent = lines[i][:len(lines[i]) - len(lines[i].lstrip())]
    if kind == "P":
        lines[i] = f'{indent}P("{js(new_text)}"),'
    elif kind == "H2":
        lines[i] = f'{indent}H2("{js(new_text)}"),'
    elif kind == "TITLE":
        lines[i] = (f'{indent}new Paragraph({{ children: [new TextRun({{ text: "{js(new_text)}", size: 32, bold: true, '
                    f'font: "Times New Roman" }})], spacing: {{ after: 240, line: 360 }} }}),')

def insert_after(prefix, new_lines):
    hits = [i for i, l in enumerate(lines) if l.strip().startswith(prefix)]
    assert len(hits) == 1, (prefix, hits)
    i = hits[0]
    indent = lines[i][:len(lines[i]) - len(lines[i].lstrip())]
    for k, (kind, txt) in enumerate(new_lines):
        lines.insert(i + 1 + k, f'{indent}{kind}("{js(txt)}"),')

# ---------------- 标题 ----------------
replace_line('new Paragraph({ children: [new TextRun({ text: "Bulk and single-cell',
    "Multi-resolution transcriptomic profiling of LONP1 in cervical cancer: epithelial-enriched expression, "
    "19p13.3 copy-number-driven variation and a mitochondrial proteostasis programme", kind="TITLE")

# ---------------- 摘要 ----------------
replace_line('P("Background: Lon peptidase 1',
    "Background: Lon peptidase 1 (LONP1) is an ATP-dependent protease of the mitochondrial matrix that maintains "
    "mitochondrial proteostasis, supports respiratory chain function and is overexpressed in several cancers. Its "
    "expression pattern, genomic determinants and transcriptional context in cervical cancer have not been "
    "systematically characterized.")
replace_line('P("Methods: We integrated RNA-seq data',
    "Methods: We integrated RNA-seq data from TCGA-CESC and GTEx (Toil recompute; 304 tumors, 13 normal cervical "
    "tissues) with three independent microarray cohorts (GSE9750, GSE7803, GSE63514; probe 209017_s_at), single-cell "
    "RNA-seq of the cervical carcinogenesis spectrum (GSE208653: 80,435 cells from 9 samples spanning HPV-negative/"
    "positive normal cervix, HSIL and cancer), and TCGA copy-number, mutation and clinical data from cBioPortal. "
    "Genome-wide co-expression (Spearman) was annotated by genomic position to separate positional (19p13.3-linked) "
    "from functional co-expression before Gene Ontology/KEGG over-representation analysis (Enrichr); immune "
    "marker-signature correlations and survival analyses (Kaplan-Meier, Cox regression) were also performed. All "
    "analysis scripts and per-sample source data are available for reproduction.")
replace_line('P("Results: LONP1 mRNA was higher',
    "Results: LONP1 mRNA was higher in cervical cancer than in normal cervix in bulk-tissue cohorts (TCGA+GTEx log2 "
    "fold change 0.47, Mann-Whitney p=7.2x10-4; GSE9750 p=2.3x10-5; GSE7803 trend p=0.052), but not in "
    "laser-capture-microdissected epithelium (GSE63514, p=0.98). Single-cell analysis resolved this discordance: "
    "LONP1 was expressed predominantly by epithelial cells (mean log1p expression 0.143 versus 0.050-0.127 in stromal "
    "and immune populations), whereas malignant epithelial cells did not show higher LONP1 than normal epithelial "
    "cells (0.113 versus 0.136), indicating that the bulk-level elevation largely reflects the higher epithelial "
    "content of tumors rather than cell-intrinsic transcriptional induction. Within tumors, LONP1 expression "
    "correlated with copy number (rho=0.41, p=3.1x10-13) and increased stepwise across GISTIC classes, while somatic "
    "mutations were rare (1/294). Genome-wide co-expression was dominated by a positional signal: 18 of the 20 "
    "strongest partners, including CLPP (rho=0.72), map to 19p13.3, and 19p13.3 genes were 82-fold over-represented "
    "among genes with rho>=0.4 (Fisher p=1.1x10-81). After excluding 19p13.3 genes, co-expressed genes remained "
    "enriched for mitochondrial inner membrane (adjusted p=2.3x10-12), mitochondrial translation and aerobic "
    "respiration, and LONP1 correlated with mitochondrial chaperones and stress-response factors located elsewhere "
    "in the genome (TRAP1 rho=0.41, PHB1/PHB2 0.40/0.39, HSPD1 0.33, ATF4 0.32, HSPE1 0.30, HSPA9 0.23); the CLPP and "
    "HSPA9 associations replicated within epithelial cells at single-cell resolution. Correlations with mitochondrial "
    "fission/fusion genes were weak (|rho|<=0.32, DNM1L rho=-0.02) and did not indicate a fission-specific "
    "programme within epithelium. LONP1 showed weak inverse correlations with M1 macrophage and dendritic cell "
    "signatures (FDR=0.0012) and was not associated with overall, disease-specific or progression-free survival "
    "(all p>0.4).")
replace_line('P("Conclusions: LONP1 is an epithelial-enriched gene',
    "Conclusions: LONP1 is an epithelial-enriched gene whose apparent overexpression in bulk cervical cancer tissue "
    "is largely attributable to tissue composition. Within tumors, its expression variation is driven by 19p13.3 "
    "copy-number gain, which also dominates its strongest co-expression partners; beyond this positional effect, "
    "LONP1 is embedded in a mitochondrial proteostasis and translation programme, without prognostic value. These "
    "findings refine the interpretation of bulk-tissue and co-expression-based biomarker claims and provide a "
    "rigorous transcriptome-level foundation for functional studies of LONP1 in cervical cancer.")
replace_line('P("Keywords:',
    "Keywords: LONP1; cervical cancer; mitochondrial proteostasis; mitochondrial unfolded protein response; copy "
    "number; co-expression; single-cell RNA-seq; TCGA")

# ---------------- Background ----------------
replace_line('P("Mitochondria are central hubs',
    "Mitochondria are central hubs of energy metabolism, redox signalling and apoptosis, and mitochondrial "
    "reprogramming is increasingly recognized as a hallmark-enabling feature of cancer [4, 5, 30]. Mitochondrial "
    "homeostasis depends on a protein quality-control network in which Lon peptidase 1 (LONP1), an AAA+ ATP-dependent "
    "protease of the mitochondrial matrix, degrades oxidized and misfolded proteins, regulates mitochondrial DNA "
    "maintenance, and supports assembly of oxidative phosphorylation (OXPHOS) complexes [6, 7, 9]. LONP1 is itself a "
    "canonical target of the mitochondrial unfolded protein response (UPRmt), together with the protease CLPP and the "
    "chaperones HSPD1, HSPE1 and HSPA9, downstream of the integrated stress response factors ATF4 and ATF5 [34-36]. "
    "LONP1 is stress-inducible and has been implicated in tumor bioenergetics and aggressive phenotypes in several "
    "cancer types [8, 10, 11, 29].")
replace_line('P("Mitochondrial morphology is governed',
    "Mitochondrial morphology is governed by a balance between fission, mediated by DNM1L/DRP1 with its receptors "
    "FIS1 and MFF, and fusion, mediated by MFN1, MFN2 and OPA1 [12-14]; shifts toward fission have been associated "
    "with proliferation and invasion in several tumors [15], and whether LONP1 expression is transcriptionally "
    "coupled to this machinery in cervical cancer is unknown. A further, often neglected, consideration for "
    "single-gene co-expression studies is genomic position: because copy-number alterations affect contiguous genes, "
    "a copy-number-driven gene co-expresses with its chromosomal neighbours irrespective of function [37-39]. A "
    "systematic description of LONP1 expression, its genomic determinants and its co-expression landscape in "
    "cervical cancer, with these confounders explicitly addressed, is lacking.")
replace_line('P("Here we address this gap',
    "Here we address this gap with a multi-cohort in silico study centred on LONP1. We quantified LONP1 expression "
    "in cervical cancer versus normal cervix across four independent datasets and at single-cell resolution, "
    "examined copy-number, mutational and histological correlates, mapped the genome-wide co-expression landscape "
    "while separating positional from functional components, characterized the mitochondrial proteostasis, "
    "biogenesis and dynamics genes associated with LONP1, assessed associations with immune cell marker signatures, "
    "and evaluated prognostic relevance. We deliberately report non-significant and discordant findings in full, "
    "including the absence of prognostic association, a negative epithelium-only comparison and the positional "
    "confounding of co-expression, to provide an unbiased reference for future functional work.")

# ---------------- Methods ----------------
replace_line('P("Filtered gene-barcode matrices of GSE208653',
    "Filtered gene-barcode matrices of GSE208653 (nine cervical tissue samples spanning HPV-negative normal (n=2), "
    "HPV-positive normal (n=2), HSIL (n=2) and cancer (n=3, two squamous and one adenocarcinoma); 10x Genomics 3' "
    "v3.1, Cell Ranger v5.0.1, GRCh38; mitochondrial genes removed and cells with >200 genes retained by the original "
    "authors) were downloaded from GEO [31]. Cells with fewer than 300 detected genes or more than 7,500 (potential "
    "doublets) were excluded, yielding 80,435 cells. Counts were normalized to 10,000 per cell and log1p-transformed "
    "in Scanpy [32]; 2,000 highly variable genes (per-sample batch key) were scaled and used for PCA (30 components), "
    "neighbourhood graph construction (k=15), Leiden clustering (resolution 0.5) [33] and UMAP. Clusters were assigned "
    "to nine major cell types by canonical markers (epithelial: EPCAM/KRT5/KRT8/KRT17/KRT18; T/NK: CD3D/CD3E/CD2/NKG7; "
    "B: MS4A1/CD79A; plasma: MZB1/IGHG1; myeloid: LYZ/CD68/AIF1; mast: TPSAB1/CPA3; endothelial: PECAM1/VWF; "
    "fibroblast: COL1A1/DCN; smooth muscle: ACTA2/MYH11). LONP1 was compared across cell types and, within epithelial "
    "cells, across disease stages per cell (Kruskal-Wallis/Mann-Whitney) and per sample (pseudobulk means), the "
    "latter reported descriptively given n=9 samples. To test whether bulk co-expression patterns are "
    "epithelium-intrinsic, epithelial cells were re-clustered (1,500 highly variable genes, 20 PCs, Leiden 0.5) and "
    "mean expression was computed per sample-by-subcluster pseudobulk unit (>=50 cells, n=63 units); Spearman "
    "correlations between LONP1 and mitochondrial proteostasis (CLPP, HSPD1, HSPA9) and dynamics genes (DNM1L, FIS1, "
    "MFF, MFN1, MFN2, OPA1) were then compared with the bulk estimates. As a bulk-level consistency check of the "
    "composition hypothesis, an epithelial marker score (mean z-score of EPCAM/KRT5/KRT8/KRT17) was computed for each "
    "TCGA+GTEx sample and LONP1 was modelled by ordinary least squares as a function of tumor status with and without "
    "adjustment for this score.")
replace_line('H2("Co-expression and enrichment")', "Co-expression, positional annotation and enrichment", kind="H2")
replace_line('P("Genome-wide co-expression with LONP1 in TCGA-CESC',
    "Genome-wide co-expression with LONP1 in TCGA-CESC (n=294) was computed server-side by cBioPortal (Spearman "
    "correlation, no threshold); of 20,337 genetic entities returned, 17,378 genes with a HUGO symbol were retained. "
    "Because LONP1 maps to 19p13.3 and its expression tracks copy number, co-expression with neighbouring genes may "
    "reflect shared copy-number variation rather than functional co-regulation [37-39]. Gene coordinates (GRCh38) were "
    "therefore obtained from the Ensembl REST API and cytogenetic bands from the UCSC hg38 cytoBand table, and every "
    "gene was flagged as located on 19p13.3 or not. Over-representation of 19p13.3 genes among co-expressed genes was "
    "tested with Fisher's exact test, and the positional structure of co-expression was visualized along chromosome "
    "19 with a 1-Mb sliding median. Genes with rho>=0.4 (n=164) were tested for over-representation against GO "
    "Biological Process 2023, GO Cellular Component 2023 and KEGG 2021 Human using Enrichr with Benjamini-Hochberg "
    "adjustment [26], both including all genes and after excluding 19p13.3 genes (n=102); the latter is presented as "
    "the primary result. Negatively correlated genes (rho<=-0.35, n=96) were tested in the same way. A curated panel "
    "of mitochondrial genes - matrix proteases (CLPP, AFG3L2, YME1L1, SPG7, PMPCB), chaperones (HSPD1, HSPE1, HSPA9, "
    "TRAP1, DNAJA3), UPRmt/integrated stress response transcription factors (ATF4, ATF5, DDIT3), prohibitins (PHB1, "
    "PHB2), biogenesis/mtDNA factors (TFAM, NRF1, PPARGC1A, POLRMT), fission (DNM1L, FIS1, MFF, MIEF1, MIEF2), fusion "
    "(MFN1, MFN2, OPA1) and mitophagy genes (PINK1, PRKN) - was examined individually using the same server-side "
    "correlations, which were verified against local computation on log2(RSEM+1) values.")
replace_line('P("All statistical analyses and figures were produced in R',
    "Statistical analyses were performed in R 4.3.3 (survival) and Python 3.10 (pandas, SciPy, lifelines 0.30, "
    "Scanpy 1.11), with key results cross-checked between the two environments; figures were drawn with plotnine "
    "0.15 (a ggplot2 implementation) and matplotlib. All tests were two-sided with alpha=0.05. Analysis scripts, "
    "per-sample source data extracted from the public repositories, and a step-by-step verification report are "
    "retained by the authors and available on request to allow full reproduction of every reported value.")

# ---------------- Results: 共表达两节重写 ----------------
replace_line('H2("LONP1 co-expression marks a fission-oriented mitochondrial programme")',
    "Genome-wide co-expression of LONP1 is dominated by a positional 19p13.3 signal", kind="H2")
replace_line('P("Across curated mitochondrial dynamics and quality-control genes (Fig. 3A)',
    "Genome-wide, 970 of 17,378 genes correlated with LONP1 at |rho|>=0.3 (628 positive, 342 negative; Fig. 3A). The "
    "strongest partners were, however, almost exclusively chromosomal neighbours of LONP1: 18 of the 20 and 35 of the "
    "50 strongest positive partners - including CLPP (rho=0.72), THOP1 (0.68), HDGFL2 (0.67), WDR18 (0.67), SGTA "
    "(0.64), NDUFA11 (0.61) and TIMM13 (0.60) - map to 19p13.3, a band that contains only 1.1% of annotated genes "
    "(188/17,324). 19p13.3 genes were 82-fold over-represented among genes with rho>=0.4 (62/164; Fisher's exact "
    "p=1.1x10-81), and 44% of all 19p13.3 genes exceeded rho=0.3 compared with 3.2% of genes elsewhere. Along "
    "chromosome 19 the correlation decayed with distance from LONP1 (median rho 0.24 for 19p13.3 genes, 0.15 for the "
    "19p arm and -0.01 for 19q, versus -0.03 for the rest of the genome; Fig. 3B). Given the copy-number dependence of "
    "LONP1 expression (Fig. 2A, B), this pattern indicates that the top co-expression partners largely reflect "
    "co-variation of the 19p13.3 region rather than functional co-regulation, and we therefore separated positional "
    "from functional components in all subsequent analyses.")
replace_line('H2("Positively co-expressed genes are enriched for OXPHOS and mitochondrial gene expression")',
    "Beyond the locus, LONP1 co-expression marks a mitochondrial proteostasis and translation programme", kind="H2")
replace_line('P("Genome-wide, 970 genes passed |rho|>=0.3',
    "After excluding 19p13.3 genes, the remaining positively co-expressed genes (rho>=0.4, n=102) were "
    "over-represented in mitochondrial membrane and inner membrane (adjusted p=2.3x10-12 for both), mitochondrial "
    "matrix (5.9x10-3) and mitochondrial ribosome (2.6x10-3), in translation (2.3x10-12), cytoplasmic translation "
    "(3.6x10-6), mitochondrial translation (8.4x10-5), mitochondrial gene expression (1.0x10-4), aerobic respiration "
    "(8.7x10-4) and cellular respiration (4.5x10-3), and in the KEGG ribosome pathway (3.3x10-8; Fig. 3E). The KEGG "
    "oxidative phosphorylation term, which was significant when 19p13.3 genes were included (adjusted p=1.0x10-3; "
    "four of its eight overlapping genes - NDUFA11, NDUFS7, UQCR11 and ATP5F1D - lie on 19p13.3), was no longer "
    "significant after their exclusion (adjusted p=0.076), indicating that the OXPHOS signal was partly positional. "
    "Among curated mitochondrial genes located outside 19p13.3 (Fig. 3C), LONP1 correlated most strongly with the "
    "chaperones TRAP1 (rho=0.41, p=2.5x10-13), HSPD1 (0.33), HSPE1 (0.30) and HSPA9 (0.23), the prohibitins PHB1 and "
    "PHB2 (0.40 and 0.39), the processing peptidase PMPCB (0.22) and the integrated stress response transcription "
    "factors ATF4 (0.32), DDIT3 (0.22) and ATF5 (0.14, p=0.018) - a set that recapitulates the composition of the "
    "mammalian UPRmt [34-36]. Biogenesis regulators (TFAM, NRF1, PPARGC1A), the inner-membrane proteases YME1L1, "
    "AFG3L2 and SPG7, and the mitophagy genes PINK1 and PRKN were not or only weakly correlated. Correlations with the "
    "fission/fusion machinery were modest: FIS1 (0.32) and MFF (0.15) positive, MFN1 (-0.23), OPA1 (-0.22) and MFN2 "
    "(-0.13) negative, and no correlation for the fission GTPase DNM1L (-0.02). Negatively correlated genes "
    "(rho<=-0.35, n=96) were enriched for Ras signalling (adjusted p=3.0x10-4), endocytosis and vesicle trafficking "
    "(Fig. 3F).")
insert_after('P("After excluding 19p13.3 genes, the remaining positively', [
    ("P", "Because bulk co-expression can also be confounded by tissue composition, we re-examined the module within "
          "the epithelial compartment of the single-cell dataset (63 sample-by-subcluster pseudobulk units; Fig. 3D). "
          "The proteostasis associations replicated within epithelium: CLPP rho=0.49 (p=4.9x10-5), HSPA9 0.46 "
          "(p=1.8x10-4) and HSPD1 0.26 (p=0.041). The dynamics genes behaved differently: within epithelium every "
          "gene examined correlated positively or not at all with LONP1 - fission genes FIS1 (0.46), MFF (0.34) and "
          "DNM1L (0.31) as well as fusion genes OPA1 (0.42), MFN1 (0.19, n.s.) and MFN2 (-0.02, n.s.) - so the "
          "inverse fusion-gene associations seen in bulk tissue did not replicate, and the pattern is more consistent "
          "with a shared mitochondrial-content signal than with a fission-specific programme. We therefore regard the "
          "mitochondrial dynamics associations as exploratory."),
])

# ---------------- Discussion ----------------
replace_line('P("This study provides a systematic transcriptome-level portrait',
    "This study provides a systematic, multi-resolution transcriptomic portrait of LONP1 in cervical cancer. Four "
    "findings are robust across analyses. First, LONP1 is elevated in bulk tumor tissue relative to normal cervix in "
    "independent cohorts, but this elevation is largely a tissue-composition effect: LONP1 is epithelial-enriched, "
    "and two orthogonal epithelium-resolved datasets showed no cell-intrinsic induction in malignant epithelium. "
    "Second, intratumoral variation in LONP1 expression is driven by 19p13.3 copy-number gain rather than mutation. "
    "Third, the same copy-number effect dominates the strongest co-expression partners of LONP1, so that "
    "co-expression must be interpreted after accounting for genomic position. Fourth, beyond this positional effect, "
    "LONP1 is embedded in a coherent mitochondrial proteostasis and translation programme - chaperones "
    "HSPD1/HSPE1/HSPA9/TRAP1, prohibitins, ATF4/DDIT3, and mitoribosomal and inner-membrane genes - that replicates "
    "within epithelial cells. These observations align with the established biology of LONP1 as a stress-responsive "
    "protease that sustains respiratory chain assembly and mitochondrial proteostasis [6-9] and with reports that Lon "
    "overexpression promotes aggressive phenotypes through mitochondrial ROS and metabolic reprogramming "
    "[8, 10, 11, 29].")
replace_line('P("The association of LONP1 with fission receptors',
    "The positional structure of the co-expression landscape deserves emphasis because it is rarely accounted for in "
    "single-gene bioinformatic studies. Eighteen of the twenty strongest LONP1 partners, including CLPP (rho=0.72) "
    "and the OXPHOS subunits NDUFA11, NDUFS7 and UQCR11, lie on 19p13.3, and 19p13.3 genes were 82-fold "
    "over-represented among strongly co-expressed genes; correspondingly, the KEGG oxidative phosphorylation "
    "enrichment disappeared once these genes were removed. Co-expression driven by shared copy-number variation is "
    "well documented [37-39], and our data show that for a copy-number-driven gene it can occupy the entire top of "
    "the partner list; we suggest that co-expression-based inferences for such genes be reported with positional "
    "annotation. The functional programme that survives this correction - mitochondrial chaperones, prohibitins, "
    "ATF4/DDIT3 and mitochondrial translation - matches the composition of the mammalian UPRmt, in which LONP1, CLPP, "
    "HSPD1, HSPE1 and HSPA9 are canonical ATF4/ATF5-dependent targets [34-36], and it generates a testable hypothesis "
    "for the functional work motivating this project: that LONP1-high cervical cancer cells operate in a state of "
    "elevated mitochondrial proteostatic demand that co-activates the UPRmt and mitochondrial translation. "
    "Mitochondrial fission-fusion genes, by contrast, showed only weak bulk correlations, no correlation for DNM1L, "
    "and uniformly positive within-epithelium correlations for both fission and fusion genes, so the transcriptomic "
    "data do not by themselves support a fission-specific programme; any effect of LONP1 on mitochondrial morphology "
    "will need to be established at the protein and organelle level rather than inferred from mRNA co-expression.")
replace_line('P("Limitations. This is an in silico study',
    "Limitations. This is an in silico study of public data without new functional experiments; causality cannot be "
    "inferred from co-expression. The normal-tissue comparator in TCGA+GTEx is small (n=13) and platform-batch "
    "effects between TCGA and GTEx, though mitigated by the uniform Toil pipeline, cannot be excluded. Microarray "
    "validation relies on a single probe set. The single-cell cohort comprises nine samples (one adenocarcinoma among "
    "three cancers), sample-level inference is therefore descriptive, and 10x data suffer dropout for moderately "
    "expressed genes such as LONP1 (about 25% of epithelial cells positive); mRNA abundance also need not track "
    "protein level or protease activity. The immune analysis used compact marker signatures rather than full "
    "deconvolution. Finally, although we separated positional from functional co-expression by excluding 19p13.3 "
    "genes, more distant copy-number co-variation (for example arm-level 19p gain, which is visible as a gradient "
    "along 19p in Fig. 3B) is not fully removed by this approach; partial correlation conditioning on genome-wide "
    "copy number would be a refinement.")

# ---------------- Conclusions ----------------
replace_line('P("LONP1 is an epithelial-enriched gene whose bulk-tissue elevation',
    "LONP1 is an epithelial-enriched gene whose bulk-tissue elevation in cervical cancer is largely a "
    "tissue-composition effect rather than cell-intrinsic transcriptional induction. Within tumors, expression "
    "variation is driven principally by 19p13.3 copy-number gain, which also accounts for the strongest co-expression "
    "partners; once this positional effect is removed, LONP1 marks a coordinated mitochondrial proteostasis and "
    "translation programme that is epithelium-intrinsic, without prognostic value in current cohorts. These results "
    "provide a reproducible, multi-resolution transcriptomic foundation and clear hypotheses for functional studies of "
    "LONP1 in cervical cancer.")

# ---------------- 图注 ----------------
replace_line('P("Figure 3. LONP1 co-expression.',
    "Figure 3. The LONP1 co-expression landscape in TCGA-CESC. (A) Genome-wide Spearman correlation of all annotated "
    "genes (n=17,378) with LONP1, ranked from highest to lowest; genes on 19p13.3 in red; dashed lines, |rho|=0.3. "
    "(B) Correlation versus genomic position along chromosome 19 (GRCh38); shaded area, 19p13.3; dashed line, LONP1; "
    "blue line, 1-Mb sliding median. (C) Curated mitochondrial proteostasis, biogenesis and dynamics genes "
    "(server-side Spearman, n=294; *p<0.05, **p<0.01, ***p<0.001; genes on 19p13.3 labelled). (D) Bulk (TCGA-CESC) "
    "versus epithelium-intrinsic (GSE208653 sample-by-subcluster pseudobulk, n=63 units) correlations for the "
    "proteostasis module and the dynamics genes. (E) Over-representation of positively co-expressed genes (rho>=0.4) "
    "after excluding 19p13.3 genes (n=102; Enrichr, Benjamini-Hochberg adjusted; dashed line, adjusted p=0.05; the "
    "KEGG oxidative phosphorylation term is shown for comparison and is not significant). (F) Enrichment of "
    "negatively correlated genes (rho<=-0.35, n=96).")

# ---------------- 参考文献 34-39 ----------------
new_refs = [
 "Zhao Q, Wang J, Levichkin IV, Stasinopoulos S, Ryan MT, Hoogenraad NJ. A mitochondrial specific stress response in mammalian cells. EMBO J. 2002;21(17):4411-9.",
 "Fiorese CJ, Schulz AM, Lin YF, Rosin N, Pellegrino MW, Haynes CM. The transcription factor ATF5 mediates a mammalian mitochondrial UPR. Curr Biol. 2016;26(15):2037-43.",
 "Quiros PM, Prado MA, Zamboni N, D'Amico D, Williams RW, Finley D, et al. Multi-omics analysis identifies ATF4 as a key regulator of the mitochondrial stress response in mammals. J Cell Biol. 2017;216(7):2027-45.",
 "Reyal F, Stransky N, Bernard-Pierrot I, Vincent-Salomon A, de Rycke Y, Elvin P, et al. Visualizing chromosomes as transcriptome correlation maps: evidence of chromosomal domains containing co-expressed genes - a study of 130 invasive ductal breast carcinomas. Cancer Res. 2005;65(4):1376-83.",
 "Fehrmann RS, Karjalainen JM, Krajewska M, Westra HJ, Maloney D, Simeonov A, et al. Gene expression analysis identifies global gene dosage sensitivity in cancer. Nat Genet. 2015;47(2):115-25.",
 "Bhattacharya A, Bense RD, Urzua-Traslavina CG, de Vries EGE, van Vugt MATM, Fehrmann RSN. Transcriptional effects of copy number alterations in a large set of human cancers. Nat Commun. 2020;11(1):715.",
]
hits = [i for i, l in enumerate(lines) if l.strip() == "];"]
assert len(hits) >= 1
i = hits[0]  # refs 数组结束
if not lines[i-1].rstrip().endswith(","): lines[i-1] = lines[i-1].rstrip() + ","
for k, r in enumerate(new_refs):
    lines.insert(i + k, f'"{js(r)}",')

out = "\n".join(lines)
assert "fission-oriented" not in out, "残留旧表述"
open(dst, "w", encoding="utf-8").write(out)
print("PATCHED ->", dst, "refs added:", len(new_refs))
