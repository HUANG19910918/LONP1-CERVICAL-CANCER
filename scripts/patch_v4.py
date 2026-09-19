# patch_v4.py —— 由 build_manuscript_v3.js 生成 build_manuscript_v4.js
# 改动：上皮内在共表达重跑（scripts/17，58个假批量单元、全基因组、扩展模块）后的 Methods / Results / Abstract /
#       Discussion / Fig3图注 同步；负相关富集面板改为 Fig. 3G。依据见 实验记录.md R16
import sys
src, dst = sys.argv[1], sys.argv[2]
lines = open(src, encoding="utf-8").read().split("\n")
def js(s): return s.replace("\\", "\\\\").replace('"', '\\"')
def find(prefix):
    hits = [i for i, l in enumerate(lines) if l.strip().startswith(prefix)]
    assert len(hits) == 1, (prefix, hits); return hits[0]
def replace_line(prefix, new_text, kind="P"):
    i = find(prefix); indent = lines[i][:len(lines[i]) - len(lines[i].lstrip())]
    lines[i] = f'{indent}{kind}("{js(new_text)}"),'
def sub_in_line(prefix, old, new, count=1):
    i = find(prefix); assert old in lines[i], (prefix, old); lines[i] = lines[i].replace(old, new, count)

# ---------- Abstract ----------
sub_in_line('P("Results: LONP1 mRNA was higher',
    "the CLPP and HSPA9 associations replicated within epithelial cells at single-cell resolution. Correlations with "
    "mitochondrial fission/fusion genes were weak (|rho|<=0.32, DNM1L rho=-0.02) and did not indicate a "
    "fission-specific programme within epithelium.",
    "within epithelial cells (58 pseudobulk units), LONP1 partners were enriched for mitochondrial respiratory chain "
    "complex assembly and mitochondrial membrane genes without any positional bias, with CLPP, HSPA9 and TRAP1 in "
    "the top quintile of genome-wide correlations, whereas ATF4 and DDIT3 did not track LONP1. Correlations with "
    "mitochondrial fission/fusion genes were weak in bulk tissue (|rho|<=0.32, DNM1L rho=-0.02), and within "
    "epithelium fission and fusion genes correlated alike with LONP1, arguing against a fission-specific programme.")

# ---------- Methods: single-cell ----------
sub_in_line('P("Filtered gene-barcode matrices of GSE208653',
    "mean expression was computed per sample-by-subcluster pseudobulk unit (>=50 cells, n=63 units); Spearman "
    "correlations between LONP1 and mitochondrial proteostasis (CLPP, HSPD1, HSPA9) and dynamics genes (DNM1L, FIS1, "
    "MFF, MFN1, MFN2, OPA1) were then compared with the bulk estimates.",
    "mean log-normalized expression was computed per sample-by-subcluster pseudobulk unit (>=50 cells; the 25 "
    "subclusters yielded 58 qualifying units). Epithelium-intrinsic co-expression was then computed genome-wide as the Spearman correlation "
    "between LONP1 and every gene detected in >=1% of epithelial cells (n=14,140) across these units. Because "
    "unit-level means share a global component (unit-to-unit differences in transcript detection shift most "
    "correlations upward), each gene's correlation is reported together with its percentile rank in the genome-wide "
    "distribution, and a sensitivity analysis centred every unit on its sample mean to remove between-sample (stage "
    "and copy-number) differences. The 300 genes most positively correlated with LONP1 within epithelium were tested "
    "for over-representation with Enrichr, and the presence of a positional (19p13.3) effect was re-examined at this "
    "level. Correlations for the curated mitochondrial genes were compared with the bulk estimates.")

# ---------- Results: epithelium-intrinsic paragraph ----------
replace_line('P("Because bulk co-expression can also be confounded by tissue composition',
    "Because bulk co-expression can also be confounded by tissue composition and by copy-number linkage, we "
    "re-examined co-expression within the epithelial compartment of the single-cell dataset (58 sample-by-subcluster "
    "pseudobulk units; Fig. 3D, F). At this level the positional signal disappeared entirely: 19p13.3 genes had a "
    "median correlation with LONP1 of 0.26 versus 0.28 for all other genes, and the two strongest bulk partners after "
    "CLPP, THOP1 and TIMM13, fell to the 43rd and 39th percentiles of the genome-wide distribution. Unit-level "
    "correlations were globally shifted upward (genome-wide median 0.28, interquartile range 0.10-0.43), so we "
    "interpret individual genes by their percentile rank. The proteostasis module remained among the strongest "
    "partners: TRAP1 (rho=0.49, 85th percentile), DNAJA3 (0.47, 82nd), PMPCB (0.46, 81st), CLPP (0.45, 79th), HSPA9 "
    "(0.45, 78th), PHB1 (0.41, 71st) and YME1L1 (0.41, 73rd), whereas HSPD1 (0.24, 44th), HSPE1 (0.31, 55th) and "
    "PHB2 (0.22, 40th) were unremarkable. The stress-response transcription factors diverged: ATF5 tracked LONP1 "
    "(0.47, 81st) but ATF4 (0.16, 32nd) and DDIT3 (0.00, 16th) did not. Biogenesis regulators that showed no bulk "
    "association were correlated within epithelium (POLRMT 0.62, 97th; PPARGC1A 0.45; TFAM 0.44), and fission and "
    "fusion genes behaved alike (FIS1 0.51, OPA1 0.47, DNM1L 0.31, MFF 0.28, MFN1 0.20, MFN2 0.05), so the inverse "
    "fusion-gene associations of bulk tissue did not replicate. Genome-wide, the 300 strongest epithelium-intrinsic "
    "partners were enriched for mitochondrial respiratory chain complex assembly (adjusted p=9.6x10-8), "
    "mitochondrial membrane (6.5x10-5), mitochondrial inner membrane (7.6x10-3) and protein transport, with no "
    "significant KEGG term; sample-centred correlations gave the same picture (mitochondrial inner membrane "
    "1.7x10-4, respiratory chain complex assembly 1.4x10-5). Within epithelium, therefore, LONP1 marks a broad "
    "mitochondrial content and biogenesis state of the cell - proteostasis, respiratory chain assembly, mtDNA "
    "transcription and both arms of the dynamics machinery together - rather than a fission-specific or "
    "ISR-driven programme, and the agreement between bulk and epithelium-intrinsic correlations across all genes "
    "was modest (Spearman 0.26, n=13,385).")
# negative enrichment now panel G
sub_in_line('P("After excluding 19p13.3 genes, the remaining positively', "(Fig. 3F).", "(Fig. 3G).")

# ---------- Discussion ----------
sub_in_line('P("This study provides a systematic, multi-resolution transcriptomic portrait',
    "Fourth, beyond this positional effect, LONP1 is embedded in a coherent mitochondrial proteostasis and "
    "translation programme - chaperones HSPD1/HSPE1/HSPA9/TRAP1, prohibitins, ATF4/DDIT3, and mitoribosomal and "
    "inner-membrane genes - that replicates within epithelial cells.",
    "Fourth, beyond this positional effect, LONP1 is embedded in a coherent mitochondrial proteostasis and "
    "translation programme - chaperones HSPD1/HSPE1/HSPA9/TRAP1, prohibitins, ATF4/DDIT3, and mitoribosomal and "
    "inner-membrane genes - whose mitochondrial core (respiratory chain assembly and membrane genes; CLPP, HSPA9, "
    "TRAP1) replicates within epithelial cells as an unbiased genome-wide signal that is free of positional bias.")
sub_in_line('P("The positional structure of the co-expression landscape deserves emphasis',
    "and it generates a testable hypothesis for the functional work motivating this project: that LONP1-high "
    "cervical cancer cells operate in a state of elevated mitochondrial proteostatic demand that co-activates the "
    "UPRmt and mitochondrial translation.",
    "and it generates a testable hypothesis for the functional work motivating this project: that LONP1-high "
    "cervical cancer cells operate in a state of elevated mitochondrial proteostatic demand that co-activates the "
    "UPRmt and mitochondrial translation. The single-cell data qualify this hypothesis in an informative way. "
    "Within epithelial cells, the proteostasis genes ranked among the strongest LONP1 partners, but so did "
    "biogenesis regulators (POLRMT, TFAM, PPARGC1A) and both fission and fusion genes, whereas ATF4 and DDIT3 did "
    "not track LONP1 at all; the epithelium-intrinsic signal therefore resembles a mitochondrial content and "
    "biogenesis state of the cell rather than a specific ISR-driven induction, and the ATF4/DDIT3 component seen in "
    "bulk tissue is more plausibly a tumour-level feature. Functional experiments should accordingly measure "
    "mitochondrial mass and respiratory chain assembly alongside UPRmt markers when LONP1 is perturbed.")
sub_in_line('P("Limitations. This is an in silico study',
    "sample-level inference is therefore descriptive,",
    "sample-level inference is therefore descriptive, the epithelium-intrinsic co-expression analysis rests on 58 "
    "pseudobulk units whose correlations are globally inflated by unit-level detection differences (hence our use "
    "of percentile ranks),")

# ---------- Figure 3 legend ----------
replace_line('P("Figure 3. The LONP1 co-expression landscape in TCGA-CESC.',
    "Figure 3. The LONP1 co-expression landscape in TCGA-CESC and within cervical epithelium. (A) Genome-wide "
    "Spearman correlation of all annotated genes (n=17,378) with LONP1, ranked from highest to lowest; genes on "
    "19p13.3 in red; dashed lines, |rho|=0.3. (B) Correlation versus genomic position along chromosome 19 (GRCh38); "
    "shaded area, 19p13.3; dashed line, LONP1; blue line, 1-Mb sliding median. (C) Curated mitochondrial "
    "proteostasis, biogenesis and dynamics genes (server-side Spearman, n=294; *p<0.05, **p<0.01, ***p<0.001; genes "
    "on 19p13.3 labelled). (D) Bulk (TCGA-CESC, n=294) versus epithelium-intrinsic (GSE208653, 58 "
    "sample-by-subcluster pseudobulk units) correlations for the curated genes; shaded band and dotted line, "
    "interquartile range and median of the epithelium-intrinsic correlations of all 14,140 detected genes (the "
    "corresponding bulk interquartile range is -0.13 to 0.08). (E) Over-representation of positively co-expressed "
    "bulk genes (rho>=0.4) after excluding 19p13.3 genes (n=102; Enrichr, Benjamini-Hochberg adjusted; dashed line, "
    "adjusted p=0.05; the KEGG oxidative phosphorylation term is shown for comparison and is not significant). (F) "
    "Over-representation of the 300 genes most positively correlated with LONP1 within epithelium (no KEGG term "
    "reached significance). (G) Enrichment of negatively correlated bulk genes (rho<=-0.35, n=96).")

out = "\n".join(lines)
assert "n=63 units" not in out and "Fig. 3F)." not in out.replace("Fig. 3D, F)", "")
open(dst, "w", encoding="utf-8").write(out)
print("PATCHED ->", dst)
