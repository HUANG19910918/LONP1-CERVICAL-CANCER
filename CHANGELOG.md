# Changelog

## V1.4.0 (2026-11-20)
Added two analyses to `supplementary_analyses/` testing whether the LONP1 co-expression program is specific to LONP1:
- `06_control_genes/`: comparison of LONP1 with matched 19p13.3 controls, mitochondrial genes on other chromosomes and random expression-matched genes (Table S16). Model B uses the same 240 tumors as `02_partial_correlation` and reproduces its LONP1 correlations.
- `07_generalization/`: the same pipeline applied to 251 seed genes, including 30 genes from single-gene cervical cancer studies, with prespecified positional, retention and specificity readouts (Tables S17A–B).
Existing analyses (V1.3.0) are unchanged.

## V1.3.0 (2026-10-14)
Added `supplementary_analyses/` for the manuscript revision:
- Multi-cohort single-cell analysis (88 samples from GSE208653, E-MTAB-11948, GSE197461, GSE236738, GSE173231 and GSE333892): QC, inferCNV malignant-cell classification with CopyKAT cross-check, pseudobulk LONP1, 19p13.3 copy-gain comparison and normal–SIL–cancer comparison (Fig. 2H–K; Figs. S5–S6; Table S10).
- Copy-number- and purity-adjusted partial Spearman co-expression with fgsea (Fig. 4H–J).
- cis-eQTL Mendelian randomization and colocalization for LONP1 and 15 partner genes, with an R MendelianRandomization cross-check (Figs. S7–S8; Tables S11–S12).
- Assembly script for the new figure panels and Additional files, and the Supplementary Methods text.
Existing analyses (V1.2.0) are unchanged.

## V1.2.0 (2026-10-02)
Enrichr custom-background sensitivity analysis; PHB1 epithelial correlations; CITATION/README updates.
