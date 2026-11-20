# Supplementary analyses (V1.4.0)

Code and derived result tables for the analyses added in revision (01–04 in V1.3.0; 06–07 in V1.4.0):

| Folder | Analysis | Manuscript |
|---|---|---|
| `01_scRNA_inferCNV/` | Multi-cohort single-cell RNA-seq (88 samples, 6 cohorts): QC, inferCNV malignant-cell calls, CopyKAT cross-check, pseudobulk LONP1, 19p13.3 gain, GSE333892 normal–SIL–cancer | Fig. 2H–K; Additional file 2: Figs. S5–S6; Table S10 |
| `02_partial_correlation/` | Copy-number- and purity-adjusted (partial Spearman) LONP1 co-expression in TCGA-CESC; fgsea | Fig. 4H–J |
| `03_MR_coloc/` | cis-eQTL Mendelian randomization and colocalization for LONP1 and 15 partner genes | Additional file 2: Figs. S7–S8; Tables S11–S12 |
| `04_assembly/` | Assembly of the new figure panels and Additional files 1–2; Supplementary Methods text | Additional files 1–2 |
| `06_control_genes/` | Specificity of the LONP1 co-expression program against control seed genes: all 19p13.3 genes (26 expression- and dosage-matched non-mitochondrial controls), 29 mitochondrial genes on other chromosomes, 300 expression-matched random genes; gene-set AUCs for raw, model A and model B (n=240) correlations | Results "The mitochondrial co-expression program is not positional…"; Table S16 |
| `07_generalization/` | The same pipeline applied to 251 seeds (LONP1, 30 literature genes, 120 dosage-sensitive genes, 100 random genes): positional dominance, retention after copy-number/purity adjustment, specificity against 30 peers; LONP1 rank on all gene sets | Same Results section; Tables S17A–B |

Scripts in 01–04 take input/output directories as command-line arguments (see the usage line at the top of each script). Comments in the scripts are partly in Chinese.

**06 and 07** locate their inputs relative to the script when run without arguments (`python -I 01_control_genes.py`, `python -I 01_generalize.py`), and write results into their own folder (`06_control_genes/`, `07_generalization/results_v2/`). They expect, under a common parent folder: `supplementary_analyses/02_partial_correlation/data/` holding the cBioPortal files `data_mrna_seq_v2_rsem.txt`, `data_log2_cna.txt`, `data_cna.txt` and `data_armlevel_cna.txt` (study `cesc_tcga_pan_can_atlas_2018`; not redistributed), and a sibling folder containing `raw_data/30_coexpression_LONP1_CESC_genomewide.csv` and `raw_data/_cache/TCGA_mastercalls.abs_tables_JSedit.fixed.txt` (both are in `data/` of this repository; copy them to such a `raw_data/` folder). `01_control_genes.py` also accepts explicit paths: `<02_partial_correlation dir> <dir containing raw_data/> <output dir> [genesets.gmt]`. Gene sets are exported from MSigDB with `00_dump_genesets.R` / `00_dump_genesets_all.R` (msigdbr 26.1.1) and are not redistributed. `literature_genes.csv` lists the 30 literature genes with their source articles.

## Data (not redistributed; download from the original sources)

* Single-cell: GEO GSE208653, GSE197461, GSE236738, GSE173231, GSE333892; ArrayExpress E-MTAB-11948. GSE333892 provides a Seurat object with log-normalized values only; integer counts are recovered by `00b_GSE333892_recover_counts.R`.
* TCGA-CESC (cBioPortal `cesc_tcga_pan_can_atlas_2018`) and ABSOLUTE purity (GDC PanCanAtlas), as for the main analysis.
* eQTL: eQTLGen cis-eQTL (2019-12-11 release; significant and full files); GTEx v8 blood, uterus, vagina and fibroblasts from the eQTL Catalogue (QTD000356, QTD000346, QTD000351, QTD000216).
* GWAS: FinnGen R12 (C3_CERVIX_UTERI_EXALLC and subtypes), Pan-UKB (ICD-10 C53; phecodes 180.1, 180.3), BioBank Japan cervical cancer (PheWeb.jp). Regional summary statistics are retrieved remotely by `03_MR_coloc/01_fetch_regions.py`.
* LD reference: 1000 Genomes phase 3 (EUR, EAS).

## Run order

```
# 01 single-cell
Rscript 00_inspect_GSE333892.R ; Rscript 00b_GSE333892_recover_counts.R
python -P 01_duplicate_check.py ; python -P 02_qc_annotate.py ; python -P 03_infercnv.py ; python -P 03b_reclassify.py
python -P 04_pseudobulk_LONP1.py ; python -P 04b_GSE333892_progression.py ; python -P 04d_SIL_grade_sensitivity.py
python -P 06a_export_for_copykat.py ; Rscript 06b_copykat.R ; python -P 06c_compare_copykat.py
python -P 05_figures.py ; python -P 04c_fig_progression.py ; python -P 04e_fig_SIL_grade.py ; python -P 06d_fig_copykat.py
# 02 partial correlation
python -I 01_partial_spearman.py ; Rscript 02_fgsea.R ; python -P 03_figures.py
# 03 Mendelian randomization / colocalization
python3 01_fetch_regions.py ; python -P 02_harmonize_mr.py ; python -P 02b_neighbour_coloc_inputs.py ; python -P 02c_eqtlgen_full_coloc_inputs.py
Rscript 03_coloc.R ; python -P 05_combine_mr.py ; Rscript 07_MR_package_crosscheck.R ; python -P 04_figures_LONP1.py ; python -P 06_figures_allgenes.py
# 04 assembly
python -P build_main_figures.py
# 06 control genes (run inside 06_control_genes/)
Rscript 00_dump_genesets.R genesets_selected.gmt ; python -I 01_control_genes.py
# 07 generalization (run inside 07_generalization/; add --quick for a short test run)
Rscript 00_dump_genesets_all.R genesets_all.gmt ; python -I 01_generalize.py
```

06 and 07 were run with Python 3.11.9 (numpy 2.4.6, pandas 3.1.1, scipy 1.17.1) and R 4.6.2 (msigdbr 26.1.1). On Windows, `RUN_*.bat` performs the same steps.

`03_infercnv.py` and `01_partial_spearman.py` additionally read gene coordinates / TCGA tables produced by the main analysis (`data/` in this repository).

## Environment

Python 3.13.5: scanpy 1.11.5, anndata 0.12.7, infercnvpy 0.6.1, pandas 2.3.2, numpy 2.3.1, scipy 1.16.1, statsmodels 0.14.6, matplotlib 3.10.5, PyMuPDF 1.27.1, openpyxl 3.1.5; pysam (remote tabix queries).
R 4.6.1: CopyKAT 1.2.5, fgsea 1.38.0, msigdbr, MendelianRandomization 0.10.0, coloc 5.2.3, data.table, Seurat (GSE333892 object only).
PLINK 1.90b7.2 (LD clumping). On Windows, `03_infercnv.py` must be run as a script (multiprocessing `__main__` guard is included).

## Result tables

`*/results/` contain the derived tables behind every reported number (e.g. `pseudobulk_LONP1_units.csv`, `infercnv_summary.csv`, `19p13.3_enrichment_by_layer.csv`, `mr_all_genes_combined.csv`, `coloc_*.csv`, `mr_package_crosscheck.csv`). Cell-level objects and regional summary statistics are not included because of their size; they are regenerated by the scripts.
