# LONP1 in cervical cancer — multi-resolution transcriptomic analysis

Code and extracted data for the manuscript *Multi-resolution transcriptomic profiling of LONP1 in cervical cancer: epithelial-enriched expression, 19p13.3 copy-number-driven variation and a mitochondrial proteostasis programme* (Wenzhou Municipal Science and Technology Bureau project).

Every number and figure in the manuscript can be regenerated from the files in this repository. All source data are public (TCGA, GTEx, NCBI GEO); the per-sample extracts used for analysis are included under `data/` so that the statistics can be recomputed without re-downloading.

## Layout

| Path | Content |
|---|---|
| `scripts/` | Numbered analysis and figure scripts (R 4.3.3 and Python 3.10/3.11), run in numeric order |
| `data/` | Per-sample extracts and derived tables (`01`–`46`; see `data/README_data_index.md`), `verification_report.txt`, `32_positional_stats.txt`, `36_epithelial_recluster_check.txt`, `40_reference_check.csv`; `_cache/` holds the cBioPortal genome-wide co-expression response, Ensembl gene coordinates and the UCSC hg38 cytoBand table |
| `figures/panels/` | Individual panels (PNG 300 dpi + PDF) |
| `figures/composite/` | Assembled Figures 1–6 |
| `docs/figure_captions.md` | Panel-level captions with the data file and script behind each panel |
| `docs/lab_notebook_zh.md` | Full analysis log (Chinese), entries R1–R16 |

## Data sources

| Dataset | Access |
|---|---|
| TCGA-CESC + GTEx cervix, RSEM log2(TPM+0.001) | UCSC Xena Toil hub, `TcgaTargetGtex_rsem_gene_tpm`, phenotype `TcgaTargetGTEX_phenotype` |
| TCGA-CESC PanCancer Atlas (expression, GISTIC2, log2 CNA, mutations, subtype, survival) | cBioPortal API, study `cesc_tcga_pan_can_atlas_2018` |
| GSE9750, GSE7803, GSE63514 (Affymetrix, probe 209017_s_at) | NCBI GEO |
| GSE208653 (10x scRNA-seq, 9 samples) | NCBI GEO, `GSE208653_RAW.tar` (md5 `e264e8e98f1005f153c3beaf1f9e29bd`) |
| GSE44001 (300 early-stage cervical cancers, DFS; Illumina DASL GPL14951) | NCBI GEO series matrix |
| ABSOLUTE tumor purity (TCGA PanCanAtlas) | GDC PanCanAtlas publications page (`TCGA_mastercalls.abs_tables_JSedit.fixed.txt`, copy in `data/_cache/`) |
| Human Protein Atlas (LONP1, antibody HPA002192) | proteinatlas.org (ENSG00000196365) |
| Gene coordinates / cytobands | Ensembl REST (GRCh38) / UCSC hg38 `cytoBand.txt.gz` |
| Enrichment | Enrichr (GO_Biological_Process_2023, GO_Cellular_Component_2023, KEGG_2021_Human) |

## Reproduction

```bash
# 1. Fetch TCGA-CESC data from cBioPortal (R) and verify core statistics
Rscript scripts/00_fetch_cbioportal.R  <analysis_dir>
Rscript scripts/01_verify_stats.R      <analysis_dir>     # -> data/verification_report.txt

# 2. Figures 1, 2, 4, 5 (plotnine); pan-cancer overview (R)
python3 scripts/08_figs_expression_plotnine.py        <analysis_dir> all
python3 scripts/09_figs_mito_surv_immune_plotnine.py  <analysis_dir> f4
python3 scripts/09_figs_mito_surv_immune_plotnine.py  <analysis_dir> f5
Rscript scripts/05_fig_pancancer.R                    <analysis_dir>
python3 scripts/12_upgrade_panels.py                  <analysis_dir> all

# 3. Single-cell (GSE208653): cell-type atlas, then epithelium-intrinsic co-expression
python3 scripts/06_scRNA_GSE208653.py <raw_dir> s0 ... s8 ; merge ; norm ; pca ; cluster ; umap ; annotate
python3 scripts/07b_fig_scRNA_plotnine.py  <analysis_dir>
python3 scripts/17_scRNA_epithelial_module.py <raw_dir> <analysis_dir> load
python3 scripts/17_scRNA_epithelial_module.py <raw_dir> <analysis_dir> recluster batch
python3 scripts/17_scRNA_epithelial_module.py <raw_dir> <analysis_dir> corr batch
python3 scripts/17_scRNA_epithelial_module.py <raw_dir> <analysis_dir> genomewide batch

# 4. Genome-wide co-expression with positional annotation, Figure 3, composites
python3 scripts/15_coexpression_positional.py  <analysis_dir>
python3 scripts/16_fig3_lonp1_centric.py       <analysis_dir> all
python3 scripts/14_supplement_panels.py        <analysis_dir> f6g
python3 scripts/10_schematic.py <analysis_dir>; python3 scripts/11_graphical_abstract.py <analysis_dir>
python3 scripts/13_assemble.py  <analysis_dir> all

# 5. Tumor purity (ABSOLUTE), external survival cohort (GSE44001), supplementary figures
python3 scripts/19_purity_adjustment.py <analysis_dir> data/_cache/TCGA_mastercalls.abs_tables_JSedit.fixed.txt
python3 scripts/20_GSE44001_survival.py <analysis_dir> GSE44001_series_matrix.txt.gz
python3 scripts/13_assemble.py <analysis_dir> FigureS1 ; python3 scripts/13_assemble.py <analysis_dir> FigureS2

# 6. Manuscript (docx) and reference check
node scripts/build_manuscript_v7.js out.docx
python3 scripts/21_full_manuscript_check.py <analysis_dir>   # -> data/46_manuscript_number_check.txt
python3 scripts/18_verify_references.py refs.txt data/40_reference_check.csv
```

`<analysis_dir>` is a directory containing `raw_data/` (= `data/` here) and `figures/`. Fonts: scripts use Times New Roman when `../fonts/times*.ttf` is present, otherwise the metric-compatible Liberation Serif.

## Environment

See `environment.txt`, which documents two environments: **A** for the original analyses (scripts 00–21) and **B** for the supplementary proportional-hazards diagnostics and time-varying models (scripts 22–30). Key versions: R 4.3.3 (survival, ggplot2, jsonlite); Python 3.10–3.11 with scanpy 1.11.5, anndata, pandas, numpy, scipy, statsmodels, lifelines 0.30, plotnine 0.15.8, matplotlib; Node 22 with `docx`.

## Notes on reproducibility

* Statistics were computed independently in R and Python (and, for the first pass, in browser JavaScript) and cross-checked; see `data/verification_report.txt`.
* Leiden clustering of the epithelial subset is deterministic for a fixed scanpy/igraph version but may differ slightly across versions; `data/36_epithelial_recluster_check.txt` records the run used in the manuscript (25 subclusters, 58 pseudobulk units).
* Negative and discordant results (no prognostic association; no difference in microdissected epithelium; positional confounding of co-expression) are reported in full.

## License

Code: MIT (see `LICENSE`). Data extracts remain subject to the terms of their original repositories (TCGA/GTEx/GEO).

## Citation

To be added on acceptance (Zenodo DOI).
