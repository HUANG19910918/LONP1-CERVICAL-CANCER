# Figure Captions（BMC Cancer格式，拼图后编号按此排）

> **文件夹结构（2026-08-31）**：`figures/子图原文件/` = 各子图（panel）的原始输出（PNG 300dpi + 矢量 PDF，78 个文件，含备用 graphical abstract）；`figures/composite/` = 按稿件编号拼好的整版图 Figure 1–6 与补充图 FigureS1–S2（投稿上传这一套）；`figures/废弃图片/` = 历史弃用版本。拼版脚本 13_assemble.py 已改为从 `子图原文件/` 读取、输出到 `composite/`。
>
> **编号映射（2026-08-30 v5起，按正文首次出现顺序重排）**：稿件 Figure 1 = 本文件 Figure 1；稿件 **Figure 2 = 本文件 Figure 6（单细胞）**；稿件 **Figure 3 = 本文件 Figure 2（基因组关联）**；稿件 **Figure 4 = 本文件 Figure 3（共表达）**；稿件 **Figure 5 = 本文件 Figure 4（生存）**；稿件 **Figure 6 = 本文件 Figure 5（免疫）**。panel 文件名保留原编号前缀（Fig6A_* 即稿件 Figure 2A）。composite/Figure*.png 已按稿件编号命名（13_assemble.py）。
> 稿件 Figure 6（免疫）面板顺序（v7）：A 签名条形、B M1 散点、C DC 散点、D 检查点基因（对应本文件 Figure 5 的 A、B、C、D）。
> 稿件 Figure 4（共表达）面板字母：A 秩图、B chr19、**C curated lollipop、D 去19p13.3富集、E bulk vs 上皮哑铃、F 上皮内在富集、G 负相关富集**（与下文"Figure 3"条目中的字母对应关系：C=C，E→D，D→E，F=F，G=G）。

## Figure 1. LONP1 expression in cervical cancer and precursor lesions.
**(A)** LONP1 mRNA expression in TCGA-CESC primary tumors (n=304) versus normal cervical tissue (GTEx n=10; TCGA adjacent normal n=3), RSEM-based log2(TPM+0.001) from the UCSC Xena Toil recompute. Violin plots with embedded boxplots (median, IQR); Mann-Whitney U test, p=7.15×10⁻⁴ (Welch t-test p=2.26×10⁻⁴).
**(B)** GSE9750 (Affymetrix U133A, probe 209017_s_at): cervical cancer tissue (n=31) versus normal cervix (n=24), log2 MAS5 signal; Mann-Whitney p=2.29×10⁻⁵.
**(C)** GSE7803 (U133A): normal cervix (n=10), HSIL (n=7), invasive squamous cell carcinoma (n=21); SCC vs normal Mann-Whitney p=0.052.
**(D)** GSE63514 (U133 Plus 2.0, laser-capture microdissected epithelium): normal (n=24), CIN1 (n=14), CIN2 (n=22), CIN3 (n=40), cancer (n=28); cancer vs normal Mann-Whitney p=0.978 (Kruskal-Wallis across groups p=0.348). Note the discordance with bulk-tissue cohorts, discussed in the text (epithelium-only versus bulk-tissue comparisons).
**(E)** Pan-cancer overview: LONP1 in TCGA tumors versus site-matched normal tissue (GTEx Normal Tissue + TCGA Solid Tissue Normal) across 32 cancer types (LAML excluded, comparator is blood). Asterisks: BH-adjusted Mann-Whitney (*FDR<0.05, **<0.01, ***<0.001); types with <5 normals shown without testing. Y-axis truncated at 2 for display; complete values in raw_data/15. CESC adjusted p=7.4×10⁻⁴; note the opposite direction in OV (tumor < normal).
Boxplots: center line=median, box=IQR, whiskers extend to 1.5×IQR; individual samples shown as dots.

## Figure 2. Genomic correlates of LONP1 expression in TCGA-CESC.
**(A)** LONP1 mRNA (log2(RSEM+1)) versus copy-number level (log2 CNA ratio), n=290; Spearman rho=0.41, p=3.1×10⁻¹³; red line = linear fit with 95% CI.
**(B)** LONP1 mRNA by GISTIC2 copy-number class (deep deletion n=6, shallow deletion n=94, diploid n=159, gain n=30, amplification n=1); Kruskal-Wallis p<0.001. LONP1 somatic mutations were rare (1/294 missense).
**(C)** LONP1 mRNA by histological subtype: squamous cell carcinoma (n=229) versus adenocarcinoma (n=43); Welch t-test p=0.022.

## Figure 3. The LONP1 co-expression landscape in TCGA-CESC（2026-08-30 LONP1中心版，取代旧Fig3）.
**(A)** Genome-wide Spearman correlation of all annotated genes (n=17,378) with LONP1, ranked; genes on 19p13.3 in red (n=188); dashed lines |rho|=0.3 (628 positive / 342 negative). Top partners CLPP 0.72, THOP1 0.68, HDGFL2 0.67, WDR18 0.67, SGTA 0.64, DUS3L 0.63 — all on 19p13.3.
**(B)** Correlation versus genomic position along chromosome 19 (GRCh38); shaded = 19p13.3 (0–6.9 Mb); dashed line = LONP1; blue line = 1-Mb sliding median. Median rho: 19p13.3 0.24, 19p arm 0.15, 19q −0.01, rest of genome −0.03.
**(C)** Curated mitochondrial proteostasis / biogenesis / dynamics genes (server-side Spearman, n=294; *p<0.05, **p<0.01, ***p<0.001; 19p13.3 genes labelled). Outside 19p13.3: TRAP1 0.41, PHB1 0.40, PHB2 0.39, HSPD1 0.33, ATF4 0.32, FIS1 0.32, HSPE1 0.30, HSPA9 0.23; DNM1L −0.02; MFN1 −0.23, OPA1 −0.22.
**(D)** Bulk (TCGA-CESC, n=294) versus epithelium-intrinsic (GSE208653, 58 sample×subcluster pseudobulk units, scripts/17) correlations for curated proteostasis / UPRmt-ISR TF / biogenesis / dynamics genes; shaded band + dotted line = IQR (0.10–0.43) and median (0.28) of the epithelium-intrinsic correlations of all 14,140 detected genes（上皮假批量相关整体上移，需按百分位解读）. Top-quintile within epithelium: POLRMT 0.62 (97th pct), FIS1 0.51, TRAP1 0.49, OPA1 0.47, DNAJA3 0.47, ATF5 0.47, PMPCB 0.46, CLPP 0.45 (79th), HSPA9 0.45 (78th); not tracking: ATF4 0.16 (32nd), DDIT3 0.00 (16th), MFN2 0.05. 19p13.3 genes show no positional effect within epithelium (median 0.26 vs 0.28).
**(E)** Enrichr over-representation of positively co-expressed genes (rho≥0.4) after excluding 19p13.3 genes (n=102; BH-adjusted; dashed line adj. p=0.05): mitochondrial (inner) membrane 2.3×10⁻¹², translation 2.3×10⁻¹², mitochondrial translation 8.4×10⁻⁵, aerobic respiration 8.7×10⁻⁴, KEGG ribosome 3.3×10⁻⁸; KEGG oxidative phosphorylation not significant after exclusion (adj. p=0.076; was 1.0×10⁻³ with 19p13.3 genes).
**(F)** Enrichr over-representation of the top-300 epithelium-intrinsic LONP1 partners (raw_data/38): mitochondrial respiratory chain complex assembly adj. p=9.6×10⁻⁸, ER membrane 8.5×10⁻⁷, mitochondrial membrane 6.5×10⁻⁵, mitochondrial inner membrane 7.6×10⁻³; KEGG none significant. Sample-centred sensitivity: mito inner membrane 1.7×10⁻⁴, complex assembly 1.4×10⁻⁵.
**(G)** Enrichment of negatively correlated bulk genes (rho≤−0.35, n=96): Ras signalling adj. p=3.0×10⁻⁴, endocytosis.

## Figure 4. LONP1 expression is not associated with survival in TCGA-CESC (negative result).
**(A-C)** Kaplan-Meier curves by median split of LONP1 expression for overall survival (n=281, 68 events; log-rank p=0.429), disease-specific survival (n=277, 51 events) and progression-free survival (n=281, 68 events). Crosses = censoring.
**(D)** Univariate Cox regression (LONP1 as continuous log2 variable): OS HR=0.84 (95%CI 0.54-1.31, p=0.451); DSS HR=0.99 (0.59-1.64, p=0.956); PFS HR=0.85 (0.55-1.32, p=0.473).

## Figure 5. Weak inverse association between LONP1 and antigen-presenting cell markers.
**(A)** Spearman correlations between LONP1 and immune marker signatures (mean z-score of member genes; gene sets listed in Methods), TCGA-CESC n=294. Colored bars: BH-FDR<0.05 (M1 macrophage rho=-0.22, dendritic cell rho=-0.22); grey: not significant.
**(B)** LONP1 versus M1 macrophage signature; line = linear fit with 95% CI.

## Figure 6. Single-cell resolution of LONP1 in the cervical carcinogenesis spectrum (GSE208653).
80,435 cells from 9 samples (HPV- normal ×2, HPV+ normal ×2, HSIL ×2, cancer ×3).
**(A)** UMAP of nine annotated major cell types (markers in Methods).
**(B)** LONP1 expression (log1p CP10K) on UMAP.
**(C)** LONP1 by cell type, ordered by mean (black dots); epithelium highest (0.143), ~2-fold over T/NK (0.061) and myeloid (0.065).
**(D)** Epithelial cells (n=20,363) by disease stage; white circles = per-sample pseudobulk means. Cancer epithelium is NOT higher than normal (0.113 vs 0.132–0.142); HSIL highest (0.196). Per-cell Kruskal-Wallis p=3.3×10⁻²³; pseudobulk cancer-vs-normal MWU p=0.057 (3 cancer samples are the lowest 3 of 9).
**(E)** Cell-type composition per stage — the tissue-composition basis for bulk-level LONP1 elevation.

## Supplementary Figure S1 (Additional file 2). Tumor purity (ABSOLUTE, n=291) — scripts/19, raw_data/42–43
**(A)** LONP1 vs purity rho=0.12 p=0.045; **(B)** epithelial score vs purity rho=0.13; **(C)** CNA vs expression purity-adjusted residuals, partial rho=0.45 p=1.2e-15.
## Supplementary Figure S2 (Additional file 2). GSE44001 DFS (n=300, 38 events) — scripts/20, raw_data/44–45
Median-split KM log-rank p=0.83; Cox HR=1.07 (0.72–1.60) p=0.73; multivariable LONP1 p=0.88, diameter HR 1.29/cm p=0.002.

---
### 检验口径速查（内部用，不投稿）
| 图 | 数据文件 | 脚本 | 关键值 |
|---|---|---|---|
| 1A | raw_data/01 | 02_fig_expression.R | MWU p=7.15e-4 |
| 1B | raw_data/04 | 同上 | p=2.29e-5 |
| 1C | raw_data/05 | 同上 | p=0.052 |
| 1D | raw_data/03 | 同上 | p=0.978 |
| 2A-B | raw_data/07+08 | 同上 | rho=0.411, p=3.1e-13 |
| 2C | raw_data/08+10 | 同上 | p=0.0216 |
| 3A-B | raw_data/30 | 15_coexpression_positional.py + 16_fig3_lonp1_centric.py | top20中18个在19p13.3；Fisher OR=82.2, p=1.1e-81 (rho≥0.4) |
| 3C | raw_data/30→33 | 16 | 见33_mito_module_spearman.csv |
| 3D | raw_data/30+37(+35) | 17+16 | 58单元；上皮内CLPP 0.45(79pct)/HSPA9 0.45/TRAP1 0.49；OPA1 +0.47；ATF4/DDIT3不相关 |
| 3E | raw_data/30→31 | 15+16 | 去19p13.3后 mito inner membrane padj=2.3e-12；OXPHOS padj=0.076 (ns) |
| 3F | raw_data/37→38 | 17+16 | 上皮内在top300：呼吸链复合体组装 padj=9.6e-8 |
| 3G | raw_data/06→21 | 16 | 负相关Ras padj=3.0e-4 |
| S1A-C | raw_data/42→43 | 19 | 纯度rho=0.12；CNA偏相关0.45 |
| S2 | raw_data/44→45 | 20 | GSE44001 HR=1.07 p=0.73 |
| (旧3B-C) | raw_data/06→11 | 03/09/12 | 含19p13.3基因的富集，仅正文对比用 |
| 4A-D | raw_data/02→13 | 04_fig_survival_immune.R | 见13_cox_forest_data.csv |
| 5A-B | raw_data/08→14 | 同上 | 见14_immune_signature_spearman.csv |
| 6A-E | raw_data/17-20 | 06_scRNA(py)+07b_plotnine | 上皮0.143 vs T/NK 0.061；癌上皮0.113<正常0.132-0.142 |
| 2D-E | raw_data/02→22 | 08_..._plotnine.py | T期p=0.493；N p=0.999 (ns) |
| (旧3E，已弃用) | raw_data/08→23 | 09 | 平衡分数rho=0.345（构造指标，v3不再展示） |
| 4E | raw_data/02→24 | 同上(lifelines) | LONP1 HR=0.98 p=0.957；N1 HR=2.91 p=0.004 |
| 5C-D | raw_data/08→25 | 同上 | 检查点基因均\|rho\|<0.15 |

**字体说明**：全图当前为Liberation Serif（TNR公制兼容替身）；将times*.ttf放入`fonts/`目录后重跑08/09/07b脚本即切换为正版Times New Roman。图注补充：Fig2新增(D)T分期(KW p=0.49)、(E)N状态(p=0.999)；Fig3已于2026-08-30整体重做为LONP1中心版（A-F，见上）；Fig4新增(E)多因素Cox森林图；Fig5新增(C)DC散点、(D)检查点基因条形图。
