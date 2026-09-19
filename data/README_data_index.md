# raw_data 数据字典（2026-08-30 整理，对应稿件 v7）

编号 = 文件前缀。"稿件图" 用稿件最终编号（Figure 1–6、S1–S2）；"脚本" 为生成/使用该文件的脚本。所有文件均可由公开数据库按脚本口径重取。

| 编号 | 文件 | 内容 | 来源/口径 | 脚本 | 稿件图/表 |
|---|---|---|---|---|---|
| 01 | TCGA_GTEx_LONP1_expression_tumor_vs_normal | 逐样本 LONP1 log2(TPM+0.001)，Tumor 304 / Normal 13 | Xena Toil `TcgaTargetGtex_rsem_gene_tpm` | 02_fig_expression.R, 08, 12(raincloud) | Fig 1B |
| 02 | TCGA_CESC_clinical_LONP1_survival | 患者级 OS/DSS/PFS、年龄、pT/pN、LONP1 log2(RSEM+1) | cBioPortal `cesc_tcga_pan_can_atlas_2018` | 04_fig_survival_immune.R, 09(f4) | Fig 3D-E, Fig 5 |
| 03 | GSE63514_LONP1_209017_s_at | LCM 上皮 128 样本，log2 归一化值 | GEO GSE63514 (GPL570) | 02, 08 | Fig 1E |
| 04 | GSE9750_LONP1_209017_s_at | 正常 24 / 癌 31 / 细胞系 9，MAS5 | GEO GSE9750 (GPL96) | 02, 08 | Fig 1C |
| 05 | GSE7803_LONP1_209017_s_at | normal 10 / HSIL 7 / SCC 21，log2 | GEO GSE7803 (GPL96) | 02, 08 | Fig 1D |
| 06 | coexpression_LONP1_CESC | 服务器端 Spearman，\|rho\|≥0.3 的 970 基因（阈值 0.3 版） | cBioPortal co-expression API | 00_fetch_cbioportal.R | （被 30 取代，保留核对） |
| 07 | CESC_LONP1_CNA | 逐样本 log2CNA + GISTIC 分类 | cBioPortal | 00 | Fig 3A-B |
| 08 | CESC_gene_panel_expression_RSEM | 58 基因（线粒体 panel + 免疫标志物）逐样本 RSEM，n=294 | cBioPortal | 00 | Fig 3, 5, 6 |
| 10 | CESC_subtype | 患者级 SUBTYPE | cBioPortal | 00 | Fig 3C |
| 11 | enrichment_enrichr_pos | rho≥0.4 全部 164 基因的 Enrichr（含 19p13.3） | Enrichr | 03_fig_mito_enrich.R | Table S2（对比用） |
| 12 | mito_panel_spearman | 19 个线粒体 panel 基因本地 Spearman | 由 08 计算 | 01_verify_stats.R | 核对用 |
| 13 | cox_forest_data | 单因素 Cox OS/DSS/PFS | 由 02 计算 | 04 | Fig 5D |
| 14 | immune_signature_spearman | 13 个免疫签名 Spearman + BH | 由 08 计算 | 04 | Fig 6A, Table S6 |
| 15 | pancancer_LONP1_per_sample | 32 癌种逐样本 | Xena Toil | 05_fig_pancancer.R | Fig 1F |
| 16 | pancancer_LONP1_stats | 各癌种 MWU + BH | 由 15 计算 | 05 | Fig 1F, Table S7 |
| 17 | GSE208653_LONP1_per_cell.csv.gz | 80,435 细胞：类型/阶段/样本/LONP1/UMAP/leiden（行序 = 06 脚本加载顺序） | GEO GSE208653_RAW.tar | 06_scRNA_GSE208653.py | Fig 2A-D |
| 18 | GSE208653_cluster_marker_scores | 28 簇 × marker 均值（注释依据） | 同上 | 06 | 核对用 |
| 19 | GSE208653_epithelial_pseudobulk_LONP1 | 上皮逐样本假批量 LONP1 | 同上 | 06/07b | Fig 2E |
| 20 | GSE208653_LONP1_by_celltype_stage | 类型 × 阶段汇总 | 同上 | 06/07b | Fig 2D-F |
| 21 | enrichment_enrichr_neg | rho≤−0.35 的 96 基因 Enrichr | Enrichr | 09 | Fig 4G, Table S2 |
| 22 | clinical_assoc_Tstage / Nstage | T/N 分组统计 | 由 02 计算 | 08 | Fig 3D-E |
| 23 | fission_fusion_balance | 分裂-融合平衡分数（v2 旧指标，v3 起不展示） | 由 08 计算 | 09 | — |
| 24 | multivariable_cox_OS | 多因素 Cox（LONP1+年龄+T+N） | 由 02 计算 | 09 | Fig 5E |
| 25 | checkpoint_gene_corr | 7 个检查点基因 Spearman + Bonferroni | 由 08 计算 | 09 | Fig 6D |
| 26 | marker_dotplot_data | 单细胞 marker 气泡图数据 | 由 17/h5ad 计算 | 12 | Fig 2C |
| 27 / 27b | epithelial_intrinsic_corr / pseudobulk_units | **旧**上皮重聚类（63 单元，R14；代码未存档） | — | 14 | 历史记录（v4 起被 34/35/37 取代） |
| 28 | TCGA_GTEx_epithelial_markers | 317 样本 EPCAM/KRT5/KRT8/KRT17 + 上皮分数 | Xena Toil | 14 | Fig 2G, S1B |
| 29 | composition_regression | LONP1 ~ tumor ± 上皮分数 OLS | 由 28 计算 | 14 | Fig 2G |
| 30 | coexpression_LONP1_CESC_genomewide | 17,378 基因 Spearman + GRCh38 坐标 + 细胞带 + 19p13.3 标记 | cBioPortal (threshold 0) + Ensembl + UCSC | 15_coexpression_positional.py | Fig 4A-C, Table S1 |
| 31 | enrichment_enrichr_pos_excl19p13.3 | rho≥0.4 且非 19p13.3 的 102 基因 Enrichr | Enrichr | 15 | Fig 4D, Table S2 |
| 32 | positional_stats.txt | 位置效应统计（Fisher、中位 rho 等） | 由 30 计算 | 15 | 正文 |
| 33 | mito_module_spearman | curated 29 基因（服务器端 rho） | 由 30 提取 | 16 | Fig 4C, Table S5 |
| 34 | epithelial_pseudobulk_units_extended | 58 个 sample×亚群单元的模块基因均值 | GSE208653 重跑 | 17_scRNA_epithelial_module.py | Table S5 |
| 35 | epithelial_intrinsic_corr_extended | 模块基因三层相关（逐细胞/假批量/逐样本） | 同上 | 17 | Table S5 |
| 36 | epithelial_recluster_check.txt | 重聚类与旧 27b 的比对、全基因组统计 | 同上 | 17 | 记录 |
| 37 | epithelial_pseudobulk_genomewide_rho | 14,140 基因上皮内在 rho + 百分位 + 样本中心化 + bulk rho | 同上 | 17 | Fig 4E, Table S3 |
| 38 | enrichment_epithelial_intrinsic_top300 | 上皮 top300 Enrichr（raw / centered） | Enrichr | 17 | Fig 4F, Table S4 |
| 39 | fig3D_bulk_vs_epithelium_module | Fig 4E 绘图数据 | 由 30+37 | 16 | Fig 4E |
| 40 | reference_check | 44 条参考文献 PubMed 核对（PMID） | E-utilities | 18_verify_references.py | 参考文献 |
| 41 | TCGA_CESC_HPV_types_cesc_tcga | 旧版 TCGA 22 例 HPV 分型（未用） | cBioPortal `cesc_tcga` | — | — |
| 42 | purity_merged | 291 样本 ABSOLUTE 纯度 + LONP1 + CNA + 上皮分数 + 免疫签名 | GDC PanCanAtlas ABSOLUTE（_cache） | 19_purity_adjustment.py | Fig S1 |
| 43 | purity_stats.txt | 纯度相关/偏相关统计 | 由 42 计算 | 19 | 正文 |
| 44 | GSE44001_LONP1_DFS | 300 例 LONP1(ILMN_1766125) + 分期/最大径/DFS | GEO GSE44001 series matrix | 20_GSE44001_survival.py | Fig S2 |
| 45 | GSE44001_survival_stats.txt | Cox/KM 结果 | 由 44 计算 | 20 | 正文 |
| 46 | manuscript_number_check.txt | 稿件全文 261 项数值程序化复核 | 全部 | 21_full_manuscript_check.py | — |
| verification_report.txt | — | R 复算报告（R1–R8，2026-08-29） | — | 01_verify_stats.R | — |
| _cache/ | coexp_full_raw.json, ensembl_positions.json, cytoBand_hg38.txt.gz, TCGA_mastercalls.abs_tables_JSedit.fixed.txt | 外部下载缓存 | — | 15, 19 | — |

## 脚本索引（scripts/）
00 取数(R) · 01 R 复算 · 02–05 R 图（参考版）· 06 单细胞全流程 · 07/07b 单细胞图 · 08/09 plotnine 主图 · 10 示意图 · 11 graphical abstract · 12 图型升级 · 13 拼版（含编号映射）· 14 补充 panel(2G/3F 旧) · 15 全基因组共表达+位置 · 16 Figure 4 · 17 单细胞上皮内在 · 18 参考文献核对 · 19 纯度 · 20 GSE44001 · 21 全文数值复核 · patch_v3…v7 + build_manuscript_v7.js 稿件生成链 · figstyle.py 字体/主题

## 稿件生成链
build_manuscript.js (v2) → patch_v3 → v3 → patch_v4 → v4 → patch_v5 + patch_v5b → v5 → patch_v6 → v6 → patch_v7 → **v7**（当前）
