# 实验记录：LONP1宫颈癌生信分析（温州市科技局结题文章）

**项目**：基于线粒体动态探索Lonp1调控宫颈癌发生发展的机制研究（2024温州市科技局）
**目标**：纯干实验（生信分析）SCI文章一篇，**目标期刊 BMC Cancer**（用户指定，格式：Background/Methods/Results/Discussion/Conclusions + Declarations，Vancouver参考文献），3个月内接收
**执行**：Claude（Cowork）+ 黄世园
**开始日期**：2026-08-29

---

## 分析环境说明

本地计算沙盒启动失败，本轮分析用**浏览器JavaScript**直接调用公开数据库API取数并现场计算。统计算法（Welch t、Mann-Whitney U、log-rank、Cox[Efron ties, Newton-Raphson]、Spearman/Pearson）为自行实现。
**待办**：沙盒恢复后，用R从相同数据源重新取数复算所有数值并绘图（ggplot2/survminer），R脚本作为最终检验口径存档于 `scripts/`。当前所有数值需与R复算结果核对一致后方可写入成稿。

## 数据来源（公开数据库，取数日期2026-08-29）

| # | 数据 | 来源与精确口径 |
|---|---|---|
| D1 | TCGA+GTEx宫颈表达 | UCSC Xena Toil hub（toil.xenahubs.net），数据集 `TcgaTargetGtex_rsem_gene_tpm`（单位log2(TPM+0.001)），表型 `TcgaTargetGTEX_phenotype.txt`。筛选：_primary_site ∈ {Cervix, Cervix Uteri}；肿瘤=TCGA Primary Tumor（n=304），正常=GTEx Normal Tissue（n=10）+ TCGA Solid Tissue Normal（n=3） |
| D2 | TCGA-CESC生存/临床/多组学 | cBioPortal API（cbioportal.org/api），研究 `cesc_tcga_pan_can_atlas_2018`，表达profile `..._rna_seq_v2_mrna`（RSEM，分析用log2(x+1)），样本列表 `..._all`（n=294有表达） |
| D3 | GEO验证集1 | GSE63514（GPL570，激光显微切割上皮：Normal 24 / CIN1 14 / CIN2 22 / CIN3 40 / Cancer 28），LONP1探针 **209017_s_at**（探针由GPL570.annot确定，全平台唯一LONP1探针） |
| D4 | GEO验证集2 | GSE9750（GPL96，MAS5值：正常宫颈24 / 宫颈癌组织31 / 细胞系9），探针209017_s_at |
| D5 | GEO验证集3 | GSE7803（GPL96，log2值：正常10 / HSIL 7 / SCC 21），探针209017_s_at |
| D6 | 富集分析 | Enrichr API（maayanlab.cloud/Enrichr），库：GO_Biological_Process_2023, GO_Cellular_Component_2023, KEGG_2021_Human |

## 结果记录（含检验口径）

### R1. 肿瘤 vs 正常差异表达（D1）✅ 显著上调 → 拟Figure 1A

- 肿瘤(n=304)：mean 5.602 / median 5.616 / SD 0.591；正常(n=13)：mean 5.133 / median 5.236 / SD 0.332（log2 TPM）
- log2FC(均值差)=0.469；Welch t=4.776, df=15.5, **p=2.26×10⁻⁴**；Mann-Whitney z=-3.384, **p=7.15×10⁻⁴**
- 原始数据：`raw_data/01_TCGA_GTEx_LONP1_expression_tumor_vs_normal.csv`（逐样本，可直接复算）
- 检验口径：对CSV分组做Welch t和MWU即可复现上述统计量

### R2. 生存分析（D2）❌ 全部不显著（如实报告）→ 拟Figure 4（KM+森林图，作为阴性结果）

单因素Cox（连续变量log2表达，Efron法）：

| 终点 | n | 事件 | HR | 95%CI | p |
|---|---|---|---|---|---|
| OS | 281 | 68 | 0.843 | 0.541–1.314 | 0.451 |
| DSS | 277 | 51 | 0.986 | 0.592–1.641 | 0.956 |
| PFS | 281 | 68 | 0.852 | 0.550–1.320 | 0.473 |

- 中位数分组log-rank（OS）：χ²=0.795, p=0.373（cutoff=10.933 log2(RSEM+1)，高/低各147）
- **决策**（用户确认）：不做最优截断值p值捕捞；文章定位改为"表达+分子关联特征"，生存作阴性结果报告
- 原始数据：`raw_data/02_TCGA_CESC_clinical_LONP1_survival.csv`（n=293，TCGA-EX-A449无随访）
- 检验口径：R `survival::coxph(Surv(OS_months,OS_event)~LONP1_log2RSEM)`、`survdiff` 复算；排除t=0样本

### R3. 线粒体动态/质控基因相关性（D2）✅ 方向清晰 → 拟Figure 3A（相关性热图/棒棒糖图）

Spearman（n=294，log2(RSEM+1)）：

| 正相关 | rho | p | 负相关 | rho | p |
|---|---|---|---|---|---|
| CLPP | 0.718 | <1e-40 | OPA1 | -0.223 | 9.1e-5 |
| HSPD1 | 0.329 | 2.6e-9 | MFN1 | -0.225 | 7.7e-5 |
| FIS1 | 0.315 | 1.4e-8 | MFN2 | -0.126 | 0.030 |
| HSPA9 | 0.228 | 6.2e-5 | YME1L1 | -0.089 | 0.128(ns) |
| NRF1 | 0.187 | 1.1e-3 | TFAM | -0.100 | 0.086(ns) |
| AFG3L2 | 0.164 | 4.6e-3 | PINK1 | -0.085 | 0.145(ns) |
| MFF | 0.152 | 8.8e-3 | DNM1L | -0.020 | 0.733(ns) |

- 解读：LONP1高表达伴随"促分裂(FIS1/MFF↑)/抑融合(OPA1/MFN1/MFN2↓)"+线粒体蛋白酶质控模块(CLPP/HSPD1/HSPA9)共上调
- 检验口径：cBioPortal取上述基因表达，log2(x+1)后Spearman；沙盒恢复后R复算并出图

### R4. 全基因组共表达+富集（D2+D6）✅ → 拟Figure 3B-D

- cBioPortal服务器端co-expression接口（molecularProfileIdA=B=rna_seq_v2_mrna, threshold=0.3, entrezGeneId=9361[LONP1], sampleList=_all）：|rho|≥0.3共970基因（正628/负342）
- 正相关top：CLPP 0.718, THOP1 0.685, HDGFL2 0.673, WDR18 0.670, SGTA 0.636, NDUFA11 0.614, TIMM13 0.597, MICOS13 0.570（注：多个top基因位于19p13.3，与LONP1同染色体区段，成稿讨论需提示CNV共变可能）
- 富集（正相关rho≥0.4，n=165，Enrichr）：
  - GO-CC：Mitochondrial Inner Membrane padj=8.1e-14；Mito Membrane 3.5e-13；Respiratory Chain Complex I 3.0e-4
  - GO-BP：Translation 5.0e-13；Aerobic Respiration 4.9e-6；Mitochondrial Gene Expression 2.4e-5；OXPHOS 8.2e-5
  - KEGG：Ribosome 9.9e-8；Oxidative phosphorylation 1.1e-3
- 负相关（rho≤-0.35，n=96）：KEGG Ras signaling padj=3.0e-4；Endocytosis 1.2e-2
- 检验口径：接口参数如上，结果确定性可重取；成稿前用R重取全表存 `raw_data/06_coexpression_LONP1_CESC.csv`（浏览器标签页意外关闭，本轮未落盘，数值以重取为准）

### R5. 免疫标志物相关性（D2）弱负相关 → 拟Figure 5（如实报告弱效应）

Spearman，签名=成员基因z-score均值（n=294）：M1巨噬细胞 rho=-0.221(p=1.1e-4)、DC rho=-0.217(p=1.4e-4)、CD4/Th -0.159、Treg -0.121、检查点配体(CD274/PDCD1LG2) -0.119(p=0.041)；CD8/NK/B/细胞毒/耗竭均ns
- 基因集定义记录：CD8(CD8A,CD8B)、CD4/Th(CD4,IL7R)、B(CD19,MS4A1,CD79A)、NK(NCAM1,KLRD1,NKG7)、M1(NOS2,IRF5,PTGS2)、M2(CD163,MRC1,MS4A4A)、单核(CD86,CSF1R)、中性粒(CEACAM8,ITGAM,CCR7)、DC(ITGAX,CD1C,NRP1)、Treg(FOXP3,CCR8,IL2RA)、耗竭(PDCD1,CTLA4,LAG3,HAVCR2,TIGIT)、配体(CD274,PDCD1LG2)、细胞毒(GZMA,GZMB,PRF1,IFNG)
- 待办：R恢复后用GSVA::ssgsea或CIBERSORTx正式版免疫浸润复核方向一致性

### R6. 表达调控：CNV/突变（D2）✅ CNV驱动 → 拟Figure 2C-D

- log2CNA vs 表达 Spearman rho=0.411, **p=2.1×10⁻¹⁴**（n=290）
- GISTIC分类表达均值：-2:10.76(n=6), -1:10.68(n=94), 0:11.02(n=159), +1:11.45(n=30), +2:11.42(n=1) → 随拷贝数递增
- 突变：仅1/294例错义突变(E654D) → LONP1上调非突变驱动，主要为19p13.3拷贝数增益
- 甲基化：cBioPortal该研究无LONP1基因级甲基化数据（两个profile均空）；可选：Xena TCGA hub 450K探针级（待R阶段决定是否补）

### R7. 亚型差异（D2）小幅差异

- 腺癌 11.119±0.527(n=43) vs 鳞癌 10.908±0.591(n=229)（log2RSEM+1），差0.21，约p~0.02（待R精确复算）

### R8. GEO独立验证（D3-D5）⚠️ 2/3支持上调，1/3不支持（如实报告）→ 拟Figure 1B-D

| 数据集 | 比较 | 结果 | p (MWU) |
|---|---|---|---|
| GSE9750 | 癌(n=31) vs 正常(n=24) | 322.3 vs 166.9（MAS5，↑~1.9倍）✅ | 2.29×10⁻⁵ |
| GSE7803 | SCC(n=21) vs 正常(n=10) | 10.83 vs 10.61（log2，↑趋势） | 0.052 |
| GSE63514 | 癌(n=28) vs 正常(n=24) | 8.17 vs 8.25（无差异）❌ | 0.978 |

- GSE63514组间：Normal 8.25 / CIN1 7.76 / CIN2 8.06 / CIN3 8.76 / Cancer 8.17
- **合理性讨论要点**：GSE63514为激光显微切割纯上皮，TCGA/GTEx为全组织匀浆（正常宫颈富含间质），提示bulk比较可能部分受组织构成影响；GSE9750为组织匀浆、显著上调。成稿Discussion必须如实讨论此不一致，不隐藏GSE63514
- 原始数据：`raw_data/03_GSE63514_...csv`、`04_GSE9750_...csv`、`05_GSE7803_...csv`（均含GSM号逐样本值，探针209017_s_at）

### R9. 泛癌表达概览（D1扩展）✅ → Figure 1E

- 口径：TcgaTargetGtex_rsem_gene_tpm，32癌种（TCGA Primary Tumor vs 同部位GTEx正常+TCGA癌旁；LAML排除因对照为血液；正常n<5不做检验），MWU+BH校正
- 结果：多数癌种肿瘤LONP1高于正常（CESC padj=7.4e-4，与Fig1A同源数据自洽✓；BRCA padj=5.1e-30、UCEC padj=2.2e-7）；**OV方向相反**（肿瘤4.68 vs 正常5.15，padj=3.8e-10）——佐证了不把卵巢癌当"验证队列"的决定，文中如实写"常见但非普遍"
- 原始数据：`raw_data/15_pancancer_LONP1_per_sample.csv`（逐样本）、`16_pancancer_LONP1_stats.csv`（各癌种统计）；脚本`scripts/05_fig_pancancer.R`
- 初稿已同步更新（Results泛癌段落+Methods口径+Figure 1E图注），docx已重新生成

### R10. 单细胞RNA-seq：LONP1细胞类型定位（GSE208653）⚠️ 重要发现，改写文章框架 → Figure 6

**数据口径**：
- GSE208653_RAW.tar，md5=e264e8e98f1005f153c3beaf1f9e29bd，2026-08-29下载自GEO
- 9样本10x 3' v3.1（CellRanger v5.0.1, GRCh38），原作者已过滤(>200基因/细胞, 线粒体UMI<20%)且**矩阵中已移除线粒体基因**（故本分析无法/无需做线粒体QC）
- 样本：NO_HPV正常×2, N_HPV(HPV+正常)×2, HSIL×2, 癌×3(SCC_4/SCC_5/ADC_6)；原文献 Guo C, et al. Clin Transl Med. 2023;13(3):e1219 (PMID 36967539)

**分析口径**（scripts/06_scRNA_GSE208653.py, scanpy 1.11.5, python 3.10）：
- QC：min_genes=300, max n_genes<7500, min_cells=3 → **80,435细胞**
- log1p(CP10K)归一化；HVG 2000(batch_key=sample)；PCA 30PC（自实现float32省内存版，randomized SVD）；neighbors k=15；Leiden res=0.5（28簇）→ marker均值argmax归并9大类
- 注释交叉验证：raw_data/18_cluster_marker_scores.csv（每簇marker分数矩阵，可人工复核）

**结果**：
1. LONP1按细胞类型（log1p均值/阳性率）：上皮0.143/24.9% > 内皮0.127 > 成纤维0.112 > … > T/NK 0.061/4.5%。**上皮富集，约为免疫细胞2倍**
2. 上皮内按阶段：正常HPV- 0.132 / 正常HPV+ 0.142 / HSIL 0.196 / **癌 0.113**。逐细胞KW p=3.3e-23；癌vs正常逐细胞MWU p=1.0e-16（方向：癌更低）；逐样本假批量：癌3样本(0.106/0.119/0.108)为9样本最低，MWU p=0.057
3. **结论：癌上皮细胞LONP1无内在上调；bulk层面"高表达"主要为组织构成效应（肿瘤上皮占比高）。与R8的GSE63514(LCM)阴性结果相互印证。瘤内CNV驱动结论(R6)不受影响。HSIL升高为新观察点**
- 原始数据：raw_data/17_LONP1_per_cell.csv.gz（8万细胞逐细胞值+UMAP坐标+注释）、18(注释依据)、19(上皮逐样本假批量)、20(细胞类型×阶段汇总)
- 检验口径：由17号文件可直接复算全部统计量；从RAW.tar+06脚本可全流程复现
- 图：Fig6A-E（**渲染引擎说明**：沙盒R环境反复损毁后改用plotnine 0.15.8[ggplot2语法Python移植]渲染，脚本07b；07号R脚本为等价标准版，留作在有R环境的机器上复现）

**文章影响（v2初稿已改写）**：标题改为"Bulk and single-cell..."；摘要/结果/讨论/结论全面改写为"上皮富集+组织构成效应+瘤内CNV驱动"框架；新增单细胞Methods小节和Figure 6图注；新增参考文献31-33

### R11. 补充panel与全图字体统一（2026-08-30）

**新增分析（检验口径）**：
- Fig2D/E 临床病理关联（数据CSV02，T按前缀归并T1/T2/T3-4，排除TX/TIS；N仅N0/N1）：T分期KW **p=0.493**、N状态MWU **p=0.999**（ns，与主框架一致，如实报告）→ raw_data/22
- Fig3D 负相关基因富集（CSV06中rho≤-0.35共96基因，Enrichr重新提交）：KEGG Ras signaling padj=3.0e-4 → raw_data/21
- Fig3E 分裂-融合平衡分数 = mean z(FIS1,MFF) − mean z(MFN1,MFN2,OPA1)，vs LONP1 Spearman **rho=0.345, p=1.2e-9** → raw_data/23
- Fig4E 多因素Cox OS（lifelines 0.30.0；LONP1+年龄+T期[T1参照]+N[N0参照]；n=174, events=30）：**LONP1 HR=0.98(0.54-1.81), p=0.957**；N1 HR=2.91(1.40-6.04), p=0.004（阳性对照成立）→ raw_data/24
- Fig5C DC签名散点；Fig5D 检查点基因逐个相关（Bonferroni校正，全部|rho|<0.15）→ raw_data/25
- 单因素Cox经lifelines复算与R/浏览器实现三方一致（OS HR=0.84/DSS 0.99/PFS 0.85）

**全图字体统一**：用户要求Times New Roman。沙盒无TNR授权字体，暂用公制兼容替身**Liberation Serif**统一重渲全部22个panel（Fig1-6）；`scripts/figstyle.py`设计为：用户将times*.ttf放入`温州市科技局结题/fonts/`后重跑脚本即自动切换正版TNR。**待办：等字体文件到位后重渲。**

**绘图引擎现状**：全部panel现由plotnine 0.15.8（ggplot2语法Python移植）渲染（脚本08/09/07b），R脚本02-05/07保留为参考版本；两套脚本统计口径一致（已交叉验证）

**图panel终数**：Fig1×5, Fig2×5, Fig3×5, Fig4×5, Fig5×4, Fig6×5 = 29 panel

### R12. 图片QC与修复（2026-08-30）

- 全29panel拼接总览目检，发现并修复：①Fig3B/3C/3D/5A/5D横向geom_col在plotnine中方向判定错误导致色块重叠错位→改为aes(类别,数值)+coord_flip标准写法，全部重渲验证通过；②KM图例残留变量名"g"→labs(color="")修复；③Fig3E/5B/5C标题补充主题行；④Fig5A的FDR文字标注改为类别名后缀星号（coord_flip下文字对齐不可靠）
- 旧版R渲染的Fig4D_forest.*移入figures/_obsolete/（沙盒无删除权限）
- **数据核对修正**：Fig5D中CD274（PD-L1）rho=-0.171、Bonferroni padj=0.0225为显著（v2初稿原文"|rho|<0.15"表述有误，已修正为如实描述CD274显著、方向与APC签名一致），检验口径见raw_data/25
- 字体：全panel统一Liberation Serif；**用户fonts/目录尚未收到times*.ttf，收到后重跑08/09/07b脚本即切换正版TNR**

### R13. 视觉升级（2026-08-30）

- **新增**：研究设计schematic（scripts/10）→ 复合图Figure1A；组成效应概念图/graphical abstract（scripts/11）→ figures/GraphicalAbstract_composition_effect（备用于投稿graphical abstract或Fig6引导panel，用户可选择手绘替换）
- **图型升级**（scripts/12）：Fig1A(旧)→raincloud；Fig2A→散点+边际直方图；Fig3B/C/D→clusterProfiler风格气泡图；新增Fig6F marker气泡矩阵（检验口径raw_data/26，同时验证细胞注释质量）
- **复合拼版**（scripts/13）：29+2 panel拼为6张双栏整版大图 → `figures/composite/Figure1-6.png/pdf`（2160px=7.2in@300dpi，panel标签A/B/C...）
- **编号变更**：Figure1: A=schematic, B=TCGA+GTEx, C=GSE9750, D=GSE7803, E=GSE63514, F=泛癌；Figure6: A=UMAP注释, B=UMAP LONP1, C=marker矩阵, D=细胞类型violin, E=上皮分期, F=构成。**手稿v2全部图引用已同步**（修复过程中曾误截断build脚本尾部，已恢复并验证重建成功）
- 散panel仍保留在figures/根目录（拼版原料，检验口径不变）

### R14. 补充分析：上皮内在自查 + 组成效应定量检验（2026-08-30）

**A. 上皮内在共表达自查**（scripts/14前置计算，epithelial.h5ad re-cluster: HVG1500/PCA20/leiden0.5→27亚群）：
- 检验单元：样本×亚群假批量（≥50细胞，n=63）；三层检验结果 → raw_data/27（逐细胞/假批量/逐样本）、27b（63单元逐值）
- **正臂复现**：CLPP rho=0.49(p<1e-4)、HSPA9 0.46、FIS1 0.46、MFF 0.34、分裂模块0.44——上皮内在✓
- **负臂不复现**：OPA1 **+0.42**（bulk为-0.22，方向翻转）、MFN1 +0.19 ns、MFN2 -0.02 ns → 抗融合信号限定为bulk现象（可能为样本间构成/CNV联动）
- **稿件修订**：摘要/Results共表达段/讨论段2/Fig3图注 已按"蛋白酶-分裂受体臂=上皮内在；抗融合臂=bulk限定"改写；新panel Fig3F（bulk vs 上皮内在哑铃图）
- 逐细胞层caveat：全基因均弱正相关(0.09-0.15)，为文库复杂度技术性共变，仅作参考

**B. 组成效应bulk内定量检验**：
- 上皮分数=mean z(EPCAM,KRT5,KRT8,KRT17)（Xena取317样本4基因 → raw_data/28）；OLS结果 → raw_data/29
- 肿瘤 vs 正常上皮分数：+0.09 vs -1.99（分离显著）；LONP1肿瘤效应：未校正0.469(p=0.0048) → 校正后0.332(**p=0.149失去显著**)
- **诚实限定**：肿瘤状态与上皮分数共线（rho=0.33, p=1.5e-9），组内相关均ns（正常n=13 rho=-0.25；肿瘤rho=0.03）→ 回归仅作一致性检验，决定性证据仍是上皮分辨率队列；此表述已写入Results与Fig6G图注
- 新panel Fig6G（散点+模型注释）；composite Figure3/6已重拼


### R15. 框架改写：从"线粒体动态"转为LONP1中心（2026-08-30）

**决策（用户确认）**：分裂/融合角度牵强——bulk中动态基因|rho|仅0.13–0.32、DNM1L无关(-0.02)、上皮内所有线粒体基因（含OPA1 +0.42、DNM1L +0.31）同向正相关，属泛线粒体信号而非分裂特异；平衡分数(旧Fig3E)属人为放大弱信号。文章改为围绕LONP1本身：多分辨率表达 + 拷贝数驱动 + 线粒体蛋白稳态/翻译程序；线粒体动态降为探索性panel。未加DepMap依赖性分析（用户暂不选）。

**新增分析（scripts/15_coexpression_positional.py，检验口径）**：
- 全基因组共表达重取：cBioPortal co-expressions接口 threshold=0.0（同profile/样本列表），20,337实体→17,378有HUGO符号基因；与CSV06（|rho|≥0.3, 970基因）逐值一致（最大差5e-16）→ raw_data/30
- 基因位置：Ensembl REST lookup/symbol (GRCh38, 17,324基因定位)；细胞带：UCSC hg38 cytoBand；缓存于raw_data/_cache/
- **位置效应（核心新发现）**：top20正相关伙伴中**18个在19p13.3**（top50中35个）；rho≥0.4的164基因中62个在19p13.3（背景1.1%），Fisher **OR=82.2, p=1.1e-81**；19p13.3基因44%达rho≥0.3 vs 其它基因3.2%；chr19沿位置梯度：19p13.3中位rho 0.24 / 19p臂0.15 / 19q -0.01 / 其余基因组-0.03 → raw_data/32_positional_stats.txt
- **去19p13.3后富集**（rho≥0.4且非19p13.3，n=102，Enrichr同三库）：线粒体膜/内膜 padj=2.3e-12、translation 2.3e-12、线粒体翻译 8.4e-5、线粒体基因表达 1.0e-4、有氧呼吸 8.7e-4、KEGG ribosome 3.3e-8；**KEGG OXPHOS降为ns（padj=0.076）**（含19p13.3时1.0e-3，8个重叠基因中NDUFA11/NDUFS7/UQCR11/ATP5F1D在19p13.3）→ raw_data/31。如实写入正文：OXPHOS信号部分为位置效应
- **线粒体模块（服务器端rho，raw_data/33）**：19p13.3外最强为TRAP1 0.41、PHB1 0.40、PHB2 0.39、HSPD1 0.33、ATF4 0.32、FIS1 0.32、HSPE1 0.30、HSPA9 0.23、DDIT3/PMPCB 0.22、ATF5 0.14(p=0.018)；TFAM/PPARGC1A/PINK1/PRKN/YME1L1/SPG7 ns；DNM1L -0.02；MFN1 -0.23/OPA1 -0.22/MFN2 -0.13。解读为UPRmt样蛋白稳态模块（LONP1/CLPP/HSPD1/HSPE1/HSPA9+ATF4/ATF5/DDIT3为经典UPRmt成分，新增参考文献34-36）
- 上皮内在复现沿用raw_data/27（未重跑scRNA；蛋白稳态臂仅CLPP/HSPA9/HSPD1三基因可用）

**图表**：Figure 3整体重做（scripts/16）：A全基因组秩图(19p13.3标红) / B chr19位置图 / C模块棒棒糖 / D bulk vs上皮分面哑铃 / E去19p13.3富集 / F负相关富集；旧3A/3B/3C/3D/3E/3F六个panel移入figures/_obsolete/。Fig1A示意图文字更新（scripts/10）并重拼Figure1；13_assemble.py字体回退路径修正
**稿件**：v3（scripts/patch_v3.py由v2脚本生成build_manuscript_v3.js）：标题改为"Multi-resolution transcriptomic profiling of LONP1 in cervical cancer: epithelial-enriched expression, 19p13.3 copy-number-driven variation and a mitochondrial proteostasis programme"；摘要/关键词/Background/Methods（共表达+位置注释+软件说明改为R+Python）/Results共表达两节重写+新增上皮内在段/Discussion段1-2/Limitations/Conclusions/Fig3图注同步；新增参考文献34-39（UPRmt：Zhao 2002, Fiorese 2016, Quirós 2017；CNV驱动共表达：Reyal 2005, Fehrmann 2015, Bhattacharya 2020）。全文无"fission-oriented"/平衡分数表述
**数值核对**：69项文中数值与CSV30/31/27/11程序化核对全部一致（修正一处：非19p13.3基因rho≥0.3比例应为3.2%，非2.7%）
**待办**：参考文献34-39建议投稿前用PubMed逐条核对页码；若后续重跑scRNA可把TRAP1/PHB1/PHB2/ATF4纳入上皮内在复现；标书IHC/siRNA预实验是否并入仍待定


### R16. 单细胞重跑：上皮内在共表达扩展到全基因组（2026-08-30，用户指令"开始跑"）

**数据与对齐（scripts/17_scRNA_epithelial_module.py）**：
- 重新下载GSE208653_RAW.tar，md5=e264e8e98f1005f153c3beaf1f9e29bd ✓；按06脚本口径加载/QC → 80,435细胞（逐样本细胞数与CSV17一致）
- 细胞注释沿用CSV17（行序=加载顺序），以逐细胞LONP1值对齐校验：最大差1.2e-3（float32舍入），99.98%细胞差<1e-4 ✓ → 上皮20,363细胞
- **上皮重聚类无法逐单元复现旧R14结果**（旧代码未存档；HVG1500/PCA20/leiden0.5 重跑得25亚群、58个≥50细胞单元，旧为27亚群/63单元；batch_key有无两种变体均不匹配）。**决策：以本次重跑（batch_key=sample变体，代码已存档）为最终检验口径**，稿件中所有上皮内在数值改用新值；旧27/27b保留作历史记录。两次结果方向与量级一致（CLPP 0.49→0.45、HSPA9 0.46→0.45、HSPD1 0.26→0.24、FIS1 0.46→0.51、OPA1 0.42→0.47、MFN2 −0.02→0.05），结论不变；无batch_key变体亦一致（raw_data/36）
- 旧值差异说明：HSPD1由p=0.041变为p=0.067（ns），稿件已如实改写

**全基因组上皮内在共表达（raw_data/37，58单元 × 上皮检出率≥1%的14,140基因）**：
- **整体上移**：假批量rho中位数0.28、IQR 0.10–0.43（单元间检出深度的全局效应）→ 文中一律附百分位解读；按样本中心化敏感性分析中位0.27（0.11–0.41），结论一致
- **上皮内无位置效应**：19p13.3基因中位rho 0.26 vs 其它0.28；THOP1/TIMM13降至43/39百分位 → 位置效应为bulk/CNV特有 ✓
- 扩展模块（假批量rho / 百分位 / p）：POLRMT 0.62/97/2.3e-7、FIS1 0.51/86/5.1e-5、TRAP1 0.49/85/8.7e-5、OPA1 0.47/82、DNAJA3 0.47/82、ATF5 0.47/81、PMPCB 0.46/81、CLPP 0.45/79/3.7e-4、HSPA9 0.45/78/4.0e-4、PPARGC1A 0.45/78、TFAM 0.44/77、PHB(=PHB1) 0.41/71、YME1L1 0.41/73；HSPE1 0.31/55、DNM1L 0.31/55、MFF 0.28/50、HSPD1 0.24/44(p=0.067)、PHB2 0.22/40、MFN1 0.20/37、ATF4 0.16/32、MFN2 0.05/19、DDIT3 0.00/16 → raw_data/35（三层）、34（单元逐值）、39（Fig3D数据）
- **解读修订**：上皮内LONP1与广义"线粒体含量/生物发生状态"共变（蛋白稳态、呼吸链组装、mtDNA转录、分裂与融合基因同向），而非分裂特异或ISR驱动；ATF4/DDIT3成分仅见于bulk（肿瘤层面特征）。稿件Results/Discussion/摘要已按此改写，不再说"蛋白稳态模块在上皮内复现"这种过强表述
- top300上皮内在伙伴富集（raw_data/38，Enrichr）：线粒体呼吸链复合体组装 padj=9.6e-8、ER膜 8.5e-7、线粒体膜 6.5e-5、线粒体内膜 7.6e-3；KEGG均ns；样本中心化版：内膜1.7e-4、复合体组装1.4e-5
- bulk与上皮内在全基因组rho一致性：Spearman 0.26（n=13,385）

**图与稿件**：Fig3D改为扩展模块哑铃图（含上皮全基因组IQR灰带+中位虚线）；新增Fig3F上皮内在富集；负相关富集改为Fig3G；composite Figure3重拼（7 panel）。稿件**v4**（patch_v4.py由v3脚本生成build_manuscript_v4.js）：摘要/Methods单细胞段/Results上皮内在段/Discussion段1-2/Limitations/Fig3图注同步。56项上皮相关数值程序化核对一致（修正一处：58单元来自25亚群中的22个）
**环境**：本次在云端沙盒完成（7GB内存，scanpy 1.11.5）；中间h5ad未回传（可由RAW.tar+脚本17全流程复现）


### R17. 投稿前收尾（2026-08-30，用户指令"先完成剩余项目"；作者信息与预实验并入两项搁置）

- **Fig(旧2)B精确p**：GISTIC五组Kruskal-Wallis H=42.96, **p=1.05e-8**（数据07+08，log2(RSEM+1)），稿件"p<0.001"改为精确值
- **HPV状态**：PanCancer Atlas研究无HPV属性；cesc_tcga（旧版）仅22例有HPV分型（raw_data/41）；Xena临床矩阵无HPV列 → 样本量不足且CESC约95%为HPV+，**不做**，记录于此
- **参考文献核对**（scripts/18_verify_references.py，PubMed E-utilities逐条esearch/esummary；6条自动检索未命中者手工核对）：**39/39全部正确**（仅页码缩写与Quirós重音为Vancouver允许差异）→ raw_data/40_reference_check.csv（含PMID）
- **补充材料**：`【3】文章初稿/supplementary/Additional_file_1_Supplementary_Tables.xlsx`（README + Table S1–S7：S1全基因组共表达+位置(30)；S2富集含/去19p13.3+负相关(11/31/21)；S3上皮内在全基因组rho(37)；S4上皮top300富集(38)；S5 curated基因bulk/上皮三层相关(33+35)；S6免疫签名定义+相关(14)；S7泛癌统计(16)）；正文相应位置已引用，新增"Additional files"节
- **代码仓库包**：`【4】代码仓库/LONP1-cervical-cancer/`（README.md英文复现说明、LICENSE MIT、CITATION.cff、environment.txt、.gitignore、SETUP_zh.md发布步骤；scripts/data/figures/docs由本机复制）。稿件Data availability改为GitHub/Zenodo占位，Software段"available on request"改为"publicly deposited"
- **语言润色**（patch_v5.py，仅表面层）：英/美拼写统一为美式（programme→program、tumour→tumor、centred/neighbour/signalling/labelled/modelled/favouring等，共23处，参考文献不动）；一处事实措辞修正（THOP1/TIMM13非"仅次于CLPP的两个"→"two of the strongest"）；Fig1图注末尾补句号；Fig1B-C→Fig. 1B, C
- **投稿前检查**（submission-checklist逐项，patch_v5b.py修正）：
  - 摘要512词超BMC Cancer上限350 → 压缩至**341词**（数值保留）
  - 图编号未按首次出现顺序（原文顺序1→6→2→3→4→5）→ **重新编号：旧Fig6→2、2→3、3→4、4→5、5→6**；面板首次出现顺序修正（Fig1E在bulk段提前提及；单细胞图"Fig. 2A, C"/"2B, D"；共表达图段落重排：curated基因先于富集，负相关富集句移至上皮段末；面板字母改为C=curated、D=去19p13.3富集、E=bulk vs上皮、F=上皮富集、G=负相关；13_assemble.py布局与输出文件名同步）
  - 参考文献未按首次引用顺序编号（原顺序1,2,3,4,5,30,6,7,9,34…）→ **按首次引用重新编号并重排列表**，程序化核验为1–39连续
  - ORCID占位行；伦理"Not applicable"✓；COI✓；Funding（立项编号待补）✓；关键词8个✓；IMRAD✓；图注完整✓
  - 校验：Fig首次出现顺序 1A-F, 2A-G, 3A-E, 4A-G, 5A-E, 6A-D ✓
- **预警期刊**：中科院《国际期刊预警名单》2020–2025各年均**不含BMC Cancer**（云南旅游职业学院/南财转载页核对）；2026年名单（37本）已发布但网页仅提供网盘下载，未能直接核对，**请投稿前在 ewl.fenqubiao.com 自查一次**
- 图形：Fig1A示意图与Fig(共表达)E面板中"programme/tumour"同步改美式；facet标签"UPRmt / ISR TF"改"ISR TF"（原被截断）
- 稿件 **v5** = v4 + 上述全部；生成链：build_manuscript_v4.js → patch_v5.py → patch_v5b.py → build_manuscript_v5.js


### R18. 三项廉价强化：肿瘤纯度校正 / 外部生存队列 / HPA蛋白证据（2026-08-30，用户确认"可以的"）

**A. TCGA-CESC 肿瘤纯度（scripts/19_purity_adjustment.py）**
- 口径：PanCanAtlas ABSOLUTE 一致性纯度 `TCGA_mastercalls.abs_tables_JSedit.fixed.txt`（GDC PanCanAtlas publications 页，api.gdc.cancer.gov/data/4f277128-…，取 call status=called），与 CSV08/07/28 按样本ID合并 → **n=291**（raw_data/42）；偏相关=两变量对纯度OLS残差后Spearman；结果汇总 raw_data/43
- LONP1 ~ 纯度 rho=**0.118, p=0.045**；上皮4基因分数 ~ 纯度 rho=0.131 (p=0.025)——切除肿瘤纯度范围窄（中位0.67, IQR 0.54–0.79），瘤内纯度梯度只能"边缘地"反映组成效应；如实写为 tumor-vs-normal 现象
- CNA–表达在校正纯度后**更强**：偏相关 rho=0.447 (p=1.2e-15)；OLS β_CNA=0.967 (p=1.1e-19)，β_purity=0.526 (p=0.002)
- 免疫签名校正纯度后仍在：M1 −0.219 (p=1.6e-4)、DC −0.192 (p=1.0e-3)、CD4/Th −0.130 (p=0.026)；注意DC/CD4签名本身与纯度强负相关（−0.54/−0.58），说明LONP1–APC负相关不是纯度伪影 → 讨论措辞相应改写
- 图：FigS1A/B/C → composite/FigureS1（Additional file 2）

**B. 外部生存队列 GSE44001（scripts/20_GSE44001_survival.py）**
- 口径：GEO series matrix（作者log2+quantile归一化），GPL14951 上唯一LONP1探针 **ILMN_1766125**（由GPL14951 family soft核对）；n=300 早期宫颈癌（IB1 217/IIA 42/IB2 28/IA2 13），DFS事件38 → raw_data/44/45
- 单因素Cox HR=**1.07 (0.72–1.60), p=0.727**（每SD HR 1.06）；中位数分组log-rank p=0.828；多因素（+分期IB2/IIA vs IA–IB1、最大径）LONP1 HR=1.03 p=0.875，最大径HR=1.29/cm p=0.002（阳性对照成立）→ TCGA阴性预后结论获独立复现
- 图：FigS2 → composite/FigureS2

**C. Human Protein Atlas 蛋白层面（raw JSON/XML: work/hpa_lonp1.xml，未入库）**
- 抗体 HPA002192（HPA002034无宫颈数据）：正常宫颈鳞状与腺上皮 staining=medium/intensity=moderate（3例）；宫颈癌肿瘤细胞 medium 8/10、not detected 2/10（9例鳞癌+1例腺癌）→ 无肿瘤细胞染色更强的证据，与"无内在上调"一致；HPA对CESC的TCGA预后判定亦为 unprognostic (p=0.205)。仅描述性引用，不放图（半定量、例数少，Limitations已注明）
- **HPV**：仍不可得（见R17）

**稿件 v6**（patch_v6.py：v5→v6）：摘要+"in TCGA or an independent 300-patient cohort"（348词）；Methods新增小节"Tumor purity, external survival cohort and protein-level data"；Results四处插入（单细胞/组成段、CNV段、生存段、免疫段）；Discussion三处；Limitations一处；Availability补GSE44001/GDC/HPA；新增Additional file 2（Fig S1–S2，PDF已生成于supplementary/）；新增参考文献5条（Carter 2012 ABSOLUTE、Taylor 2018 PanCanAtlas aneuploidy、Lee 2013 GSE44001、Uhlén 2015、Uhlen 2017）并按首次引用重排为1–44，PubMed逐条核对通过（raw_data/40更新为44条）
**环境**：lifelines 0.30.3（autograd-gamma需手工安装）、statsmodels


### R20. 全面回顾与复核（2026-08-30，用户指令"回顾、检查、整理全部数据"）

**A. 全文数值程序化复核（scripts/21_full_manuscript_check.py → raw_data/46）**
- 覆盖摘要、Results 全部小节、6 张图注、Additional file 1/2 图注中的 **261 项数值**（样本量、均值/SD、rho、p、HR/CI、富集 padj、百分位、纯度、GSE44001 等），从 raw_data 逐项重算比对
- 结果：**261/261 一致**。过程中发现一处表述不够明确：单细胞逐样本假批量 MWU p=0.057 是"癌 3 样本 vs 正常 4 样本"（若与全部 6 个非癌样本比较则 p=0.024），稿件已改为写明比较对象
- 参考文献 44 条 PubMed 核对（raw_data/40）全部通过

**B. 一致性核对（图注 ↔ 版面 ↔ 正文）**
- 发现 **Figure 6 面板顺序与图注不符**（版面 B=检查点基因、C=M1、D=DC；图注与正文为 B=M1、C=DC、D=检查点）→ 修改 13_assemble.py 版面为 A 签名条形 / B M1 / C DC / D 检查点，重拼 Figure6（composite/Figure6.png/pdf 已更新）
- Figure 1–5、S1–S2 图注与版面逐面板核对一致；图与参考文献首次出现顺序均为连续（1A–F, 2A–G, 3A–E, 4A–G, 5A–E, 6A–D；引用 1–44）
- Methods 与实际口径核对：软件版本改为 "Python 3.10–3.11（加 statsmodels）"（纯度分析用 statsmodels OLS；本轮环境 3.11）

**C. 语言细节（patch_v7.py，v6→v7）**
- 语法错误 "were with the chaperones" → "were the chaperones"（v5b 段落重排遗留）
- Fig 4C 图注 "server-side Spearman" → "Spearman correlation computed by cBioPortal"；Fig 5E 图注 p=0.957 → 0.96（与正文一致）
- 单细胞注释句改写（避免循环表述）；摘要 Methods 提及独立生存队列（348 词）

**D. 数据整理**
- 新建 `raw_data/README_data_index.md`：46 个数据文件的数据字典（内容/来源口径/生成脚本/对应稿件图表）+ 脚本索引 + 稿件生成链；标注 06/23/27/27b/41 为历史或未用文件（保留以便复现旧版）
- figures/ 已于 R19 清理；composite 为 Figure1–6 + S1–S2
- 仓库包同步：scripts（含 21、patch_v7、build_manuscript_v7.js）、data（46 + README）、figures/composite/Figure6、docs

**当前状态**：稿件 **v7**（内容与 v6 相同，仅细节修正）。待用户：作者信息/立项编号/ORCID；GitHub+Zenodo 发布；字体文件；预警名单自查


## 拟定图表规划（R绘图，单panel输出+caption，用户手动拼图）

- Fig1 表达：A TCGA+GTEx小提琴/箱线图；B-D GEO三队列箱线图（含GSE63514阴性结果）
- Fig2 调控：A 泛癌表达概览(待做)；B 亚型比较；C CNA相关散点；D GISTIC分组箱线
- Fig3 线粒体动态关联：A 关键基因相关性棒棒糖/热图；B 共表达火山式散点；C-D GO/KEGG富集条形图
- Fig4 生存（阴性结果）：A-C OS/DSS/PFS KM曲线；D 森林图
- Fig5 免疫：签名相关性条形图 + 代表散点
- 每张图caption注明：数据来源、n、统计方法、精确p值

## 待办清单

- [x] 泛癌LONP1表达概览（R9）
- [x] R复算全部数值（R复算与图表节）
- [x] 共表达全表落盘（raw_data/06→30）
- [x] 标书核对（预期成果=SCI 1篇）
- [x] BMC Cancer格式稿件（v7）；基金号待用户提供
- [x] 预警名单（2020–2025不含BMC Cancer；2026名单请自查）+ submission checklist（R17）
- [ ] 作者/单位/邮箱/ORCID/立项编号（用户）
- [ ] GitHub+Zenodo发布并回填链接（用户）
- [ ] Times New Roman字体文件到位后重渲全部panel

## 诚信声明

所有数据来自真实公开数据库，按上述口径可完全复现。不显著/不一致结果（R2生存、R8-GSE63514）如实记录并将写入成稿。不进行选择性报告、最优截断值捕捞或数据操纵。

## R复算与图表（2026-08-29晚）

- 计算沙盒恢复；micromamba安装R 4.3.3 + ggplot2/survival/svglite（conda-forge）
- `scripts/00_fetch_cbioportal.R`：重取cBioPortal数据→raw_data/06(共表达970基因,与浏览器一致)/07(CNA)/08(基因panel表达)/10(亚型)
- `scripts/01_verify_stats.R`→`raw_data/verification_report.txt`：**R复算与浏览器计算逐项一致**（R1 t=4.776/p=2.262e-4；R2 OS HR=0.843/p=0.4508等）。差异说明：中位数分组log-rank在R中排除t=0样本后chi2=0.626/p=0.429（原0.795/0.373），结论不变；亚型差异精确p=0.0216。成稿一律采用R复算值
- GSE63514补充：Kruskal-Wallis 5组 p=0.348；CIN3 vs Normal p=0.243（均ns）
- `scripts/02-04`：13个图panel（PDF+300dpi PNG）→ `figures/`，图注见 `figures/figure_captions.md`（含每图检验口径速查表）
- 衍生数据落盘：raw_data/11(富集)/12(线粒体panel相关)/13(Cox森林数据)/14(免疫签名相关)

## 时间线

- 2026-08-29：环境故障→改浏览器方案；完成R1-R8全部核心分析；原始数据01-05落盘
- 2026-08-29晚：R环境搭建；R复算全部通过；13个图panel完成；标书核对完成（预期成果=SCI 1篇；标书无立项编号，待用户提供）
- 2026-08-29晚：BMC Cancer格式英文初稿完成→`【3】文章初稿/LONP1_cervical_cancer_manuscript_v1.docx`（正文数值均取自R复算/落盘CSV；含阴性结果与GSE63514不一致的完整讨论；生成脚本`scripts/build_manuscript.js`）。待补：作者/单位/邮箱、基金立项编号、预实验数据（如有则增加验证章节）

---

## R22 2026-09-03 用户端 v7 手改稿状态核对（参考文献迁移中）
- 背景：用户在本机对 `【3】文章初稿/LONP1_cervical_cancer_manuscript_v7.docx` 做了手动修改（mtime 2026-09-02，42,405 字节；容器内构建版为 30,807 字节）。暂存后与容器 v7 逐行比对。
- 已完成的手改：
  1. 标题改为无冒号版 "Multi-resolution transcriptomic analysis of LONP1 in cervical cancer reveals an epithelial-enriched, copy-number-driven mitochondrial proteostasis program"。
  2. 摘要内 4 处科学计数法改为真上标（7.2×10⁻⁴、2.3×10⁻⁵、3.1×10⁻¹³、1.1×10⁻⁸¹）。
  3. Background 全部引用改由 Zotero 域代码插入（document.xml 含 7 处 ADDIN ZOTERO_ITEM、6 条 CSL_CITATION），并重编号为新 1–9；其中旧 2/3 换成 Schiffman 2007、Tewari 2025，旧 5/6/9/10/12/15/16 未保留，正文相应改为只提 ATF5。
- 尚未完成 / 当前不一致处（程序化核对）：
  - 参考文献表仅 9 条，而正文最大引用号仍为 44；Methods/Results/Discussion 共 22 处引用位点、涉及 34 条旧文献仍是旧编号。
  - Discussion 仍引用 [10-12]（ATF4/ATF5 系列）与 [13-16]（LONP1 肿瘤系列），但这些文献已在 Background 精简时移出文献表，需在 Zotero 中重新加入或改写该处句子。
  - 正文余部仍有 31 处 `x10-`、4 处 `+/-`、7 处 `>=`、3 处 `<=` 未转符号；`rho=` 19 处未转 ρ。
- 产出：`【3】文章初稿/参考文献补插对照表.md`（按 Methods/Results/Discussion 顺序列出每个待插引用位点、上下文与对应旧文献全称，另附去重旧文献清单）。
- 处理原则：因稿件已含 Zotero 活动域代码，**不再对该 docx 做程序化改写**（会破坏域），后续编号与符号统一由用户在 Word/Zotero 内完成；容器内 v7 仅作比对基线保留。

---

## R23 2026-09-06 补回 Fig 5 阴性结果声明句；确认参考文献迁移已完成
- 背景：用户已在本机通过 Zotero 完成 Background 之后全部引用的重新插入（对照 `参考文献补插对照表.md` 手工操作）。程序化核对：正文引用号与文献表已 1–30 完全对齐，无孤儿引用，document.xml 中 ZOTERO 域代码 21 处、CSL_CITATION 20 条，均为活动域。
- 发现 Fig 5E 段落中原有的诚信声明句在手动编辑过程中被意外删除：
  "We report this negative result in full; no alternative cut-points were explored."
  该句对应课题组"如实报告、不做 optimal-cutoff p-hacking"的原则，在 Results 中的落地表述。
- 处理：仅对 `word/document.xml` 中该处纯文本 `<w:t>` run 做字符串替换插回原句，不触碰任何 `w:fldChar`/`w:instrText`（Zotero 域代码），插回后 pandoc 转纯文本核对，diff 确认仅新增该句，其余内容逐字符一致；ZOTERO/CSL_CITATION 计数不变，域代码未受影响。已写回 OneDrive 原文件（同路径覆盖，expectedMtimeMs 校验通过，无并发冲突）。
- 记录在案、尚待用户后续处理的格式一致性问题（手动编辑引入，未改动）：
  - "microdissected" / "micro dissected" 混用（3 vs 3）
  - "pseudobulk" / "pseudo bulk" 混用（3 vs 5）
  - "proteostasis" / "protein homeostasis" 混用（4 vs 5）
  - "TCGA+GTEx" / "TCGA + GTEx" 间距不统一（3 vs 4）
  - "in silico" / "in-silico" 连字符不统一（各 1 处）

---

## R24 2026-09-20 Methods 逐节同行评议式复核（Single-cell → Statistics），含比例风险假设补充分析

**工作方式**：用户逐条提出【修改建议/理由/建议改为】，每条先核对代码与原始数据再判断同意与否，获用户确认后以 Word 审阅模式（作者 Claude）改写稿件；全程只改纯文本 `<w:t>`，不触碰 `w:fldChar`/`w:instrText`。**改写后核验：ZOTERO 域 29 处、CSL_CITATION 29 条、fldChar 30 组，与改写前逐项一致，Zotero 活动域未受影响**（R23 禁令未被违反）。每轮均通过 docx 往返校验（去除修订后与原文逐字符一致，段落数 132 不变）并渲染 PDF 目视确认。

**A. 术语与口径更正（据代码核对，非推定）**
- `pseudobulk` → 按实际算法改为 sample-level / mean-expression profile（方法为对 log 归一化值取均值，非计数汇总）
- Enrichr 背景：脚本用标准 `addList`+`enrich`，未上传自定义背景 → 正文如实写为各文库默认背景
- 百分位排名：原文把全局正偏移归因于 transcript detection，脚本未做该定量验证 → 改为可能性表述，并注明百分位不校正技术混杂
- 样本中心化：原文称去除 stage/copy-number 差异，实际仅去除样本均值 → 改为如实表述
- 线粒体面板显著性星号用的是名义 p（`res$sig <- cut(res$p, ...)`，无 p.adjust）→ 注明为名义 p、探索性
- `partial Spearman`：脚本对**原值**回归取残差再算 Spearman，非先取秩的标准偏相关 → Results 两处改名，Fig S1C 标题重画（脚本 30）。两法差异极小（CNA 0.447 vs 0.430；M1 −0.219 vs −0.209）
- HPA：实验记录中 staining=medium / intensity=moderate 两字段，正文误将等级值 medium 安在 intensity 上 → 改为 protein expression level，并注明正常仅 3 例、非配对描述性比较
- `All tests were two-sided`：`fisher_exact(..., alternative="greater")` 为单侧、Enrichr 上尾、KW 为总体检验 → 按实际分述

**B. 比例风险假设（PH）诊断 —— 本轮唯一实质性发现**
- 此前从未检验过 PH（全目录检索 Schoenfeld/cox.zph 零命中）。补做后：**OS p=0.021、PFS p=0.0041 违反 PH**；DSS p=0.093 未检出；多因素 OS 全局 p=0.35 未检出
- 时变系数模型 log HR(t)=b0+b1(log t − c)：三个 TCGA 终点 HR 均随时间单调下降、约在中位事件时间穿过 1（OS 6 月 1.55 → 60 月 0.45，p=0.013；PFS 1.42 → 0.24，p=0.0064；DSS p=0.089）；按中位事件时间分段（前后各 34 事件）方向一致
- 线性假设无问题（样条 LRT p≥0.77）；**外部队列 GSE44001 未见 PH 违反（p=0.75/0.88），时变项 p=0.56 且方向相反 → TCGA 的时间依赖模式未获独立复现**
- 首版时变模型因 LONP1（均值 10.6）与自身交互项共线产出 HR(3月)=160 的伪结果，中心化后重算方得稳定估计（已在脚本 24 注释中记录该错误）
- 定位：事后诊断、晚期风险集小（60 月时 OS/PFS 仅 41/36 例）、DSS 未显著、三终点同批病例 → 作探索性报告；正文全面改为 **time-constant** 口径（摘要、Results 小标题与段落、Discussion、Conclusions、图 5 图注），并说明原 HR 为时间平均效应，不等于"无关联"
- Limitations 补两条：多因素模型排除 107/281 例（其中 105 例因 path_N 为 NX 或缺失，与是否手术分期相关，非随机缺失）；60 月后估计依赖不足 45 例

**C. 其他核对与补充**
- 多因素 Cox 变量编码与排除规则补入 Methods（T1 亚期合并、N 仅取 N0/N1、n=174/30 事件）
- GSE44001 直径为连续变量（cm）；发现 14 例直径记录为 0（IB1 9、IA2 5，事件数 0），剔除后 HR 1.294→1.274（p=0.0019→0.0047），结论不变（脚本 29）
- 17,378 与 17,324 之差为 54 个缺基因组定位的基因，已在 Methods 说明
- LONP1 自身未进共表达排名（cBioPortal 接口不返回查询基因；单细胞侧 `.drop("LONP1")`），已在正文注明
- HPA 原始 XML（`work/hpa_lonp1.xml`）已丢失，按实验记录转录为 `raw_data/54_HPA_LONP1_staining.csv` 入库，注明来源性质

**D. 产出**
- 脚本：`22_PH_assumption_schoenfeld.py`、`23_PH_step1_data_checks.py`、`24_PH_step2_timevarying.py`、`25_PH_schoenfeld_plots.py`、`26_PH_GSE44001.py`、`27_figS3_PH_diagnostics.py`、`28_tableS8_PH_results.py`、`29_GSE44001_diameter_sensitivity.py`、`30_regen_figS1C.py`
- 数据：`raw_data/47`–`54`
- 图表：**Fig S3**（A 缩放 Schoenfeld 残差；B HR(t) 及 95%CI，含 GSE44001 对照）；**Table S8**（68 行，已写入 Additional file 1，原表备份为 `_backup_before_TableS8.xlsx`）；Fig S1C 重画
- `environment.txt` 改为 A（原始分析）/ B（补充分析）两段；lifelines 记为 0.30.0→0.30.3（R18/R19 两轮），statsmodels 原始版本未记录，但 Environment B 以 0.15.0 复算得到相同的纯度校正估计（0.4472）

**E. 待办**
- Fig S3 尚未并入 `Additional_file_2_Supplementary_Figures.pdf`
- 代码仓库需推送 GitHub 并发 Release 生成新 Zenodo 版本；稿件 Availability 现引 DOI 10.5281/zenodo.22840982 对应旧内容，发版后需更新（或改引 concept DOI）
