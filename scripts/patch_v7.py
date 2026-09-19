# patch_v7.py —— 由 build_manuscript_v6.js 生成 build_manuscript_v7.js：全文复核后的细节修正（依据 实验记录.md R20）
#  (1) 语法错误 "were with the chaperones" (2) 假批量 MWU 比较对象写明（vs 四个正常样本）(3) 单细胞注释句改写
#  (4) 软件版本与 statsmodels (5) 图注: Fig4C 术语、Fig5E p 值统一、Fig6 面板顺序与版面一致 (6) 摘要 Methods 提及独立生存队列
import sys
src, dst = sys.argv[1], sys.argv[2]
lines = open(src, encoding="utf-8").read().split("\n")
def find(prefix):
    hits = [i for i, l in enumerate(lines) if l.strip().startswith(prefix)]
    assert len(hits) == 1, (prefix, hits); return hits[0]
def sub(prefix, old, new):
    i = find(prefix); assert old in lines[i], (prefix, old[:70]); lines[i] = lines[i].replace(old, new, 1)

sub('P("Among curated mitochondrial genes (Fig. 4C',
    "those located outside 19p13.3 that correlated most strongly with LONP1 were with the chaperones TRAP1",
    "those located outside 19p13.3 that correlated most strongly with LONP1 were the chaperones TRAP1")
sub('P("We analyzed 80,435 single cells',
    "per-sample pseudobulk means placed the three cancer samples lowest among the nine samples, Mann-Whitney p=0.057)",
    "per-sample pseudobulk means placed the three cancer samples lowest among the nine samples, Mann-Whitney versus the four normal samples p=0.057)")
sub('P("We analyzed 80,435 single cells',
    "identified nine major cell types whose canonical marker profiles confirmed the annotation (Fig. 2A, C),",
    "identified nine major cell types, with canonical marker profiles supporting the annotation (Fig. 2A, C),")
sub('P("Statistical analyses were performed in R 4.3.3',
    "Python 3.10 (pandas, SciPy, lifelines 0.30, Scanpy 1.11)",
    "Python 3.10-3.11 (pandas, SciPy, statsmodels, lifelines 0.30, Scanpy 1.11)")
sub('P("Figure 4. The LONP1 co-expression landscape', "(server-side Spearman, n=294;", "(Spearman correlation computed by cBioPortal, n=294;")
sub('P("Figure 5. Survival (negative result).', "LONP1 HR=0.98 (0.54-1.81, p=0.957)", "LONP1 HR=0.98 (0.54-1.81, p=0.96)")
sub('P("Figure 6. Immune marker signatures.',
    "(B) LONP1 vs M1 macrophage signature. (C) LONP1 vs dendritic cell signature. (D) Individual immune checkpoint genes;",
    "(B) LONP1 vs M1 macrophage signature. (C) LONP1 vs dendritic cell signature. (D) Individual immune checkpoint genes;")  # 版面已按此顺序重拼（13_assemble.py）
sub('P("Methods: We integrated TCGA-CESC and GTEx RNA-seq data',
    "single-cell RNA-seq of the cervical carcinogenesis spectrum (GSE208653; 80,435 cells)", "single-cell RNA-seq (GSE208653; 80,435 cells)")
sub('P("Methods: We integrated TCGA-CESC and GTEx RNA-seq data',
    "immune-signature correlations and survival analyses were also performed.",
    "immune-signature correlations and survival analyses (TCGA and an independent cohort) were also performed.")
open(dst, "w", encoding="utf-8").write("\n".join(lines)); print("PATCHED(v7) ->", dst)
