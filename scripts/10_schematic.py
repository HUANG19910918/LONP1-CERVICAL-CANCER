# 10_schematic.py —— Figure 1A 研究设计示意图（matplotlib手绘）
import sys, os, glob
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

base = sys.argv[1]; fig_dir = os.path.join(base, "figures", "子图原文件"); os.makedirs(fig_dir, exist_ok=True)
for p in glob.glob(os.path.join(os.path.dirname(base), "fonts", "*.ttf")):
    font_manager.fontManager.addfont(p)
names = {f.name for f in font_manager.fontManager.ttflist}
FAM = "Times New Roman" if "Times New Roman" in names else "Liberation Serif"
matplotlib.rcParams["font.family"] = FAM; matplotlib.rcParams["pdf.fonttype"] = 42

C = dict(bulk="#4DBBD5", epi="#00A087", omic="#8491B4", mod="#F5F0EB",
         modb="#B09C85", res="#E64B35", txt="#333333")

fig, ax = plt.subplots(figsize=(11.5, 4.6), dpi=300)
ax.set_xlim(0, 100); ax.set_ylim(0, 46); ax.axis("off")

def box(x, y, w, h, fc, ec, lw=1.0, r=1.2, alpha=1.0):
    b = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                       fc=fc, ec=ec, lw=lw, alpha=alpha, zorder=2)
    ax.add_patch(b); return b

def txt(x, y, s, size=8.4, weight="normal", color=C["txt"], ha="center", va="center", style="normal"):
    ax.text(x, y, s, fontsize=size, fontweight=weight, color=color, ha=ha, va=va,
            family=FAM, style=style, zorder=3, linespacing=1.35)

def arrow(x1, y1, x2, y2, color="#999999", lw=1.4):
    a = FancyArrowPatch((x1,y1),(x2,y2), arrowstyle="-|>", mutation_scale=11,
                        color=color, lw=lw, zorder=1, shrinkA=2, shrinkB=2)
    ax.add_patch(a)

# ===== 左列：数据层 =====
txt(14, 44.3, "DATA  (all public)", 9.5, "bold")
# bulk cohorts
box(2, 33.5, 24, 8.6, "#EAF6FA", C["bulk"], 1.2)
txt(14, 40.2, "Bulk tissue - 4 cohorts", 8.8, "bold")
txt(14, 36.6, "TCGA-CESC + GTEx  (304 T / 13 N)\nGSE9750 (31 T / 24 N)   GSE7803 (21 T / 10 N)", 7.4)
# epithelium-resolved
box(2, 22.6, 24, 9.2, "#E8F5F1", C["epi"], 1.2)
txt(14, 29.9, "Epithelium-resolved - 2 cohorts", 8.8, "bold")
txt(14, 25.8, "GSE63514  LCM epithelium (128 samples)\nGSE208653  scRNA-seq, 80,435 cells\n(normal / HPV+ / HSIL / cancer)", 7.4)
# multi-omics
box(2, 13.6, 24, 7.2, "#EFF0F6", C["omic"], 1.2)
txt(14, 19.2, "TCGA multi-omics", 8.8, "bold")
txt(14, 16.2, "copy number - mutation - subtype\nclinical follow-up (OS / DSS / PFS)", 7.4)

# ===== 中列：分析模块 =====
txt(50, 44.3, "ANALYSIS", 9.5, "bold")
mods = [
    ("Expression across resolutions", "tumor vs normal at bulk, LCM and single-cell level", 37.4),
    ("Genomic & clinical correlates", "CNA / GISTIC / mutation - subtype - T & N stage", 30.6),
    ("Co-expression landscape", "genome-wide Spearman - 19p13.3 positional effect - GO / KEGG", 23.8),
    ("Survival modelling", "KM - univariate & multivariable Cox (age, T, N)", 17.0),
    ("Immune landscape", "13 marker signatures - checkpoint genes", 10.2),
]
for title, sub, y in mods:
    box(37, y, 26, 5.6, C["mod"], C["modb"], 1.0)
    txt(50, y+3.9, title, 8.3, "bold")
    txt(50, y+1.7, sub, 6.9, color="#666666")

# ===== 右列：结论 =====
txt(85.5, 44.3, "KEY FINDINGS", 9.5, "bold")
res = [
    ("LONP1 is epithelial-enriched;\nbulk elevation = tissue composition", 36.4, True),
    ("No cell-intrinsic induction in\nmalignant epithelium (LCM + scRNA)", 29.6, True),
    ("Intratumoral variation driven by\n19p13.3 copy-number gain", 22.8, False),
    ("Mitochondrial proteostasis / translation\nprogram beyond the 19p13.3 locus", 16.0, False),
    ("Not prognostic; weak inverse link\nto APC signatures", 9.2, False),
]
for s, y, hot in res:
    box(72, y, 27, 5.8, "#FDECEA" if hot else "#FAFAFA", C["res"] if hot else "#BBBBBB", 1.2 if hot else 0.9)
    txt(85.5, y+2.9, s, 7.4, "bold" if hot else "normal")

# ===== 箭头 =====
for y in (37.8, 27.2, 17.2):
    arrow(26.5, y, 36.5, y)
for _,_,y in mods:
    arrow(63.5, y+2.8, 71.5, y+2.8)

plt.tight_layout(pad=0.4)
for ext in ("png","pdf"):
    fig.savefig(os.path.join(fig_dir, f"Fig1_schematic.{ext}"), bbox_inches="tight",
                facecolor="white", dpi=300)
print("SCHEMATIC_DONE font:", FAM)
