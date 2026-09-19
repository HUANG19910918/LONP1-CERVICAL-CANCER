# 11_graphical_abstract.py —— "组成效应"概念图 / graphical abstract
# 核心信息：癌上皮细胞LONP1并未内在上调；bulk层面的"高表达"来自肿瘤组织上皮占比更高
import sys, os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch, Rectangle

base = sys.argv[1]; fig_dir = os.path.join(base, "figures", "子图原文件"); os.makedirs(fig_dir, exist_ok=True)
for p in glob.glob(os.path.join(os.path.dirname(base), "fonts", "*.ttf")):
    font_manager.fontManager.addfont(p)
names = {f.name for f in font_manager.fontManager.ttflist}
FAM = "Times New Roman" if "Times New Roman" in names else "Liberation Serif"
matplotlib.rcParams["font.family"] = FAM; matplotlib.rcParams["pdf.fonttype"] = 42

EPI = "#E64B35"; STROMA = "#9DBAD1"; IMMUNE = "#C7CBE0"; TXT="#333333"
fig, ax = plt.subplots(figsize=(10.6, 5.0), dpi=300)
ax.set_xlim(0, 106); ax.set_ylim(0, 50); ax.axis("off")

def txt(x,y,s,size=8.5,weight="normal",color=TXT,ha="center",va="center",style="normal"):
    ax.text(x,y,s,fontsize=size,fontweight=weight,color=color,ha=ha,va=va,family=FAM,
            style=style,linespacing=1.35,zorder=6)

rng = np.random.default_rng(7)
def tissue(x0, y0, w, h, frac_epi, title, n=90):
    ax.add_patch(FancyBboxPatch((x0,y0), w, h, boxstyle="round,pad=0,rounding_size=1.6",
                 fc="#FBF9F7", ec="#8C7B6B", lw=1.3, zorder=1))
    txt(x0+w/2, y0+h+2.2, title, 9.5, "bold")
    # 网格+抖动布点
    cols, rows = 10, 9
    pts = []
    for i in range(cols):
        for j in range(rows):
            px = x0+2.2 + (w-4.4)*i/(cols-1) + rng.uniform(-0.7,0.7)
            py = y0+2.2 + (h-4.4)*j/(rows-1) + rng.uniform(-0.7,0.7)
            pts.append((px,py))
    rng.shuffle(pts)
    n_epi = int(len(pts)*frac_epi)
    for k,(px,py) in enumerate(pts):
        if k < n_epi: fc, r = EPI, 1.05
        elif k < n_epi + (len(pts)-n_epi)//2: fc, r = STROMA, 0.9
        else: fc, r = IMMUNE, 0.8
        ax.add_patch(Circle((px,py), r, fc=fc, ec="white", lw=0.5, alpha=0.95, zorder=2))
    return n_epi, len(pts)

nN, tot = tissue(3, 11, 26, 27, 0.28, "Normal cervix")
nT, _   = tissue(33, 11, 26, 27, 0.72, "Cervical cancer")
txt(16, 8.3, "epithelium ~30%", 8, color="#7A6A5B")
txt(46, 8.3, "epithelium ~70%", 8, color="#7A6A5B")

# 图例
ax.add_patch(Circle((6.5, 4.2), 1.0, fc=EPI, ec="white", lw=0.5))
txt(8.3, 4.2, "epithelial cell - LONP1-high", 7.8, ha="left")
ax.add_patch(Circle((34, 4.2), 0.9, fc=STROMA, ec="white", lw=0.5))
ax.add_patch(Circle((36, 4.2), 0.8, fc=IMMUNE, ec="white", lw=0.5))
txt(37.8, 4.2, "stromal / immune cell - LONP1-low", 7.8, ha="left")

# 中间等式说明
txt(31, 44.5, "Per-cell LONP1 in epithelium:  cancer = normal   (LCM microarray + scRNA-seq, 80,435 cells)", 8.8, "bold")

# 右侧：bulk测量的箭头与柱状示意
arrow = FancyArrowPatch((60.5, 24.5), (67.5, 24.5), arrowstyle="-|>", mutation_scale=14,
                        color="#8C7B6B", lw=1.8, zorder=3)
ax.add_patch(arrow)
txt(64, 27.3, "bulk\nRNA-seq", 7.6, color="#7A6A5B")

# bar示意
bx = 70; bw = 7.5
ax.add_patch(Rectangle((bx, 13), bw, 10, fc=STROMA, ec="none", alpha=0.85, zorder=2))
ax.add_patch(Rectangle((bx, 23), bw, 4.5, fc=EPI, ec="none", zorder=2))
ax.add_patch(Rectangle((bx+10, 13), bw, 4.5, fc=STROMA, ec="none", alpha=0.85, zorder=2))
ax.add_patch(Rectangle((bx+10, 17.5), bw, 12.5, fc=EPI, ec="none", zorder=2))
ax.plot([bx-1.5, bx+2*10+bw-1], [13,13], color="#555555", lw=1.1)
txt(bx+bw/2, 10.8, "Normal", 8)
txt(bx+10+bw/2, 10.8, "Tumor", 8)
txt(bx+10+bw/2+9.5, 21.5, '"up-regulated"\nin bulk', 8.2, style="italic", color=EPI)
txt(bx+9, 33.5, "Measured bulk LONP1", 8.6, "bold")
txt(bx+9, 6.0, "signal = per-cell level x epithelial fraction", 7.8, color="#7A6A5B")

# 底部结论条
ax.add_patch(FancyBboxPatch((3, 0.2), 100, 1.2, boxstyle="round,pad=0.35,rounding_size=0.8",
             fc="#FDECEA", ec=EPI, lw=1.2, zorder=4))
txt(53, 0.85, "Apparent LONP1 “overexpression” in cervical cancer bulk tissue is a tissue-composition effect, not cell-intrinsic induction",
    8.8, "bold")

plt.tight_layout(pad=0.3)
for ext in ("png","pdf"):
    fig.savefig(os.path.join(fig_dir, f"GraphicalAbstract_composition_effect.{ext}"),
                bbox_inches="tight", facecolor="white", dpi=300)
print("ABSTRACT_DONE font:", FAM)
