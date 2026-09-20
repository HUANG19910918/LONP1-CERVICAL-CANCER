# 13_assemble.py —— 把panel拼成整版复合大图（BMC双栏宽度7.2in @300dpi）
# 输入 figures/子图原文件/Fig*.png（散panel）  输出 figures/composite/Figure1..6 + FigureS1..S2 (.png/.pdf)
import sys, os
from PIL import Image, ImageDraw, ImageFont

base = sys.argv[1]
figroot = os.path.join(base, "figures")
# 散panel源目录：2026-09 目录重组后移到 【1】图片源文件/(1) 子图原文件，保留旧路径作回退
fig = os.path.join(figroot, "【1】图片源文件", "(1) 子图原文件")
if not os.path.isdir(fig):
    fig = os.path.join(figroot, "子图原文件")
out = os.path.join(figroot, "composite")
os.makedirs(out, exist_ok=True)
W = 2160  # 7.2in * 300dpi
MARG, GAP = 26, 22
FB = None
for cand in [os.path.join(os.path.dirname(base), "fonts", "timesbd.ttf"),
             "/usr/share/fonts/truetype/liberation2/LiberationSerif-Bold.ttf",
             "/usr/share/fonts/truetype/liberation/LiberationSerif-Bold.ttf"]:
    try:
        FB = ImageFont.truetype(cand, 64); break
    except Exception:
        continue

def load(name): return Image.open(os.path.join(fig, name))

# 编号说明（2026-08-30 v5）：复合图按稿件正文首次出现顺序编号。panel 文件名保留原始编号前缀（Fig6*=稿件Figure 2 等），
# 映射：稿件Figure1=Fig1*/schematic；Figure2=Fig6*（单细胞）；Figure3=Fig2*（基因组关联）；Figure4=Fig3*（共表达）；
#       Figure5=Fig4*（生存）；Figure6=Fig5*（免疫）
LAYOUTS = {
 "Figure1": [
   [("Fig1_schematic.png", 1.0, "A")],
   [("Fig1A_TCGA_GTEx_expression.png", 0.27, "B"), ("Fig1B_GSE9750_expression.png", 0.24, "C"),
    ("Fig1C_GSE7803_expression.png", 0.24, "D"), ("Fig1D_GSE63514_expression.png", 0.25, "E")],
   [("Fig1E_pancancer_LONP1.png", 1.0, "F")],
 ],
 "Figure2": [   # 单细胞（原Figure6）
   [("Fig6A_UMAP_celltype.png", 0.55, "A"), ("Fig6B_UMAP_LONP1.png", 0.45, "B")],
   [("Fig6F_marker_dotplot.png", 1.0, "C")],
   [("Fig6C_violin_celltype.png", 0.40, "D"), ("Fig6D_epithelial_stage.png", 0.32, "E"),
    ("Fig6E_composition.png", 0.28, "F")],
   [("Fig6G_composition_check.png", 0.55, "G")],
 ],
 "Figure3": [   # 基因组/临床关联（原Figure2）
   [("Fig2A_CNA_expression_scatter.png", 0.36, "A"), ("Fig2B_expression_by_GISTIC.png", 0.34, "B"),
    ("Fig2C_subtype.png", 0.30, "C")],
   [("Fig2D_T_stage.png", 0.52, "D"), ("Fig2E_N_stage.png", 0.48, "E")],
 ],
 "Figure4": [   # 共表达全景（原Figure3；字母按正文首次出现顺序：C=curated, D=去19p13.3富集, E=bulk vs上皮, F=上皮富集, G=负相关）
   [("Fig3A_coexp_rank.png", 0.50, "A"), ("Fig3B_chr19_positional.png", 0.50, "B")],
   [("Fig3C_mito_module_lollipop.png", 0.54, "C"), ("Fig3E_enrichment_excl19p13.png", 0.46, "D")],
   [("Fig3D_bulk_vs_intrinsic.png", 0.46, "E"), ("Fig3F_epithelium_enrichment.png", 0.54, "F")],
   [("Fig3G_neg_enrichment.png", 0.46, "G")],
 ],
 "Figure5": [   # 生存（原Figure4）
   [("Fig4A_KM_OS.png", 0.34, "A"), ("Fig4B_KM_DSS.png", 0.33, "B"), ("Fig4C_KM_PFS.png", 0.33, "C")],
   [("Fig4D_forest_univariate.png", 0.50, "D"), ("Fig4E_forest_multivariable.png", 0.50, "E")],
 ],
 "Figure6": [   # 免疫（原Figure5）；v7: 面板顺序与图注/正文一致 A=签名条形, B=M1, C=DC, D=检查点基因
   [("Fig5A_immune_bar.png", 0.55, "A"), ("Fig5B_M1_scatter.png", 0.45, "B")],
   [("Fig5C_DC_scatter.png", 0.50, "C"), ("Fig5D_checkpoint_genes.png", 0.50, "D")],
 ],
 "FigureS1": [   # 补充图S1：肿瘤纯度（scripts/19）
   [("FigS1A_purity_LONP1.png", 0.34, "A"), ("FigS1B_purity_episcore.png", 0.33, "B"), ("FigS1C_CNA_purity_adjusted.png", 0.33, "C")],
 ],
 "FigureS2": [   # 补充图S2：GSE44001 外部生存队列（scripts/20）
   [("FigS2_GSE44001_KM.png", 0.55, "")],
 ],
 "FigureS3": [   # 补充图S3：比例风险假设诊断与时变系数模型（scripts/27）
   [("FigS3A_schoenfeld.png", 1.0, "A")],
   [("FigS3B_HR_over_time.png", 1.0, "B")],
 ],
}

which = sys.argv[2] if len(sys.argv)>2 else "all"
for name, rows in LAYOUTS.items():
    if which != "all" and which != name: continue
    rendered = []
    for row in rows:
        usable = W - 2*MARG - GAP*(len(row)-1)
        imgs = []
        for fn, wt, lab in row:
            im = load(fn)
            tw = int(usable*wt)
            im = im.resize((tw, int(im.height*tw/im.width)), Image.LANCZOS)
            imgs.append((im, lab))
        rh = max(im.height for im,_ in imgs)
        rendered.append((imgs, rh))
    H = 2*MARG + sum(rh for _,rh in rendered) + GAP*(len(rendered)-1)
    canvas = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(canvas)
    y = MARG
    for imgs, rh in rendered:
        x = MARG
        for im, lab in imgs:
            canvas.paste(im, (x, y + (rh-im.height)//2))
            if lab: d.text((x+2, y-6), lab, font=FB, fill="black")
            x += im.width + GAP
        y += rh + GAP
    canvas.save(os.path.join(out, f"{name}.png"), dpi=(300,300))
    canvas.save(os.path.join(out, f"{name}.pdf"), "PDF", resolution=300)
    print(name, canvas.size)
print("ASSEMBLE_DONE")
