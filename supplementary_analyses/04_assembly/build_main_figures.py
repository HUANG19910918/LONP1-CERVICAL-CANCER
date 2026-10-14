# build_main_figures.py —— 组装投稿用新增图与补充材料（只生成新文件，不改动【5】最终版本）
# 1 Figure 2 新增面板 H–K 拼图（H 重绘：修正 Adjacent/HSIL 组标签重叠）
# 2 Figure 4 新增面板 H–J 拼图
# 3 Additional file 2：原 Figure S1–S4 + 新 Figure S5–S8（版式同原 PDF：A4，Times Bold 12 标题，图宽 415 pt）
# 4 Additional file 1：原 Table S1–S9 + 新 Table S10–S12
# 输出均为矢量 PDF + 600 dpi PNG/TIFF（拼图）
# 用法: python -P build_main_figures.py <supplementary-analyses 目录> <【5】最终版本 目录>
import sys, os, io, shutil
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import fitz, openpyxl, docx
from openpyxl.styles import Font
from PIL import Image

sa, final = sys.argv[1:3]
out = os.path.join(sa, "main_figures"); os.makedirs(out, exist_ok=True)
sc = os.path.join(sa, "01_scRNA_inferCNV"); pc = os.path.join(sa, "02_partial_correlation"); mr = os.path.join(sa, "03_MR_coloc")
pub = os.path.join(sc, "figures", "publication")

# ---------- 1a 重绘 H（同 05_figures.py 面板 A，只改组标签位置） ----------
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "font.size": 7, "axes.linewidth": 0.7, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.labelsize": 8, "legend.fontsize": 6.5, "legend.frameon": False,
                     "xtick.major.width": 0.7, "ytick.major.width": 0.7})
NPG = {"GSE208653": "#E64B35", "E-MTAB-11948": "#3C5488", "GSE197461": "#00A087", "GSE236738": "#F39B7F", "GSE333892": "#8491B4"}
CLS = {"malignant": "#DC0000", "unresolved": "#D9D9D9", "nonmalignant": "#4DBBD5"}
short = lambda s: s.split("_", 1)[1] if s.startswith("GSM") else s.replace("_patient", " P")
S = pd.read_csv(os.path.join(sc, "results", "infercnv_summary.csv"))
s95 = S[S.rule == "q95"].copy(); s99 = S[S.rule == "q99"].set_index(["dataset", "sample"])
s95["og"] = s95.group.map({"Normal": 0, "AdjacentNormal": 1, "HSIL": 2, "SIL": 2, "Cancer": 3})
s95 = s95.sort_values(["og", "dataset", "sample"]).reset_index(drop=True)
f, ax = plt.subplots(figsize=(7.0, 2.5))
x = np.arange(len(s95)); bottom = np.zeros(len(s95))
for c in ["malignant", "unresolved", "nonmalignant"]:
    v = 100 * s95[c] / s95.n_epi; ax.bar(x, v, bottom=bottom, color=CLS[c], width=0.8, lw=0, label=c.capitalize()); bottom += v
ax.scatter(x, [s99.loc[(d, s), "pct_malignant"] for d, s in zip(s95.dataset, s95["sample"])], s=6, c="black", zorder=3,
           label="Malignant (q99 rule)")
ax.set_xticks(x); ax.set_xticklabels([short(s) for s in s95["sample"]], rotation=90, fontsize=5.5)
for t, d in zip(ax.get_xticklabels(), s95.dataset): t.set_color(NPG[d])
ax.set_ylabel("Epithelial cells (%)"); ax.set_ylim(0, 100); ax.set_xlim(-0.6, len(s95) - 0.4)
gpos = s95.groupby("group", sort=False).apply(lambda g: (g.index.min(), g.index.max()))
for g, (a, b) in gpos.items():
    lab = {"AdjacentNormal": "Adjacent"}.get(g, g); small = b - a < 3
    ax.annotate(lab, ((a + b) / 2, 111 if g == "AdjacentNormal" else (102 if small else 103)), ha="center", va="bottom",
                fontsize=5.5 if small else 7, annotation_clip=False)
    if a > 0: ax.axvline(a - 0.5, color="black", lw=0.5)
ax.legend(ncol=4, loc="upper left", bbox_to_anchor=(0, 1.36))
for i, d in enumerate(NPG):
    if d in set(s95.dataset): ax.text(1.005, 0.9 - i * 0.1, d, color=NPG[d], transform=ax.transAxes, fontsize=6)
f.tight_layout(); H_pdf = os.path.join(out, "panel_Fig2H_infercnv_classification.pdf"); f.savefig(H_pdf); plt.close(f)

# ---------- 拼图工具 ----------
def place(page, src, rect, clip=None, label=None):
    s = fitz.open(src); r = s[0].rect if clip is None else clip(s[0].rect)
    sc_ = min(rect.width / r.width, rect.height / r.height)
    tgt = fitz.Rect(rect.x0, rect.y0, rect.x0 + r.width * sc_, rect.y0 + r.height * sc_)
    page.show_pdf_page(tgt, s, 0, clip=r)
    if label: page.insert_text((rect.x0 - 2, rect.y0 + 9), label, fontname="hebo", fontsize=11)
    return tgt
cut = lambda f0, f1: (lambda r: fitz.Rect(r.x0 + r.width * f0, r.y0, r.x0 + r.width * f1, r.y1))
left, right = cut(0, 0.5), cut(0.5, 1)
def export(doc, name):
    p = os.path.join(out, name + ".pdf"); doc.save(p)
    pix = doc[0].get_pixmap(dpi=600); png = os.path.join(out, name + ".png"); pix.save(png)
    Image.open(png).convert("RGB").save(os.path.join(out, name + ".tif"), compression="tiff_lzw", dpi=(600, 600))
    print("saved", name)

# ---------- 1 Figure 2 H–K（宽 504 pt = 7 in） ----------
F1 = os.path.join(pub, "Fig_S1F_LONP1_progression_pooled.pdf"); F_D = os.path.join(pub, "Fig_S1D_LONP1_19p13.3_gain.pdf")
d = fitz.open(); pg = d.new_page(width=520, height=390)
place(pg, H_pdf, fitz.Rect(12, 4, 516, 184), label="H")
w = 168
place(pg, F1, fitz.Rect(12, 196, 12 + w, 386), clip=cut(0.524, 1), label="I")
place(pg, F_D, fitz.Rect(12 + w + 4, 196, 12 + 2 * w + 4, 386), clip=left, label="J")
place(pg, F1, fitz.Rect(12 + 2 * w + 8, 196, 516, 386), clip=cut(0, 0.523), label="K")
export(d, "Figure2_new_panels_H-K")

# ---------- 2 Figure 4 H–J ----------
d = fitz.open(); pg = d.new_page(width=520, height=470)
place(pg, os.path.join(pc, "figures", "Fig_P2a_19p13.3_enrichment_OR.pdf"), fitz.Rect(12, 4, 270, 210), label="H")
place(pg, os.path.join(pc, "figures", "Fig_P2b_raw_vs_partial_rho.pdf"), fitz.Rect(282, 4, 516, 220), label="I")
place(pg, os.path.join(pc, "figures", "Fig_P2c_fgsea_mito_pathways.pdf"), fitz.Rect(12, 226, 516, 466), label="J")
export(d, "Figure4_new_panels_H-J")

# ---------- 3 Additional file 2 ----------
add2 = fitz.open(os.path.join(final, "Supplementary Figures.pdf"))
def sfig(title, items):  # items: [(pdf, label or None)], 纵向排列，总宽 415 pt
    pg = add2.new_page(width=595, height=842)
    pg.insert_text((90, 83), title, fontname="tibo", fontsize=12)
    y = 92
    for src, lab in items:
        r = fitz.open(src)[0].rect; h = min(415 * r.height / r.width, 742 - y)
        t = place(pg, src, fitz.Rect(90 + (12 if lab else 0), y, 505, y + h), label=lab); y = t.y1 + 8
sfig("Figure S5", [(os.path.join(pub, "Fig_S1E_copykat_crossvalidation.pdf"), None)])
sfig("Figure S6", [(os.path.join(pub, "Fig_S1G_SIL_grade_sensitivity.pdf"), None)])
# S7：A 竖长，放左；B、C 叠放右
pg = add2.new_page(width=595, height=842); pg.insert_text((90, 83), "Figure S7", fontname="tibo", fontsize=12)
place(pg, os.path.join(mr, "figures", "Fig_MR_A_LONP1_forest.pdf"), fitz.Rect(102, 92, 338, 760), label="A")
place(pg, os.path.join(mr, "figures", "Fig_MR_B_LONP1_regional.pdf"), fitz.Rect(352, 92, 520, 254), label="B")
place(pg, os.path.join(mr, "figures", "Fig_MR_C_coloc_heatmap.pdf"), fitz.Rect(352, 262, 520, 560), label="C")
sfig("Figure S8", [(os.path.join(mr, "figures", "Fig_MR_D_allgenes_forest.pdf"), "A"),
                   (os.path.join(mr, "figures", "Fig_MR_E_allgenes_coloc.pdf"), "B")])
add2.save(os.path.join(out, "Additional file 2_Supplementary Figures S1-S8.pdf")); print("saved Additional file 2 (figures only)")
# Supplementary Methods（A4，Times 10 pt，版心同图页 90–505 pt）置于 S1 之前
css = ("* {font-family: serif;} h1 {font-size: 12pt; font-weight: bold; margin: 0 0 2pt 0;} "
       ".sub {font-size: 9pt; font-style: italic; margin: 0 0 6pt 0;} h2 {font-size: 10pt; font-weight: bold; margin: 7pt 0 2pt 0;} "
       "p {font-size: 9.5pt; text-align: justify; margin: 0 0 4pt 0; line-height: 1.25;}")
story = fitz.Story(open(os.path.join(out, "supplementary_methods.html"), encoding="utf-8").read(), user_css=css)
buf = io.BytesIO(); wr = fitz.DocumentWriter(buf); more = 1
while more:
    dev = wr.begin_page(fitz.paper_rect("a4")); more, _ = story.place(fitz.Rect(90, 71, 505, 780)); story.draw(dev); wr.end_page()
wr.close()
meth = fitz.open("pdf", buf.getvalue()); meth.insert_pdf(add2)
meth.save(os.path.join(out, "Additional file 2_Supplementary Methods and Figures S1-S8.pdf"))
print("saved Additional file 2 with methods: pages", meth.page_count, "(methods", meth.page_count - add2.page_count, ")")

# ---------- 4 Additional file 1 ----------
dst = os.path.join(out, "Additional file 1_Supplementary Tables S1-S12.xlsx")
shutil.copy(os.path.join(final, "Supplementary Tables.xlsx"), dst)
wb = openpyxl.load_workbook(dst)
def sheet(name, title, header, rows, note=None):
    ws = wb.create_sheet(name); ws.append([title]); ws["A1"].font = Font(bold=True)
    ws.append(header)
    for c in ws[2]: c.font = Font(bold=True)
    for r in rows: ws.append(list(r))
    if note: ws.append([]); ws.append([note])
sheet("Table S10", "Table S10. Single-cell datasets used in the supplementary analyses",
      ["Dataset", "Samples used", "Platform", "Cells after QC", "Use"],
      [["GSE208653", "4 normal, 2 HSIL, 3 cancer", "10x", 74892, "Normal reference; HSIL; tumors"],
       ["E-MTAB-11948", "3 cancer + 3 adjacent normal", "10x 3' v3", 64592, "Normal reference (2 adjacent); tumors"],
       ["GSE197461", "8 cancer (3 SCC, 5 ADC; 6 HPV+)", "10x 5'", 29217, "Tumors"],
       ["GSE236738", "3 treatment-naive cancer", "10x", 18977, "Tumors"],
       ["GSE173231", "10 normal biopsies", "10x 5'", 20771, "Not used as epithelial reference (723 epithelial cells)"],
       ["GSE333892", "23 cancer, 10 normal, 19 SIL", "10x / BD (per README)", 352881, "Normal-SIL-cancer comparison; tumors; pooled normal"],
       ["Excluded", "GSE168652 (duplicate of E-MTAB-11948 P3); GSE171894 (normalized only)", "", "", ""]],
      "Cells after QC are summed over the samples used. GSE333892 counts were reconstructed from the deposited normalized object (see Supplementary Methods).")
t = docx.Document(os.path.join(sa, "Table_S11_STROBE-MR_checklist_20261012.docx")).tables[0]
R = [[c.text for c in r.cells] for r in t.rows]
sheet("Table S11", "Table S11. STROBE-MR checklist", R[0], R[1:])
w2 = openpyxl.load_workbook(os.path.join(sa, "Table_S12_MR_coloc_full_results_20261012.xlsx"))
for s in w2.sheetnames:
    src = w2[s]; nm = "Table S12 README" if s == "README" else "Table " + s.replace("_", " ", 1)[:28]
    ws = wb.create_sheet(nm[:31])
    for row in src.iter_rows(values_only=True): ws.append(list(row))
if "README" in wb.sheetnames:
    rd = wb["README"]; rd.append([]); rd.append(["Tables S10-S12 were added in the revision (supplementary single-cell, partial-correlation and Mendelian randomization analyses)."])
wb.save(dst); print("saved Additional file 1:", wb.sheetnames)
