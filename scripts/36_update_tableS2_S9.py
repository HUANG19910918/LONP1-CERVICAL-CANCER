# 36_update_tableS2_S9.py —— 补全 Table S2 "all genes" 组 overlapping_gene_symbols；新增 Table S9（自定义背景敏感性）并写入 README
# 输入: raw_data/60_, 61_（由 35 脚本生成）；输出: 直接更新 【5】最终版本/Additional_file_1_Supplementary_Tables.xlsx
import sys, pandas as pd, openpyxl, copy
from openpyxl.styles import Font, PatternFill, Alignment
base = sys.argv[1]; xl = base + "/../【5】最终版本/Additional_file_1_Supplementary_Tables.xlsx"
d = pd.read_csv(base + "/raw_data/61_enrichment_bulk_pos_all_default_with_genes.csv").set_index(["library", "term"])
c = pd.read_csv(base + "/raw_data/60_enrichment_custom_background.csv")
wb = openpyxl.load_workbook(xl)
ws = wb["Table S2"]; n = 0
for row in ws.iter_rows(min_row=4):
    if row[0].value and "all genes" in str(row[0].value) and row[6].value in (None, ""):
        x = d.loc[(row[1].value, row[2].value)]
        assert x.n_overlap == row[5].value and abs(x.p / row[3].value - 1) < 1e-6
        row[6].value = x.genes; n += 1
print("filled", n)
ws["A1"].value = ws["A1"].value.replace(" Overlapping gene symbols were not exported for the all-genes set.", "")
ws["A1"].value += " Enrichr default background (all genes in each library). Sensitivity analyses with custom backgrounds are in Table S9."
# Table S9
if "Table S9" in wb.sheetnames: del wb["Table S9"]
w = wb.create_sheet("Table S9")
w["A1"].value = ("Table S9. Enrichr sensitivity analysis with custom backgrounds (Speedrichr), restricted to genes that could enter each ranking. "
                 "Bulk sets: background = 17,378 genes ranked in the TCGA-CESC genome-wide co-expression (Table S1); for the set excluding 19p13.3 genes, both foreground and background exclude 19p13.3 genes (17,190 genes). "
                 "Epithelium-intrinsic sets: background = 14,140 genes with detection >=1% in epithelial cells that entered the ranking (Table S3). Enrichr intersects the background with each library. "
                 "All returned terms are listed; adjusted p values are Benjamini-Hochberg values computed by Enrichr within each library.")
w["A1"].font = copy.copy(ws["A1"].font); w["A1"].alignment = copy.copy(ws["A1"].alignment)
hdr = ["gene_set", "background_genes", "library", "term", "p_value", "adjusted_p_BH", "overlap_genes", "overlapping_gene_symbols"]
bgn = {"bulk_pos_excl19p13.3_n102": "17,190 ranked genes outside 19p13.3", "epi_top300_raw": "14,140 epithelial ranked genes", "epi_top300_sample_centered": "14,140 epithelial ranked genes"}
for j, h in enumerate(hdr, 1):
    cell = w.cell(3, j, h); cell.font = copy.copy(ws["A3"].font); cell.fill = copy.copy(ws["A3"].fill); cell.alignment = copy.copy(ws["A3"].alignment)
for i, r in enumerate(c.itertuples(), 4):
    for j, v in enumerate([r.gene_set, bgn.get(r.gene_set, "17,378 TCGA-CESC ranked genes"), r.library, r.term, r.p, r.padj, r.n_overlap, r.genes], 1): w.cell(i, j, v)
for col, wd in zip("ABCDEFGH", [34, 32, 28, 60, 14, 14, 12, 70]): w.column_dimensions[col].width = wd
w.freeze_panes = "A4"
rd = wb["README"]
for row in rd.iter_rows():
    if row[0].value and str(row[0].value).startswith("Table S8"):
        r8 = row[0].row
rd.insert_rows(r8 + 1)
rd.cell(r8 + 1, 1, "Table S9  Enrichr sensitivity analysis with custom (rankable-gene) backgrounds for the bulk and epithelium-intrinsic gene sets (Fig. 4D, F, G)")
rd.cell(r8 + 1, 1)._style = copy.copy(rd.cell(r8, 1)._style)
wb.save(xl); print("saved")
