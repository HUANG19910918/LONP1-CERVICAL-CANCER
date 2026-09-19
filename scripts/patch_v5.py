# patch_v5.py —— 由 build_manuscript_v4.js 生成 build_manuscript_v5.js
# 改动：(1) 语言润色（表面层：英式/美式拼写统一为美式；不改句式与措辞强度）
#       (2) Fig2B Kruskal-Wallis 精确 p (1.1x10-8)  (3) 一处事实性措辞修正（THOP1/TIMM13 非"仅次于CLPP的两个"）
#       (4) 补充材料 Additional file 1 (Table S1-S7) 正文引用 + Additional files 节
#       (5) Data availability 改为仓库/DOI 占位  依据见 实验记录.md R17
import sys, re
src, dst = sys.argv[1], sys.argv[2]
lines = open(src, encoding="utf-8").read().split("\n")
def js(s): return s.replace("\\", "\\\\").replace('"', '\\"')
def find(prefix):
    hits = [i for i, l in enumerate(lines) if l.strip().startswith(prefix)]
    assert len(hits) == 1, (prefix, hits); return hits[0]
def sub_in_line(prefix, old, new):
    i = find(prefix); assert old in lines[i], (prefix, old); lines[i] = lines[i].replace(old, new, 1)
def insert_after(prefix, new_lines):
    i = find(prefix); indent = lines[i][:len(lines[i]) - len(lines[i].lstrip())]
    for k, (kind, txt) in enumerate(new_lines): lines.insert(i + 1 + k, f'{indent}{kind}("{js(txt)}"),')

# ---------- (1) 拼写统一：仅作用于正文 P()/H1()/H2()/标题行，不动参考文献数组 ----------
ref_start = find("const refs = ["); ref_end = next(i for i in range(ref_start, len(lines)) if lines[i].strip() == "];")
SPELL = [(r"\bprogrammes\b", "programs"), (r"\bprogramme\b", "program"), (r"\btumour\b", "tumor"), (r"\btumours\b", "tumors"),
         (r"\bcentred\b", "centered"), (r"\bcentre\b", "center"), (r"\bneighbours\b", "neighbors"), (r"\bneighbouring\b", "neighboring"),
         (r"\bneighbourhood\b", "neighborhood"), (r"\bsignalling\b", "signaling"), (r"\bmodelled\b", "modeled"), (r"\blabelled\b", "labeled"),
         (r"\bfavouring\b", "favoring"), (r"\bcolour\b", "color"), (r"\bbehaviour\b", "behavior"), (r"\bcharacterised\b", "characterized"),
         (r"\banalysed\b", "analyzed"), (r"\bvisualised\b", "visualized"), (r"\bnormalised\b", "normalized")]
n_spell = 0
for i, l in enumerate(lines):
    if ref_start <= i <= ref_end: continue
    s = l.strip()
    if not (s.startswith("P(") or s.startswith("H1(") or s.startswith("H2(") or s.startswith("new Paragraph({ children: [new TextRun({ text:")): continue
    new = l
    for pat, rep in SPELL:
        new, k = re.subn(pat, rep, new); n_spell += k
    lines[i] = new
print("spelling replacements:", n_spell)

# ---------- (2) Fig2B 精确 p ----------
sub_in_line('P("LONP1 expression correlated with its copy-number ratio', "Kruskal-Wallis p<0.001; Fig. 2B)", "Kruskal-Wallis p=1.1x10-8; Fig. 2B)")
sub_in_line('P("Figure 2. Genomic and clinicopathological correlates', "(Kruskal-Wallis p<0.001)", "(Kruskal-Wallis p=1.1x10-8)")

# ---------- (3) 措辞修正 ----------
sub_in_line('P("Because bulk co-expression can also be confounded by tissue composition',
    "and the two strongest bulk partners after CLPP, THOP1 and TIMM13, fell to the 43rd and 39th percentiles",
    "and two of the strongest bulk partners, THOP1 and TIMM13, fell to the 43rd and 39th percentiles")
sub_in_line('P("Figure 1. Study design and LONP1 expression', "whiskers 1.5xIQR", "whiskers 1.5xIQR.")
sub_in_line('P("In contrast, in GSE63514, where expression was profiled', "Fig. 1B-C may reflect", "Fig. 1B, C may reflect")

# ---------- (4) 补充材料引用 ----------
sub_in_line('P("To place the cervical finding in a pan-cancer context', "across 32 cancer types (Fig. 1F).", "across 32 cancer types (Fig. 1F; Additional file 1: Table S7).")
sub_in_line('P("Genome-wide, 970 of 17,378 genes correlated with LONP1', "(628 positive, 342 negative; Fig. 3A).", "(628 positive, 342 negative; Fig. 3A; Additional file 1: Table S1).")
sub_in_line('P("After excluding 19p13.3 genes, the remaining positively', "in the KEGG ribosome pathway (3.3x10-8; Fig. 3E).", "in the KEGG ribosome pathway (3.3x10-8; Fig. 3E; Additional file 1: Table S2).")
sub_in_line('P("After excluding 19p13.3 genes, the remaining positively', "Among curated mitochondrial genes located outside 19p13.3 (Fig. 3C),", "Among curated mitochondrial genes located outside 19p13.3 (Fig. 3C; Additional file 1: Table S5),")
sub_in_line('P("Because bulk co-expression can also be confounded by tissue composition', "pseudobulk units; Fig. 3D, F).", "pseudobulk units; Fig. 3D, F; Additional file 1: Tables S3-S5).")
sub_in_line('P("Because bulk co-expression can also be confounded by tissue composition', "and protein transport, with no significant KEGG term;", "and protein transport, with no significant KEGG term (Additional file 1: Table S4);")
sub_in_line('P("Among thirteen immune marker signatures', "signatures showed no association (Fig. 5A-C).", "signatures showed no association (Fig. 5A-C; Additional file 1: Table S6).")
sub_in_line('P("Thirteen immune signatures were defined', "following marker conventions of published immunome compendia [28].", "following marker conventions of published immunome compendia [28] (Additional file 1: Table S6).")

# ---------- (5) Data availability + software statement ----------
sub_in_line('P("All data analyzed are publicly available',
    "Analysis scripts and extracted per-sample data are available from the corresponding author on reasonable request.",
    "All analysis scripts, the per-sample data extracted from these repositories, derived tables and a step-by-step verification report are deposited at GitHub ([GitHub URL]) and archived at Zenodo ([Zenodo DOI]). Supplementary Tables S1-S7 are provided as Additional file 1.")
sub_in_line('P("Statistical analyses were performed in R 4.3.3',
    "Analysis scripts, per-sample source data extracted from the public repositories, and a step-by-step verification report are retained by the authors and available on request to allow full reproduction of every reported value.",
    "Analysis scripts, per-sample source data extracted from the public repositories, and a step-by-step verification report are publicly deposited (see Availability of data and materials) to allow full reproduction of every reported value.")
# Additional files 节（Figure legends 之后）
i = find('P("Figure 6. Single-cell resolution of LONP1')
indent = lines[i][:len(lines[i]) - len(lines[i].lstrip())]
if not lines[i].rstrip().endswith(","): lines[i] = lines[i].rstrip() + ","
lines.insert(i + 1, f'{indent}H1("Additional files"),')
lines.insert(i + 2, f'{indent}P("{js("Additional file 1: Supplementary Tables S1-S7 (XLSX). Table S1, genome-wide LONP1 co-expression in TCGA-CESC with GRCh38 coordinates and cytobands; Table S2, Enrichr results for bulk co-expressed genes with and without 19p13.3 genes and for negatively correlated genes; Table S3, epithelium-intrinsic genome-wide co-expression (GSE208653 pseudobulk); Table S4, Enrichr results for epithelium-intrinsic LONP1 partners; Table S5, curated mitochondrial genes with bulk and epithelium-intrinsic correlations; Table S6, immune marker signature definitions and correlations; Table S7, pan-cancer tumor-versus-normal statistics.")}"),')

out = "\n".join(lines)
body = "\n".join(l for k, l in enumerate(lines) if not (ref_start <= k <= ref_end))
for w in ["programme", "tumour", "centred", "neighbour", "signalling", "labelled", "modelled", "p<0.001; Fig. 2B"]:
    assert w not in body, w
open(dst, "w", encoding="utf-8").write(out); print("PATCHED ->", dst)
