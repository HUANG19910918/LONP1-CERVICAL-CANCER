# patch_v5b.py —— 投稿前合规修正（在 patch_v5 输出上再处理，输出仍为 build_manuscript_v5.js）
#  (1) 摘要压缩至 BMC Cancer 上限(350词)以内  (2) 图按正文首次出现顺序重新编号：旧Fig6→2, 2→3, 3→4, 4→5, 5→6
#  (3) 面板首次出现顺序修正（Fig1E 提前提及；单细胞图 A,C / B,D；共表达图段落重排 + 面板字母 C/D/E 调整）
#  (4) 参考文献按首次引用顺序重新编号并重排列表  (5) ORCID 占位   依据见 实验记录.md R17
import sys, re
src, dst = sys.argv[1], sys.argv[2]
lines = open(src, encoding="utf-8").read().split("\n")
def js(s): return s.replace("\\", "\\\\").replace('"', '\\"')
def find(prefix):
    hits = [i for i, l in enumerate(lines) if l.strip().startswith(prefix)]
    assert len(hits) == 1, (prefix, hits); return hits[0]
def replace_line(prefix, new_text, kind="P"):
    i = find(prefix); indent = lines[i][:len(lines[i]) - len(lines[i].lstrip())]
    lines[i] = f'{indent}{kind}("{js(new_text)}"),'
def sub_in_line(prefix, old, new):
    i = find(prefix); assert old in lines[i], (prefix, old); lines[i] = lines[i].replace(old, new, 1)

# ---------- (1) 摘要 ----------
replace_line('P("Background: Lon peptidase 1',
    "Background: Lon peptidase 1 (LONP1) is a mitochondrial matrix protease that maintains mitochondrial proteostasis and is "
    "overexpressed in several cancers, but its expression pattern, genomic determinants and transcriptional context in cervical "
    "cancer have not been characterized.")
replace_line('P("Methods: We integrated RNA-seq data',
    "Methods: We integrated TCGA-CESC and GTEx RNA-seq data (304 tumors, 13 normal cervical tissues), three microarray cohorts "
    "(GSE9750, GSE7803, GSE63514), single-cell RNA-seq of the cervical carcinogenesis spectrum (GSE208653; 80,435 cells) and "
    "TCGA copy-number, mutation and clinical data. Genome-wide co-expression was annotated by genomic position to separate "
    "positional (19p13.3-linked) from functional co-expression before enrichment analysis; immune-signature correlations and "
    "survival analyses were also performed.")
replace_line('P("Results: LONP1 mRNA was higher',
    "Results: LONP1 mRNA was higher in tumors than in normal cervix in bulk-tissue cohorts (TCGA+GTEx log2 fold change 0.47, "
    "p=7.2x10-4; GSE9750 p=2.3x10-5) but not in laser-capture-microdissected epithelium (GSE63514, p=0.98). Single-cell "
    "analysis resolved this discordance: LONP1 was epithelial-enriched, yet "
    "malignant epithelial cells did not express more LONP1 than normal epithelial cells, indicating that the bulk elevation "
    "largely reflects tissue composition. Within tumors, expression tracked copy number (rho=0.41, p=3.1x10-13) and mutations "
    "were rare (1/294). Genome-wide co-expression was dominated by a positional signal: 18 of the 20 strongest partners, "
    "including CLPP (rho=0.72), map to 19p13.3, and 19p13.3 genes were 82-fold over-represented among genes with rho>=0.4 "
    "(p=1.1x10-81). After excluding 19p13.3 genes, partners remained enriched for mitochondrial inner membrane, mitochondrial "
    "translation and aerobic respiration, and LONP1 correlated with chaperones and stress-response factors elsewhere in the "
    "genome (TRAP1 rho=0.41, PHB1 0.40, HSPD1 0.33, ATF4 0.32). Within epithelial cells, LONP1 partners were enriched for "
    "respiratory chain complex assembly without positional bias (CLPP, HSPA9 and TRAP1 in the top quintile); fission and fusion genes correlated alike with LONP1, arguing against a fission-specific program. LONP1 was "
    "weakly inversely correlated with M1 macrophage and dendritic cell signatures and was not prognostic (all p>0.4).")
replace_line('P("Conclusions: LONP1 is an epithelial-enriched gene',
    "Conclusions: LONP1 is an epithelial-enriched gene whose apparent overexpression in bulk cervical cancer tissue is largely "
    "a tissue-composition effect. Its intratumoral variation is driven by 19p13.3 copy-number gain, which also dominates its "
    "strongest co-expression partners; beyond this positional effect, LONP1 marks a mitochondrial proteostasis and translation "
    "program without prognostic value.")

# ---------- (3a) 面板顺序：Fig1E 提前提及；单细胞图 ----------
sub_in_line('P("The overall design - four bulk-tissue cohorts',
    "(SCC versus normal p=0.052; Fig. 1D).",
    "(SCC versus normal p=0.052; Fig. 1D); the fourth cohort, microdissected epithelium, is considered separately below (Fig. 1E).")
sub_in_line('P("We analyzed 80,435 single cells', "identified nine major cell types (Fig. 6A), whose canonical marker profiles confirmed the annotation (Fig. 6C),",
            "identified nine major cell types whose canonical marker profiles confirmed the annotation (Fig. 6A, C),")
sub_in_line('P("We analyzed 80,435 single cells', "with intermediate levels in endothelium and fibroblasts (Fig. 6B-D).", "with intermediate levels in endothelium and fibroblasts (Fig. 6B, D).")

# ---------- (3b) 共表达图：段落重排（curated genes 先于富集；负相关富集移至上皮段末）；字母 C=lollipop, D=enrichment excl, E=dumbbell, F=epi enrichment, G=negative ----------
i = find('P("After excluding 19p13.3 genes, the remaining positively'); old = lines[i]
m = re.search(r'P\("(.*)"\),\s*$', old); txt = m.group(1)
a = txt.index("After excluding 19p13.3 genes"); b = txt.index("Among curated mitochondrial genes"); c = txt.index("Negatively correlated genes")
enrich = txt[a:b].strip(); curated = txt[b:c].strip(); neg = txt[c:].strip()
curated = curated.replace("Among curated mitochondrial genes located outside 19p13.3 (Fig. 3C; Additional file 1: Table S5), LONP1 correlated most strongly",
                          "Among curated mitochondrial genes (Fig. 3C; Additional file 1: Table S5), those located outside 19p13.3 that correlated most strongly with LONP1 were")
enrich = enrich.replace("(3.3x10-8; Fig. 3E; Additional file 1: Table S2)", "(3.3x10-8; Fig. 3D; Additional file 1: Table S2)")
assert "Fig. 3D" in enrich and "Fig. 3C" in curated
lines[i] = old[:old.index('P("')] + f'P("{js(curated + " " + enrich)}"),'
j = find('P("Because bulk co-expression can also be confounded by tissue composition')
lines[j] = lines[j].replace("pseudobulk units; Fig. 3D, F; Additional file 1: Tables S3-S5).", "pseudobulk units; Fig. 3E, F; Additional file 1: Tables S3-S5).")
assert "Fig. 3E, F" in lines[j]
neg = neg.replace("(Fig. 3G).", "(Fig. 3G).")
lines[j] = lines[j].rstrip()[:-3] + " " + js(neg) + '"),'  # 追加到该段末尾
# 图注 3 重写字母
replace_line('P("Figure 3. The LONP1 co-expression landscape',
    "Figure 3. The LONP1 co-expression landscape in TCGA-CESC and within cervical epithelium. (A) Genome-wide Spearman correlation "
    "of all annotated genes (n=17,378) with LONP1, ranked from highest to lowest; genes on 19p13.3 in red; dashed lines, |rho|=0.3. "
    "(B) Correlation versus genomic position along chromosome 19 (GRCh38); shaded area, 19p13.3; dashed line, LONP1; blue line, "
    "1-Mb sliding median. (C) Curated mitochondrial proteostasis, biogenesis and dynamics genes (server-side Spearman, n=294; "
    "*p<0.05, **p<0.01, ***p<0.001; genes on 19p13.3 labeled). (D) Over-representation of positively co-expressed bulk genes "
    "(rho>=0.4) after excluding 19p13.3 genes (n=102; Enrichr, Benjamini-Hochberg adjusted; dashed line, adjusted p=0.05; the KEGG "
    "oxidative phosphorylation term is shown for comparison and is not significant). (E) Bulk (TCGA-CESC, n=294) versus "
    "epithelium-intrinsic (GSE208653, 58 sample-by-subcluster pseudobulk units) correlations for the curated genes; shaded band "
    "and dotted line, interquartile range and median of the epithelium-intrinsic correlations of all 14,140 detected genes (the "
    "corresponding bulk interquartile range is -0.13 to 0.08). (F) Over-representation of the 300 genes most positively correlated "
    "with LONP1 within epithelium (no KEGG term reached significance). (G) Enrichment of negatively correlated bulk genes "
    "(rho<=-0.35, n=96).")

# ---------- (2) 图重新编号 ----------
FIGMAP = {"6": "2", "2": "3", "3": "4", "4": "5", "5": "6", "1": "1"}
ref_start = find("const refs = ["); ref_end = next(k for k in range(ref_start, len(lines)) if lines[k].strip() == "];")
def is_body(k):
    s = lines[k].strip()
    return not (ref_start <= k <= ref_end) and (s.startswith("P(") or s.startswith("H2("))
for k in range(len(lines)):
    if not is_body(k): continue
    lines[k] = re.sub(r"Fig\. (\d)", lambda m: "Fig. §" + FIGMAP[m.group(1)], lines[k])
    lines[k] = re.sub(r"Figure (\d)\.", lambda m: "Figure §" + FIGMAP[m.group(1)] + ".", lines[k])
    lines[k] = lines[k].replace("§", "")
# 图注段落按新编号排序
leg_idx = [k for k in range(len(lines)) if lines[k].strip().startswith('P("Figure ')]
legs = sorted([lines[k] for k in leg_idx], key=lambda l: int(re.search(r'P\("Figure (\d)', l).group(1)))
for k, l in zip(leg_idx, legs): lines[k] = l

# ---------- (4) 参考文献按首次引用重新编号 ----------
body_text = "\n".join(lines[k] for k in range(len(lines)) if is_body(k))
order = []
for grp in re.findall(r"\[([\d,\s\-]+)\]", body_text):
    for part in grp.split(","):
        part = part.strip()
        if "-" in part:
            x, y = part.split("-"); nums = range(int(x), int(y) + 1)
        elif part.isdigit(): nums = [int(part)]
        else: continue
        for n in nums:
            if n not in order: order.append(n)
refs = [lines[k] for k in range(ref_start + 1, ref_end)]
assert len(order) == len(refs), (len(order), len(refs))
newnum = {old: new for new, old in enumerate(order, 1)}
def remap(m):
    nums = []
    for part in m.group(1).split(","):
        part = part.strip()
        if "-" in part:
            x, y = part.split("-"); nums += list(range(int(x), int(y) + 1))
        elif part.isdigit(): nums.append(int(part))
    nums = sorted(newnum[n] for n in nums)
    out, s = [], 0
    while s < len(nums):
        e = s
        while e + 1 < len(nums) and nums[e + 1] == nums[e] + 1: e += 1
        out.append(f"{nums[s]}-{nums[e]}" if e - s >= 2 else (f"{nums[s]}, {nums[e]}" if e > s else f"{nums[s]}"))
        s = e + 1
    return "[" + ", ".join(out) + "]"
for k in range(len(lines)):
    if is_body(k): lines[k] = re.sub(r"\[([\d,\s\-]+)\]", remap, lines[k])
new_refs = [refs[old - 1] for old in order]
for idx, l in enumerate(new_refs):
    if not l.rstrip().endswith(","): new_refs[idx] = l.rstrip() + ","
lines[ref_start + 1:ref_end] = new_refs

# ---------- (5) ORCID ----------
i = find('P("* Correspondence: [email]'); indent = lines[i][:len(lines[i]) - len(lines[i].lstrip())]
lines.insert(i + 1, f'{indent}P("ORCID: [corresponding author ORCID iD]"),')

open(dst, "w", encoding="utf-8").write("\n".join(lines)); print("PATCHED(v5b) ->", dst, "| ref order:", order)
