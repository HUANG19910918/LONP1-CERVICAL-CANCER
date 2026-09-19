# patch_v6.py —— 由 build_manuscript_v5.js 生成 build_manuscript_v6.js
# 新增：(1) TCGA 肿瘤纯度（ABSOLUTE）校正  (2) 外部生存队列 GSE44001  (3) Human Protein Atlas 蛋白层面证据
#       (4) Additional file 2（Fig. S1–S2）  (5) 新参考文献并按首次引用重新编号   依据见 实验记录.md R18
import sys, re
src, dst = sys.argv[1], sys.argv[2]
lines = open(src, encoding="utf-8").read().split("\n")
def js(s): return s.replace("\\", "\\\\").replace('"', '\\"')
def find(prefix):
    hits = [i for i, l in enumerate(lines) if l.strip().startswith(prefix)]
    assert len(hits) == 1, (prefix, hits); return hits[0]
def sub_in_line(prefix, old, new):
    i = find(prefix); assert old in lines[i], (prefix, old[:60]); lines[i] = lines[i].replace(old, new, 1)
def insert_after(prefix, items):
    i = find(prefix); indent = lines[i][:len(lines[i]) - len(lines[i].lstrip())]
    if not lines[i].rstrip().endswith(","): lines[i] = lines[i].rstrip() + ","
    for k, (kind, txt) in enumerate(items): lines.insert(i + 1 + k, f'{indent}{kind}("{js(txt)}"),')

# ---- 新参考文献（临时编号 40-44，稍后按首次引用重排）----
ref_start = find("const refs = ["); ref_end = next(k for k in range(ref_start, len(lines)) if lines[k].strip() == "];")
new_refs = [
 "Carter SL, Cibulskis K, Helman E, McKenna A, Shen H, Zack T, et al. Absolute quantification of somatic DNA alterations in human cancer. Nat Biotechnol. 2012;30(5):413-21.",
 "Taylor AM, Shih J, Ha G, Gao GF, Zhang X, Berger AC, et al. Genomic and functional approaches to understanding cancer aneuploidy. Cancer Cell. 2018;33(4):676-689.e3.",
 "Lee YY, Kim TJ, Kim JY, Choi CH, Do IG, Song SY, et al. Genetic profiling to predict recurrence of early cervical cancer. Gynecol Oncol. 2013;131(3):650-4.",
 "Uhlen M, Fagerberg L, Hallstrom BM, Lindskog C, Oksvold P, Mardinoglu A, et al. Tissue-based map of the human proteome. Science. 2015;347(6220):1260419.",
 "Uhlen M, Zhang C, Lee S, Sjostedt E, Fagerberg L, Bidkhori G, et al. A pathology atlas of the human cancer transcriptome. Science. 2017;357(6352):eaan2507.",
]
if not lines[ref_end - 1].rstrip().endswith(","): lines[ref_end - 1] = lines[ref_end - 1].rstrip() + ","
for k, r in enumerate(new_refs): lines.insert(ref_end + k, f'"{js(r)}",')
ref_end += len(new_refs)
n_old = ref_end - ref_start - 1 - len(new_refs)  # 39
R = {"carter": n_old + 1, "taylor": n_old + 2, "lee": n_old + 3, "uhlen15": n_old + 4, "uhlen17": n_old + 5}

# ---- Abstract ----
sub_in_line('P("Results: LONP1 mRNA was higher', "and was not prognostic (all p>0.4).", "and was not prognostic in TCGA or an independent 300-patient cohort (all p>0.4).")

# ---- Methods: 新小节 ----
insert_after('P("Patients with follow-up time >0 were included', [
 ("H2", "Tumor purity, external survival cohort and protein-level data"),
 ("P", f"ABSOLUTE tumor purity estimates for TCGA-CESC were taken from the PanCanAtlas consensus calls [{R['carter']}, {R['taylor']}] (n=291 with expression data). "
       "Associations of LONP1 with purity, and of copy number and immune signatures with LONP1 after adjustment for purity (Spearman correlation of "
       "residuals from linear regression on purity; ordinary least squares with purity as covariate), were computed. As an external survival cohort, "
       f"GSE44001 [{R['lee']}] (300 patients with early-stage cervical cancer treated by radical hysterectomy at a single center; Illumina HumanHT-12 DASL, "
       "GPL14951; log2 quantile-normalized values as deposited; LONP1 probe ILMN_1766125, the only LONP1 probe on the platform) was analyzed for "
       "disease-free survival by univariate Cox regression, median-split Kaplan-Meier analysis and a multivariable Cox model including FIGO stage "
       f"(IA-IB1 versus IB2-IIA) and largest tumor diameter. Protein-level information was obtained from the Human Protein Atlas [{R['uhlen15']}, {R['uhlen17']}] "
       "(antibody HPA002192; immunohistochemistry of normal cervix and cervical cancer; accessed August 2026)."),
])

# ---- Results ----
sub_in_line('P("We analyzed 80,435 single cells',
    "which is precisely why epithelium-resolved data were required for a decisive answer.",
    "which is precisely why epithelium-resolved data were required for a decisive answer. Within tumors, LONP1 correlated only weakly with ABSOLUTE "
    "purity (rho=0.12, p=0.045, n=291; Additional file 2: Fig. S1A) and the epithelial marker score itself tracked purity weakly (rho=0.13, p=0.025; "
    "Fig. S1B), as expected given the narrow purity range of resected tumors (median 0.67, interquartile range 0.54-0.79); the composition effect is "
    "thus a tumor-versus-normal phenomenon that the within-tumor purity gradient captures only marginally. Protein-level data from the Human Protein "
    "Atlas were consistent with the absence of cell-intrinsic induction: antibody HPA002192 stained squamous and glandular epithelium of normal cervix "
    "at medium intensity and tumor cells at medium intensity in 8 of 10 cervical cancers (not detected in the remaining 2), with no indication of "
    f"stronger staining in tumor cells than in normal epithelium [{R['uhlen15']}, {R['uhlen17']}].")
sub_in_line('P("LONP1 expression correlated with its copy-number ratio',
    "Somatic LONP1 mutations were nearly absent",
    "The copy-number association was independent of tumor purity (purity-adjusted partial rho=0.45, p=1.2x10-15; linear model coefficient for the "
    "log2 copy-number ratio 0.97, p=1.1x10-19, with purity as covariate; Additional file 2: Fig. S1C). Somatic LONP1 mutations were nearly absent")
sub_in_line('P("LONP1 expression was not associated with overall survival',
    "We report this negative result in full; no alternative cut-points were explored.",
    "We report this negative result in full; no alternative cut-points were explored. The absence of prognostic association was replicated in an "
    "independent cohort of 300 early-stage cervical cancers (GSE44001; 38 recurrences): disease-free survival was unrelated to LONP1 as a continuous "
    "variable (HR=1.07 per log2 unit, 95% CI 0.72-1.60, p=0.73) or by median split (log-rank p=0.83), whereas tumor diameter carried the expected "
    "prognostic weight in the multivariable model (HR=1.29 per cm, 1.10-1.52, p=0.002; Additional file 2: Fig. S2).")
sub_in_line('P("Among thirteen immune marker signatures',
    "arguing against a strong direct link between LONP1 and checkpoint regulation in bulk tissue.",
    "arguing against a strong direct link between LONP1 and checkpoint regulation in bulk tissue. These associations were not explained by tumor "
    "purity: after adjustment for ABSOLUTE purity the partial correlations with the M1 macrophage (rho=-0.22, p=1.6x10-4), dendritic cell (-0.19, "
    "p=1.0x10-3) and CD4/Th (-0.13, p=0.026) signatures persisted, even though the dendritic cell and CD4/Th signatures themselves correlated strongly "
    "with purity (rho=-0.54 and -0.58).")

# ---- Discussion ----
sub_in_line('P("Equally important are the results that did not support',
    "whether analyzed continuously or by median split. We chose not to search",
    "whether analyzed continuously or by median split, nor for disease-free survival in an independent cohort of 300 early-stage patients. We chose not to search")
sub_in_line('P("Equally important are the results that did not support',
    "Whether LONP1 protein, activity, or stress-induced regulation differs between malignant and normal epithelium remains open, and protein-level validation in situ (immunohistochemistry comparing tumor cells with normal epithelium) is the direct next step.",
    "The available protein-level evidence points the same way: Human Protein Atlas immunohistochemistry shows medium LONP1 staining in normal cervical "
    "epithelium and medium or absent staining in tumor cells of cervical cancers, although these data are semi-quantitative and based on few cases. "
    "Whether LONP1 protein turnover, activity or stress-induced regulation differs between malignant and normal epithelium remains open, and quantitative "
    "in situ comparison of tumor cells with adjacent normal epithelium is the direct next step.")
sub_in_line('P("The weak inverse association between LONP1 and M1 macrophage',
    "but the effect sizes here are small (|rho|~0.22) and marker-signature approaches in bulk RNA-seq cannot resolve cell-intrinsic from compositional effects.",
    "but the effect sizes here are small (|rho|~0.22) and, although they persisted after adjustment for tumor purity, marker-signature approaches in "
    "bulk RNA-seq cannot resolve cell-intrinsic from compositional effects.")
sub_in_line('P("Limitations. This is an in silico study',
    "The immune analysis used compact marker signatures rather than full deconvolution.",
    "The immune analysis used compact marker signatures rather than full deconvolution, and purity adjustment relies on DNA-based ABSOLUTE estimates "
    "that capture only the within-tumor range of tissue composition.")

# ---- Availability / Additional files ----
sub_in_line('P("All data analyzed are publicly available',
    "NCBI GEO (accessions GSE9750, GSE7803, GSE63514, GSE208653).",
    "NCBI GEO (accessions GSE9750, GSE7803, GSE63514, GSE208653, GSE44001), the GDC PanCanAtlas publication page (ABSOLUTE purity calls) and the Human Protein Atlas (https://www.proteinatlas.org).")
sub_in_line('P("All data analyzed are publicly available', "Supplementary Tables S1-S7 are provided as Additional file 1.",
            "Supplementary Tables S1-S7 are provided as Additional file 1 and Supplementary Figures S1-S2 as Additional file 2.")
insert_after('P("Additional file 1: Supplementary Tables S1-S7', [
 ("P", "Additional file 2: Supplementary Figures S1-S2 (PDF). Figure S1, tumor purity (ABSOLUTE) in TCGA-CESC: (A) LONP1 expression versus purity; "
       "(B) epithelial marker score versus purity; (C) copy number versus expression after adjustment for purity. Figure S2, disease-free survival by "
       "LONP1 (median split) in the independent GSE44001 cohort of 300 early-stage cervical cancers."),
])

# ---- 参考文献按首次引用重新编号 ----
def is_body(k):
    s = lines[k].strip(); return not (ref_start <= k <= ref_end) and (s.startswith("P(") or s.startswith("H2("))
ref_end = next(k for k in range(ref_start, len(lines)) if lines[k].strip() == "];")
body_text = "\n".join(lines[k] for k in range(len(lines)) if is_body(k))
order = []
for grp in re.findall(r"\[([\d,\s\-]+)\]", body_text):
    for part in grp.split(","):
        part = part.strip()
        if "-" in part: x, y = part.split("-"); nums = range(int(x), int(y) + 1)
        elif part.isdigit(): nums = [int(part)]
        else: continue
        for n in nums:
            if n not in order: order.append(n)
refs = [lines[k] for k in range(ref_start + 1, ref_end)]
assert len(order) == len(refs), (len(order), len(refs), sorted(set(range(1, len(refs) + 1)) - set(order)))
newnum = {old: new for new, old in enumerate(order, 1)}
def remap(m):
    nums = []
    for part in m.group(1).split(","):
        part = part.strip()
        if "-" in part: x, y = part.split("-"); nums += list(range(int(x), int(y) + 1))
        elif part.isdigit(): nums.append(int(part))
    nums = sorted(newnum[n] for n in nums); out, s = [], 0
    while s < len(nums):
        e = s
        while e + 1 < len(nums) and nums[e + 1] == nums[e] + 1: e += 1
        out.append(f"{nums[s]}-{nums[e]}" if e - s >= 2 else (f"{nums[s]}, {nums[e]}" if e > s else f"{nums[s]}")); s = e + 1
    return "[" + ", ".join(out) + "]"
for k in range(len(lines)):
    if is_body(k): lines[k] = re.sub(r"\[([\d,\s\-]+)\]", remap, lines[k])
new_list = [refs[old - 1] for old in order]
new_list = [l if l.rstrip().endswith(",") else l.rstrip() + "," for l in new_list]
lines[ref_start + 1:ref_end] = new_list
open(dst, "w", encoding="utf-8").write("\n".join(lines)); print("PATCHED(v6) ->", dst, "| refs:", len(new_list))
