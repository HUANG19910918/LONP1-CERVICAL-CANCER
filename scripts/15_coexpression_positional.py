# 15_coexpression_positional.py —— LONP1全基因组共表达（threshold=0）+ 基因组位置注释 + 去19p13.3后富集
# 复现口径：
#   共表达: cBioPortal /molecular-profiles/co-expressions/fetch, profileA=B=cesc_tcga_pan_can_atlas_2018_rna_seq_v2_mrna,
#           threshold=0.0, entrezGeneId=9361 (LONP1), sampleListId=cesc_tcga_pan_can_atlas_2018_all (服务器端Spearman, n=294)
#   基因位置: Ensembl REST /lookup/symbol/homo_sapiens (GRCh38)；细胞带: UCSC hg38 cytoBand.txt.gz
#   富集: Enrichr (GO_Biological_Process_2023, GO_Cellular_Component_2023, KEGG_2021_Human)，输入=rho>=0.4且不在19p13.3的基因
# 输出: raw_data/30_coexpression_LONP1_CESC_genomewide.csv, 31_enrichment_enrichr_pos_excl19p13.3.csv,
#       32_positional_stats.txt
# 用法: python3 15_coexpression_positional.py <分析根目录>
import sys, os, json, gzip, time, io, urllib.request
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

base = sys.argv[1]; rd = os.path.join(base, "raw_data"); tmp = os.path.join(base, "raw_data", "_cache")
os.makedirs(tmp, exist_ok=True)
API = "https://www.cbioportal.org/api"; ENS = "https://rest.ensembl.org"
PROF = "cesc_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"; SL = "cesc_tcga_pan_can_atlas_2018_all"

def jpost(url, body, hdr=None):
    req = urllib.request.Request(url, data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json", "Accept": "application/json", **(hdr or {})})
    return json.load(urllib.request.urlopen(req, timeout=600))

# ---------- 1. 全基因组共表达 ----------
cache = os.path.join(tmp, "coexp_full_raw.json")
if os.path.exists(cache):
    co = json.load(open(cache))
else:
    co = jpost(f"{API}/molecular-profiles/co-expressions/fetch?molecularProfileIdA={PROF}&molecularProfileIdB={PROF}&threshold=0.0",
               {"entrezGeneId": 9361, "sampleListId": SL})
    json.dump(co, open(cache, "w"))
co = pd.DataFrame(co)
co = co[co.geneticEntityType == "GENE"].copy()
co["entrez"] = co.geneticEntityId.astype(int)
co = co.rename(columns={"spearmansCorrelation": "rho", "pValue": "p"})[["entrez", "rho", "p"]].dropna(subset=["rho"])
print("co-expression rows:", len(co))

# entrez -> symbol
sym = {}
ids = co.entrez.tolist()
for i in range(0, len(ids), 1000):
    g = jpost(f"{API}/genes/fetch?geneIdType=ENTREZ_GENE_ID", [str(x) for x in ids[i:i+1000]])
    for x in g: sym[int(x["entrezGeneId"])] = x["hugoGeneSymbol"]
co["gene"] = co.entrez.map(sym)
co = co.dropna(subset=["gene"])
print("with symbol:", len(co))

# ---------- 2. 基因位置 (Ensembl, GRCh38) ----------
poscache = os.path.join(tmp, "ensembl_positions.json")
pos = json.load(open(poscache)) if os.path.exists(poscache) else {}
todo = [g for g in co.gene.unique() if g not in pos]
for i in range(0, len(todo), 900):
    chunk = todo[i:i+900]
    for attempt in range(3):
        try:
            r = jpost(f"{ENS}/lookup/symbol/homo_sapiens", {"symbols": chunk}); break
        except Exception as e:
            print("ensembl retry", attempt, e); time.sleep(3); r = {}
    for g in chunk:
        v = r.get(g)
        pos[g] = [v["seq_region_name"], v["start"], v["end"]] if v else None
    json.dump(pos, open(poscache, "w"))
    print("positions", min(i+900, len(todo)), "/", len(todo))
co["chr"] = co.gene.map(lambda g: pos.get(g)[0] if pos.get(g) else None)
co["start"] = co.gene.map(lambda g: pos.get(g)[1] if pos.get(g) else np.nan)
co["end"] = co.gene.map(lambda g: pos.get(g)[2] if pos.get(g) else np.nan)
main = [str(c) for c in range(1, 23)] + ["X", "Y"]
co.loc[~co.chr.isin(main), ["chr", "start", "end"]] = [None, np.nan, np.nan]
print("mapped to primary chromosome:", co.chr.notna().sum())

# ---------- 3. cytoband (UCSC hg38) ----------
cb_path = os.path.join(tmp, "cytoBand_hg38.txt.gz")
if not os.path.exists(cb_path):
    open(cb_path, "wb").write(urllib.request.urlopen("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz", timeout=60).read())
cb = pd.read_csv(cb_path, sep="\t", header=None, names=["chrom", "s", "e", "band", "stain"])
cb["chr"] = cb.chrom.str.replace("chr", "")
def band_of(chrom, start):
    if chrom is None or pd.isna(start): return None
    x = cb[(cb.chr == chrom) & (cb.s <= start) & (cb.e > start)]
    return (chrom + x.band.iloc[0]) if len(x) else None
co["cytoband"] = [band_of(c, s) for c, s in zip(co.chr, co.start)]
co["on_19p13_3"] = co.cytoband.eq("19p13.3")
lon = co[co.gene == "LONP1"]
LON_START = int(lon.start.iloc[0]) if len(lon) else 5691845
co["dist_to_LONP1_Mb"] = np.where(co.chr == "19", (co.start - LON_START) / 1e6, np.nan)
co = co.sort_values("rho", ascending=False).reset_index(drop=True)
co["rank"] = np.arange(1, len(co) + 1)
co[["gene", "entrez", "rho", "p", "chr", "start", "end", "cytoband", "on_19p13_3", "dist_to_LONP1_Mb", "rank"]].to_csv(
    os.path.join(rd, "30_coexpression_LONP1_CESC_genomewide.csv"), index=False)

# ---------- 4. 位置统计 ----------
ann = co[co.chr.notna()].copy()
lines = []
def L(s): print(s); lines.append(s)
L(f"genome-wide co-expression: {len(co)} genes; {len(ann)} with GRCh38 position")
L(f"|rho|>=0.3: {(co.rho.abs()>=0.3).sum()} (pos {(co.rho>=0.3).sum()}, neg {(co.rho<=-0.3).sum()})  [CSV06 recorded 970]")
bg = ann.on_19p13_3.mean()
L(f"background fraction of 19p13.3 genes: {ann.on_19p13_3.sum()}/{len(ann)} = {bg:.4f}")
for N in (20, 30, 50, 100):
    top = ann.head(N)
    L(f"top {N} positive partners on 19p13.3: {top.on_19p13_3.sum()}/{N} = {top.on_19p13_3.mean():.2f}")
for thr in (0.3, 0.4, 0.5):
    sel = ann[ann.rho >= thr]; rest = ann[ann.rho < thr]
    tab = [[sel.on_19p13_3.sum(), (~sel.on_19p13_3).sum()], [rest.on_19p13_3.sum(), (~rest.on_19p13_3).sum()]]
    orr, pv = fisher_exact(tab, alternative="greater")
    L(f"rho>={thr}: n={len(sel)}, 19p13.3={sel.on_19p13_3.sum()} ({sel.on_19p13_3.mean():.3f}); Fisher OR={orr:.1f}, p={pv:.2e}")
chr19 = ann[ann.chr == "19"]
near = chr19[chr19.dist_to_LONP1_Mb.abs() <= 2]; far = chr19[chr19.dist_to_LONP1_Mb.abs() > 2]
L(f"chr19 genes within 2 Mb of LONP1: n={len(near)}, median rho={near.rho.median():.3f}; >2 Mb: n={len(far)}, median rho={far.rho.median():.3f}; MWU p={mannwhitneyu(near.rho, far.rho).pvalue:.2e}")
L(f"19p13.3 genes: n={ann[ann.on_19p13_3].shape[0]}, median rho={ann[ann.on_19p13_3].rho.median():.3f}; other genes median rho={ann[~ann.on_19p13_3].rho.median():.3f}")
# 保留的功能模块基因（非19p13.3）示例
for g in ["CLPP", "HSPD1", "HSPA9", "HSPE1", "TRAP1", "PHB1", "PHB2", "DNAJA3", "ATF4", "ATF5", "DDIT3", "SHMT2", "MTHFD2",
          "PMPCB", "YME1L1", "AFG3L2", "SPG7", "TIMM13", "TIMM44", "NDUFA11", "MRPL54", "POLRMT", "TFAM", "NRF1", "PPARGC1A",
          "FIS1", "MFF", "DNM1L", "MFN1", "MFN2", "OPA1"]:
    x = co[co.gene == g]
    if len(x): L(f"  {g}: rho={x.rho.iloc[0]:.3f} p={x.p.iloc[0]:.2e} band={x.cytoband.iloc[0]} rank={x['rank'].iloc[0]}")

# ---------- 5. 去19p13.3后的富集 (Enrichr) ----------
def enrichr(genes, libs, tag):
    data = urllib.parse.urlencode({"list": "\n".join(genes), "description": tag}).encode()
    import urllib.parse as up
    body = ("--XX\r\nContent-Disposition: form-data; name=\"list\"\r\n\r\n" + "\n".join(genes) +
            "\r\n--XX\r\nContent-Disposition: form-data; name=\"description\"\r\n\r\n" + tag + "\r\n--XX--\r\n").encode()
    req = urllib.request.Request("https://maayanlab.cloud/Enrichr/addList", data=body,
                                 headers={"Content-Type": "multipart/form-data; boundary=XX"})
    uid = json.load(urllib.request.urlopen(req, timeout=120))["userListId"]
    rows = []
    for lib in libs:
        r = json.load(urllib.request.urlopen(f"https://maayanlab.cloud/Enrichr/enrich?userListId={uid}&backgroundType={lib}", timeout=120))
        for e in r[lib]:
            rows.append({"lib": lib, "term": e[1], "p": e[2], "padj": e[6], "n_overlap": len(e[5]), "genes": ";".join(e[5])})
    return pd.DataFrame(rows)
import urllib.parse
libs = ["GO_Biological_Process_2023", "GO_Cellular_Component_2023", "KEGG_2021_Human"]
pos_all = co[(co.rho >= 0.4)].gene.tolist()
pos_ex = co[(co.rho >= 0.4) & (~co.on_19p13_3)].gene.tolist()
L(f"enrichment input: rho>=0.4 all n={len(pos_all)}; excluding 19p13.3 n={len(pos_ex)}")
e = enrichr(pos_ex, libs, "LONP1_pos_rho0.4_excl19p13.3")
e = e.sort_values(["lib", "padj"])
e.to_csv(os.path.join(rd, "31_enrichment_enrichr_pos_excl19p13.3.csv"), index=False)
for lib in libs:
    L(f"[{lib}] top5 (excl. 19p13.3):")
    for _, r in e[e.lib == lib].head(5).iterrows():
        L(f"   {r.term}  padj={r.padj:.1e} n={r.n_overlap}")
open(os.path.join(rd, "32_positional_stats.txt"), "w").write("\n".join(lines) + "\n")
print("DONE15")
