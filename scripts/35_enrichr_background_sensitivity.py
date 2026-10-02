# 35_enrichr_background_sensitivity.py —— Enrichr 自定义背景敏感性分析（审阅意见第12项）+ 补全 Table S2 "all genes" 交集基因
# 复现口径：
#   接口: Speedrichr (https://maayanlab.cloud/speedrichr/api: addList -> addbackground -> backgroundenrich)，背景与文库基因取交集由服务器完成
#   前景/背景:
#     bulk pos all (rho>=0.4, n=164)            背景 = 30号表全部可排名基因
#     bulk pos excl 19p13.3 (n=102)             前景与背景均剔除19p13.3基因
#     bulk neg (rho<=-0.35, n=96)               背景 = 30号表全部可排名基因
#     epithelium top300 raw / sample-centered   背景 = 37号表全部14,140个可排名基因
#   另: 默认背景重跑 bulk pos all，补全 Table S2 "all genes" 组的 overlapping_gene_symbols
# 输出: raw_data/60_enrichment_custom_background.csv, 61_enrichment_bulk_pos_all_default_with_genes.csv, 62_enrichment_background_sensitivity.txt
# 用法: python3 35_enrichr_background_sensitivity.py <分析根目录>
import sys, os, json, time, urllib.request, urllib.parse
import pandas as pd

rd = os.path.join(sys.argv[1], "raw_data")
SR = "https://maayanlab.cloud/speedrichr/api"
LIBS = ["GO_Biological_Process_2023", "GO_Cellular_Component_2023", "KEGG_2021_Human"]

def post_form(url, fields):
    b = "".join(f"--XX\r\nContent-Disposition: form-data; name=\"{k}\"\r\n\r\n{v}\r\n" for k, v in fields.items()) + "--XX--\r\n"
    req = urllib.request.Request(url, data=b.encode(), headers={"Content-Type": "multipart/form-data; boundary=XX"})
    for i in range(4):
        try: return json.load(urllib.request.urlopen(req, timeout=180))
        except Exception as ex:
            print("retry", url, ex); time.sleep(5)
    raise RuntimeError(url)

def enrich(genes, libs, background=None):
    uid = post_form(f"{SR if background is not None else "https://maayanlab.cloud/Enrichr"}/addList", {"list": "\n".join(genes), "description": "LONP1"})["userListId"]
    if background is not None:
        bid = post_form(f"{SR}/addbackground", {"background": "\n".join(background)})["backgroundid"]
    rows = []
    for lib in libs:
        if background is None:
            r = json.load(urllib.request.urlopen(f"https://maayanlab.cloud/Enrichr/enrich?userListId={uid}&backgroundType={lib}", timeout=180))
        else:
            r = post_form(f"{SR}/backgroundenrich", {"userListId": uid, "backgroundType": lib, "backgroundid": bid})
        for e in r[lib]:
            rows.append({"library": lib, "term": e[1], "p": e[2], "padj": e[6], "n_overlap": len(e[5]), "genes": ";".join(e[5])})
    return pd.DataFrame(rows)

co = pd.read_csv(os.path.join(rd, "30_coexpression_LONP1_CESC_genomewide.csv"))
ep = pd.read_csv(os.path.join(rd, "37_epithelial_pseudobulk_genomewide_rho.csv"))
bg_bulk = co.gene.tolist(); bg_bulk_ex = co[~co.on_19p13_3].gene.tolist(); bg_epi = ep.gene.tolist()
sets = [
 ("bulk_pos_all_n164", co[co.rho >= 0.4].gene.tolist(), bg_bulk, LIBS),
 ("bulk_pos_excl19p13.3_n102", co[(co.rho >= 0.4) & ~co.on_19p13_3].gene.tolist(), bg_bulk_ex, LIBS),
 ("bulk_neg_rho<=-0.35_n96", co[co.rho <= -0.35].gene.tolist(), bg_bulk, LIBS),
 ("epi_top300_raw", ep.sort_values("rho_pb", ascending=False).head(300).gene.tolist(), bg_epi, LIBS),
 ("epi_top300_sample_centered", ep.sort_values("rho_pb_sample_centered", ascending=False).head(300).gene.tolist(), bg_epi, LIBS),
]
out, L = [], []
def P(s): print(s, flush=True); L.append(s)
P(f"background sizes: bulk={len(bg_bulk)}, bulk excl 19p13.3={len(bg_bulk_ex)}, epithelium={len(bg_epi)}")
for name, fg, bg, libs in sets:
    P(f"{name}: foreground n={len(fg)} (all in background: {set(fg) <= set(bg)})")
    e = enrich(fg, libs, bg).sort_values(["library", "padj"]); e.insert(0, "gene_set", name); out.append(e)
    for lib in libs:
        sub = e[e.library == lib].head(5)
        P(f"  [{lib}] n_terms={int((e.library == lib).sum())}; sig(padj<0.05)={int(((e.library == lib) & (e.padj < 0.05)).sum())}; top: " +
          "; ".join(f"{r.term} padj={r.padj:.1e} n={r.n_overlap}" for r in sub.itertuples()))
pd.concat(out).to_csv(os.path.join(rd, "60_enrichment_custom_background.csv"), index=False)
d = enrich(sets[0][1], LIBS).sort_values(["library", "padj"])
d.to_csv(os.path.join(rd, "61_enrichment_bulk_pos_all_default_with_genes.csv"), index=False)
P(f"default-background rerun of bulk_pos_all: {len(d)} terms")
open(os.path.join(rd, "62_enrichment_background_sensitivity.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
print("DONE35")
