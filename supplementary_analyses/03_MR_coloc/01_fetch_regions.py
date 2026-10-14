# 01_fetch_regions.py —— 远程 tabix / API 按基因 ±1 Mb 窗口抓取暴露与结局汇总统计，以及 1000G 参考基因型（在 WSL python3 + pysam 下运行）
# 结局：FinnGen R12（GRCh38）C3_CERVIX_UTERI_EXALLC（主）、CD2_INSITU_CERVIX_UTERI_EXALLC、C3_CERVIX_SQUAM_EXALLC、C3_CERVIX_ADENO_EXALLC；
#       Pan-UKB（GRCh37）icd10-C53、phecode-180.1、phecode-180.3（女性）；BBJ CeC（PheWeb.jp region API，GRCh37）。
# 暴露：GTEx v8（eQTL Catalogue，GRCh38）血液 QTD000356、子宫 QTD000346、阴道 QTD000351、成纤维细胞 QTD000216，只保留目标基因行；
#       eQTLGen 由 02 号脚本从本地文件提取。
# LD：1000 Genomes phase3（GRCh37）区段 VCF，供 plink 计算 EUR/EAS LD。
# 输出：<out>/<gene>/<source>.tsv（原始行 + 表头）
# 用法: python3 01_fetch_regions.py <genes_coordinates.csv> <out 目录> [基因逗号列表] [窗口 bp，默认 1000000]
import sys, os, csv, json, time, urllib.request, pysam
coords, outd = sys.argv[1], sys.argv[2]
only = set(sys.argv[3].split(",")) if len(sys.argv) > 3 and sys.argv[3] else None
W = int(sys.argv[4]) if len(sys.argv) > 4 else 1000000
FG = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_{}.gz"
FG_EP = ["C3_CERVIX_UTERI_EXALLC", "CD2_INSITU_CERVIX_UTERI_EXALLC", "C3_CERVIX_SQUAM_EXALLC", "C3_CERVIX_ADENO_EXALLC"]
UKB = {"UKB_C53": "icd10-C53-both_sexes", "UKB_180.1": "phecode-180.1-females", "UKB_180.3": "phecode-180.3-females"}
UKB_URL = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/{}.tsv.bgz"
UKB_TBI = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files_tabix/{}.tsv.bgz.tbi"
GTEX = {"GTEx_blood": "QTD000356", "GTEx_uterus": "QTD000346", "GTEx_vagina": "QTD000351", "GTEx_fibroblast": "QTD000216"}
GTEX_URL = "http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/{0}/{0}.all.tsv.gz"
KG = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
def retry(fn, n=4):
    for i in range(n):
        try: return fn()
        except Exception as e:
            print("  retry", i, e, flush=True); time.sleep(2 ** (i + 1))
    raise RuntimeError("failed after retries")
def tab(url, chrom, a, b, path, header=None, index=None, keep=None):
    if os.path.exists(path): return
    def go():
        f = pysam.TabixFile(url, index=index) if index else pysam.TabixFile(url)
        h = header or (f.header[-1] if f.header else None)
        rows = [r for r in f.fetch(chrom, max(0, a), b) if (keep is None or keep(r))]
        with open(path + ".tmp", "w") as o:
            if h: o.write(h.lstrip("#") + "\n")
            o.write("\n".join(rows) + ("\n" if rows else ""))
        os.replace(path + ".tmp", path); return len(rows)
    n = retry(go); print(f"  {os.path.basename(path)}: {n} rows", flush=True)
genes = list(csv.DictReader(open(coords)))
for g in genes:
    name = g["gene"]
    if only and name not in only: continue
    d = os.path.join(outd, name); os.makedirs(d, exist_ok=True); print(name, flush=True)
    c38, s38, e38 = g["chr_GRCh38"].split(".")[0], int(float(g["start_GRCh38"])), int(float(g["end_GRCh38"]))
    c37, s37, e37 = g["chr_GRCh37"].split(".")[0], int(float(g["start_GRCh37"])), int(float(g["end_GRCh37"]))
    for ep in FG_EP:
        tab(FG.format(ep), c38, s38 - W, e38 + W, os.path.join(d, f"FG_{ep}.tsv"))
    for k, v in UKB.items():
        tab(UKB_URL.format(v), c37, s37 - W, e37 + W, os.path.join(d, f"{k}.tsv"), index=UKB_TBI.format(v))
    ensg = g["ensg_GRCh38"]
    for k, q in GTEX.items():
        tab(GTEX_URL.format(q), c38, s38 - W, e38 + W, os.path.join(d, f"{k}.tsv"), keep=lambda r, e=ensg: e in r)
    p = os.path.join(d, "BBJ_CeC.tsv")
    if not os.path.exists(p):
        def bbj():
            url = ("https://pheweb.jp/api/region/CeC/lz-results/?filter=analysis%20in%203%20and%20chromosome%20in%20%27{}%27"
                   "%20and%20position%20ge%20{}%20and%20position%20le%20{}").format(c37, max(1, s37 - W), e37 + W)
            return json.load(urllib.request.urlopen(url, timeout=600))["data"]
        dat = retry(bbj); keys = list(dat.keys())
        with open(p, "w") as o:
            o.write("\t".join(keys) + "\n")
            for i in range(len(dat[keys[0]])): o.write("\t".join(str(dat[k][i]) for k in keys) + "\n")
        print(f"  BBJ_CeC.tsv: {len(dat[keys[0]])} rows", flush=True)
    vcf = os.path.join(d, "1kg_phase3.vcf")
    if os.path.exists(vcf) and os.path.getsize(vcf) < 100000: os.remove(vcf)  # 旧版本（要求 rs 号）写出的空文件
    if not os.path.exists(vcf):
        def kg():
            # 1000G phase3 VCF 的 ID 列为 "."（无 rs 号），改用 GRCh37 "chr:pos" 作为 ID；只保留双等位 SNP
            f = pysam.TabixFile(KG.format(c37)); hdr = list(f.header)
            with open(vcf + ".tmp", "w") as o:
                o.write("\n".join(hdr) + "\n")
                for r in f.fetch(c37, max(0, s37 - W), e37 + W):
                    x = r.split("\t", 8)
                    if len(x[3]) == 1 and len(x[4]) == 1 and x[4] in "ACGT":
                        x[2] = f"{x[0]}:{x[1]}"; o.write("\t".join(x) + "\n")
            os.replace(vcf + ".tmp", vcf)
        retry(kg); print("  1kg_phase3.vcf written", flush=True)
print("FETCH_DONE")
