# 02_harmonize_mr.py —— cis-eQTL 双样本 MR：数据协调、工具变量选择（plink clump）、Wald/IVW/Egger/加权中位数、Steiger；并导出 coloc 输入
# 统一坐标：GRCh37 "chr:pos"；GRCh38 来源（FinnGen、GTEx）用 UCSC hg38ToHg19 链文件（pyliftover）转换。
# 等位基因对齐到 1000G phase3 的 REF/ALT（效应等位 = ALT）；链互补自动翻转；回文 SNP 若 EUR AF 在 0.42–0.58 之间剔除。
# 暴露：eQTLGen 血液（显著 cis 文件，Zscore→beta/se：beta=z/sqrt(2f(1-f)(N+z²))，f 取 1000G EUR ALT 频率）；
#       GTEx v8 血液/子宫/阴道/成纤维细胞（eQTL Catalogue，beta/se 为 ALT 效应）。
# 结局：FinnGen R12 四个终点；Pan-UKB C53、180.1、180.3（EUR 列）；BBJ CeC（东亚，LD 用 1000G EAS）。
# 工具变量：P<5e-8，plink1.9 --clump r²<0.1、1000 kb（欧洲结局用 EUR 参考，BBJ 用 EAS 参考）；F=(b/se)²，剔除 F<10。
# MR：1 个 SNP 用 Wald ratio（一阶 SE）；>=2 个用 IVW（多于 3 个且 Cochran Q 显著时用乘性随机效应）；>=3 个加 MR-Egger 与加权中位数（bootstrap 1000）。
# Steiger：r²=t²/(t²+N)，比较暴露与结局 r²（结局 N 用有效样本量 4/(1/Ncase+1/Nctrl)）。
# 输出：results/mr_<gene>.csv（每个 暴露×结局 一行）、results/instruments_<gene>.csv、coloc_input/<gene>__<exposure>__<outcome>.tsv
# 用法: python -P 02_harmonize_mr.py <03_MR_coloc 目录> <regions 目录> <ref 目录> <eQTLGen 显著文件 或 -> <基因逗号列表>
import sys, os, gzip, subprocess, warnings
import numpy as np, pandas as pd
from scipy.stats import norm
from pyliftover import LiftOver
warnings.filterwarnings("ignore")
base, regd, refd, eqtlgen_f, genes_arg = sys.argv[1:6]
res = os.path.join(base, "results"); cid = os.path.join(base, "coloc_input"); os.makedirs(cid, exist_ok=True)
coords = pd.read_csv(os.path.join(res, "genes_coordinates.csv")).set_index("gene")
lo = LiftOver(os.path.join(refd, "hg38ToHg19.over.chain.gz"))
OUT = {  # 名称: (文件, 解析类型, 病例, 对照, 人群)
    "FinnGen_CC": ("FG_C3_CERVIX_UTERI_EXALLC.tsv", "fg", 878, 222078, "EUR"),
    "FinnGen_CIS": ("FG_CD2_INSITU_CERVIX_UTERI_EXALLC.tsv", "fg", 3701, 219058, "EUR"),
    "FinnGen_SCC": ("FG_C3_CERVIX_SQUAM_EXALLC.tsv", "fg", 222, 222078, "EUR"),
    "FinnGen_ADC": ("FG_C3_CERVIX_ADENO_EXALLC.tsv", "fg", 142, 222078, "EUR"),
    "UKB_C53": ("UKB_C53.tsv", "ukb", 273, 420258, "EUR"),
    "UKB_180.1": ("UKB_180.1.tsv", "ukb", 299, 198725, "EUR"),
    "UKB_180.3": ("UKB_180.3.tsv", "ukb", 2186, 198725, "EUR"),
    "BBJ_CeC": ("BBJ_CeC.tsv", "bbj", 967, 60614, "EAS"),
}
EXP = {"eQTLGen_blood": (None, 31684), "GTEx_blood": ("GTEx_blood.tsv", 670), "GTEx_uterus": ("GTEx_uterus.tsv", 129),
       "GTEx_vagina": ("GTEx_vagina.tsv", 141), "GTEx_fibroblast": ("GTEx_fibroblast.tsv", 483)}
GTEX_COLS = ["molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples", "maf", "pvalue", "beta", "se",
             "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id", "median_tpm", "rsid"]
COMP = str.maketrans("ACGT", "TGCA")
def wsl(cmd):
    return subprocess.run(["wsl", "-d", "Ubuntu-24.04", "--", "bash", "-c", cmd], capture_output=True, text=True,
                          env={**os.environ, "MSYS_NO_PATHCONV": "1"})
def winpath_to_wsl(p):
    p = os.path.abspath(p).replace("\\", "/"); return "/mnt/" + p[0].lower() + p[2:]

def kg_table(gdir):
    rows = []
    with open(os.path.join(gdir, "1kg_phase3.vcf")) as f:
        for l in f:
            if l.startswith("#"): continue
            x = l.split("\t", 8)
            info = dict(kv.split("=", 1) for kv in x[7].split(";") if "=" in kv)
            rows.append((x[2], x[3], x[4], float(info.get("EUR_AF", "nan").split(",")[0]), float(info.get("EAS_AF", "nan").split(",")[0])))
    return pd.DataFrame(rows, columns=["key", "kref", "kalt", "eur_af", "eas_af"]).drop_duplicates("key", keep=False).set_index("key")

def bfiles(gdir):
    out = {}
    pan = pd.read_csv(os.path.join(refd, "1kg.panel"), sep="\t")
    for pop in ("EUR", "EAS"):
        pre = os.path.join(gdir, f"1kg_{pop}")
        if not os.path.exists(pre + ".bed"):
            ids = pan[pan.super_pop == pop]["sample"]
            pd.DataFrame({"a": ids, "b": ids}).to_csv(pre + ".keep", sep="\t", header=False, index=False)
            r = wsl(f"cd '{winpath_to_wsl(gdir)}' && plink1.9 --vcf 1kg_phase3.vcf --double-id --keep 1kg_{pop}.keep --snps-only just-acgt "
                    f"--maf 0.005 --make-bed --out 1kg_{pop} > /dev/null")
        out[pop] = pre
    return out

def lift38(chrom, pos):
    r = lo.convert_coordinate("chr" + str(chrom), int(pos) - 1)
    return f"{str(r[0][0]).replace('chr', '')}:{r[0][1] + 1}" if r else None

def align(df, kg, pop):
    # df: key, ea(效应等位), oa, beta, se ；返回对齐到 kalt 的 beta
    d = df.join(kg, on="key", how="inner")
    same = (d.ea == d.kalt) & (d.oa == d.kref); swap = (d.ea == d.kref) & (d.oa == d.kalt)
    ea_c = d.ea.str.translate(COMP); oa_c = d.oa.str.translate(COMP)
    same_c = (ea_c == d.kalt) & (oa_c == d.kref); swap_c = (ea_c == d.kref) & (oa_c == d.kalt)
    pal = ((d.kref + d.kalt).isin(["AT", "TA", "CG", "GC"]))
    af = d.eur_af if pop == "EUR" else d.eas_af
    amb = pal & af.between(0.42, 0.58)
    sign = np.where(same | (same_c & ~pal), 1, np.where(swap | (swap_c & ~pal), -1, 0))
    d["beta"] = d.beta * sign; d = d[(sign != 0) & ~amb]
    return d

def load_exposure(name, gdir, gene, kg):
    if name == "eQTLGen_blood":
        if EQ is None: return None
        e = EQ[EQ.GeneSymbol == coords.loc[gene, "alias_GRCh37"]].copy()
        if not len(e): e = EQ[EQ.Gene == coords.loc[gene, "ensg_GRCh38"]].copy()
        if not len(e): return None
        e["key"] = e.SNPChr.astype(str) + ":" + e.SNPPos.astype(str)
        e = e.join(kg[["eur_af", "kalt", "kref"]], on="key", how="inner")
        f = np.where(e.AssessedAllele == e.kalt, e.eur_af, 1 - e.eur_af); n = e.NrSamples.values
        den = np.sqrt(2 * f * (1 - f) * (n + e.Zscore ** 2))
        e = pd.DataFrame({"key": e.key.values, "rsid": e.SNP.values, "ea": e.AssessedAllele.values, "oa": e.OtherAllele.values,
                          "beta": e.Zscore.values / den, "se": 1 / den, "p": e.Pvalue.values, "n": n})
        return e.dropna()
    fn, N = EXP[name]; p = os.path.join(gdir, fn)
    if not os.path.exists(p) or os.path.getsize(p) == 0: return None
    g = pd.read_csv(p, sep="\t", header=None, names=GTEX_COLS, low_memory=False)
    g = g[(g.type == "SNP") & (g.molecular_trait_id == coords.loc[gene, "ensg_GRCh38"])]
    g["key"] = [lift38(c, x) for c, x in zip(g.chromosome, g.position)]
    return pd.DataFrame({"key": g.key, "rsid": g.rsid, "ea": g.alt, "oa": g.ref, "beta": g.beta, "se": g.se, "p": g.pvalue,
                         "n": N}).dropna(subset=["key"])

def load_outcome(name, gdir):
    fn, kind, nca, nco, pop = OUT[name]; p = os.path.join(gdir, fn)
    if not os.path.exists(p): return None
    if kind == "fg":
        o = pd.read_csv(p, sep="\t", low_memory=False)
        o["key"] = [lift38(c, x) for c, x in zip(o.chrom, o.pos)]
        o = pd.DataFrame({"key": o.key, "ea": o.alt, "oa": o.ref, "beta": o.beta, "se": o.sebeta, "p": o.pval, "af_out": o.af_alt})
    elif kind == "ukb":
        h = open(os.path.join(base, "_meta", name + ".header")).read().strip().split("\t")
        o = pd.read_csv(p, sep="\t", header=None, names=h, na_values=["NA"], low_memory=False)
        sfx = "EUR"
        o = pd.DataFrame({"key": o.chr.astype(str) + ":" + o.pos.astype(str), "ea": o.alt, "oa": o.ref, "beta": o[f"beta_{sfx}"],
                          "se": o[f"se_{sfx}"], "p": 10 ** (-o[f"neglog10_pval_{sfx}"]), "af_out": o[f"af_controls_{sfx}"]})
    else:
        o = pd.read_csv(p, sep="\t", low_memory=False)
        o = pd.DataFrame({"key": o.chr.astype(str) + ":" + o.position.astype(str), "ea": o.alt, "oa": o.ref, "beta": o.beta,
                          "se": o.sebeta, "p": o.pvalue, "af_out": o.af})
    return o.dropna(subset=["key", "beta", "se"]).drop_duplicates("key", keep=False)

def clump(ins, bfile, gdir, tag):
    if len(ins) <= 1: return ins
    f = os.path.join(gdir, f"clump_{tag}.txt"); ins[["key", "p_x"]].rename(columns={"key": "SNP", "p_x": "P"}).to_csv(f, sep="\t", index=False)
    r = wsl(f"cd '{winpath_to_wsl(gdir)}' && plink1.9 --bfile {os.path.basename(bfile)} --clump clump_{tag}.txt --clump-p1 5e-8 "
            f"--clump-p2 1 --clump-r2 0.1 --clump-kb 1000 --out clump_{tag} > /dev/null; cat clump_{tag}.clumped 2>/dev/null")
    lead = [l.split()[2] for l in r.stdout.strip().split("\n")[1:] if l.strip()]
    if not lead:  # 参考面板缺失这些 SNP 时退回：只取最显著的一个
        return ins.nsmallest(1, "p_x")
    return ins[ins.key.isin(lead)]

def mr(h):
    bx, by, sx, sy = h.bx.values, h.by.values, h.sx.values, h.sy.values
    out = {"n_snp": len(h)}
    if len(h) == 1:
        b = by[0] / bx[0]; se = sy[0] / abs(bx[0]); out.update(method="Wald ratio", b=b, se=se)
    else:
        w = bx ** 2 / sy ** 2; b = np.sum(bx * by / sy ** 2) / np.sum(w); se = 1 / np.sqrt(np.sum(w))
        res_ = (by - b * bx) / sy; Q = float(np.sum(res_ ** 2)); df = len(h) - 1
        from scipy.stats import chi2
        out.update(Q=Q, Q_p=float(chi2.sf(Q, df)))
        if len(h) > 3 and out["Q_p"] < 0.05: se = se * np.sqrt(Q / df); out["method"] = "IVW (multiplicative RE)"
        else: out["method"] = "IVW (fixed)"
        out.update(b=b, se=se)
        if len(h) >= 3:
            s = np.sign(bx); X = np.column_stack([np.ones(len(h)), bx * s]); W = np.diag(1 / sy ** 2)
            beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ (by * s)); resid = by * s - X @ beta
            sig2 = max(1, float(resid @ W @ resid) / (len(h) - 2)); cov = np.linalg.inv(X.T @ W @ X) * sig2
            out.update(egger_b=beta[1], egger_se=np.sqrt(cov[1, 1]), egger_int=beta[0], egger_int_p=2 * norm.sf(abs(beta[0] / np.sqrt(cov[0, 0]))))
            def wm(bx_, by_, sy_):
                r = by_ / bx_; wr = (bx_ / sy_) ** 2; o = np.argsort(r); r, wr = r[o], wr[o]
                c = (np.cumsum(wr) - 0.5 * wr) / wr.sum(); return np.interp(0.5, c, r)
            rng = np.random.default_rng(1)
            boots = [wm(rng.normal(bx, sx), rng.normal(by, sy), sy) for _ in range(1000)]
            out.update(wm_b=wm(bx, by, sy), wm_se=float(np.std(boots)))
    out["p"] = 2 * norm.sf(abs(out["b"] / out["se"]))
    out["OR"] = np.exp(out["b"]); out["OR_lo"] = np.exp(out["b"] - 1.96 * out["se"]); out["OR_hi"] = np.exp(out["b"] + 1.96 * out["se"])
    return out

EQ = None
if eqtlgen_f != "-" and os.path.exists(eqtlgen_f):
    EQ = pd.read_csv(eqtlgen_f, sep="\t", usecols=["Pvalue", "SNP", "SNPChr", "SNPPos", "AssessedAllele", "OtherAllele", "Zscore",
                                                     "Gene", "GeneSymbol", "NrSamples"], low_memory=False)
for gene in genes_arg.split(","):
    gdir = os.path.join(regd, gene)
    if not os.path.exists(os.path.join(gdir, "1kg_phase3.vcf")): print("skip", gene, "(no region data)"); continue
    kg = kg_table(gdir); bf = bfiles(gdir)
    outs = {k: load_outcome(k, gdir) for k in OUT}
    rows, inst_rows = [], []
    for en in EXP:
        ex = load_exposure(en, gdir, gene, kg)
        if ex is None or not len(ex): rows.append({"gene": gene, "exposure": en, "outcome": "-", "note": "no exposure data"}); continue
        exa = align(ex, kg, "EUR")
        for on, o in outs.items():
            pop = OUT[on][4]; nca, nco = OUT[on][2], OUT[on][3]
            if o is None or not len(o): continue
            oa = align(o, kg, pop)
            m = exa[["key", "rsid", "beta", "se", "p", "n", "eur_af", "eas_af"]].merge(oa[["key", "beta", "se", "p", "af_out"]], on="key",
                                                                                    suffixes=("_x", "_y"))
            # coloc 输入（全区段）
            m.rename(columns={"beta_x": "bx", "se_x": "sx", "p_x": "px", "beta_y": "by", "se_y": "sy", "p_y": "py"}).to_csv(
                os.path.join(cid, f"{gene}__{en}__{on}.tsv"), sep="\t", index=False)
            ins = m[m.p_x < 5e-8].copy(); ins["F"] = (ins.beta_x / ins.se_x) ** 2; ins = ins[ins.F > 10]
            if not len(ins):
                rows.append({"gene": gene, "exposure": en, "outcome": on, "note": "no instrument at P<5e-8"}); continue
            ins = clump(ins, bf[pop], gdir, f"{en}_{pop}")
            h = ins.rename(columns={"beta_x": "bx", "se_x": "sx", "beta_y": "by", "se_y": "sy"})
            r = mr(h)
            neff = 4 / (1 / nca + 1 / nco)
            r2x = float(np.sum((h.bx / h.sx) ** 2 / ((h.bx / h.sx) ** 2 + h.n))); r2y = float(np.sum((h.by / h.sy) ** 2 / ((h.by / h.sy) ** 2 + neff)))
            rows.append({"gene": gene, "exposure": en, "outcome": on, "cases": nca, "controls": nco, "pop": pop, **r,
                         "mean_F": float(h.F.mean()), "r2_exposure": r2x, "r2_outcome": r2y, "steiger_correct_direction": r2x > r2y,
                         "snps": ";".join(h.rsid.astype(str))})
            for x in h.itertuples(): inst_rows.append({"gene": gene, "exposure": en, "outcome": on, "rsid": x.rsid, "key": x.key, "bx": x.bx,
                                                       "sx": x.sx, "px": x.p_x, "F": x.F, "by": x.by, "sy": x.sy, "py": x.p_y})
    R = pd.DataFrame(rows); R.to_csv(os.path.join(res, f"mr_{gene}.csv"), index=False)
    pd.DataFrame(inst_rows).to_csv(os.path.join(res, f"instruments_{gene}.csv"), index=False)
    cols = [c for c in ["exposure", "outcome", "n_snp", "method", "OR", "OR_lo", "OR_hi", "p", "mean_F", "steiger_correct_direction", "note"] if c in R]
    print(f"== {gene}"); print(R[cols].round(4).to_string(index=False))
