# 01_generalize.py -- step 2 of the generalization analysis (cervical cancer, TCGA-CESC PanCancer Atlas)
#
# For every seed gene, run the same pipeline used for LONP1 and report three readouts:
#   1. positional dominance : fraction of the seed's top-50 raw co-expression partners on the seed's own chromosome arm
#   2. retention            : of the 10 gene sets most enriched among the seed's raw partners (what a GSEA in a
#                             single-gene paper would report), the fraction still in the top 10% of sets after
#                             deconfounding (partial Spearman adjusted for copy number [+ arm CN + purity] and with
#                             the seed's own arm removed from the background)
#   3. specificity          : for the seed's top deconfounded gene set T*, the seed's AUC compared with the AUCs that
#                             up to 15 other members of T* (off the seed's arm) obtain for the same set when used as
#                             seeds themselves. A seed at the median of its peers is a typical member of a shared
#                             module; a seed above the 80th percentile is more tightly linked to T* than its members.
# Seeds: literature genes (literature_genes.csv), a systematic sample of dosage-sensitive genes in recurrently altered
# regions (GISTIC), random expressed genes, and LONP1.
# Partial correlation: identical to 02_partial_correlation/01_partial_spearman.py (ranks -> OLS residuals on ranked
# covariates -> Pearson of residuals). Model A = seed CN + partner CN; model B = A + seed-arm mean CN + ABSOLUTE purity.
# Gene-set enrichment readout: AUC (probability that a member gene ranks above a non-member; 0.5 = no enrichment).
#
# usage (no arguments): run from supplementary-analyses/07_generalization/ ; inputs are located relative to it.
# optional: python -I 01_generalize.py --quick   (literature + LONP1 + 20 random, 5 peers; for a fast test run)
import sys, os, glob, time
import numpy as np, pandas as pd
from scipy.stats import rankdata, hypergeom, binom, fisher_exact, mannwhitneyu
from scipy import sparse

QUICK = "--quick" in sys.argv
here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(os.path.dirname(here))
base = os.path.join(os.path.dirname(here), "02_partial_correlation")
dd = os.path.join(base, "data")
hits = glob.glob(os.path.join(root, "*", "raw_data", "30_coexpression_LONP1_CESC_genomewide.csv"))
assert hits, "cannot find raw_data/30_coexpression_LONP1_CESC_genomewide.csv"
v8 = os.path.dirname(os.path.dirname(hits[0]))
gmt = os.path.join(here, "genesets_all.gmt")
out = os.path.join(here, "results_v2_quick" if QUICK else "results_v2")
os.makedirs(out, exist_ok=True)
log = []
def L(s=""):
    print(s, flush=True); log.append(str(s))
def save_log():
    open(os.path.join(out, "01_generalize_log.txt"), "w", encoding="utf-8").write("\n".join(log) + "\n")
t0 = time.time()

# ---------------- data ----------------
def read_cbio(path):
    d = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    d = d.dropna(subset=["Entrez_Gene_Id"]); d["Entrez_Gene_Id"] = d["Entrez_Gene_Id"].astype(int)
    d = d.drop_duplicates("Entrez_Gene_Id").set_index("Entrez_Gene_Id")
    sym = d["Hugo_Symbol"]
    m = d.drop(columns=[c for c in ("Hugo_Symbol", "Cytoband") if c in d.columns]).apply(pd.to_numeric, errors="coerce")
    return m, sym

expr, esym = read_cbio(os.path.join(dd, "data_mrna_seq_v2_rsem.txt"))
cnl, _ = read_cbio(os.path.join(dd, "data_log2_cna.txt"))
gis, _ = read_cbio(os.path.join(dd, "data_cna.txt"))
pos = pd.read_csv(os.path.join(v8, "raw_data", "30_coexpression_LONP1_CESC_genomewide.csv")).drop_duplicates("entrez").set_index("entrez")
pur = pd.read_csv(os.path.join(v8, "raw_data", "_cache", "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"), sep="\t")
pur = pur[pur["call status"] == "called"]; pur["sid"] = pur["array"].str[:15]
pur = pur.drop_duplicates("sid").set_index("sid")["purity"]
LON = 9361

E = expr.dropna(axis=0, how="any"); E = E[E.std(axis=1) > 0]
S = [s for s in E.columns if s in cnl.columns]
C = cnl.loc[:, S]
genes = [g for g in E.index if g in C.index and not C.loc[g].isna().any()]
E = E.loc[genes, S]; C = C.loc[genes]
sym = esym.reindex(genes).astype(str).values
cyto = pos.cytoband.reindex(genes).astype(str)
cyto.loc[LON] = "19p13.3"
cyto = cyto.values
def arm_of(c):
    c = str(c)
    for k, ch in enumerate(c):
        if ch in "pq":
            return c[:k + 1]
    return "NA"
arm = np.array([arm_of(c) for c in cyto])
has_pos = arm != "NA"
S_pur = [s for s in S if s in pur.index and not np.isnan(pur[s])]
ip = np.array([S.index(s) for s in S_pur])
arm_cn = {}
for a in np.unique(arm[has_pos]):
    arm_cn[a] = C.values[arm == a].mean(axis=0)
gi = {g: i for i, g in enumerate(genes)}
L(f"genes with expression+CN: {len(genes)} ({has_pos.sum()} with cytoband); samples n={len(S)} (model A), n={len(S_pur)} (model B)")

def rk(M):
    return np.apply_along_axis(rankdata, -1, M).astype(np.float64)
def zs(M):
    M = M - M.mean(-1, keepdims=True); sd = np.sqrt((M ** 2).sum(-1, keepdims=True)); sd[sd == 0] = np.nan
    return M / sd
Ev, Cv = E.values, C.values
RE_A, RC_A = rk(Ev), rk(Cv)
RE_B, RC_B = rk(Ev[:, ip]), rk(Cv[:, ip])
ZE_A = zs(RE_A)
mean_expr = np.log2(Ev + 1).mean(1)
cn_expr_rho = np.nansum(zs(RE_A) * zs(RC_A), axis=1)
purv = pur[S_pur].values

def resid(M, Z):
    Q, _ = np.linalg.qr(Z)
    return M - (M @ Q) @ Q.T

def partial(i, RE, RC, extra):
    """partial Spearman of seed i with every gene, adjusting for [seed CN, extra..., partner CN]."""
    n = RE.shape[1]
    Z = np.column_stack([np.ones(n), RC[i]] + [rankdata(e) for e in extra])
    x = resid(RE[i][None, :], Z)[0]
    Y = resid(RE, Z); K = resid(RC, Z)
    kn = np.sqrt((K ** 2).sum(1)); ok = kn > 1e-6 * np.sqrt(n) * n
    xn = x / np.sqrt((x ** 2).sum())
    Yn = Y / np.sqrt((Y ** 2).sum(1, keepdims=True))
    Kn = np.zeros_like(K); Kn[ok] = K[ok] / kn[ok, None]
    a = Yn @ xn; b = Kn @ xn; c = (Yn * Kn).sum(1)
    r = np.where(ok, (a - b * c) / np.sqrt(np.clip((1 - b ** 2) * (1 - c ** 2), 1e-12, None)), a)
    r[i] = np.nan
    return r

def layers_for(i):
    raw = ZE_A @ ZE_A[i]; raw[i] = np.nan
    A = partial(i, RE_A, RC_A, [])
    extra = [purv] + ([arm_cn[arm[i]][ip]] if arm[i] in arm_cn else [])
    B = partial(i, RE_B, RC_B, extra)
    return {"raw": raw, "A": A, "B": B}

# ---------------- gene sets ----------------
if os.path.exists(gmt):
    sets = {}
    for line in open(gmt, encoding="utf8"):
        f = line.rstrip("\n").split("\t"); sets[f[0]] = set(f[2:])
    L(f"gene sets from {os.path.basename(gmt)}: {len(sets)}")
else:
    sys.exit("genesets_all.gmt not found: run  Rscript 00_dump_genesets_all.R genesets_all.gmt  first")
sym_idx = {}
for k, s in enumerate(sym):
    sym_idx.setdefault(s, k)
names, rows, cols = [], [], []
for nm, members in sets.items():
    idx = sorted({sym_idx[g] for g in members if g in sym_idx})
    if 15 <= len(idx) <= 500:
        names += [nm]; rows += [len(names) - 1] * len(idx); cols += idx
M = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(names), len(genes)))
names = np.array(names)
L(f"gene sets with 15-500 measured genes: {len(names)}")

def auc_vec(r, universe, return_z=False):
    """AUC of every gene set for correlation vector r within universe (bool). NaN if <10 members in universe.
    With return_z, also the standardized score z = (AUC - 0.5) / SD_null(n1, n0) (Mann-Whitney normal approximation),
    used to rank sets so that small sets are not favoured by chance."""
    u = universe & ~np.isnan(r)
    R = np.zeros(len(r)); R[u] = rankdata(r[u])
    U = sparse.diags(u.astype(float))
    Mu = M @ U
    n1 = np.asarray(Mu.sum(1)).ravel(); n0 = u.sum() - n1
    s = Mu @ R
    a = (s - n1 * (n1 + 1) / 2) / (n1 * n0)
    a[n1 < 10] = np.nan
    if not return_z:
        return a
    z = (a - 0.5) / np.sqrt((n0 + n1 + 1) / (12.0 * n1 * n0))
    z[n1 < MIN_SET] = np.nan
    return a, z

# ---------------- seeds ----------------
rng = np.random.default_rng(20261115)
seeds = {}  # index -> group
lit = pd.read_csv(os.path.join(here, "literature_genes.csv"))
missing = []
ALIAS = {"MPC1": ["BRP44L"]}   # older symbols used in the cBioPortal PanCancer Atlas matrix
lit_name = {}
for g in lit.gene:
    hit = next((a for a in [g] + ALIAS.get(g, []) if a in sym_idx), None)
    if hit is None: missing.append(g); continue
    seeds[sym_idx[hit]] = "literature" if g != "LONP1" else "LONP1"; lit_name[sym_idx[hit]] = g
seeds[gi[LON]] = "LONP1"
L(f"literature seeds found: {sum(v == 'literature' for v in seeds.values())}; not found: {missing}")

# systematic: dosage-sensitive genes in recurrently altered regions (GISTIC discrete calls)
G = gis.reindex(genes)[S]
gain = (G >= 1).mean(axis=1).values; loss = (G <= -1).mean(axis=1).values
alt = np.maximum(gain, loss)
ALT_FREQ = 0.25   # gene gained or lost (GISTIC) in >= 25% of tumours
cand = np.where(has_pos & (alt >= ALT_FREQ) & (cn_expr_rho >= 0.3) & (mean_expr >= np.quantile(mean_expr, 0.25)))[0]
cand = [i for i in cand if i not in seeds]
L(f"systematic candidates (alteration freq >= {ALT_FREQ}, CN-expr rho >= 0.3, expressed): {len(cand)}; "
  f"arms: {pd.Series(arm[cand]).value_counts().to_dict()}")
pick = []
for a, grp in pd.Series(cand, index=arm[cand]).groupby(level=0):
    if len(grp) >= 20:      # arm must carry a recurrent region, not a handful of genes
        pick += list(rng.choice(grp.values, size=min(8, len(grp)), replace=False))
n_sys = 0 if QUICK else 120
if len(pick) > n_sys: pick = list(rng.choice(pick, size=n_sys, replace=False))
for i in pick: seeds[int(i)] = "systematic"
rand_pool = np.where(has_pos & (mean_expr >= np.quantile(mean_expr, 0.25)))[0]
rand_pool = [i for i in rand_pool if i not in seeds]
for i in rng.choice(rand_pool, size=20 if QUICK else 100, replace=False):
    seeds[int(i)] = "random"
L("seeds: " + str(pd.Series(seeds).value_counts().to_dict()))
N_PEERS = 5 if QUICK else 30
MIN_SET = 25          # sets ranked/selected only if >= 25 members in the universe
PROGRAM_Z = 5.0       # a deconfounded program is considered present if its z >= 5 and AUC >= 0.60
PROGRAM_AUC = 0.60
SPEC_PCT, SPEC_DELTA = 90, 0.05   # seed-specific: above 90% of peers AND >= 0.05 AUC above peer median
PRIMARY = "B"
TOP_N = 100           # partners used for the over-representation (ORA) readout, as in "top co-expressed genes" analyses
POS_DOMINANT = 0.50   # positional-dominant: >= 50% of top-50 raw partners on the seed's own arm
POS_P = 1e-3          # positional-enriched: binomial P < 1e-3 against the arm's share of genes

# ---------------- run ----------------
cache = {}
def get_layers(i):
    if i not in cache: cache[i] = layers_for(i)
    return cache[i]
peer_cache = {}
def peer_layer(j, lay):
    """single layer for a peer gene, not cached (keeps memory low)."""
    if j in cache: return cache[j][lay]
    if lay == "A": return partial(j, RE_A, RC_A, [])
    extra = [purv] + ([arm_cn[arm[j]][ip]] if arm[j] in arm_cn else [])
    return partial(j, RE_B, RC_B, extra)
def peer_auc(j, set_k, lay):
    key = (j, set_k, lay)
    if key not in peer_cache:
        r = peer_layer(j, lay)
        u = np.ones(len(genes), bool) & (arm != arm[j]); u[j] = False
        peer_cache[key] = auc_vec(r, u)[set_k]
    return peer_cache[key]

def bh(p):
    p = np.asarray(p, float); q = np.full(p.shape, np.nan); ok = ~np.isnan(p)
    o = np.argsort(p[ok]); m = ok.sum(); v = p[ok][o] * m / (np.arange(m) + 1)
    v = np.minimum.accumulate(v[::-1])[::-1]; t = np.empty(m); t[o] = np.minimum(v, 1); q[ok] = t
    return q
def ora(r, universe, n=None):
    """hypergeometric over-representation of every gene set among the top-n partners of r within universe."""
    n = n or TOP_N
    u = universe & ~np.isnan(r)
    idx = np.where(u)[0]; top = idx[np.argsort(-r[idx])[:n]]
    ind = np.zeros(len(r)); ind[top] = 1
    Mu = M @ sparse.diags(u.astype(float))
    n1 = np.asarray(Mu.sum(1)).ravel(); k = Mu @ ind
    p = hypergeom.sf(k - 1, u.sum(), n1, len(top)); p[n1 < 10] = np.nan
    return p, bh(p), top

rows_out, top_out, peer_out = [], [], []
lon_check = None
for n_done, (i, grp) in enumerate(seeds.items()):
    lay = get_layers(i)
    u_all = np.ones(len(genes), bool); u_all[i] = False
    u_adj = u_all & (arm != arm[i])
    r = lay["raw"]
    f_arm = float(np.mean(arm[has_pos] == arm[i]))
    row = dict(entrez=genes[i], gene=lit_name.get(i, sym[i]), group=grp, cytoband=cyto[i], arm=arm[i], mean_log2_expr=mean_expr[i],
               cn_expr_rho=cn_expr_rho[i], gistic_alt_freq=alt[i], arm_gene_fraction=f_arm)
    # ---- 1. positional dominance: top-50 partners on the seed's own arm, raw and after CN adjustment (all genes kept) ----
    for lname in ("raw", "A", "B"):
        t50 = np.argsort(-np.nan_to_num(lay[lname], nan=-9))[:50]
        row[f"top50_same_arm_{lname}"] = float(np.mean(arm[t50] == arm[i]))
        if lname == "raw":
            row["top50_same_band_raw"] = float(np.mean(cyto[t50] == cyto[i]))
            k = int(np.sum(arm[t50] == arm[i]))
            row["pos_binom_p"] = float(binom.sf(k - 1, 50, f_arm)); row["pos_fold"] = float(k / 50 / f_arm)
    row["positional"] = ("dominant" if row["top50_same_arm_raw"] >= POS_DOMINANT else
                         "enriched" if row["pos_binom_p"] < POS_P else "none")
    # ---- 2. retention of what a paper would report ----
    #  (a) top-N partner ORA (sensitive to positional partners)   (b) genome-wide rank-based enrichment (GSEA-like)
    p_raw, q_raw, top_raw = ora(r, u_all)
    sig_raw = np.where(q_raw < 0.05)[0]
    rep_raw = sig_raw[np.argsort(p_raw[sig_raw])[:10]]
    row["ora_top_same_arm"] = float(np.mean(arm[top_raw] == arm[i]))
    row["ora_n_sig_raw"] = int(len(sig_raw))
    a_raw, z_raw = auc_vec(r, u_all, return_z=True)
    top10_raw = np.argsort(-np.nan_to_num(z_raw, nan=-1e9))[:10]
    for lname in ("A", "B"):
        p_adj, q_adj, _ = ora(lay[lname], u_adj)
        row[f"ora_n_sig_{lname}"] = int(np.nansum(q_adj < 0.05))
        row[f"ora_retention_{lname}"] = float(np.mean(q_adj[rep_raw] < 0.05)) if len(rep_raw) else np.nan
        a_adj, z_adj = auc_vec(lay[lname], u_adj, return_z=True)
        top10_adj = np.argsort(-np.nan_to_num(z_adj, nan=-1e9))[:10]
        row[f"gsea_top10_overlap_{lname}"] = len(set(top10_raw) & set(top10_adj)) / 10
        row[f"gsea_top10_meanAUC_raw"] = float(np.nanmean(a_raw[top10_raw]))
        row[f"gsea_top10_meanAUC_{lname}"] = float(np.nanmean(a_adj[top10_raw]))
        if lname == PRIMARY:
            sig_adj = np.where(q_adj < 0.05)[0]
            for lst, ks, pv, qv in (("ORA_raw", rep_raw, p_raw, q_raw), ("ORA_B", sig_adj[np.argsort(p_adj[sig_adj])[:10]], p_adj, q_adj)):
                for rank_, k in enumerate(ks):
                    top_out.append(dict(seed=row["gene"], group=grp, list=lst, rank=rank_ + 1, set=names[k], p=pv[k], q=qv[k],
                                        q_raw=q_raw[k], q_adj=q_adj[k]))
            # ---- 3. specificity of the top deconfounded program ----
            kstar = int(np.nanargmax(z_adj))
            row["Tstar"] = names[kstar]; row["Tstar_AUC_seed"] = float(a_adj[kstar]); row["Tstar_z_seed"] = float(z_adj[kstar])
            row["program_present"] = bool(z_adj[kstar] >= PROGRAM_Z and a_adj[kstar] >= PROGRAM_AUC)
            members = np.asarray(M[kstar].todense()).ravel() > 0
            pool = np.where(members & (arm != arm[i]) & has_pos)[0]
            pool = pool[pool != i]
            peers = rng.choice(pool, size=min(N_PEERS, len(pool)), replace=False) if len(pool) else []
            pa = np.array([peer_auc(int(j), kstar, lname) for j in peers])
            for j, v in zip(peers, pa):
                peer_out.append(dict(seed=row["gene"], Tstar=names[kstar], peer=sym[j], peer_AUC=v))
            row["n_peers"] = len(pa); row["peer_median_AUC"] = float(np.nanmedian(pa)) if len(pa) else np.nan
            row["specificity_pct"] = float(100 * np.mean(pa < a_adj[kstar])) if len(pa) else np.nan
            row["seed_in_Tstar"] = bool(members[i])
            for lst, ks in (("GSEA_raw", top10_raw), ("GSEA_B", top10_adj)):
                for rank_, k in enumerate(ks):
                    top_out.append(dict(seed=row["gene"], group=grp, list=lst, rank=rank_ + 1, set=names[k], AUC_raw=a_raw[k], AUC_adj=a_adj[k],
                                        z_raw=z_raw[k], z_adj=z_adj[k]))
            if i == gi[LON]:
                lon_check = lay
                d = pd.DataFrame(dict(set=names, AUC_raw=a_raw, z_raw=z_raw, AUC_B=a_adj, z_B=z_adj, ORA_q_raw=q_raw, ORA_q_B=q_adj))
                d["rank_raw"] = d.z_raw.rank(ascending=False); d["rank_B"] = d.z_B.rank(ascending=False)
                d.sort_values("rank_B").to_csv(os.path.join(out, "LONP1_all_sets.csv"), index=False)
    rows_out.append(row)
    if n_done % 10 == 0:
        L(f"seed {n_done + 1}/{len(seeds)} done ({time.time() - t0:.0f} s)"); save_log()
T = pd.DataFrame(rows_out)

# ---------------- classification (prespecified thresholds, set before the run) ----------------
def classify(r):
    if not r.program_present: return "no robust program"
    if pd.isna(r.specificity_pct): return "unclassified"
    spec = r.specificity_pct >= SPEC_PCT and (r.Tstar_AUC_seed - r.peer_median_AUC) >= SPEC_DELTA
    return "seed-specific" if spec else "shared module"
T["class"] = T.apply(classify, axis=1)
T.to_csv(os.path.join(out, "seed_metrics.csv"), index=False)
pd.DataFrame(top_out).to_csv(os.path.join(out, "seed_top_sets.csv"), index=False)
pd.DataFrame(peer_out).to_csv(os.path.join(out, "seed_peers.csv"), index=False)

# ---------------- validation against the published LONP1 table ----------------
old_p = os.path.join(base, "results", "LONP1_partial_spearman_genomewide.csv")
if os.path.exists(old_p) and lon_check is not None:
    old = pd.read_csv(old_p).set_index("entrez")
    for lname, col in [("raw", "rho_raw"), ("A", "rho_partial_A")]:
        new = pd.Series(lon_check[lname], index=genes).reindex(old.index)
        m = new.notna() & old[col].notna()
        L(f"[validation] LONP1 {lname}: n={m.sum()} r={np.corrcoef(new[m], old.loc[m, col])[0, 1]:.5f}")

# ---------------- summary ----------------
def show(title, df):
    L("\n" + title); L(df.to_string())
show("1. POSITIONAL: class by group (dominant = >=50% of top-50 raw partners on own arm; enriched = binomial P<1e-3)",
     pd.crosstab(T.group, T.positional, margins=True))
show("   median same-arm fraction of top-50 partners: raw -> model A -> model B (chance = arm_gene_fraction)",
     T.groupby("group")[["arm_gene_fraction", "top50_same_arm_raw", "top50_same_arm_A", "top50_same_arm_B", "pos_fold"]].median().round(3))
show("2. RETENTION by group (medians). ORA = top-%d partners; GSEA = genome-wide ranking" % TOP_N,
     T.groupby("group")[["ora_n_sig_raw", "ora_n_sig_B", "ora_retention_A", "ora_retention_B", "gsea_top10_overlap_A", "gsea_top10_overlap_B"]].median().round(2))
show("   same, split by positional class",
     T.groupby("positional")[["ora_top_same_arm", "ora_n_sig_raw", "ora_n_sig_B", "ora_retention_B", "gsea_top10_overlap_B"]].median().round(2))
T["no_raw_ORA"] = T.ora_n_sig_raw == 0
show("   seeds whose raw top-%d partners give NO significant gene set (FDR<0.05), by positional class" % TOP_N,
     pd.crosstab(T.positional, T.no_raw_ORA, margins=True))
show("3. SPECIFICITY: class by group", pd.crosstab(T.group, T["class"], margins=True))
show("   medians by group", T.groupby("group")[["Tstar_z_seed", "Tstar_AUC_seed", "peer_median_AUC", "specificity_pct"]].median().round(3))
L(f"   chance rate of 'seed-specific' if a seed were exchangeable with its peers: about {100 - SPEC_PCT}% "
  f"(percentile >= {SPEC_PCT} among {N_PEERS} peers), before the AUC-difference requirement")
lit_, rnd_ = T[T.group == "literature"], T[T.group == "random"]
if len(lit_) and len(rnd_):
    a, b = (lit_["class"] == "seed-specific").sum(), (rnd_["class"] == "seed-specific").sum()
    L(f"   literature vs random, seed-specific: {a}/{len(lit_)} vs {b}/{len(rnd_)}; Fisher P={fisher_exact([[a, len(lit_) - a], [b, len(rnd_) - b]])[1]:.3f}")
    for c in ("specificity_pct", "Tstar_z_seed", "top50_same_arm_raw", "ora_retention_B"):
        x, y = lit_[c].dropna(), rnd_[c].dropna()
        if len(x) and len(y): L(f"   literature vs random, {c}: medians {x.median():.3f} vs {y.median():.3f}; Mann-Whitney P={mannwhitneyu(x, y).pvalue:.3f}")
show("most frequent top deconfounded program (T*) across all seeds", T.Tstar.value_counts().head(12))
show("literature seeds + LONP1", T[T.group.isin(["literature", "LONP1"])][["gene", "cytoband", "top50_same_arm_raw", "top50_same_arm_B", "positional",
     "ora_n_sig_raw", "ora_retention_B", "Tstar", "Tstar_AUC_seed", "peer_median_AUC", "specificity_pct", "class"]].round(3).set_index("gene"))
L(f"\ntotal time {time.time() - t0:.0f} s")
save_log()
print("DONE")
