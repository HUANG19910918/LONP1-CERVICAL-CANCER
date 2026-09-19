# 18_verify_references.py —— 用 PubMed E-utilities 逐条核对参考文献（作者/刊名/年卷期页），输出核对表
# 用法: python3 18_verify_references.py <refs.txt(每行"n. 文献")> <out.csv>
import sys, re, json, time, urllib.request, urllib.parse
src, out = sys.argv[1], sys.argv[2]
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(u):
    for k in range(4):
        try: return json.load(urllib.request.urlopen(u, timeout=30))
        except Exception as e: time.sleep(1.5)
    return None
def parse(ref):
    # "Authors. Title. Journal. Year;Vol(Issue):Pages."
    m = re.match(r"^(?P<auth>[^.]+(?:\.[^.]*?)?)\. (?P<title>.+?)\. (?P<jour>[^.]+?)\. (?P<year>\d{4});(?P<vol>[^(:]+)?(?:\((?P<iss>[^)]+)\))?:?(?P<pages>[^.]+)?\.?$", ref)
    return m.groupdict() if m else None
rows = []
for line in open(src, encoding="utf-8"):
    line = line.strip()
    if not line: continue
    n, ref = re.match(r"^\s*(\d+)\.\s*(.*)$", line).groups()
    p = parse(ref)
    title = p["title"] if p else ref[:80]
    fa = ref.split(",")[0].split(" ")[0]
    ids = get(E + "esearch.fcgi?db=pubmed&retmode=json&term=" + urllib.parse.quote(f'"{title}"[Title]'))
    idl = ids["esearchresult"]["idlist"] if ids else []
    if not idl:  # 退化：首作者+年份+标题关键词
        kw = " ".join(re.findall(r"[A-Za-z]{5,}", title)[:5])
        ids = get(E + "esearch.fcgi?db=pubmed&retmode=json&term=" + urllib.parse.quote(f"{fa}[Author - First] AND {p['year'] if p else ''}[DP] AND ({kw})"))
        idl = ids["esearchresult"]["idlist"] if ids else []
    rec = {"n": int(n), "original": ref, "pmid": "", "status": "NOT_FOUND", "pubmed_authors": "", "pubmed_title": "", "pubmed_source": "", "issues": ""}
    if idl:
        s = get(E + "esummary.fcgi?db=pubmed&retmode=json&id=" + ",".join(idl[:5]))
        best = None
        for pid in idl[:5]:
            d = s["result"][pid]
            if best is None or (p and p["year"] in d.get("pubdate", "")): best = d
        d = best
        auth = d.get("authors", []); al = [a["name"] for a in auth]
        pm_auth = ", ".join(al[:6]) + (", et al." if len(al) > 6 else ".")
        vol, iss, pages = d.get("volume", ""), d.get("issue", ""), d.get("pages", "")
        year = d.get("pubdate", "")[:4]
        src_str = f"{d.get('source','')}. {year};{vol}" + (f"({iss})" if iss else "") + (f":{pages}" if pages else "")
        issues = []
        if p:
            if p["year"] != year: issues.append(f"year {p['year']}->{year}")
            if (p["vol"] or "").strip() != vol: issues.append(f"vol {p['vol']}->{vol}")
            if (p["iss"] or "") != (iss or ""): issues.append(f"issue {p['iss']}->{iss}")
            op = (p["pages"] or "").strip()
            if op and pages and op.replace("–", "-") != pages: issues.append(f"pages {op}->{pages}")
            oj = p["jour"].replace(".", "").strip().lower(); nj = d.get("source", "").replace(".", "").strip().lower()
            if oj != nj: issues.append(f"journal '{p['jour']}'->'{d.get('source','')}'")
            ofa = ref.split(" ")[0].strip(",")
            if al and al[0].split(" ")[0].lower() != ofa.lower(): issues.append(f"first author {ofa}->{al[0]}")
            ot = re.sub(r"[^a-z0-9]", "", p["title"].lower()); nt = re.sub(r"[^a-z0-9]", "", d.get("title", "").lower())
            if ot != nt: issues.append("title differs")
        rec.update(pmid=best["uid"], status="OK" if not issues else "CHECK", pubmed_authors=pm_auth,
                   pubmed_title=d.get("title", ""), pubmed_source=src_str, issues="; ".join(issues))
    rows.append(rec); print(rec["n"], rec["status"], rec["pmid"], rec["issues"], flush=True)
    time.sleep(0.4)
import csv
with open(out, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("DONE", sum(r["status"] == "OK" for r in rows), "OK /", len(rows))
