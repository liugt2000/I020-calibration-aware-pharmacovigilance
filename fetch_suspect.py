"""
Suspect-only recomputation of the I020 PV pipeline.

Single-variable sensitivity analysis: exposure is redefined from
    generic_name:"DRUG"  (any reporting role)
to
    generic_name:"DRUG" AND drugcharacterization:"1"  (suspect role only)
Everything else (event PT sets, controls, comparators, N) is unchanged,
so any difference isolates the effect of the role restriction.

Outputs (parallel to the primary pipeline, original files untouched):
    counts_suspect.csv
    counts_active_suspect.csv
    counts_cache_suspect.json
"""
import json, time, hashlib, csv
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests

HERE = Path(__file__).parent
CFG = json.loads((HERE / "config.json").read_text(encoding="utf-8"))
CACHE = HERE / "counts_cache_suspect.json"
cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}
BASE = "https://api.fda.gov/drug/event.json"
SUSPECT = 'patient.drug.drugcharacterization:"1"'

G = CFG["drug_field"]
R = CFG["reaction_field"]


def query(search, tries=5):
    k = hashlib.md5(("suspect||" + search).encode()).hexdigest()[:12]
    if k in cache:
        return cache[k]
    s = requests.Session()
    s.headers.update({"User-Agent": "I020-PV-suspect/0.1 (academic)"})
    for i in range(tries):
        try:
            r = s.get(BASE, params={"search": search}, timeout=60)
            if r.status_code == 404:
                cache[k] = 0; return 0
            d = r.json()
            if "error" in d:
                time.sleep(1.5 * (i + 1)); continue
            v = d["meta"]["results"]["total"]
            cache[k] = v; return v
        except Exception:
            time.sleep(1.5 * (i + 1))
    raise RuntimeError(f"failed: {search}")


def drug_q(d):
    return f'({G}:"{CFG["drugs"][d]}" AND {SUSPECT})'


def pt_q(p):
    return f'{R}:"{p}"'


def or_q(terms):
    return "(" + " OR ".join(terms) + ")"


# ---------- build all unique query strings ----------
queries = {}   # key -> search string


def register(key, search):
    queries[key] = search


# drug totals (suspect)
for d in CFG["drugs"]:
    register(f"dtot::{d}", drug_q(d))
# index pairs
for name, spec in CFG["index_outcomes"].items():
    pq = or_q([pt_q(p) for p in spec["pts"]])
    register(f"idx_a::{name}", f'{drug_q(spec["drug"])} AND {pq}')
    register(f"et::{name}", pq)  # event total stays all-role (single-variable change)
# calibration panels
for d in CFG["index_drugs_for_calibration"]:
    for panel, pts in [("N1", CFG["negative_controls_N1"]),
                       ("N2", CFG["negative_controls_N2"]),
                       ("PC", CFG["positive_controls"])]:
        for pt in pts:
            ck = f"CAL::{d}::{panel}::{pt}"
            register(f"a::{ck}", f'{drug_q(d)} AND {pt_q(pt)}')
            register(f"et::{ck}", pt_q(pt))
# active comparators
for name, spec in CFG["index_outcomes"].items():
    d = spec["drug"]
    comps = [c for c in CFG["active_comparators"].get(d, []) if c in CFG["drugs"]]
    if not comps:
        continue
    E = or_q([pt_q(p) for p in spec["pts"]])
    C = or_q([drug_q(c) for c in comps])
    register(f"act_a::{name}", f'{drug_q(d)} AND {E}')
    register(f"act_c::{name}", f'{C} AND {E}')
    register(f"act_ctot::{name}", C)

print(f"suspect queries to fetch: {len(queries)} (cached: {len([k for k in queries if hashlib.md5(('suspect||'+queries[k]).encode()).hexdigest()[:12] in cache])})", flush=True)

# ---------- fetch in parallel ----------
new = 0
with ThreadPoolExecutor(max_workers=4) as ex:
    futs = {ex.submit(query, q): k for k, q in queries.items()}
    for i, fut in enumerate(as_completed(futs), 1):
        fut.result(); new += 1
        if i % 25 == 0:
            CACHE.write_text(json.dumps(cache), encoding="utf-8")
            print(f"  {i}/{len(queries)} done", flush=True)
CACHE.write_text(json.dumps(cache), encoding="utf-8")

# ---------- assemble counts_suspect.csv ----------
rows = []
for name, spec in CFG["index_outcomes"].items():
    rows.append(dict(label=f"IDX::{name}", drug=spec["drug"], outcome_key=name,
                     pts="|".join(spec["pts"]),
                     a=cache[hashlib.md5(("suspect||"+queries[f"idx_a::{name}"]).encode()).hexdigest()[:12]],
                     drug_total=cache[hashlib.md5(("suspect||"+queries[f"dtot::{spec['drug']}"]).encode()).hexdigest()[:12]],
                     event_total=cache[hashlib.md5(("suspect||"+queries[f"et::{name}"]).encode()).hexdigest()[:12]],
                     N=CFG["database_total"]))
for d in CFG["index_drugs_for_calibration"]:
    for panel, pts in [("N1", CFG["negative_controls_N1"]),
                       ("N2", CFG["negative_controls_N2"]),
                       ("PC", CFG["positive_controls"])]:
        for pt in pts:
            ck = f"CAL::{d}::{panel}::{pt}"
            rows.append(dict(label=ck, drug=d, outcome_key=f"{panel}::{pt}", pts=pt,
                             a=cache[hashlib.md5(("suspect||"+queries[f"a::{ck}"]).encode()).hexdigest()[:12]],
                             drug_total=cache[hashlib.md5(("suspect||"+queries[f"dtot::{d}"]).encode()).hexdigest()[:12]],
                             event_total=cache[hashlib.md5(("suspect||"+queries[f"et::{ck}"]).encode()).hexdigest()[:12]],
                             N=CFG["database_total"]))
with (HERE / "counts_suspect.csv").open("w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=["label", "drug", "outcome_key", "pts", "a", "drug_total", "event_total", "N"])
    w.writeheader(); w.writerows(rows)
print(f"wrote counts_suspect.csv ({len(rows)} rows)", flush=True)

# ---------- assemble counts_active_suspect.csv ----------
act = []
for name, spec in CFG["index_outcomes"].items():
    d = spec["drug"]
    comps = [c for c in CFG["active_comparators"].get(d, []) if c in CFG["drugs"]]
    if not comps:
        continue
    g = lambda k: cache[hashlib.md5(("suspect||"+queries[k]).encode()).hexdigest()[:12]]
    act.append(dict(pair=name, index=d, comparators="|".join(comps),
                    a=g(f"act_a::{name}"), index_total=g(f"dtot::{d}"),
                    c=g(f"act_c::{name}"), comparator_total=g(f"act_ctot::{name}")))
with (HERE / "counts_active_suspect.csv").open("w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=["pair", "index", "comparators", "a", "index_total", "c", "comparator_total"])
    w.writeheader(); w.writerows(act)
print(f"wrote counts_active_suspect.csv ({len(act)} rows)", flush=True)
print(f"cache: {len(cache)} queries", flush=True)
