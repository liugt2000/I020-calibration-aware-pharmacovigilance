"""
Fetch active-comparator 2x2 counts.
For each index outcome (drug d, event E) and its comparator pool C:
  a = reports with d AND E
  b = reports with d AND NOT E   = drug_total(d) - a
  c = reports with ANY drug in C AND E
  d = reports with ANY drug in C AND NOT E = comparator_pool_total - c
Then ROR_active = (a/b)/(c/d), which holds indication roughly constant.
"""
import json, hashlib, time
from pathlib import Path
import requests

HERE = Path(__file__).parent
CFG = json.loads((HERE / "config.json").read_text(encoding="utf-8"))
CACHE = HERE / "counts_cache.json"
cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}
BASE = "https://api.fda.gov/drug/event.json"
sess = requests.Session()
sess.headers.update({"User-Agent": "I020-PV-research/0.1"})


def api_total(search, tries=4):
    k = hashlib.md5(("total||" + search).encode()).hexdigest()[:12]
    if k in cache:
        return cache[k]
    for i in range(tries):
        try:
            r = sess.get(BASE, params={"search": search}, timeout=60)
            if r.status_code == 404:
                cache[k] = 0; return 0
            d = r.json()
            if "error" in d:
                time.sleep(2*(i+1)); continue
            v = d["meta"]["results"]["total"]
            cache[k] = v; return v
        except Exception:
            time.sleep(2*(i+1))
    raise RuntimeError(search)


def drug_q(d): return f'{CFG["drug_field"]}:"{CFG["drugs"][d]}"'
def pt_q(p): return f'{CFG["reaction_field"]}:"{p}"'
def or_q(ts): return "(" + " OR ".join(ts) + ")"


rows = []
for name, spec in CFG["index_outcomes"].items():
    d = spec["drug"]
    comps = [c for c in CFG["active_comparators"].get(d, []) if c in CFG["drugs"]]
    if not comps:
        continue
    E = or_q([pt_q(p) for p in spec["pts"]])
    C = or_q([drug_q(c) for c in comps])
    a = api_total(f"{drug_q(d)} AND {E}")
    dtot = api_total(drug_q(d))
    c = api_total(f"{C} AND {E}")
    ctot = api_total(C)
    rows.append(dict(pair=name, index=d, comparators="|".join(comps),
                     a=a, index_total=dtot, c=c, comparator_total=ctot))
    print(f"{name:24s} a={a:5d} c={c:5d} cmp_n={ctot}", flush=True)
    CACHE.write_text(json.dumps(cache), encoding="utf-8")

import csv
out = HERE / "counts_active.csv"
with out.open("w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=["pair", "index", "comparators", "a", "index_total", "c", "comparator_total"])
    w.writeheader(); w.writerows(rows)
print("wrote", out, flush=True)
