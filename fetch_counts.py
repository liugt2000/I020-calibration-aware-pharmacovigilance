"""
Fetch 2x2 cell counts from openFDA drug/event API for I020 PV pipeline.
For each (drug, outcome-PT-set) we need:
  a = reports with drug AND any PT in set
  a+b = reports with drug
  a+c = reports with any PT in set
  N = total reports (constant)
Caching to JSON to avoid repeated calls. Rate-limited.
"""
import json, time, os, hashlib
from pathlib import Path
import requests

HERE = Path(__file__).parent
CFG = json.loads((HERE / "config.json").read_text(encoding="utf-8"))
CACHE = HERE / "counts_cache.json"
BASE = "https://api.fda.gov/drug/event.json"

sess = requests.Session()
sess.headers.update({"User-Agent": "I020-PV-research/0.1 (academic pharmacovigilance analysis)"})

cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}


def key(*parts):
    return hashlib.md5("||".join(map(str, parts)).encode()).hexdigest()[:12]


def api_total(search, tries=4):
    k = key("total", search)
    if k in cache:
        return cache[k]
    for attempt in range(tries):
        try:
            r = sess.get(BASE, params={"search": search}, timeout=60)
            if r.status_code == 404:
                cache[k] = 0
                return 0
            d = r.json()
            if "error" in d:
                if "No match" in str(d["error"].get("message", "")):
                    cache[k] = 0
                    return 0
                time.sleep(2 * (attempt + 1))
                continue
            total = d["meta"]["results"]["total"]
            cache[k] = total
            return total
        except Exception:
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"failed: {search}")


G = CFG["drug_field"]
R = CFG["reaction_field"]


def drug_q(drug):
    return f'{G}:"{CFG["drugs"][drug]}"'


def pt_q(pt):
    return f'{R}:"{pt}"'


def or_q(terms):
    return "(" + " OR ".join(terms) + ")"


def build_targets():
    """List of (label, drug_key, outcome_key, pts) tasks."""
    tasks = []
    for name, spec in CFG["index_outcomes"].items():
        tasks.append((f"IDX::{name}", spec["drug"], name, spec["pts"]))
    for drug in CFG["index_drugs_for_calibration"]:
        for panel, pts in [("N1", CFG["negative_controls_N1"]),
                           ("N2", CFG["negative_controls_N2"]),
                           ("PC", CFG["positive_controls"])]:
            for pt in pts:
                tasks.append((f"CAL::{drug}::{panel}::{pt}", drug, f"{panel}::{pt}", [pt]))
    return tasks


def main():
    drug_totals = {d: api_total(drug_q(d)) for d in CFG["drugs"]}
    print("drug totals:", json.dumps(drug_totals, ensure_ascii=False), flush=True)

    tasks = build_targets()
    print(f"total tasks: {len(tasks)}", flush=True)
    results = []
    pt_total_cache = {}
    for i, (label, drug, outcome_key, pts) in enumerate(tasks):
        pq = or_q([pt_q(p) for p in pts])
        a = api_total(f"{drug_q(drug)} AND {pq}")
        if outcome_key not in pt_total_cache:
            pt_total_cache[outcome_key] = api_total(pq)
        ac = pt_total_cache[outcome_key]
        results.append({
            "label": label,
            "drug": drug,
            "outcome_key": outcome_key,
            "pts": pts,
            "a": a,
            "drug_total": drug_totals[drug],
            "event_total": ac,
            "N": CFG["database_total"],
        })
        if (i + 1) % 20 == 0:
            print(f"  fetched {i+1}/{len(tasks)}", flush=True)
            CACHE.write_text(json.dumps(cache), encoding="utf-8")

    CACHE.write_text(json.dumps(cache), encoding="utf-8")
    out = HERE / "counts.csv"
    import csv
    with out.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["label", "drug", "outcome_key", "pts", "a", "drug_total", "event_total", "N"])
        w.writeheader()
        for r in results:
            rr = dict(r); rr["pts"] = "|".join(rr["pts"])
            w.writerow(rr)
    print(f"wrote {out} ({len(results)} rows)", flush=True)
    print(f"cache size: {len(cache)} queries", flush=True)


if __name__ == "__main__":
    main()
