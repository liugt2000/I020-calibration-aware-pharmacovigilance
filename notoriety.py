#!/usr/bin/env python3
"""I020 C5 — Notoriety / temporal stability analysis.

For each index pair we build a QUARTERLY series of the drug-event reporting count
as a fraction of the drug's total reports in that quarter, minus the concurrent
fraction for that drug's pooled negative-control outcomes (N1). A notoriety pattern
is a step increase in the index fraction that is NOT mirrored in the control fraction.

Regulatory anchor dates (verified, from FDA labels / DailyMed archives):
  osimertinib: initial approval 2015-Q4; 1st-line indication 2018-Q2 (s018 label 2020-Q3)
  tepotinib:   accelerated approval 2021-Q2 (MET exon14 NSCLC)
  lazertinib:  not FDA-approved as single agent (IIT context) -> no regulatory anchor
  amivantamab: accelerated approval 2021-Q2 (EGFR exon20 insertion NSCLC)

Because ILD/QT warnings were present in the initial osimertinib label, there is no
single "sentinel communication" for these pairs; we therefore use data-driven
change-point detection AND report alignment with the indication-expansion quarters.

Data: openFDA receivedate counts (public). No fabrication of dates or counts.
"""
import requests, json, time, os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
S = requests.Session(); S.headers.update({"User-Agent":"I020-pv-notoriety/0.1"})
B = "https://api.fda.gov/drug/event.json"
G = "patient.drug.openfda.generic_name"
R = "patient.reaction.reactionmeddrapt"

DRUGS = {
 "osimertinib":"OSIMERTINIB","tepotinib":"TEPOTINIB",
 "lazertinib":"LAZERTINIB","amivantamab":"AMIVANTAMAB",
}
EVENTS = {
 "osimertinib":{
   "ILD":["Interstitial lung disease","Pneumonitis","Interstitial pneumonitis","Pulmonary fibrosis","Acute interstitial pneumonia"],
   "QT":["Electrocardiogram QT prolonged","Long QT syndrome"]},
 "tepotinib":{"AKI":["Acute kidney injury","Renal impairment","Blood creatinine increased","Renal failure acute"]},
 "lazertinib":{"ILD":["Interstitial lung disease","Pneumonitis","Interstitial pneumonitis","Pulmonary fibrosis","Acute interstitial pneumonia"]},
 "amivantamab":{
   "IRR":["Infusion related reaction","Infusion related reactions"],
   "ILD":["Interstitial lung disease","Pneumonitis","Interstitial pneumonitis","Pulmonary fibrosis","Acute interstitial pneumonia"]},
}
# N1 inert controls per drug (subset used in calibration) — concurrent reference
N1 = ["Cataract","Osteoarthritis","Carpal tunnel syndrome","Inguinal hernia",
      "Diverticulitis","Cholelithiasis","Benign prostatic hyperplasia"]

ANCHORS = {  # regulatory anchor quarter (for interpretation only, not imposed)
 "osimertinib":"2018Q2","tepotinib":"2021Q2","lazertinib":None,"amivantamab":"2021Q2",
}

def qkey(y,q): return f"{y}Q{q}"
def date_counts(search):
    """return dict date_str(YYYYMMDD)->count by counting receivedate"""
    for attempt in range(3):
        try:
            r=S.get(B,params={"search":search,"count":"receivedate"},timeout=90)
            d=r.json()
            if "error" in d:
                if "No matches" in str(d["error"].get("message","")): return {}
                return None
            return {x["time"]:int(x["count"]) for x in d.get("results",[])}
        except Exception:
            time.sleep(1.5)
    return None

def to_quarter(dc):
    out={}
    for dt,n in dc.items():
        y=int(dt[:4]); m=int(dt[4:6]); q=(m-1)//3+1
        k=qkey(y,q); out[k]=out.get(k,0)+n
    return out

def drug_search(drug):
    return f'{G}:"{drug}"'

def event_search(drug, pts):
    terms=" OR ".join(f'{R}:"{p}"' for p in pts)
    return f'{drug_search(drug)} AND ({terms})'

def main():
    series={}
    for dname,dterm in DRUGS.items():
        print(f"\n=== {dname} ===",flush=True)
        # drug total by date (this is a large query; count receivedate over all drug reports)
        total_dc=date_counts(drug_search(dterm))
        if total_dc is None: print("  total query failed",flush=True); continue
        total_q=to_quarter(total_dc)
        print(f"  total quarters={len(total_q)} span={min(total_q)}..{max(total_q)}",flush=True)
        # control events
        ctrl_dc=date_counts(event_search(dterm,N1))
        ctrl_q=to_quarter(ctrl_dc or {})
        # each outcome
        for outcome,pts in EVENTS[dname].items():
            idx_dc=date_counts(event_search(dterm,pts))
            idx_q=to_quarter(idx_dc or {})
            quarters=sorted(set(total_q)|set(idx_q)|set(ctrl_q))
            rows=[]
            for k in quarters:
                t=total_q.get(k,0)
                if t<20: continue   # require minimum drug exposure per quarter
                a=idx_q.get(k,0); c=ctrl_q.get(k,0)
                rows.append({"quarter":k,"drug_total":t,"idx_n":a,"ctrl_n":c,
                             "idx_frac":a/t,"ctrl_frac":c/t,"diff_frac":(a-c)/t})
            series[f"{dname}_{outcome}"]={"anchor":ANCHORS[dname],"rows":rows}
            print(f"  {outcome}: {len(rows)} usable quarters, "
                  f"idx total={sum(r['idx_n'] for r in rows)}, ctrl total={sum(r['ctrl_n'] for r in rows)}",flush=True)
            time.sleep(0.5)
        time.sleep(0.5)

    # change-point: largest |mean(diff after)-mean(diff before)| over all split points
    cps={}
    for name,obj in series.items():
        rows=obj["rows"]
        vals=np.array([r["diff_frac"] for r in rows])
        if len(vals)<8: 
            cps[name]={"change_quarter":None,"delta":None,"anchor":obj["anchor"]}; continue
        best=None
        for i in range(4,len(vals)-4):
            d=vals[i:].mean()-vals[:i].mean()
            if best is None or abs(d)>abs(best[1]): best=(rows[i]["quarter"],d,i)
        cps[name]={"change_quarter":best[0],"delta":round(float(best[1]),5),
                   "split_index":best[2],"anchor":obj["anchor"],
                   "max_diff_frac":round(float(vals.max()),5),
                   "min_diff_frac":round(float(vals.min()),5)}
    print("\n=== CHANGE POINTS ===",flush=True)
    for name,c in cps.items():
        print(f"  {name:22s} cp={c['change_quarter']} delta={c['delta']} anchor={c['anchor']}",flush=True)

    with open(os.path.join(HERE,"notoriety_results.json"),"w",encoding="utf-8") as f:
        json.dump({"meta":{"generated":time.strftime("%Y-%m-%d %H:%M:%S"),
                           "method":"quarterly index-vs-control reporting fraction; max change-point",
                           "anchors":ANCHORS,"min_drug_reports_per_quarter":20},
                   "series":series,"change_points":cps},f,indent=2)
    print("\nWROTE notoriety_results.json",flush=True)

if __name__=="__main__":
    main()
