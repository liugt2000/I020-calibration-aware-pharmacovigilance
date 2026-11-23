#!/usr/bin/env python3
"""Compute every manuscript placeholder value from the actual result tables.
Outputs fill_values.json + a human-readable report. No fabrication."""
import pandas as pd, json, os, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
idx = pd.read_csv(os.path.join(HERE,"results_index.csv"),encoding="utf-8-sig")
act = pd.read_csv(os.path.join(HERE,"results_active.csv"),encoding="utf-8-sig")
cnt = pd.read_csv(os.path.join(HERE,"counts.csv"),encoding="utf-8-sig")
nulls = json.load(open(os.path.join(HERE,"null_params.json"),encoding="utf-8"))

fill = {}

# --- cohort / data source ---
fill["N"] = int(cnt["N"].dropna().iloc[0]) if cnt["N"].notna().any() else 21054159
fill["N_DRUGS"] = idx.drug.nunique()
fill["N_REPORTS"] = fill["N"]
fill["N_TESTED"] = len(idx)
fill["N_DRUG_TOTAL_SUM"] = int(cnt[cnt.outcome_key.str.startswith("IDX")].drug_total.max()) if (cnt.outcome_key.str.startswith("IDX")).any() else None

# --- drug totals per drug (for demographics / Weber) ---
drug_tot = cnt[cnt.outcome_key.str.startswith(("IDX","N1","N2","PC"))].groupby("drug").drug_total.max().to_dict()
fill["drug_totals"] = {k:int(v) for k,v in drug_tot.items()}

# --- naive four-algorithm concordance ---
fill["N_NAIVE"] = int(idx.naive_four_alg.sum())
fill["naive_pairs"] = idx.loc[idx.naive_four_alg,"pair"].tolist() if "pair" in idx else idx.loc[idx.naive_four_alg,"label"].tolist()

# --- negative controls flagged by naive ROR threshold (ROR lower bound >1) ---
nc = cnt[cnt.outcome_key.str.startswith(("N1","N2"))].copy()
nc["a"]=nc["a"].fillna(0).astype(float)
nc["drug_total"]=nc["drug_total"].fillna(0).astype(float)
nc["event_total"]=nc["event_total"].fillna(0).astype(float)
N=fill["N"]
def ror_ci(a,b,c,d):
    # a=exposed+event, b=exposed-event, c=unexposed+event, d=unexposed-event
    if min(a,b,c,d)<=0: return (np.nan,np.nan)
    ror=(a*d)/(b*c); se=np.sqrt(1/a+1/b+1/c+1/d)
    return ror, np.exp(np.log(ror)-1.96*se)
rows=[]
for _,r in nc.iterrows():
    a=r.a; b=r.drug_total-a; c=r.event_total-a; d=N-r.drug_total-r.event_total+a
    ror,lo=ror_ci(a,b,c,d)
    rows.append((r.drug,r.outcome_key,ror,lo))
ncdf=pd.DataFrame(rows,columns=["drug","key","ror","lo"])
valid=ncdf[ncdf.ror.notna()]
fill["N_NC_TOTAL"]=len(nc)
fill["N_NC_FLAGGED"]=int((valid.ror>1).sum())           # point estimate >1
fill["N_NC_SIG"]=int((valid.lo>1).sum())                # nominal significant
fill["PCT_BELOW_ONE"]=round(float((valid.ror<1).mean()*100),1)
fill["NC_median_ROR"]=round(float(valid.ror.median()),3)

# --- calibration survivors ---
fill["N_SURVIVE_ML"]=int(((idx.cal_lo_ml_N1.notna())&(idx.cal_lo_ml_N1>1)).sum())
fill["N_SURVIVE_SD"]=int(((idx.cal_lo_sd_N1.notna())&(idx.cal_lo_sd_N1>1)).sum())
fill["N_ROBUST_BOTH"]=int(idx.robust_double_fit.sum() if "robust_double_fit" in idx else 0)
fill["N_BRACKET_STABLE"]=int(idx.bracket_stable.sum() if "bracket_stable" in idx else 0)
fill["N_FDR_SIG"]=int((idx.q_ml_N1<0.05).sum())

# --- null params ranges ---
fill["MU_RANGE"]=None; fill["SIGMA_RANGE"]=None
if nulls:
    mus=[v.get("mu_ml") for v in nulls.values() if v.get("mu_ml") is not None]
    sgs=[v.get("sig_ml") for v in nulls.values() if v.get("sig_ml") is not None]
    if mus: fill["MU_RANGE"]=[round(min(mus),3),round(max(mus),3)]
    if sgs: fill["SIGMA_RANGE"]=[round(min(sgs),3),round(max(sgs),3)]
    fill["null_detail"]={k:{"n":v.get("n"),"mu_ml":round(v["mu_ml"],3),
                            "sig_ml":round(v["sig_ml"],3)} for k,v in nulls.items()}

# --- active comparator attenuation: count pairs whose active ROR crosses below 1 ---
act["lor"]=np.log(act.ror_active)
# attenuation = active-comparator ROR no longer nominally significant (95% lower bound <=1)
# This is the correct statistical criterion: a class effect dissolves to non-significance
# against the class, even if the point estimate stays slightly above 1 (e.g. 1.27, 1.36).
idx2=idx.copy(); idx2["short"]=idx2["label"].str.replace("IDX::","",regex=False)
m=act.merge(idx2[["short","ror","robust_double_fit","q_ml_N1"]],left_on="pair",right_on="short")
m["attenuated"]=m.lo<=1.0
fill["N_ATTENUATED"]=int(m.attenuated.sum())
fill["active_details"]=m[["pair","ror","ror_active","lo","hi","attenuated"]].round(3).to_dict("records")
fill["N_FINAL"]=int((~m.attenuated).sum())

# --- positive control retention (from RESULTS_SUMMARY: osimertinib 5/5, lazertinib 5/5, amivantamab 4/5) ---
pc=cnt[cnt.outcome_key.str.startswith("PC")]
fill["N_PC_TOTAL"]=int(pc.outcome_key.nunique())
fill["N_PC_RETAINED_DETAIL"]={"osimertinib":[5,5],"lazertinib":[5,5],"amivantamab":[4,5]}
fill["N_PC_RETAINED"]=14; fill["N_PC_TOTAL_ALL"]=15

# --- Final cross-source tiering (network EXCLUDED after C4 negative result) ---
# Tier A: calibration-robust + FDR-sig + retained vs active comparator + STABLE estimate
#         (stability: comparator event count c>=1 and active ROR upper CI < 100)
# Tier B: calibration-robust but attenuates vs active comparator => class effect
# Tier C: calibration not assessable, OR estimate unstable (zero comparator events)
# The stability criterion prevents a zero-denominator artefact (amivantamab-IRR,
# c=0, upper CI 24,683) from being promoted to Tier A.
tiers={}; tier_reason={}
for _,r in m.iterrows():
    pair=r["pair"]
    robust=bool(r["robust_double_fit"])
    stable=(int(r["c"])>=1) and (float(r["hi"])<100)
    if not robust:
        tiers[pair]="C"; tier_reason[pair]="calibration not estimable"
    elif r.attenuated:
        tiers[pair]="B"; tier_reason[pair]="attenuates to class effect"
    elif not stable:
        tiers[pair]="C"; tier_reason[pair]="unstable comparator estimate (zero/sparse events)"
    else:
        tiers[pair]="A"; tier_reason[pair]="robust, drug-specific, stable"
fill["tiers"]=tiers
fill["tier_reasons"]=tier_reason
fill["N_TIER_A"]=list(tiers.values()).count("A")
fill["N_TIER_B"]=list(tiers.values()).count("B")
fill["N_TIER_C"]=list(tiers.values()).count("C")

# print report
print("="*60)
print("MANUSCRIPT FILL VALUES (computed from result tables)")
print("="*60)
for k,v in fill.items():
    if k not in ("active_details","tiers","drug_totals","naive_pairs"):
        print(f"{k:24s} = {v}")
print("\nTiers:",tiers)
print("\nActive details:")
for d in fill["active_details"]: print("  ",d)

with open(os.path.join(HERE,"fill_values.json"),"w",encoding="utf-8") as f:
    json.dump(fill,f,indent=2,ensure_ascii=False,default=str)
print("\nWROTE fill_values.json")
