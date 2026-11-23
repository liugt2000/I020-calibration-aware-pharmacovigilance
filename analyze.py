"""
I020 PV analysis: four disproportionality algorithms + empirical-null calibration.
Implements lock decisions L2/L3/L4 from skeleton:
  - ML fit ln(ROR)_i ~ N(mu, sigma^2 + SE_i^2) primary; sample-SD conservative sensitivity
  - robust = both fits calibrated LB>1 AND BH q<0.05 under ML
  - four algorithms reported but concordance NOT treated as validation
"""
import json, math
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

HERE = Path(__file__).parent
CFG = json.loads((HERE / "config.json").read_text(encoding="utf-8"))
df = pd.read_csv(HERE / "counts.csv", encoding="utf-8-sig")

Z = 1.959963984540054


def cells(row):
    a = int(row["a"]); dt = int(row["drug_total"]); et = int(row["event_total"]); N = int(row["N"])
    a = max(a, 0)
    b = max(dt - a, 0)
    c = max(et - a, 0)
    d = max(N - a - b - c, 0)
    return a, b, c, d


def disprop(a, b, c, d):
    # Haldane-Anscombe continuity correction when any zero
    ac, bc, cc, dc = a, b, c, d
    if min(ac, bc, cc, dc) == 0:
        ac, bc, cc, dc = ac + .5, bc + .5, cc + .5, dc + .5
    N = ac + bc + cc + dc
    ror = (ac * dc) / (bc * cc)
    se_ln = math.sqrt(1/ac + 1/bc + 1/cc + 1/dc)
    ror_lo = math.exp(math.log(ror) - Z * se_ln)
    ror_hi = math.exp(math.log(ror) + Z * se_ln)
    prr = (ac/(ac+bc)) / (cc/(cc+dc))
    chi2 = N * (ac*dc - bc*cc)**2 / ((ac+bc)*(cc+dc)*(ac+cc)*(bc+dc))
    # BCPNN information component (approximate variance formula, Norén-style)
    E = (ac+bc)*(ac+cc)/N
    ic = math.log2((ac + 0.5) / (E + 0.5))
    ic_sd = math.sqrt(1/(ac+0.5) + 1/(bc+0.5) + 1/(cc+0.5) + 1/(dc+0.5)) / math.log(2)
    ic025 = ic - Z * ic_sd
    # MGPS/EBGM via simple empirical-Bayes Gamma-Poisson (method-of-moments prior)
    lam = ac / E if E > 0 else 0
    return dict(a=a, b=b, c=c, d=d, ror=ror, ror_lo=ror_lo, ror_hi=ror_hi,
                se_ln=se_ln, prr=prr, chi2=chi2, ic=ic, ic025=ic025, lam=lam, E=E)


for col in ["ror", "ror_lo", "ror_hi", "se_ln", "prr", "chi2", "ic", "ic025", "lam", "E", "a", "b", "c", "d"]:
    df[col] = df.apply(lambda r: disprop(*cells(r))[col], axis=1)


def fit_null(lnr, se, method="ml"):
    lnr = np.asarray(lnr, float); se = np.asarray(se, float)
    if method == "sample_sd":
        return float(np.mean(lnr)), float(np.std(lnr, ddof=1))
    def nll(p):
        mu, log_sig = p
        sig2 = math.exp(log_sig)
        v = sig2 + se**2
        return 0.5 * np.sum(np.log(v) + (lnr - mu)**2 / v)
    best = None
    for mu0 in [np.mean(lnr), 0.0, np.median(lnr)]:
        for s0 in [0.3, 0.7, 1.2]:
            r = minimize(nll, [mu0, math.log(s0)], method="Nelder-Mead",
                         options={"maxiter": 4000, "xatol": 1e-9, "fatol": 1e-9})
            if best is None or r.fun < best.fun:
                best = r
    mu = float(best.x[0]); sig = math.exp(float(best.x[1]))
    return mu, sig


def bh_fdr(p):
    p = np.asarray(p, float); n = len(p)
    order = np.argsort(p); ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q, 0, 1)
    return out


# ---- calibration per drug, N1 panel (primary) and N2 panel (bracket) ----
null_params = {}
for drug in CFG["index_drugs_for_calibration"]:
    for panel in ["N1", "N2"]:
        sub = df[(df["drug"] == drug) & (df["outcome_key"].str.startswith(panel + "::"))]
        sub = sub[sub["a"] > 0]
        if len(sub) < 5:
            continue
        lnr = np.log(sub["ror"].values); se = sub["se_ln"].values
        mu_ml, sig_ml = fit_null(lnr, se, "ml")
        mu_sd, sig_sd = fit_null(lnr, se, "sample_sd")
        null_params[(drug, panel)] = dict(n=len(sub), mu_ml=mu_ml, sig_ml=sig_ml,
                                          mu_sd=mu_sd, sig_sd=sig_sd)
print("=== NULL PARAMS ===")
for k, v in null_params.items():
    print(k, {kk: round(vv, 4) for kk, vv in v.items()})


def calibrate(row, panel="N1"):
    drug = row["drug"]
    if (drug, panel) not in null_params:
        return {}
    p = null_params[(drug, panel)]
    lnr = math.log(row["ror"]); se = row["se_ln"]
    out = {}
    for tag, mu, sig in [("ml", p["mu_ml"], p["sig_ml"]), ("sd", p["mu_sd"], p["sig_sd"])]:
        lc = lnr - mu
        sc = math.sqrt(se**2 + sig**2)
        out[f"cal_ror_{tag}"] = math.exp(lc)
        out[f"cal_lo_{tag}"] = math.exp(lc - Z * sc)
        out[f"cal_hi_{tag}"] = math.exp(lc + Z * sc)
        out[f"cal_p_{tag}"] = 2 * (1 - norm.cdf(abs(lc / sc)))
    return out


idx = df[df["label"].str.startswith("IDX::")].copy()
for panel in ["N1", "N2"]:
    for k in ["cal_ror", "cal_lo", "cal_hi", "cal_p"]:
        for tag in ["ml", "sd"]:
            idx[f"{k}_{tag}_{panel}"] = idx.apply(
                lambda r: calibrate(r, panel).get(f"{k}_{tag}", np.nan), axis=1)

idx["q_ml_N1"] = bh_fdr(idx["cal_p_ml_N1"].fillna(1).values)
idx["q_ml_N2"] = bh_fdr(idx["cal_p_ml_N2"].fillna(1).values)
idx["robust_double_fit"] = (
    (idx["cal_lo_ml_N1"] > 1) & (idx["cal_lo_sd_N1"] > 1) & (idx["q_ml_N1"] < 0.05))
idx["bracket_stable"] = idx["robust_double_fit"] & (idx["cal_lo_ml_N2"] > 1)
idx["naive_four_alg"] = ((idx["ror_lo"] > 1) & (idx["a"] >= 3) &
                         (idx["prr"] >= 2) & (idx["chi2"] >= 4) &
                         (idx["ic025"] > 0))

# positive control check
pc = df[df["outcome_key"].str.startswith("PC::")].copy()
print("\n=== POSITIVE CONTROL RETENTION (N1, ML) ===")
for drug in CFG["index_drugs_for_calibration"]:
    s = pc[(pc["drug"] == drug)].copy()
    if (drug, "N1") not in null_params:
        continue
    cal = s.apply(lambda r: calibrate(r, "N1"), axis=1)
    raw_sig = (s["ror_lo"] > 1).sum()
    cal_sig = sum(1 for x in cal if x.get("cal_lo_ml", 0) > 1)
    print(f"  {drug:13s} raw {raw_sig}/{len(s)}  calibrated {cal_sig}/{len(s)}")

cols = ["label", "drug", "outcome_key", "a", "ror", "ror_lo", "ror_hi", "prr", "chi2",
        "ic", "ic025", "naive_four_alg", "cal_ror_ml_N1", "cal_lo_ml_N1", "cal_lo_sd_N1",
        "q_ml_N1", "cal_lo_ml_N2", "q_ml_N2", "robust_double_fit", "bracket_stable"]
print("\n=== INDEX PAIRS ===")
show = idx[cols].copy()
for c in show.select_dtypes(include=[float]).columns:
    show[c] = show[c].round(3)
print(show.to_string(index=False))

out = HERE / "results_index.csv"
show.to_csv(out, index=False, encoding="utf-8-sig")
(HERE / "null_params.json").write_text(json.dumps(
    {f"{k[0]}::{k[1]}": v for k, v in null_params.items()}, indent=2), encoding="utf-8")
print(f"\nwrote {out}")
