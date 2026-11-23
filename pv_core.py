"""
Shared core for I020 PV analysis: disproportionality, empirical-null calibration,
active-comparator ROR. Parameterised by input CSV so the primary and suspect-only
pipelines run through identical logic.
"""
import json, math
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

HERE = Path(__file__).parent
CFG = json.loads((HERE / "config.json").read_text(encoding="utf-8"))
Z = 1.959963984540054


def cells(row):
    a = max(int(row["a"]), 0); dt = int(row["drug_total"])
    et = int(row["event_total"]); N = int(row["N"])
    b = max(dt - a, 0); c = max(et - a, 0); d = max(N - a - b - c, 0)
    return a, b, c, d


def disprop(a, b, c, d):
    ac, bc, cc, dc = a, b, c, d
    if min(ac, bc, cc, dc) == 0:
        ac, bc, cc, dc = ac + .5, bc + .5, cc + .5, dc + .5
    N = ac + bc + cc + dc
    ror = (ac * dc) / (bc * cc)
    se = math.sqrt(1/ac + 1/bc + 1/cc + 1/dc)
    prr = (ac/(ac+bc)) / (cc/(cc+dc))
    chi2 = N * (ac*dc - bc*cc)**2 / ((ac+bc)*(cc+dc)*(ac+cc)*(bc+dc))
    E = (ac+bc)*(ac+cc)/N
    ic = math.log2((ac + 0.5) / (E + 0.5))
    ic_sd = math.sqrt(1/(ac+0.5) + 1/(bc+0.5) + 1/(cc+0.5) + 1/(dc+0.5)) / math.log(2)
    return dict(a=a, b=b, c=c, d=d, ror=ror,
                ror_lo=math.exp(math.log(ror) - Z*se),
                ror_hi=math.exp(math.log(ror) + Z*se), se_ln=se,
                prr=prr, chi2=chi2, ic=ic, ic025=ic - Z*ic_sd,
                lam=ac/E if E > 0 else 0, E=E)


def fit_null(lnr, se, method="ml"):
    lnr = np.asarray(lnr, float); se = np.asarray(se, float)
    if method == "sample_sd":
        return float(np.mean(lnr)), float(np.std(lnr, ddof=1))
    def nll(p):
        sig2 = math.exp(p[1]); v = sig2 + se**2
        return 0.5 * np.sum(np.log(v) + (lnr - p[0])**2 / v)
    best = None
    for mu0 in [np.mean(lnr), 0.0, np.median(lnr)]:
        for s0 in [0.3, 0.7, 1.2]:
            r = minimize(nll, [mu0, math.log(s0)], method="Nelder-Mead",
                         options={"maxiter": 4000, "xatol": 1e-9, "fatol": 1e-9})
            if best is None or r.fun < best.fun:
                best = r
    return float(best.x[0]), math.exp(float(best.x[1]))


def bh_fdr(p):
    p = np.asarray(p, float); n = len(p)
    order = np.argsort(p); ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q, 0, 1)
    return out


def run_index_pipeline(counts_csv, verbose=False):
    """Return (idx DataFrame with calibration columns, null_params dict, pc_retention)."""
    df = pd.read_csv(HERE / counts_csv, encoding="utf-8-sig")
    for col in ["ror", "ror_lo", "ror_hi", "se_ln", "prr", "chi2", "ic", "ic025", "lam", "E", "a", "b", "c", "d"]:
        df[col] = df.apply(lambda r: disprop(*cells(r))[col], axis=1)

    nulls = {}
    for drug in CFG["index_drugs_for_calibration"]:
        for panel in ["N1", "N2"]:
            sub = df[(df["drug"] == drug) & (df["outcome_key"].str.startswith(panel + "::"))]
            sub = sub[sub["a"] > 0]
            if len(sub) < 5:
                continue
            mu_ml, sig_ml = fit_null(np.log(sub["ror"].values), sub["se_ln"].values, "ml")
            mu_sd, sig_sd = fit_null(np.log(sub["ror"].values), sub["se_ln"].values, "sample_sd")
            nulls[(drug, panel)] = dict(n=len(sub), mu_ml=mu_ml, sig_ml=sig_ml,
                                        mu_sd=mu_sd, sig_sd=sig_sd)

    def calibrate(row, panel="N1"):
        drug = row["drug"]
        if (drug, panel) not in nulls:
            return {}
        p = nulls[(drug, panel)]
        lnr = math.log(row["ror"]); se = row["se_ln"]
        out = {}
        for tag, mu, sig in [("ml", p["mu_ml"], p["sig_ml"]), ("sd", p["mu_sd"], p["sig_sd"])]:
            lc = lnr - mu; sc = math.sqrt(se**2 + sig**2)
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
    idx["robust_double_fit"] = ((idx["cal_lo_ml_N1"] > 1) & (idx["cal_lo_sd_N1"] > 1) & (idx["q_ml_N1"] < 0.05))
    idx["bracket_stable"] = idx["robust_double_fit"] & (idx["cal_lo_ml_N2"] > 1)
    idx["naive_four_alg"] = ((idx["ror_lo"] > 1) & (idx["a"] >= 3) & (idx["prr"] >= 2) &
                             (idx["chi2"] >= 4) & (idx["ic025"] > 0))

    pc = df[df["outcome_key"].str.startswith("PC::")].copy()
    pc_ret = {}
    for drug in CFG["index_drugs_for_calibration"]:
        s = pc[pc["drug"] == drug]
        if (drug, "N1") not in nulls or len(s) == 0:
            continue
        cal = s.apply(lambda r: calibrate(r, "N1"), axis=1)
        pc_ret[drug] = dict(n=len(s), raw=int((s["ror_lo"] > 1).sum()),
                            calibrated=int(sum(1 for x in cal if x.get("cal_lo_ml", 0) > 1)))
    return idx, nulls, pc_ret


def run_active_pipeline(counts_csv):
    df = pd.read_csv(HERE / counts_csv, encoding="utf-8-sig")
    rows = []
    for _, r in df.iterrows():
        a = int(r["a"]); it = int(r["index_total"]); c = int(r["c"]); ct = int(r["comparator_total"])
        b = it - a; d = ct - c
        if min(a, b, c, d) == 0:
            a, b, c, d = a+.5, b+.5, c+.5, d+.5
        ror = (a*d)/(b*c)
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
        rows.append(dict(pair=r["pair"], index=r["index"], comparators=r["comparators"],
                         a=int(r["a"]), c=int(r["c"]), ror_active=ror,
                         lo=np.exp(np.log(ror) - Z*se), hi=np.exp(np.log(ror) + Z*se),
                         lnr_active=np.log(ror)))
    return pd.DataFrame(rows)
