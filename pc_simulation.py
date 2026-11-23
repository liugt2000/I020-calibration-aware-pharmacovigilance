#!/usr/bin/env python3
"""
pc_simulation.py -- Calibration simulation to test the positive-control claim.

Motivation
----------
The manuscript's "14/15 positive controls retained" argument is weak: label-established
positive controls have very large true effects and are retained almost by construction,
so their survival does not demonstrate that calibration avoids over-shrinkage.

This simulation injects KNOWN true effects into the fitted reporting-bias process and asks:
  (A) Type-I error: under a true null (logROR_true = 0), does calibration control the
      false-positive rate at the nominal 5%, whereas the naive (uncalibrated) test does not?
  (B) Power curve: for true effects spanning borderline-to-strong (ROR 1.2 .. 6), what
      fraction does calibration detect at alpha 0.05 / BH q 0.05? Where is the 80% point?
  (C) Borderline positive-control panel: 15 synthetic controls with deliberately modest
      true effects (ROR ~ 1.4-2.5). If calibration retains most of them, the "no
      over-shrinkage" claim is supported by effects far weaker than label controls.

The reporting-bias process uses the fitted empirical nulls: obs ~ N(mu + delta, sigma^2 + SE^2).
SEs are resampled (with replacement) from the empirical SE distribution of real index pairs.

Fixed seed 20261121. Output: simulation_results.json (console table also printed).
"""
import json, math
import numpy as np
from scipy.stats import norm

RNG = np.random.default_rng(20261121)
N_SIM = 20000
Z = 1.959963984540054
ALPHA = 0.05

NULLS = json.load(open("null_params.json", encoding="utf-8"))

# Empirical SE values from the 6 real index pairs (2x2 derived).
# se = sqrt(1/a + 1/b + 1/c + 1/d); a=index events, b=drug_total-a,
# c=event_total-a, d=N-drug_total-event_total+a
INDEX_SE = {
    "osimertinib_ILD": (33148, 1091, 100416),
    "osimertinib_QT":  (33148, 227, 38342),
    "tepotinib_AKI":   (691, 84, 330458),
    "lazertinib_ILD":  (865, None, None),      # filled from counts below
    "amivantamab_IRR": (3858, None, None),
    "amivantamab_ILD": (3858, None, None),
}

def se_from_counts(a, drug_total, event_total, N=21054159):
    b = drug_total - a; c = event_total - a; d = N - drug_total - event_total + a
    if min(a, b, c, d) <= 0: return None
    return math.sqrt(1/a + 1/b + 1/c + 1/d)

# Pull all real index SEs directly from counts.csv for accuracy.
import csv
se_pool = []
with open("counts.csv", encoding="utf-8-sig") as f:
    for row in csv.DictReader(f):
        if row["label"].startswith("IDX::"):
            s = se_from_counts(int(row["a"]), int(row["drug_total"]), int(row["event_total"]))
            if s: se_pool.append(s)
se_pool = np.array(sorted(set(round(x, 6) for x in se_pool)))
print("Empirical index-pair SEs:", np.round(se_pool, 4).tolist())

def bh_q(p):
    """Benjamini-Hochberg adjusted q for a 1-D array of p-values."""
    p = np.asarray(p, float); n = len(p)
    order = np.argsort(p); ranked = p[order]
    q = ranked * n / (np.arange(1, n+1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[order] = np.minimum(q, 1.0)
    return out

def run_null(mu, sigma, label):
    """(A) Type-I error under true null. Returns naive vs calibrated FPR."""
    se = RNG.choice(se_pool, size=N_SIM, replace=True)
    z = RNG.standard_normal(N_SIM)
    obs = mu + np.sqrt(sigma**2 + se**2) * z          # true logROR = 0
    p_naive = 2 * norm.sf(np.abs(obs) / se)            # assumes null=0, var=SE^2
    p_cal   = 2 * norm.sf(np.abs(obs - mu) / np.sqrt(se**2 + sigma**2))
    return {
        "null": label, "mu": round(mu, 4), "sigma": round(sigma, 4),
        "naive_fpr": round(float(np.mean(p_naive < ALPHA)), 4),
        "cal_fpr":   round(float(np.mean(p_cal < ALPHA)), 4),
        "naive_fpr_fdr": round(float(np.mean(bh_q(p_naive) < ALPHA)), 4),
        "cal_fpr_fdr":   round(float(np.mean(bh_q(p_cal) < ALPHA)), 4),
    }

def run_power(mu, sigma, label, rors=(1.2, 1.5, 2, 3, 4, 6)):
    """(B) Detection power across true effect sizes."""
    se = RNG.choice(se_pool, size=N_SIM, replace=True)
    z = RNG.standard_normal(N_SIM)
    rows = []
    for ror in rors:
        delta = math.log(ror)
        obs = mu + delta + np.sqrt(sigma**2 + se**2) * z
        p_cal = 2 * norm.sf(np.abs(obs - mu) / np.sqrt(se**2 + sigma**2))
        rows.append({
            "true_ROR": ror,
            "power_p05": round(float(np.mean(p_cal < ALPHA)), 4),
            "power_q05": round(float(np.mean(bh_q(p_cal) < ALPHA)), 4),
        })
    return {"null": label, "mu": round(mu,4), "sigma": round(sigma,4), "curve": rows}

def run_borderline_panel(mu, sigma, label, n_ctrl=15, seed_offset=0):
    """(C) 15 synthetic controls with modest true effects (ROR 1.4-2.5)."""
    rg = np.random.default_rng(20261121 + seed_offset)
    true_rors = np.round(np.linspace(1.4, 2.5, n_ctrl), 3)
    delta = np.log(true_rors)
    se = rg.choice(se_pool, size=n_ctrl, replace=True)
    z = rg.standard_normal(n_ctrl)
    obs = mu + delta + np.sqrt(sigma**2 + se**2) * z
    p_cal = 2 * norm.sf(np.abs(obs - mu) / np.sqrt(se**2 + sigma**2))
    q_cal = bh_q(p_cal)
    # naive comparison (what a label-free naive screen would call)
    p_naive = 2 * norm.sf(np.abs(obs) / se)
    return {
        "null": label, "mu": round(mu,4), "sigma": round(sigma,4),
        "true_rors": true_rors.tolist(),
        "p_cal": np.round(p_cal, 6).tolist(),
        "q_cal": np.round(q_cal, 6).tolist(),
        "retained_p05": int(np.sum(p_cal < ALPHA)),
        "retained_q05": int(np.sum(q_cal < ALPHA)),
        "naive_retained_p05": int(np.sum(p_naive < ALPHA)),
        "n_ctrl": n_ctrl,
    }

def main():
    out = {"n_sim": N_SIM, "alpha": ALPHA, "seed": 20261121,
           "se_pool": np.round(se_pool, 6).tolist()}

    # Two null families spanning the bracket: osimertinib N1 (primary, mu<0)
    # and osimertinib N2 (upper-bound confounded sentinel, sigma large).
    null_families = {
        "N1_primary": (NULLS["osimertinib::N1"]["mu_ml"], NULLS["osimertinib::N1"]["sig_ml"]),
        "N2_upper":   (NULLS["osimertinib::N2"]["mu_ml"], NULLS["osimertinib::N2"]["sig_ml"]),
    }

    out["typeI_error"] = []
    out["power"] = []
    out["borderline_panel"] = []
    for label, (mu, sigma) in null_families.items():
        out["typeI_error"].append(run_null(mu, sigma, label))
        out["power"].append(run_power(mu, sigma, label))
        out["borderline_panel"].append(run_borderline_panel(mu, sigma, label))

    # Also a zero-bias reference (perfect world, mu=0 sigma=0) to show the ceiling.
    out["power"].append(run_power(0.0, 0.0, "ideal_no_bias"))

    with open("simulation_results.json", "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

    # ---- console report ----
    print("\n=== (A) TYPE-I ERROR (true null) ===")
    for r in out["typeI_error"]:
        print(f"  {r['null']:12s} mu={r['mu']:+.3f} sig={r['sigma']:.3f} | "
              f"naive FPR={r['naive_fpr']:.3f}  cal FPR={r['cal_fpr']:.3f}  "
              f"(FDR: naive={r['naive_fpr_fdr']:.3f} cal={r['cal_fpr_fdr']:.3f})")
    print("\n=== (B) POWER CURVE (calibrated, p<0.05) ===")
    hdr = "  ROR true | " + " ".join(f"{k:>12s}" for k in null_families) + "   ideal"
    print(hdr)
    rors = [c["true_ROR"] for c in out["power"][0]["curve"]]
    for i, ror in enumerate(rors):
        vals = [next(c for c in out["power"] if c["null"]==k)["curve"][i]["power_p05"]
                for k in null_families]
        ideal = out["power"][-1]["curve"][i]["power_p05"]
        print(f"  {ror:8.1f} | " + " ".join(f"{v:12.3f}" for v in vals) + f"   {ideal:6.3f}")
    print("\n=== (C) BORDERLINE PANEL (15 controls, ROR 1.4-2.5) ===")
    for r in out["borderline_panel"]:
        print(f"  {r['null']:12s} mu={r['mu']:+.3f} sig={r['sigma']:.3f} | "
              f"cal retained {r['retained_p05']}/{r['n_ctrl']} (p<.05), "
              f"{r['retained_q05']}/{r['n_ctrl']} (BH q<.05); "
              f"naive {r['naive_retained_p05']}/{r['n_ctrl']}")
    print("\nSaved simulation_results.json")

if __name__ == "__main__":
    main()
