#!/usr/bin/env python3
"""
pc_misspec.py -- Model-misspecification stress test for empirical calibration.

Motivation (reviewer R2, major concern)
---------------------------------------
The base simulation (pc_simulation.py) generates data from the SAME normal
null that calibration fits, so it can only show "if the model is right,
calibration works" -- a self-referential demonstration. This companion
script deliberately VIOLATES the fitted model in three ways and asks whether
calibration (still applied with the fitted mu, sigma) keeps type-I error
near nominal, degrades gracefully, or fails:

  M0  correctly specified: bias ~ N(mu_hat, sigma_hat^2 + SE^2)       [reference]
  M1  variance misspecification: true sigma = k * sigma_hat (k = 1.5, 2.0)
      -- the fitted null UNDERSTATES dispersion (e.g. sparse control panel).
  M2  heavy-tailed bias: bias drawn from a t distribution scaled to the same
      sigma_hat, with df = 5 (leptokurtic reporting bias / outlier reports).
  M3  location drift: each pair's true null location is mu_hat + eps,
      eps ~ N(0, tau^2), tau = 0.15, 0.30 -- the control panel mis-centres the null.

For every scenario we report naive FPR (uncalibrated, null=0), calibrated FPR
using the FITTED parameters (what an analyst would actually do), and an ORACLE
calibrated FPR using the TRUE generating parameters (ceiling / diagnostic).

Interpretation logic: if calibrated FPR stays near 5% under mild misspecification
and degrades gradually (not catastrophically) under severe misspecification, the
base simulation's conclusion is robust to its own assumptions. If it fails
catastrophically even under mild misspecification, that is an applicability
boundary that must be stated in the manuscript. Either way this removes the
circularity objection by testing the procedure off-model.

Fixed seed 20261121. Output: misspec_results.json.
"""
import json, math
import numpy as np
from scipy.stats import norm, t as tdist

RNG = np.random.default_rng(20261121)
N_SIM = 20000
ALPHA = 0.05

NULLS = json.load(open("null_params.json", encoding="utf-8"))
SE_POOL = np.array(json.load(open("simulation_results.json", encoding="utf-8"))["se_pool"])

MU_HAT = NULLS["osimertinib::N1"]["mu_ml"]     # -0.836
SIG_HAT = NULLS["osimertinib::N1"]["sig_ml"]   # 0.214


def fpr_summary(bias, se, mu_fit, sig_fit):
    """FPR for naive and calibrated tests; bias is the true null draw (delta=0)."""
    obs = bias + np.sqrt(se**2) * 0.0  # bias already includes dispersion scale
    # reconstruct observed logROR: true logROR=0, so obs = bias_draw + SE*Z_meas
    # bias_draw already carries sigma-scale dispersion; add pure measurement error se*Z
    z_meas = RNG.standard_normal(len(se))
    obs = bias + se * z_meas
    p_naive = 2 * norm.sf(np.abs(obs) / se)
    p_cal = 2 * norm.sf(np.abs(obs - mu_fit) / np.sqrt(se**2 + sig_fit**2))
    return {
        "naive_fpr": round(float(np.mean(p_naive < ALPHA)), 4),
        "cal_fpr": round(float(np.mean(p_cal < ALPHA)), 4),
    }


def run():
    se = RNG.choice(SE_POOL, size=N_SIM, replace=True)
    out = {
        "n_sim": N_SIM, "seed": 20261121, "alpha": ALPHA,
        "fitted_null": {"mu_hat": round(MU_HAT, 4), "sigma_hat": round(SIG_HAT, 4)},
        "scenarios": [],
    }

    def add(name, desc, bias, mu_fit, sig_fit, mu_true, sig_true):
        res = fpr_summary(bias, se, mu_fit, sig_fit)
        # oracle uses true generating params
        z_meas = RNG.standard_normal(N_SIM)
        obs = bias + se * z_meas
        p_oracle = 2 * norm.sf(np.abs(obs - mu_true) / np.sqrt(se**2 + sig_true**2))
        res.update({
            "scenario": name, "description": desc,
            "oracle_cal_fpr": round(float(np.mean(p_oracle < ALPHA)), 4),
            "true_params": {"mu": round(float(mu_true), 4) if np.isscalar(mu_true) else "pairwise",
                            "sigma": round(float(sig_true), 4)},
            "fit_params": {"mu": round(float(mu_fit), 4) if np.isscalar(mu_fit) else "pairwise",
                           "sigma": round(float(sig_fit), 4)},
        })
        out["scenarios"].append(res)

    # M0 correctly specified
    bias0 = MU_HAT + SIG_HAT * RNG.standard_normal(N_SIM)
    add("M0_correct", "bias ~ N(mu_hat, sigma_hat^2); model correct",
        bias0, MU_HAT, SIG_HAT, MU_HAT, SIG_HAT)

    # M1 variance misspecification: true dispersion inflated, fit unchanged
    for k in (1.5, 2.0):
        sig_true = k * SIG_HAT
        biask = MU_HAT + sig_true * RNG.standard_normal(N_SIM)
        add(f"M1_sigma_x{k}", f"true sigma = {k}x fitted (panel underestimates dispersion)",
            biask, MU_HAT, SIG_HAT, MU_HAT, sig_true)

    # M2 heavy-tailed bias: t(df=5) scaled to sigma_hat (unit variance)
    df = 5
    tv = tdist.rvs(df, size=N_SIM, random_state=RNG) / math.sqrt(df / (df - 2))
    bias2 = MU_HAT + SIG_HAT * tv
    add("M2_heavy_t5", "bias ~ t(5) scaled to sigma_hat; normal fit mis-specified tail",
        bias2, MU_HAT, SIG_HAT, MU_HAT, SIG_HAT)

    # M3 location drift: pairwise eps ~ N(0, tau^2)
    for tau in (0.15, 0.30):
        eps = RNG.normal(0, tau, N_SIM)
        bias3 = MU_HAT + eps + SIG_HAT * RNG.standard_normal(N_SIM)
        # fitted still uses scalar mu_hat; oracle uses pairwise true mu = mu_hat+eps
        se_s = se
        z_meas = RNG.standard_normal(N_SIM)
        obs = bias3 + se_s * z_meas
        p_naive = 2 * norm.sf(np.abs(obs) / se_s)
        p_cal = 2 * norm.sf(np.abs(obs - MU_HAT) / np.sqrt(se_s**2 + SIG_HAT**2))
        mu_true_vec = MU_HAT + eps
        p_oracle = 2 * norm.sf(np.abs(obs - mu_true_vec) / np.sqrt(se_s**2 + SIG_HAT**2))
        out["scenarios"].append({
            "scenario": f"M3_drift_tau{tau}",
            "description": f"pairwise null location mu_hat + eps, eps~N(0,{tau}^2)",
            "naive_fpr": round(float(np.mean(p_naive < ALPHA)), 4),
            "cal_fpr": round(float(np.mean(p_cal < ALPHA)), 4),
            "oracle_cal_fpr": round(float(np.mean(p_oracle < ALPHA)), 4),
            "true_params": {"mu": "pairwise", "sigma": round(SIG_HAT, 4)},
            "fit_params": {"mu": round(MU_HAT, 4), "sigma": round(SIG_HAT, 4)},
        })

    with open("misspec_results.json", "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

    print("\n=== MODEL-MISSPECIFICATION STRESS TEST (FPR at alpha=0.05) ===")
    print(f"  fitted null: mu={MU_HAT:+.3f}, sigma={SIG_HAT:.3f}")
    print(f"  {'scenario':18s} {'naive':>8s} {'calibrated':>11s} {'oracle':>8s}")
    for r in out["scenarios"]:
        print(f"  {r['scenario']:18s} {r['naive_fpr']:8.3f} {r['cal_fpr']:11.3f} {r['oracle_cal_fpr']:8.3f}")
    print("\nSaved misspec_results.json")


if __name__ == "__main__":
    run()
