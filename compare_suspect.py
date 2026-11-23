"""
Compare primary (all-role) vs suspect-only pipelines. Quantifies the effect of
restricting exposure to suspect reports. Outputs a machine-readable comparison
and a markdown table for the manuscript / supplement.
"""
import json
from pathlib import Path
import pandas as pd
from pv_core import run_index_pipeline, run_active_pipeline

HERE = Path(__file__).parent

idx_p, null_p, pc_p = run_index_pipeline("counts.csv")
idx_s, null_s, pc_s = run_index_pipeline("counts_suspect.csv")
act_p = run_active_pipeline("counts_active.csv")
act_s = run_active_pipeline("counts_active_suspect.csv")

# ---- exposure counts ----
cnt_p = pd.read_csv(HERE / "counts.csv", encoding="utf-8-sig")
cnt_s = pd.read_csv(HERE / "counts_suspect.csv", encoding="utf-8-sig")
drugs = ["osimertinib", "tepotinib", "lazertinib", "amivantamab"]
exposure = []
for d in drugs:
    a = cnt_p[cnt_p.drug == d]["drug_total"].iloc[0]
    b = cnt_s[cnt_s.drug == d]["drug_total"].iloc[0]
    exposure.append(dict(drug=d, all_roles=int(a), suspect=int(b),
                         delta_pct=round(100*(a-b)/a, 2) if a else None))
exp_df = pd.DataFrame(exposure)
print("=== EXPOSURE COUNTS (drug_total) ===")
print(exp_df.to_string(index=False))

# ---- index pair ROR comparison ----
cols = ["label", "a", "ror", "ror_lo", "ror_hi", "robust_double_fit", "bracket_stable", "naive_four_alg"]
cmp = idx_p[cols].merge(idx_s[cols], on="label", suffixes=("_all", "_sus"))
cmp["ror_delta_pct"] = (100*(cmp.ror_sus - cmp.ror_all)/cmp.ror_all).round(2)
print("\n=== INDEX PAIRS: all vs suspect ===")
print(cmp[["label", "a_all", "a_sus", "ror_all", "ror_sus", "ror_delta_pct",
           "robust_double_fit_all", "robust_double_fit_sus"]].round(3).to_string(index=False))

# ---- active comparator comparison ----
ac = act_p[["pair", "a", "c", "ror_active", "lo", "hi"]].merge(
    act_s[["pair", "a", "c", "ror_active", "lo", "hi"]], on="pair", suffixes=("_all", "_sus"))
ac["ror_delta_pct"] = (100*(ac.ror_active_sus - ac.ror_active_all)/ac.ror_active_all).round(2)
print("\n=== ACTIVE COMPARATOR ROR: all vs suspect ===")
print(ac[["pair", "a_all", "a_sus", "c_all", "c_sus", "ror_active_all", "ror_active_sus", "ror_delta_pct"]].round(3).to_string(index=False))

# ---- null parameter comparison ----
null_rows = []
for k in sorted(set(null_p) | set(null_s)):
    p = null_p.get(k, {}); s = null_s.get(k, {})
    null_rows.append(dict(panel="::".join(k), n_all=p.get("n"), n_sus=s.get("n"),
                          mu_all=round(p.get("mu_ml", float('nan')), 4),
                          mu_sus=round(s.get("mu_ml", float('nan')), 4),
                          sig_all=round(p.get("sig_ml", float('nan')), 4),
                          sig_sus=round(s.get("sig_ml", float('nan')), 4)))
ndf = pd.DataFrame(null_rows)
print("\n=== NULL PARAMS: all vs suspect ===")
print(ndf.to_string(index=False))

# ---- summary ----
summary = dict(
    max_exposure_delta_pct=float(exp_df.delta_pct.max()),
    max_index_ror_delta_pct=float(cmp.ror_delta_pct.abs().max()),
    tier_changes=[],
    attenuation_calls_all=[],
    attenuation_calls_sus=[],
)
# check robustness calls unchanged
for _, r in cmp.iterrows():
    if r.robust_double_fit_all != r.robust_double_fit_sus:
        summary["tier_changes"].append(str(r.label))
# positive controls
summary["positive_controls"] = dict(all=pc_p, suspect=pc_s)

out = HERE / "suspect_comparison.json"
out.write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
exp_df.to_csv(HERE / "suspect_exposure.csv", index=False, encoding="utf-8-sig")
cmp.to_csv(HERE / "suspect_index_compare.csv", index=False, encoding="utf-8-sig")
ac.to_csv(HERE / "suspect_active_compare.csv", index=False, encoding="utf-8-sig")
ndf.to_csv(HERE / "suspect_null_compare.csv", index=False, encoding="utf-8-sig")
print(f"\n=== SUMMARY ===")
print(json.dumps(summary, indent=2, default=str))
print(f"\nwrote suspect_*.csv and {out.name}")
