#!/usr/bin/env python3
"""decision_impact.py -- Quantify the decision consequences of naive vs
calibration-aware signal triage, using the real tier calls and simulation FPRs.

Pipeline stages (which pairs remain eligible for an AGENT-SPECIFIC safety action):
  Stage 0  naive four-algorithm concordance        -> all 6 flagged
  Stage 1  + empirical calibration (robust both fits, FDR) -> 5 (tepotinib not estimable)
  Stage 2  + active-comparator re-estimation        -> 3 (three ILD attenuate to class effects)
  Stage 3  + estimate-stability requirement (Tier A)-> 1 (osimertinib-QT; IRR unstable)

Outputs decision_impact.json and console report.
"""
import json

f = json.load(open("fill_values.json", encoding="utf-8"))
tiers = f["tiers"]  # pair -> A/B/C
sim = json.load(open("simulation_results.json", encoding="utf-8"))

# Stage membership from the frozen tier calls
all_pairs = list(tiers.keys())
stage0 = all_pairs                                             # 6, all naive-flagged
stage1 = [p for p in all_pairs if tiers[p] in ("A", "B", "C")
          and p != "P3_tepotinib_AKI"]                        # calibration-estimable
stage2 = [p for p in stage1 if tiers[p] == "A" or tiers[p] == "C"]  # not class-effect (ILD=B)
stage3 = [p for p in stage2 if tiers[p] == "A"]               # stable estimate

# False-positive rate from simulation (true null)
fpr = {r["null"]: r["naive_fpr"] for r in sim["typeI_error"]}
cal_fpr = {r["null"]: r["cal_fpr"] for r in sim["typeI_error"]}

# Expected false-positive agent-specific actions if a naive screen of N pairs
# were adopted at the observed reporting-process null: N * FPR (binomial mean).
n = len(stage0)
expected_fp_naive = n * sum(fpr.values())/len(fpr)
expected_fp_cal   = n * sum(cal_fpr.values())/len(cal_fpr)

out = {
    "stages": [
        {"stage": "0 Naive four-algorithm concordance", "n": len(stage0),
         "pairs": stage0, "action": "all pairs enter agent-specific label/signal review"},
        {"stage": "1 + Empirical calibration (robust both fits + FDR)", "n": len(stage1),
         "pairs": stage1, "removed": ["P3_tepotinib_AKI"],
         "action": "tepotinib-AKI excluded: null not estimable at current volume"},
        {"stage": "2 + Active-comparator re-estimation", "n": len(stage2),
         "pairs": stage2,
         "removed": [p for p in stage1 if tiers[p] == "B"],
         "action": "three ILD pairs reclassified as class effects (class-level monitoring, not agent-specific)"},
        {"stage": "3 + Estimate stability (final Tier A)", "n": len(stage3),
         "pairs": stage3,
         "removed": [p for p in stage2 if tiers[p] == "C"],
         "action": "amivantamab-IRR excluded: comparator zero events, UB=24,683"},
    ],
    "naive_to_final_reduction": {"from": len(stage0), "to": len(stage3),
                                 "pct": round(100*(1-len(stage3)/len(stage0)), 1)},
    "false_positive_actions_if_naive_adopted": {
        "fpr_naive_mean": round(sum(fpr.values())/len(fpr), 4),
        "fpr_calibrated_mean": round(sum(cal_fpr.values())/len(cal_fpr), 4),
        "expected_false_actions_per_6_naive": round(expected_fp_naive, 2),
        "expected_false_actions_per_6_calibrated": round(expected_fp_cal, 2),
        "interpretation": "Under the fitted reporting-process null, a naive screen of 6 pairs would be expected to generate ~5.6 false agent-specific actions; calibration reduces this to ~0.3.",
    },
    "decision_rule_note": "An 'agent-specific action' means a label-change discussion or targeted product-specific monitoring decision. Class-effect reclassification is not counted as a loss: it still triggers class-level monitoring.",
}

json.dump(out, open("decision_impact.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)

print("=== Decision funnel (agent-specific actions) ===")
for s in out["stages"]:
    print(f"  Stage {s['stage'][:1]}: {s['n']} pairs  | {s['action']}")
print(f"\nReduction: {out['naive_to_final_reduction']['from']} -> "
      f"{out['naive_to_final_reduction']['to']} "
      f"({out['naive_to_final_reduction']['pct']}% fewer agent-specific triggers)")
print(f"\nExpected false actions / 6 pairs: naive {expected_fp_naive:.2f}, "
      f"calibrated {expected_fp_cal:.2f}")
print("Saved decision_impact.json")
