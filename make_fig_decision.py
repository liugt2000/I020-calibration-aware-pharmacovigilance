#!/usr/bin/env python3
"""Fig. 8 -- Decision-impact funnel: naive vs calibration-aware signal triage."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

D = json.load(open("decision_impact.json", encoding="utf-8"))
stages = D["stages"]

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 600})

fig, ax = plt.subplots(figsize=(13, 4.8))
ax.set_xlim(0, 13); ax.set_ylim(0, 5); ax.axis("off")

labels = ["Naive four-algorithm\nconcordance", "+ Empirical\ncalibration",
          "+ Active-comparator\nre-estimation", "+ Estimate stability\n(final Tier A)"]
counts = [s["n"] for s in stages]
colors = ["#B23A48", "#D98324", "#2A7F9E", "#1B6B4A"]
xs = [1.0, 4.0, 7.0, 10.2]

for i, (x, lab, cnt, col) in enumerate(zip(xs, labels, counts, colors)):
    # box height scales with count
    h = 0.9 + 1.35 * cnt / 6
    y = 3.1 - h/2
    box = FancyBboxPatch((x, y), 1.9, h, boxstyle="round,pad=0.08,rounding_size=0.12",
                         facecolor=col, edgecolor="white", linewidth=2)
    ax.add_patch(box)
    ax.text(x+0.95, y+h/2, str(cnt), ha="center", va="center", fontsize=26,
            color="white", fontweight="bold")
    ax.text(x+0.95, 0.75, lab, ha="center", va="center", fontsize=9.5,
            color="#222", linespacing=1.3)
    if i < len(xs)-1:
        ax.annotate("", xy=(xs[i+1]-0.15, 3.1), xytext=(x+2.0, 3.1),
                    arrowprops=dict(arrowstyle="-|>", lw=1.6, color="#666"))

# removal annotations above arrows
rems = ["−1 not estimable\n(tepotinib-AKI)", "−3 class effects\n(all ILD)", "−1 unstable\n(amivantamab-IRR)"]
for i, r in enumerate(rems):
    ax.text((xs[i]+2.0+xs[i+1]-0.15)/2, 4.35, r, ha="center", va="bottom",
            fontsize=8.5, color="#555")

ax.text(6.5, 4.85, "Decision-impact funnel: pairs eligible for an agent-specific safety action",
        ha="center", va="center", fontsize=12.5, fontweight="bold")

fp = D["false_positive_actions_if_naive_adopted"]
ax.text(6.5, 0.12,
        f"Expected false agent-specific actions per 6 screened pairs: "
        f"naive {fp['expected_false_actions_per_6_naive']:.1f}  vs  calibrated {fp['expected_false_actions_per_6_calibrated']:.1f}"
        f"    |    net reduction 6 → 1 ({D['naive_to_final_reduction']['pct']}% fewer triggers)",
        ha="center", va="bottom", fontsize=9, color="#333",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#F2F5F7", edgecolor="#CBD5DC"))

fig.savefig("figures/Fig7_decision.png", dpi=600, bbox_inches="tight")
print("saved figures/Fig7_decision.png")
