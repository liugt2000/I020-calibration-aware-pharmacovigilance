#!/usr/bin/env python3
"""make_fig_misspec.py -- Figure for the model-misspecification stress test."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

res = json.load(open("misspec_results.json", encoding="utf-8"))
sc = res["scenarios"]
names = [s["scenario"] for s in sc]
naive = [s["naive_fpr"] * 100 for s in sc]
cal = [s["cal_fpr"] * 100 for s in sc]
oracle = [s["oracle_cal_fpr"] * 100 for s in sc]

labels = ["M0\ncorrect", "M1\n$\\sigma\\times1.5$", "M1\n$\\sigma\\times2$",
          "M2\nheavy tail", "M3\ndrift .15", "M3\ndrift .30"]
x = np.arange(len(names))
w = 0.26

fig, ax = plt.subplots(figsize=(8.6, 4.2), dpi=600)
b1 = ax.bar(x - w, naive, w, label="Naive (null = 0)", color="#c0392b")
b2 = ax.bar(x, cal, w, label="Calibrated (fitted $\\mu,\\sigma$)", color="#2c6e8f")
b3 = ax.bar(x + w, oracle, w, label="Oracle (true params)", color="#7f8c8d")
ax.axhline(5, color="#2c3e50", ls="--", lw=1.1, label="Nominal 5%")
for bars in (b1, b2, b3):
    for bar in bars:
        h = bar.get_height()
        ax.annotate(f"{h:.0f}", xy=(bar.get_x() + bar.get_width()/2, h),
                    xytext=(0, 2), textcoords="offset points", ha="center",
                    va="bottom", fontsize=6.4)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("False-positive rate at $\\alpha$=0.05 (%)", fontsize=9)
ax.set_ylim(0, 112)
ax.set_title("Empirical calibration under model misspecification", fontsize=10.5, pad=10)
ax.legend(fontsize=7.6, loc="upper left", ncol=2, frameon=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="y", labelsize=8)
fig.text(0.5, -0.02,
         "M0 model correct; M1 fitted null underestimates dispersion; M2 heavy-tailed (t5) bias; "
         "M3 pairwise null-centre drift.",
         ha="center", fontsize=7, color="#555")
fig.tight_layout()
fig.savefig("figures/SuppFig_S2_misspec.png", bbox_inches="tight", dpi=600)
print("Saved figures/SuppFig_S2_misspec.png")
