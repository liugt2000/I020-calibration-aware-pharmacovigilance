#!/usr/bin/env python3
"""Fig. 7 -- Simulation-based validation of empirical calibration.
Panel A: Type-I error, naive vs calibrated (two null families).
Panel B: Power curve across true ROR 1.2-6 (N1 primary, N2 upper, ideal reference).
Panel C: Borderline positive-control panel retention (ROR 1.4-2.5).
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = json.load(open("simulation_results.json", encoding="utf-8"))

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.9, "figure.dpi": 600,
})

fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
C_NAIVE = "#B23A48"; C_CAL = "#2A7F9E"; C_N2 = "#D98324"; C_IDEAL = "#7A8B99"

# ---- Panel A: Type-I error ----
ax = axes[0]
t1 = D["typeI_error"]
labels = [r["null"].replace("_", " ") for r in t1]
naive = [r["naive_fpr"]*100 for r in t1]
cal = [r["cal_fpr"]*100 for r in t1]
x = np.arange(len(labels)); w = 0.34
b1 = ax.bar(x-w/2, naive, w, color=C_NAIVE, label="Naive (uncalibrated)", edgecolor="white")
b2 = ax.bar(x+w/2, cal, w, color=C_CAL, label="Calibrated", edgecolor="white")
for bars in (b1, b2):
    for b in bars:
        ax.text(b.get_x()+b.get_width()/2, b.get_height()+2, f"{b.get_height():.1f}%",
                ha="center", va="bottom", fontsize=9)
ax.axhline(5, color="#333", ls="--", lw=1, label="Nominal 5%")
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel("False-positive rate at α=0.05 (%)")
ax.set_ylim(0, 112)
ax.set_title("A  Type-I error under a true null", loc="left", fontweight="bold", fontsize=11)
ax.legend(frameon=False, fontsize=8, loc="upper center")

# ---- Panel B: Power curve ----
ax = axes[1]
rors = [c["true_ROR"] for c in D["power"][0]["curve"]]
series = {
    "N1 primary":  [c for c in D["power"] if c["null"]=="N1_primary"][0]["curve"],
    "N2 upper":    [c for c in D["power"] if c["null"]=="N2_upper"][0]["curve"],
    "Ideal (no bias)": D["power"][-1]["curve"],
}
for (name, curve), color, ls, mk in zip(
        series.items(), [C_CAL, C_N2, C_IDEAL], ["-", "-", ":"], ["o", "s", "^"]):
    ax.plot(rors, [r["power_p05"] for r in curve], color=color, ls=ls, marker=mk,
            ms=4.5, lw=1.8, label=name)
ax.axhline(0.8, color="#888", ls="--", lw=1)
ax.text(1.2, 0.82, "80% power", fontsize=8, color="#666")
ax.set_xlabel("True ROR (log scale)"); ax.set_ylabel("Calibrated detection power")
ax.set_xscale("log"); ax.set_xticks(rors)
ax.set_xticklabels([str(r) for r in rors])
ax.set_ylim(0, 1.05)
ax.set_title("B  Power to detect a true effect", loc="left", fontweight="bold", fontsize=11)
ax.legend(frameon=False, fontsize=8, loc="lower right")

# ---- Panel C: Borderline panel ----
ax = axes[2]
bp = D["borderline_panel"]
labels = [r["null"].replace("_", " ") for r in bp]
cal_ret = [r["retained_q05"] for r in bp]
naive_ret = [r["naive_retained_p05"] for r in bp]
x = np.arange(len(labels)); w = 0.34
b1 = ax.bar(x-w/2, naive_ret, w, color=C_NAIVE, label="Naive retained", edgecolor="white")
b2 = ax.bar(x+w/2, cal_ret, w, color=C_CAL, label="Calibrated retained (BH q<0.05)", edgecolor="white")
for bars in (b1, b2):
    for b in bars:
        ax.text(b.get_x()+b.get_width()/2, b.get_height()+0.25,
                f"{int(b.get_height())}/15", ha="center", va="bottom", fontsize=9)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel("Synthetic controls retained (of 15)")
ax.set_ylim(0, 17)
ax.set_title("C  Borderline controls, true ROR 1.4–2.5", loc="left", fontweight="bold", fontsize=11)
ax.legend(frameon=False, fontsize=7.5, loc="upper center")

fig.suptitle("Simulation-based validation of empirical calibration (20,000 replicates, seed 20261121)",
             fontsize=12, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig("figures/Fig6_simulation.png", dpi=600, bbox_inches="tight")
print("saved figures/Fig6_simulation.png")
