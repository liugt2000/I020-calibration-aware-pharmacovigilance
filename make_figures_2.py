#!/usr/bin/env python3
"""Fig.1 study-flow and Fig.6 READUS-PV compliance map."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE=os.path.dirname(os.path.abspath(__file__)); FIG=os.path.join(HERE,"figures")
plt.rcParams.update({"font.family":"DejaVu Sans"})

# ---------- Fig 1: flow ----------
fig,ax=plt.subplots(figsize=(9,10)); ax.set_xlim(0,10); ax.set_ylim(0,18); ax.axis("off")
ax.set_title("Figure 1. Study flow and negative-control definition",fontsize=12,fontweight="bold",loc="left")

def box(x,y,w,h,text,fc="#EAF1FB",ec="#2E6F8E",fs=8.5):
    p=FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.02,rounding_size=0.06",
                     linewidth=1,edgecolor=ec,facecolor=fc)
    ax.add_patch(p); ax.text(x+w/2,y+h/2,text,ha="center",va="center",fontsize=fs)
    return (x+w/2,y)

def arrow(p1,p2):
    ax.add_patch(FancyArrowPatch(p1,p2,arrowstyle="-|>",mutation_scale=12,
                                 linewidth=1.1,color="#5F6B72"))

c1=box(2.5,16.4,5,1.0,"openFDA FAERS interface\n21,054,159 reports (snapshot 2026-11-21)")
c2=box(2.5,14.6,5,1.0,"Count extraction\n4 index drugs · 6 index pairs")
arrow(box(5,16.4,0,0),box(5,15.6,0,0)) if False else None
arrow((5,16.4),(5,15.6))
# branch: negative controls
nc=box(0.3,12.6,4.3,1.3,"Negative-control adjudication\nN1 inert / N2 contextual / N3 temporal\nlabel + RCT + meta decontamination",fc="#F1EFE8",ec="#888780",fs=8)
idxb=box(5.4,12.6,4.3,1.3,"Index pairs (prespecified)\nosimertinib ILD/QT · tepotinib AKI\nlazertinib ILD · amivantamab IRR/ILD",fc="#EAF1FB")
arrow((3.6,14.6),(2.5,13.9)); arrow((7.5,14.6),(7.5,13.9))
cal=box(2.0,10.4,6,1.2,"Empirical calibration  ln(ROR)_i ~ N(mu, sigma^2 + SE_i^2)\nML primary + sample-SD sensitivity · BH-FDR · min 10 usable controls",fc="#E1F5EE",ec="#2E7D5B",fs=8.3)
arrow((2.45,12.6),(4.0,11.6)); arrow((7.55,12.6),(6.0,11.6))
comp=box(2.0,8.6,6,1.0,"Active-comparator re-estimation\n(same indication, different agent)",fc="#E1F5EE",ec="#2E7D5B")
arrow((5,10.4),(5,9.6))
# parallel validations
t1=box(0.3,6.6,4.3,1.2,"Temporal / notoriety analysis\nquarterly index-minus-control fraction\nchange-point vs regulatory anchors",fc="#FAEEDA",ec="#C77716",fs=8)
t2=box(5.4,6.6,4.3,1.2,"Network discriminating test\nChEMBL targets + Open Targets tox sets\nSTRING interactome · 10k degree-matched perms",fc="#FAEEDA",ec="#C77716",fs=8)
arrow((3.5,8.6),(2.45,7.8)); arrow((6.5,8.6),(7.55,7.8))
tier=box(2.0,4.4,6,1.2,"Cross-source tiering (network excluded after negative result)\nTier A: retained (n=2) · Tier B: class effects (n=3) · Tier C: uncalibratable (n=1)",fc="#EEEDFE",ec="#534AB7",fs=8.2)
arrow((2.45,6.6),(4.0,5.6)); arrow((7.55,6.6),(6.0,5.6))
rpt=box(2.0,2.4,6,0.9,"READUS-PV reporting · frozen definitions + hash chain",fc="#F1EFE8",ec="#888780",fs=8.5)
arrow((5,4.4),(5,3.3))
ax.text(5,1.2,"Prespecified plan frozen before outcome access (timestamped deposit + SHA-256 chain).",
        ha="center",fontsize=7.5,style="italic",color="#666")
ax.text(5,0.6,"C4 (network discrimination) and C5 (notoriety) were NOT confirmed and are reported as negative results.",
        ha="center",fontsize=7.5,color="#B23A48")
fig.savefig(os.path.join(FIG,"Fig1_flow.png"),dpi=600,bbox_inches="tight"); plt.close(fig)

# ---------- Fig 6: READUS-PV compliance map ----------
fig,ax=plt.subplots(figsize=(8,9)); ax.set_xlim(0,10); ax.set_ylim(0,15); ax.axis("off")
ax.set_title("Figure 6. READUS-PV compliance map",fontsize=12,fontweight="bold",loc="left")
items=[
 ("Data source, country, period, version","§2.2"),
 ("Deduplication method","§2.2 / note"),
 ("Drug coding (PS vs all roles) & dictionary","§2.2 / §3.2 / S2b"),
 ("Event definition (MedDRA ver., PT list)","§2.2 / S2"),
 ("Four algorithm formulae & thresholds","§2.4"),
 ("Comparator definition (DB + active)","§2.3"),
 ("Multiplicity handling","§2.5"),
 ("Negative-control rationale & cleaning","§2.5 / Fig.1"),
 ("Sensitivity / robustness analyses","§2.5, §2.6, §3.6"),
 ("Time-at-onset & notoriety handling","§2.6 / §3.5"),
 ("Limitations (no causality/incidence)","Discussion"),
 ("Code / data availability","§2.8"),
]
cols=[0.4,6.4]; cw=[5.9,3.2]
ax.text(cols[0]+0.1,13.9,"READUS-PV item",fontsize=9,fontweight="bold")
ax.text(cols[1]+0.1,13.9,"Manuscript location",fontsize=9,fontweight="bold")
y=13.3
for i,(item,loc) in enumerate(items):
    fc="#F7F9FB" if i%2==0 else "white"
    ax.add_patch(plt.Rectangle((cols[0],y-0.45),9.3,0.9,facecolor=fc,edgecolor="none"))
    ax.text(cols[0]+0.15,y,item,fontsize=8.3,va="center")
    ax.add_patch(plt.Rectangle((cols[1],y-0.3),cw[1],0.6,facecolor="#E1F5EE",edgecolor="#2E7D5B",linewidth=0.6))
    ax.text(cols[1]+cw[1]/2,y,loc,fontsize=8,ha="center",va="center",color="#1B5E40")
    y-=0.93
ax.text(5,0.9,"Full item-by-item checklist: Supplementary Table S10.",ha="center",fontsize=8,style="italic",color="#666")
ax.text(5,0.35,"Note: suspect-role restriction tested as sensitivity analysis (§3.2, <2.4% change); case-level deduplication not possible with aggregate counts — stated as a limitation.",
        ha="center",fontsize=6.8,color="#B23A48")
fig.savefig(os.path.join(FIG,"SuppFig_S1_readuspv.png"),dpi=600,bbox_inches="tight"); plt.close(fig)

print("Wrote Fig1_flow.png, SuppFig_S1_readuspv.png")
for f in ["Fig1_flow.png","SuppFig_S1_readuspv.png"]:
    print(" ",f,os.path.getsize(os.path.join(FIG,f)),"bytes")
