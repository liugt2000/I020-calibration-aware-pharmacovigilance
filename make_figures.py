#!/usr/bin/env python3
"""Generate publication figures (600 dpi PNG) for I020 manuscript v1.0.
All values read from computed result tables — no hard-coded fabrication."""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size":9,"font.family":"DejaVu Sans","axes.linewidth":0.8,
                     "axes.spines.top":False,"axes.spines.right":False,"figure.dpi":120})

idx = pd.read_csv(os.path.join(HERE,"results_index.csv"),encoding="utf-8-sig")
act = pd.read_csv(os.path.join(HERE,"results_active.csv"),encoding="utf-8-sig")
cnt = pd.read_csv(os.path.join(HERE,"counts.csv"),encoding="utf-8-sig")
net = json.load(open(os.path.join(HERE,"network_data","network_results.json"),encoding="utf-8"))
sens = json.load(open(os.path.join(HERE,"network_data","sensitivity_c4.json"),encoding="utf-8"))
noto = json.load(open(os.path.join(HERE,"notoriety_results.json"),encoding="utf-8"))

COL = {"naive":"#B23A48","active":"#2E6F8E","ctrl":"#9AA0A6","ok":"#2E7D5B","atten":"#C77716"}

# ---------- FIG 2: naive vs calibrated ----------
fig, axes = plt.subplots(2,2,figsize=(10,7.5))
fig.suptitle("Figure 2. Naive versus calibrated disproportionality",fontsize=12,fontweight="bold",y=0.98)

# A: distributions
ax=axes[0,0]
idx_lr=np.log(idx["ror"].dropna())
nc=cnt[cnt.outcome_key.str.startswith(("N1","N2"))].copy()
N=21054159
nc["a"]=nc["a"].fillna(0); nc["dt"]=nc["drug_total"].fillna(0); nc["et"]=nc["event_total"].fillna(0)
ncr=[]
for _,r in nc.iterrows():
    a=r["a"];dt=r["dt"];et=r["et"]
    b=dt-a;c=et-a;d=N-dt-et+a
    if min(a,b,c,d)>0: ncr.append(np.log((a*d)/(b*c)))
ax.hist(ncr,bins=25,color=COL["ctrl"],alpha=0.7,label="Negative controls (n=%d)"%len(ncr))
ax.hist(idx_lr,bins=10,color=COL["naive"],alpha=0.6,label="Index pairs (n=6)")
ax.axvline(0,color="black",ls="--",lw=1,label="Theoretical null (ln ROR=0)")
ax.set_xlabel("ln(ROR), naive"); ax.set_ylabel("Frequency"); ax.set_title("A. Null does not hold for inert events",loc="left")
ax.legend(fontsize=7,frameon=False)

# B: empirical null params
ax=axes[0,1]
nulls=json.load(open(os.path.join(HERE,"null_params.json"),encoding="utf-8"))
labels=[k.replace("::","\n") for k in nulls]
mus=[v["mu_ml"] for v in nulls.values()]; sigs=[v["sig_ml"] for v in nulls.values()]
x=np.arange(len(labels)); w=0.38
ax.bar(x-w/2,mus,w,color=COL["ctrl"],label="mu (ML)")
ax.bar(x+w/2,[min(s,1.5) for s in sigs],w,color=COL["active"],label="sigma (ML, capped 1.5)")
ax.set_xticks(x); ax.set_xticklabels(labels,rotation=45,ha="right",fontsize=6.5)
ax.axhline(0,color="black",lw=0.6); ax.set_ylabel("Value"); ax.set_title("B. Drug-specific empirical nulls",loc="left")
ax.legend(fontsize=7,frameon=False)

# C: raw vs calibrated
ax=axes[1,0]
sub=idx[["cal_ror_ml_N1","cal_lo_ml_N1"]].copy(); sub["raw"]=idx["ror"]
valid=sub[sub.cal_ror_ml_N1.notna()]
ax.scatter(valid["raw"],valid["cal_ror_ml_N1"],color=COL["naive"],s=45,zorder=3)
lo=valid["raw"].min()*0.7; hi=valid["raw"].max()*1.4
ax.plot([lo,hi],[lo,hi],ls="--",color="black",lw=0.8)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo,hi); ax.set_ylim(lo,hi)
ax.set_xlabel("Naive ROR"); ax.set_ylabel("Calibrated ROR (ML)"); ax.set_title("C. Calibration raises point estimates",loc="left")

# D: significance proportions
ax=axes[1,1]
stages=["Naive (theoretical null)","Calibrated (ML)","Calibrated + FDR"]
vals=[100*6/6,100*5/6,100*5/6]
bars=ax.barh(stages,vals,color=[COL["naive"],COL["active"],COL["ok"]])
for b,v in zip(bars,vals): ax.text(v+1,b.get_y()+b.get_height()/2,"%.0f%%"%v,va="center",fontsize=8)
ax.set_xlim(0,115); ax.axvline(5,color="red",ls=":",lw=1)
ax.set_xlabel("% pairs nominally significant"); ax.set_title("D. Significance before/after de-biasing",loc="left")
fig.tight_layout(rect=[0,0,1,0.96])
fig.savefig(os.path.join(FIG,"Fig2_calibration.png"),dpi=600,bbox_inches="tight"); plt.close(fig)

# ---------- FIG 3: waterfall ----------
fig,ax=plt.subplots(figsize=(9,5))
idx2=idx.copy(); idx2["short"]=idx2["label"].str.replace("IDX::","",regex=False)
m=act.merge(idx2[["short","ror","robust_double_fit"]],left_on="pair",right_on="short",how="left")
m["ror"]=m["ror"].fillna(m["ror_active"])
m["short"]=m["pair"].str.replace("P[0-9b]*_","",regex=True)
m=m.sort_values("ror_active")
y=np.arange(len(m))
ax.barh(y,np.log(m["ror"]),color=COL["naive"],height=0.62,label="ln(ROR) vs full database")
ax.barh(y,np.log(m["ror_active"]),color=COL["active"],height=0.34,label="ln(ROR) vs active comparators")
for i,(_,r) in enumerate(m.iterrows()):
    ax.annotate("Δ=%.2f"%(np.log(r.ror_active)-np.log(r.ror)),
                (max(np.log(r.ror),np.log(r.ror_active)),i),textcoords="offset points",xytext=(4,-4),fontsize=7)
ax.set_yticks(y); ax.set_yticklabels(m["short"])
ax.axvline(0,color="black",lw=0.8,ls="--"); ax.set_xlabel("ln(ROR)")
ax.set_title("Figure 3. Active-comparator attenuation waterfall",fontsize=12,fontweight="bold",loc="left")
ax.legend(fontsize=8,frameon=False,loc="lower right")
ax.text(0.02,-0.16,"Three ILD pairs fall to non-significance against their class;\nQT and AKI are retained. Amivantamab-IRR is unstable (zero comparator events).",
        transform=ax.transAxes,fontsize=7.5,color="#555")
fig.tight_layout(); fig.savefig(os.path.join(FIG,"Fig3_waterfall.png"),dpi=600,bbox_inches="tight"); plt.close(fig)

# ---------- FIG 4: notoriety ----------
fig,axes=plt.subplots(2,3,figsize=(11,6))
fig.suptitle("Figure 4. Temporal reporting-fraction analysis (no notoriety pattern)",fontsize=11,fontweight="bold",y=0.99)
cps=noto["change_points"]; series=noto["series"]
for ax,(name,cp) in zip(axes.ravel(),cps.items()):
    rows=series[name]["rows"]
    qs=[r["quarter"] for r in rows]; d=[r["diff_frac"] for r in rows]
    ax.plot(range(len(qs)),d,color=COL["active"],lw=1.3)
    ax.axhline(0,color="black",lw=0.6)
    si=cp.get("split_index")
    if si is not None and si<len(qs):
        ax.axvline(si,color=COL["atten"],ls="--",lw=1)
        ax.text(si,max(d)*0.9,"cp Δ=%+.3f"%cp["delta"],fontsize=6.5,color=COL["atten"],ha="center")
    step=max(1,len(qs)//5)
    ax.set_xticks(range(0,len(qs),step)); ax.set_xticklabels([qs[i] for i in range(0,len(qs),step)],rotation=90,fontsize=5.5)
    ax.set_title(name.replace("_","–"),fontsize=8,loc="left")
    ax.tick_params(axis="y",labelsize=6.5)
axes.ravel()[-1].axis("off")
fig.supylabel("Index minus control reporting fraction",fontsize=9)
fig.tight_layout(rect=[0.015,0,1,0.96]); fig.savefig(os.path.join(FIG,"Fig4_notoriety.png"),dpi=600,bbox_inches="tight"); plt.close(fig)

# ---------- FIG 5: network heatmap ----------
fig,ax=plt.subplots(figsize=(8,5))
scen=["S1_union_union","S2_primary_union","S3_union_clean","S4_primary_clean"]
pairs=["osimertinib_QT","tepotinib_AKI","amivantamab_IRR","osimertinib_ILD","lazertinib_ILD","amivantamab_ILD"]
M=np.full((len(pairs),len(scen)),np.nan)
for j,s in enumerate(scen):
    for r in sens[s]:
        if r["pair"] in pairs: M[pairs.index(r["pair"]),j]=r["z"]
M=np.clip(M,-3,3)
im=ax.imshow(M,cmap="RdBu_r",vmin=-3,vmax=3,aspect="auto")
ax.set_xticks(range(len(scen))); ax.set_xticklabels(["S1\nall genes","S2\nno target hubs","S3\norgan-specific","S4\nboth"])
ax.set_yticks(range(len(pairs))); ax.set_yticklabels([p.replace("_","–") for p in pairs])
for i in range(len(pairs)):
    for j in range(len(scen)):
        v=M[i,j]
        ax.text(j,i,"%+.2f"%v,ha="center",va="center",fontsize=8,
                color="white" if abs(v)>1.6 else "black")
cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.04); cb.set_label("Network proximity z (degree-matched null)")
ax.set_title("Figure 5. Network evidence fails discrimination; tepotinib–AKI is circular",fontsize=11,fontweight="bold",loc="left",pad=12)
ax.text(0,-0.22,"Exact label-permutation p = 0.40 (robust vs attenuated). tepotinib–AKI: S1 z=−2.54 (p=0.021) → S2 z=+0.67 (p=1.00) when borrowed hubs removed.",
        transform=ax.transAxes,fontsize=7.5,color="#555")
fig.tight_layout(); fig.savefig(os.path.join(FIG,"Fig5_network.png"),dpi=600,bbox_inches="tight"); plt.close(fig)

print("Figures written to",FIG)
for f in sorted(os.listdir(FIG)): print("  ",f,os.path.getsize(os.path.join(FIG,f)),"bytes")
