#!/usr/bin/env python3
"""
C4 sensitivity: test whether the single significant proximity (tepotinib-AKI)
is an artifact of CIRCULAR hub inclusion.

Scenarios:
  S1 union_targets  + union_tox     (baseline, contains published hub overlap)
  S2 primary_targets + union_tox    (remove target-side circular hubs)
  S3 union_targets  + clean_tox     (remove toxicity-side promiscuous hubs)
  S4 primary_targets + clean_tox    (double-clean, most stringent)

If tepotinib-AKI loses significance in S2/S3/S4 -> proximity result is circular.
"""
import json, os, math
import numpy as np
import networkx as nx
from network_proximity import (
    fetch_string_network, largest_component, build_shortest_paths,
    proximity, degree_bins, degree_matched_null, PAIRS, N_PERM,
)

HERE = os.path.dirname(os.path.abspath(__file__))

# clean organ-specific toxicity cores: remove promiscuous signalling hubs that were
# borrowed from the very network-toxicology papers we are trying to test.
CLEAN_TOX = {
    "ILD": ["TGFB1", "MMP1", "MMP7", "COL1A1", "COL3A1", "PDGFA", "CTGF",
            "LOX", "TIMP1", "SERPINE1", "SMAD3"],
    "QT":  ["KCNH2", "KCNQ1", "SCN5A", "KCNE1", "KCNE2", "KCNJ2", "CACNA1C",
            "KCND3", "SCN1B"],
    "AKI": ["HAVCR1", "LCN2", "SLC34A1", "LRP2", "CUBN", "AQP1", "UMOD",
            "TLR4", "NOD1", "TRPM2"],
    "IRR": ["FCGR1A", "FCGR2A", "FCGR3A", "MS4A2", "TPSAB1", "CPA3", "HDC",
            "C5AR1", "HRH1", "PTGDR2"],
}

def main():
    gs = json.load(open(os.path.join(HERE, "gene_sets.json"), encoding="utf-8"))
    allgenes = set()
    for d in gs["drug_targets"].values(): allgenes |= set(d["union"])
    for t in gs["toxicity_sets"].values(): allgenes |= set(t["union"])
    for v in CLEAN_TOX.values(): allgenes |= set(v)
    allgenes = sorted(allgenes)

    G, n_edges = fetch_string_network(allgenes, 0.4)
    G = largest_component(G)
    inf_val = nx.diameter(G) + 1
    sp = build_shortest_paths(G)
    binof, deg = degree_bins(G, G.nodes())
    print(f"network nodes={len(G)} edges={n_edges} diameter={inf_val-1}", flush=True)

    scenarios = {
        "S1_union_union":    ("union", "union"),
        "S2_primary_union":  ("primary", "union"),
        "S3_union_clean":    ("union", "clean"),
        "S4_primary_clean":  ("primary", "clean"),
    }
    out = {}
    for sname, (tmode, dmode) in scenarios.items():
        print(f"\n--- {sname}: targets={tmode}, tox={dmode} ---", flush=True)
        rows = []
        for label, drug, outcome, group in PAIRS:
            tgt = (gs["drug_targets"][drug]["primary_chembl"] if tmode == "primary"
                   else gs["drug_targets"][drug]["union"])
            dis = (CLEAN_TOX[outcome] if dmode == "clean"
                   else gs["toxicity_sets"][outcome]["union"])
            d_obs, nt, nd = proximity(tgt, dis, sp, inf_val)
            null = degree_matched_null(tgt, dis, G, sp, inf_val, binof, N_PERM)
            z = (d_obs - null.mean()) / null.std(ddof=1)
            p = float((null <= d_obs).mean())
            rows.append({"pair": label, "group": group, "d": d_obs, "z": float(z), "p": p})
            print(f"  {label:18s} d={d_obs:.3f} z={z:6.2f} p={p:.4f} [{group}]", flush=True)
        out[sname] = rows

    path = os.path.join(HERE, "sensitivity_c4.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nWROTE {path}", flush=True)

if __name__ == "__main__":
    main()
