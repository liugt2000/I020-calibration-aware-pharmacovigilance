#!/usr/bin/env python3
"""
I020 C4 — Network evidence as a DISCRIMINATING test.

Method (Zhou et al., Nat Commun 2020 network proximity; Guney et al. 2016):
  For drug target set T and toxicity gene set D on a PPI graph G:
    d(s,D) = min_{t in D} shortest_path(s,t)   for each s in T
    proximity d_c(T,D) = mean_s d(s,D)          (smaller = closer)
  Significance vs degree-matched random gene sets (10,000 perms):
    z = (d_obs - mean_null)/sd_null             (more negative = closer than chance)
    p_emp = P(d_null <= d_obs)

Discriminating test (prespecified C4):
  robust pairs (survive calibration + active comparator) vs attenuated pairs.
  With n=3 vs 3, Mann-Whitney two-sided min p = 0.1 -> UNDERPOWERED.
  Use EXACT label permutation over C(6,3)=20 assignments (one-sided) instead.

Network: STRING induced subnetwork on union of all gene sets (combined score >= 0.4,
primary; >=0.7 sensitivity), fetched live via STRING API for reproducibility.
"""
import requests, json, time, itertools, math, os
import numpy as np
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20261121)  # fixed seed, reproducible
S = requests.Session(); S.headers.update({"User-Agent": "I020-pv-network/0.1"})

# ---- index pairs: (label, drug, outcome_key, group) — fixed by calibration results ----
PAIRS = [
    ("osimertinib_ILD", "osimertinib", "ILD", "attenuated"),
    ("osimertinib_QT",  "osimertinib", "QT",  "robust"),
    ("tepotinib_AKI",   "tepotinib",   "AKI", "robust"),
    ("lazertinib_ILD",  "lazertinib",  "ILD", "attenuated"),
    ("amivantamab_IRR", "amivantamab", "IRR", "robust"),
    ("amivantamab_ILD", "amivantamab", "ILD", "attenuated"),
]
N_PERM = 10000


def fetch_string_network(genes, min_score):
    """Induced STRING network on gene symbols; returns nx.Graph weighted by score."""
    r = S.get("https://string-db.org/api/json/network",
              params={"identifiers": "\r".join(genes), "species": 9606,
                      "required_score": int(min_score * 1000)}, timeout=120)
    r.raise_for_status()
    edges = r.json()
    G = nx.Graph()
    for e in edges:
        a, b = e["preferredName_A"], e["preferredName_B"]
        G.add_edge(a, b, score=float(e["score"]))
    return G, len(edges)


def largest_component(G):
    if len(G) == 0:
        return G
    return G.subgraph(max(nx.connected_components(G), key=len)).copy()


def build_shortest_paths(G):
    return dict(nx.all_pairs_shortest_path_length(G))


def proximity(T, D, sp, inf_val):
    """d_c = mean over s in T of min shortest path s->D."""
    T = [g for g in T if g in sp]; D = [g for g in D if g in sp]
    if not T or not D:
        return None, 0, 0
    ds = []
    for s in T:
        dmin = min((sp[s].get(t, inf_val) for t in D), default=inf_val)
        ds.append(dmin)
    return float(np.mean(ds)), len(T), len(D)


def degree_bins(G, nodes):
    """Map node -> degree bin (quintiles of log-degree) for degree matching."""
    deg = dict(G.degree())
    vals = np.array([math.log1p(deg[n]) for n in G.nodes()])
    qs = np.quantile(vals, np.linspace(0, 1, 11))
    def binof(n):
        v = math.log1p(deg.get(n, 0))
        return int(np.searchsorted(qs, v, side="right") - 1)
    return binof, deg


def degree_matched_null(T, D, G, sp, inf_val, binof, n_perm):
    """Sample toxicity-like sets matched to D's degree distribution; compute d_c null."""
    nodes = list(G.nodes())
    by_bin = {}
    for n in nodes:
        by_bin.setdefault(binof(n), []).append(n)
    d_bins = [binof(g) for g in D if g in G]
    if not d_bins:
        d_bins = [binof(g) for g in nodes[:len(D)]]
    null = []
    for _ in range(n_perm):
        sample = set()
        for b in d_bins:
            pool = by_bin.get(b) or nodes
            sample.add(pool[rng.integers(len(pool))])
        d, _, _ = proximity(T, list(sample), sp, inf_val)
        if d is not None:
            null.append(d)
    return np.array(null)


def exact_label_perm_test(z_by_pair, group_by_pair):
    """Exact one-sided permutation p for mean(z)_robust < mean(z)_attenuated."""
    pairs = list(z_by_pair.keys())
    z = np.array([z_by_pair[p] for p in pairs])
    grp = np.array([group_by_pair[p] for p in pairs])
    obs = z[grp == "robust"].mean() - z[grp == "attenuated"].mean()
    k = (grp == "robust").sum()
    counts = 0; total = 0
    for combo in itertools.combinations(range(len(pairs)), k):
        mask = np.zeros(len(pairs), dtype=bool); mask[list(combo)] = True
        stat = z[mask].mean() - z[~mask].mean()
        counts += int(stat <= obs + 1e-12)
        total += 1
    return obs, counts / total, total


def main():
    gs = json.load(open(os.path.join(HERE, "gene_sets.json"), encoding="utf-8"))
    allgenes = set()
    for d in gs["drug_targets"].values(): allgenes |= set(d["union"])
    for t in gs["toxicity_sets"].values(): allgenes |= set(t["union"])
    allgenes = sorted(allgenes)

    results = {}
    for score_thr in [0.4, 0.7]:
        tag = "score040" if score_thr == 0.4 else "score070"
        G, n_edges = fetch_string_network(allgenes, score_thr)
        G = largest_component(G)
        inf_val = nx.diameter(G) + 1 if len(G) > 1 else 10
        sp = build_shortest_paths(G)
        binof, deg = degree_bins(G, G.nodes())
        print(f"\n[{tag}] nodes={len(G)} edges={n_edges} diameter={inf_val-1}", flush=True)

        recs = []
        for label, drug, outcome, group in PAIRS:
            T = gs["drug_targets"][drug]["union"]
            D = gs["toxicity_sets"][outcome]["union"]
            d_obs, nt, nd = proximity(T, D, sp, inf_val)
            null = degree_matched_null(T, D, G, sp, inf_val, binof, N_PERM)
            if len(null) < 30:
                z = None; p_emp = None
            else:
                z = (d_obs - null.mean()) / null.std(ddof=1)
                p_emp = float((null <= d_obs).mean())
            rec = {"pair": label, "drug": drug, "outcome": outcome, "group": group,
                   "n_targets": nt, "n_tox_genes": nd, "d_obs": d_obs,
                   "null_mean": float(null.mean()) if len(null) else None,
                   "null_sd": float(null.std(ddof=1)) if len(null) else None,
                   "z": None if z is None else float(z), "p_emp": p_emp,
                   "n_perm_valid": len(null)}
            recs.append(rec)
            print(f"  {label:18s} d={d_obs:.3f} null={null.mean():.3f}±{null.std(ddof=1):.3f} "
                  f"z={rec['z']:.2f} p={p_emp:.4f} [{group}]", flush=True)

        zmap = {r["pair"]: r["z"] for r in recs if r["z"] is not None}
        gmap = {r["pair"]: r["group"] for r in recs if r["z"] is not None}
        if len(zmap) == 6:
            obs_diff, p_perm, n_perm_label = exact_label_perm_test(zmap, gmap)
            print(f"  DISCRIMINATING robust-attenuated mean z diff = {obs_diff:.3f}, "
                  f"exact perm p(one-sided) = {p_perm:.4f} ({n_perm_label} perms)", flush=True)
        else:
            obs_diff, p_perm, n_perm_label = None, None, None
        results[tag] = {"network": {"nodes": len(G), "edges": n_edges,
                                    "score_threshold": score_thr, "diameter": inf_val-1},
                        "records": recs,
                        "discriminating": {"mean_z_diff_robust_minus_atten": obs_diff,
                                           "p_one_sided_exact_perm": p_perm,
                                           "n_label_perms": n_perm_label}}

    out = os.path.join(HERE, "network_results.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump({"meta": {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
                            "method": "Zhou 2020 network proximity + degree-matched perm (10k)",
                            "seed": 20261121, "n_perm": N_PERM},
                   "results": results}, f, indent=2)
    print(f"\nWROTE {out}", flush=True)


if __name__ == "__main__":
    main()
