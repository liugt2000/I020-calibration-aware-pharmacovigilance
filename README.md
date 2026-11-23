# One agent-specific signal, three class effects, two unresolved — analysis deposit

**Version 1.0 · 2026-11-23**
Companion deposit to: *One agent-specific signal, three class effects, two unresolved:
calibration-aware pharmacovigilance in NSCLC targeted therapy* (Journal of Translational
Medicine, Special Issue on Computational Modelling and Network Medicine in Drug Toxicology
and Clinical Pharmacovigilance).

This deposit contains the complete, version-locked analysis pipeline behind every number in
the manuscript: openFDA/FAERS aggregate queries, the four disproportionality algorithms,
empirical-null calibration, active-comparator re-estimation, the Monte-Carlo and
misspecification simulations, the notoriety time series, the network-proximity
discrimination test, and the decision-impact model.

## Data sources (all public; no patient-level data)

| Source | Use | Locked version / date |
|---|---|---|
| openFDA drug-event API (`api.fda.gov/drug/event.json`) | FAERS aggregate 2x2 counts | snapshot 2026-11-21; 21,054,159 reports |
| ChEMBL API | drug targets (molecular_mechanism = 1) | ChEMBL_37 (2026-05-01) |
| Open Targets Platform GraphQL API | disease-target associations | v26.9.0, data release 2026-09 |
| STRING REST API | induced interactome | v12.5; combined score >= 0.4 (primary), >= 0.7 (sensitivity) |

## Environment

Python 3.13.14; pandas 2.3.3, numpy 2.5.4, scipy 1.18.1, networkx 3.7, matplotlib 3.11.2,
requests 2.34.2. See `environment_lock.txt`. Fixed random seeds: 20261121 (calibration
simulation) and 20261121 (degree-matched network permutation). Two-sided alpha = 0.05.

## Reproduction

Run from the repository root. The `fetch_*` scripts query live APIs and require network
access; all downstream steps run on the deposited CSVs and are deterministic.

```
python fetch_counts.py             # -> counts.csv            (index pairs, whole database)
python fetch_active_comparators.py # -> counts_active.csv     (active comparators)
python fetch_suspect.py            # -> counts_suspect.csv, counts_active_suspect.csv
python analyze.py                  # -> results_index.csv, null_params.json, fill_values.json
python analyze_active.py           # -> results_active.csv
python compare_suspect.py          # -> suspect_*.csv, suspect_comparison.json
python pc_simulation.py            # -> simulation_results.json   (Fig. 6)
python pc_misspec.py               # -> misspec_results.json      (Supplementary Fig. S2)
python notoriety.py                # -> notoriety_results.json    (Fig. 4)
cd network_data && python build_gene_sets.py && python network_proximity.py     && python sensitivity_c4.py     # -> network_results.json, sensitivity_c4.json (Fig. 5)
cd .. && python decision_impact.py # -> decision_impact.json      (Fig. 7)
python make_figures_2.py && python make_fig_sim.py     && python make_fig_decision.py && python make_fig_misspec.py  # 600 dpi figures
```

Because openFDA exposes no fixed closing quarter and no MedDRA version, re-running the
`fetch_*` scripts today will not reproduce the deposited counts exactly. The deposited
`counts*.csv` files are the authoritative snapshot for the manuscript; cite them, not a
fresh query.

## Directory layout

```
*.py                     analysis code, flat so that relative paths work as in the paper
config.json              drug, outcome and negative-control definitions
counts*.csv              openFDA aggregate counts (2x2 cells)
results_*.csv            disproportionality, calibration and active-comparator results
null_params.json         fitted empirical-null parameters per drug and control panel
simulation_results.json  Monte-Carlo type-I error / power results
misspec_results.json     off-model misspecification stress test
notoriety_results.json   quarterly index-minus-control reporting fractions
network_data/            gene sets, network proximity, scenario sensitivity
figures/                 Fig. 1-7 and Supplementary Figs. S1-S2 (600 dpi PNG)
manuscript_admin/        reference-verification tooling and audit logs
RESULTS_SUMMARY.md       human-readable summary of every headline number
```

## Licence

Code: MIT (`LICENSE`). Derived data and figures: CC BY 4.0. The underlying FAERS records
remain subject to the FDA openFDA terms of use.

## Contact

Corresponding author (details in the manuscript).
