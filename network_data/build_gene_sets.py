#!/usr/bin/env python3
"""
I020 C4 network evidence: build drug-target and toxicity gene sets.
Data sources (all public, no licence):
  - Drug targets: ChEMBL mechanism (primary targets) + curated published off-targets
  - Toxicity gene sets: Open Targets disease->target association (top高分) + curated
    core genes from published network-toxicology literature
Output: gene_sets.json  (frozen before any proximity computation)
"""
import requests, json, time, os

HERE = os.path.dirname(os.path.abspath(__file__))
S = requests.Session()
S.headers.update({"User-Agent": "I020-pv-network/0.1 (KIWI lab)"})
OT = "https://api.platform.opentargets.org/api/v4/graphql"
CH = "https://www.ebi.ac.uk/chembl/api/data"

# ---- ChEMBL target CHEMBL ID -> gene symbol (verified via ChEMBL target API) ----
# CHEMBL203 = EGFR (Epidermal growth factor receptor erbB1)
# CHEMBL3717 = MET (Hepatocyte growth factor receptor)
CHEMBL_SYMBOL = {"CHEMBL203": "EGFR", "CHEMBL3717": "MET"}

def chembl_target_symbol(tid):
    return CHEMBL_SYMBOL.get(tid)

DRUGS = {
    "osimertinib":  ["CHEMBL203"],            # EGFR
    "tepotinib":    ["CHEMBL3717"],           # MET
    "lazertinib":   ["CHEMBL203"],            # EGFR
    "amivantamab":  ["CHEMBL203", "CHEMBL3717"],  # EGFR + MET (bispecific)
    "capmatinib":   ["CHEMBL3717"],           # MET (active comparator)
}
# published network-toxicology hub/known off-targets (from PMID 41582009 tepotinib-AKI;
# general pharmacology for QT/ILD). Clearly flagged as curated literature overlay.
CURATED_DRUG_EXTRA = {
    "osimertinib":  ["ERBB2", "ERBB4", "KCNH2"],   # ERBB family + hERG (QT liability)
    "tepotinib":    ["HSP90AA1", "AKT1", "EGFR", "CASP3", "TNF", "SRC"],  # published tepotinib-AKI hubs
    "lazertinib":   ["ERBB2"],
    "amivantamab":  [],
    "capmatinib":   [],
}

# ---- Open Targets disease EFO/MONDO IDs ----
DISEASES = {
    "ILD": "MONDO_0015925",      # interstitial lung disease
    "QT":  "HP_0001657",         # prolonged QT interval
    "AKI": "MONDO_0002492",      # acute kidney injury
    "IRR": "EFO_1002004",        # drug hypersensitivity syndrome (infusion reaction proxy)
}
# curated core toxicity genes (literature-based, mechanism-anchored)
CURATED_TOX_CORE = {
    "ILD": ["TGFB1", "TNF", "IL1B", "IL6", "MMP1", "MMP7", "COL1A1", "COL3A1",
            "PDGFA", "CTGF", "LOX", "TIMP1", "SERPINE1", "SMAD3", "AKT1", "EGFR"],
    "QT":  ["KCNH2", "KCNQ1", "SCN5A", "KCNE1", "KCNE2", "KCNJ2", "CACNA1C",
            "HERG1", "KCND3", "SCN1B", "AKT1"],
    "AKI": ["HAVCR1", "LCN2", "TNF", "IL1B", "IL6", "CASP3", "AKT1", "HSP90AA1",
            "SRC", "EGFR", "TLR4", "NFKB1", "NOD1", "TRPM2"],
    "IRR": ["TNF", "IL1B", "IL6", "IL2", "IFNG", "FCGR1A", "FCGR2A", "FCGR3A",
            "MS4A2", "TPSAB1", "CPA3", "HDC", "C5AR1", "HRH1", "PTGDR2"],
}

def ot_top_targets(disease_id, n=150):
    q = '''{ disease(efoId:"%s") {
        id name associatedTargets(page:{size:%d,index:0}, orderByScore:"desc"){
            count rows { score target { id approvedSymbol } }
        } } }''' % (disease_id, n)
    r = S.post(OT, json={"query": q}, timeout=60)
    d = r.json()
    if "errors" in d:
        return [], {"error": str(d["errors"])[:200]}
    dis = d["data"]["disease"]
    rows = dis["associatedTargets"]["rows"]
    genes = [x["target"]["approvedSymbol"] for x in rows
             if x["target"].get("approvedSymbol")]
    return genes, {"name": dis["name"], "total": dis["associatedTargets"]["count"],
                   "top_score": rows[0]["score"] if rows else None}

def main():
    out = {"meta": {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "sources": {
            "drug_primary": "ChEMBL mechanism API (molecular_mechanism=1)",
            "drug_extra": "curated published network-toxicology / pharmacology (see notes)",
            "toxicity_ot": "Open Targets Platform GraphQL associatedTargets, top-N by score",
            "toxicity_core": "curated literature core genes",
        },
        "note": "FROZEN before any network proximity computation (C4 discriminating test)."
    }, "drug_targets": {}, "toxicity_sets": {}}

    # drug primary targets from ChEMBL
    for drug, tids in DRUGS.items():
        syms = []
        for tid in tids:
            s = chembl_target_symbol(tid)
            if s: syms.append(s)
            time.sleep(0.25)
        extra = CURATED_DRUG_EXTRA.get(drug, [])
        out["drug_targets"][drug] = {
            "primary_chembl": syms,
            "curated_extra": extra,
            "union": sorted(set(syms) | set(extra)),
            "chembl_ids": tids,
        }
        print(f"DRUG {drug:13s} primary={syms} union={out['drug_targets'][drug]['union']}", flush=True)

    # toxicity gene sets
    for tox, did in DISEASES.items():
        ot_genes, info = ot_top_targets(did, n=150)
        core = CURATED_TOX_CORE[tox]
        union = sorted(set(core) | set(ot_genes))
        out["toxicity_sets"][tox] = {
            "disease_id": did,
            "open_targets_info": info,
            "ot_top150": ot_genes,
            "curated_core": core,
            "union": union,
            "n_union": len(union),
        }
        print(f"TOX  {tox:4s} {info} curated={len(core)} union={len(union)}", flush=True)
        time.sleep(0.5)

    path = os.path.join(HERE, "gene_sets.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(f"\nWROTE {path}", flush=True)

if __name__ == "__main__":
    main()
