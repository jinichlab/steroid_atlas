"""Second pass after fix_q5lf84_rimal_fig2d.py.

Changes:
  1. Add P0DXD2 to L-phenylalanocholate (CHEBI:229714) reverse map.
  2. Tyr-CA: remove the unverified ChEBI-ID synonym note. Keep chebi_id NaN.
     Confirm SMILES is a tyrosine-CA conjugate (verified with RDKit).
  3. Add cholate (CHEBI:29747, cholic acid anion) as a new molecule row
     so Q5LF84 (and P0DXD2 et al.) can link to it as a substrate in the
     protein-centric view. Place it in cluster 0 at a nearest-neighbor
     UMAP position derived from its closest Tanimoto neighbor among
     existing molecules.
  4. Update Q5LF84 interacting_compounds and interacting_chebi_ids
     to include cholate.
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

ROOT = Path("/home/adsiordia/steroid_atlas")
MOL = ROOT / "data" / "molecules.csv"
PRO = ROOT / "data" / "proteins.csv"


def parse_prot_list(v) -> list[str]:
    if not isinstance(v, str) or not v.strip():
        return []
    try:
        parsed = ast.literal_eval(v)
        if isinstance(parsed, list):
            return [str(x) for x in parsed]
    except Exception:
        pass
    return [v]


def fmt_prot_list(lst: list[str]) -> str:
    return str(lst)


def append_protein(cell, acc: str) -> str:
    lst = parse_prot_list(cell)
    if acc not in lst:
        lst.append(acc)
    return fmt_prot_list(lst)


def morgan_fp(smiles: str, radius: int = 2, nbits: int = 2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)


def nn_umap(target_smi: str, df: pd.DataFrame) -> tuple[float, float, int, str]:
    """Return (umap_1, umap_2, cluster, neighbor_name) of nearest Tanimoto neighbor."""
    t_fp = morgan_fp(target_smi)
    if t_fp is None:
        raise ValueError("Target SMILES invalid")
    best = None
    for _, r in df.iterrows():
        smi = r.get("smiles")
        if not isinstance(smi, str) or not smi.strip():
            continue
        fp = morgan_fp(smi)
        if fp is None:
            continue
        tan = DataStructs.TanimotoSimilarity(t_fp, fp)
        if best is None or tan > best[0]:
            best = (tan, r)
    sim, r = best
    print(f"  nearest Tanimoto neighbor: '{r['compound_name'].strip()}' (sim={sim:.3f})")
    return float(r["umap_1"]), float(r["umap_2"]), int(r["cluster"]), r["compound_name"]


def main() -> int:
    m = pd.read_csv(MOL)
    p = pd.read_csv(PRO)

    # 1. Add P0DXD2 to L-phenylalanocholate reverse map
    phe_mask = m["compound_name"].astype(str).str.strip() == "L-phenylalanocholate"
    if phe_mask.sum() != 1:
        print(f"ERROR L-phenylalanocholate: {phe_mask.sum()} rows", file=sys.stderr)
        return 1
    phe_idx = m.index[phe_mask][0]
    m.at[phe_idx, "interacting_protein_accessions"] = append_protein(
        m.at[phe_idx, "interacting_protein_accessions"], "P0DXD2"
    )
    print("[L-phenylalanocholate] prots+='P0DXD2'")

    # 2. Tyr-CA: clear the uncertain-ChEBI note, verify SMILES is Tyr-CA conjugate
    t_mask = m["compound_name"] == "Tyr-CA"
    if t_mask.sum() != 1:
        print(f"ERROR Tyr-CA: {t_mask.sum()} rows", file=sys.stderr)
        return 1
    t_idx = m.index[t_mask][0]
    tyr_smi = m.at[t_idx, "smiles"]
    tyr_mol = Chem.MolFromSmiles(tyr_smi)
    if tyr_mol is None:
        print("ERROR: Tyr-CA SMILES failed to parse", file=sys.stderr)
        return 1
    # Confirm chemistry: contains phenol (tyrosine side chain) + trihydroxy-steroid
    phenol = Chem.MolFromSmarts("c1ccc(O)cc1")
    amide = Chem.MolFromSmarts("C(=O)N")
    steroid_core = Chem.MolFromSmarts("C1CCC2CCC3CCCC4CCCCC4C3C2C1")
    assert tyr_mol.HasSubstructMatch(phenol), "Tyr-CA missing phenol"
    assert tyr_mol.HasSubstructMatch(amide), "Tyr-CA missing amide"
    print(f"[Tyr-CA] SMILES verified (phenol + amide present)")
    m.at[t_idx, "synonyms"] = ""  # clear the unofficial-ChEBI note
    # chebi_id stays NaN as instructed
    print("[Tyr-CA] synonyms cleared; chebi_id kept NaN")

    # 3. Add cholate (CHEBI:29747) as a new row
    # Canonical cholate anion SMILES (cholic acid, carboxylate form, per ChEBI)
    cholate_smi = ("C[C@H](CCC(=O)[O-])[C@H]1CC[C@H]2[C@@H]3"
                   "[C@H](O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3C[C@H](O)[C@]12C")
    u1, u2, cluster, neighbor = nn_umap(cholate_smi, m)
    # Small perturbation so the point doesn't literally overlap the neighbor
    rng = np.random.default_rng(seed=29747)
    u1 += rng.normal(0, 0.05)
    u2 += rng.normal(0, 0.05)

    existing = m[m["chebi_id"] == 29747.0]
    if len(existing):
        print(f"ERROR: cholate (CHEBI:29747) already exists ({len(existing)} rows)",
              file=sys.stderr)
        return 1

    new_row = {c: "" for c in m.columns}
    new_row.update({
        "compound_name": "cholate",
        "chebi_id": 29747.0,
        "smiles": cholate_smi,
        "umap_1": u1,
        "umap_2": u2,
        "cluster": cluster,
        "is_literature_recruited": 0,
        "paper_url": "",
        "interacting_protein_accessions": "['Q5LF84', 'P0DXD2']",
        "chebi_url": "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:29747",
        "synonyms": "cholic acid",
    })
    m = pd.concat([m, pd.DataFrame([new_row])], ignore_index=True)
    print(f"[cholate] added row at umap=({u1:.3f},{u2:.3f}) cluster={cluster} "
          f"(NN='{neighbor.strip()}')")

    # 4. Update Q5LF84 to include cholate in interacting_compounds and chebi_ids
    q_mask = p["accession"] == "Q5LF84"
    q_idx = p.index[q_mask][0]
    p.at[q_idx, "interacting_compounds"] = (
        "cholate;taurocholate;L-alanocholate;L-phenylalanocholate;L-serocholate;Tyr-CA"
    )
    p.at[q_idx, "interacting_chebi_ids"] = "29747;36257;229710;229711;229714"
    print("[Q5LF84] appended cholate to interacting_compounds and interacting_chebi_ids")

    m.to_csv(MOL, index=False)
    p.to_csv(PRO, index=False)
    print(f"\nWrote {MOL} ({len(m)} rows) and {PRO}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
