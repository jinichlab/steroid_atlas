"""Stereochemistry audit across all molecules in `data/molecules.csv`.

Non-destructive. Reports:
  1. True duplicates     — same canonical InChIKey, different compound_name
  2. Name conflicts      — same compound_name, different InChIKey
  3. Stereochem pairs    — known distinct pairs (CDCA/UDCA, LCA/isoLCA,
                           CA/alloCA, α/β/γ-MCA) must have different InChIKeys
  4. Ionization dupes    — same non-stereo InChIKey, different protonation
                           (anion/acid pair for the same scaffold)
  5. L/D amino-acid      — any amidate stem present in both L- and D- variants
  6. Undefined stereo    — molecules with ≥1 unassigned stereocenter

Output: `literature/results/stereo_audit.tsv` (per-molecule InChIKey +
formula + undefined-stereo count) and a console report.
"""
from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.inchi import MolToInchiKey

from chemistry_utils import MOLECULES_CSV, ensure_results_dir  # type: ignore[import-not-found]

csv.field_size_limit(sys.maxsize)
RDLogger.DisableLog("rdApp.*")


def analyze_row(name: str, smiles: str) -> dict[str, object] | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    canon_smi = Chem.MolToSmiles(mol, isomericSmiles=True)
    inchikey = MolToInchiKey(mol)
    flat = Chem.MolFromSmiles(canon_smi)
    Chem.RemoveStereochemistry(flat)
    flat_inchikey = MolToInchiKey(flat)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    stereo = Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                       useLegacyImplementation=False)
    defined = sum(1 for _, flag in stereo if flag != "?")
    undefined = sum(1 for _, flag in stereo if flag == "?")
    return {
        "name": name,
        "canon_smi": canon_smi,
        "inchikey": inchikey,
        "flat_inchikey": flat_inchikey,
        "formula": formula,
        "stereo_defined": defined,
        "stereo_undefined": undefined,
    }


def main() -> int:
    rows = list(csv.DictReader(MOLECULES_CSV.open()))
    print(f"Loaded {len(rows)} molecules")

    analyzed: list[dict[str, object]] = []
    for r in rows:
        smi = (r.get("smiles") or "").strip()
        name = (r.get("compound_name") or "").strip()
        if not smi:
            continue
        rec = analyze_row(name, smi)
        if rec is not None:
            analyzed.append(rec)
    print(f"Analyzed successfully: {len(analyzed)}/{len(rows)}\n")

    # 1. True duplicates
    by_inchi: dict[str, list[dict]] = defaultdict(list)
    for a in analyzed:
        by_inchi[a["inchikey"]].append(a)
    true_dupes = {k: v for k, v in by_inchi.items()
                  if len({x["name"] for x in v}) > 1}
    print(f"[1] True duplicates (same InChIKey, different names): {len(true_dupes)} groups")

    # 2. Name conflicts
    by_name: dict[str, list[dict]] = defaultdict(list)
    for a in analyzed:
        by_name[a["name"]].append(a)
    name_conflicts = {n: v for n, v in by_name.items()
                      if len({x["inchikey"] for x in v}) > 1}
    print(f"[2] Name conflicts (same name, different InChIKey):   {len(name_conflicts)}")

    # 3. Known stereochem pairs
    pair_specs: list[tuple[str, list[str], list[str]]] = [
        ("CDCA (7α) vs UDCA (7β)",
         ["chenodeoxycholate", "chenodeoxycholic"],
         ["ursodeoxycholate", "ursodeoxycholic"]),
        ("LCA (3α) vs iso-LCA (3β)",
         ["lithocholate", "lithocholic acid"],
         ["iso-lca", "isolithocholate", "isolithocholic"]),
        ("CA (5β) vs allo-CA (5α)",
         ["cholate", "cholic acid"],
         ["allocholate", "allocholic", "allo-ca"]),
        ("α-MCA vs β-MCA",
         ["α-muricholic", "alpha-muricholic", "αmca", "alphamca"],
         ["β-muricholic", "beta-muricholic", "βmca", "betamca"]),
        ("α-MCA vs γ-MCA",
         ["α-muricholic", "alpha-muricholic", "αmca", "alphamca"],
         ["γ-muricholic", "gamma-muricholic", "γmca", "gammamca"]),
        ("β-MCA vs γ-MCA",
         ["β-muricholic", "beta-muricholic", "βmca", "betamca"],
         ["γ-muricholic", "gamma-muricholic", "γmca", "gammamca"]),
        ("HDCA (6α-OH) vs DCA (12α-OH)",
         ["hyodeoxycholate", "hyodeoxycholic"],
         ["deoxycholate", "deoxycholic"]),
    ]

    def find_first(patterns: list[str]) -> dict | None:
        excludes = ("glyco", "tauro", "tauryl",
                    "-cys", "-ala", "-asn", "-arg")
        for a in analyzed:
            nlo = a["name"].lower()
            if any(p in nlo for p in patterns) and not any(x in nlo for x in excludes):
                return a
        return None

    pair_failures = 0
    for desc, pa, pb in pair_specs:
        a = find_first(pa)
        b = find_first(pb)
        if a is None or b is None:
            continue
        if a["inchikey"] == b["inchikey"]:
            print(f"    ✗ FAIL  {desc} share InChIKey — "
                  f"{a['name']!r} vs {b['name']!r}")
            pair_failures += 1
    print(f"[3] Known stereochem-pair failures:                   {pair_failures}")

    # 4. Ionization dupes
    by_flat: dict[str, list[dict]] = defaultdict(list)
    for a in analyzed:
        by_flat[a["flat_inchikey"]].append(a)
    ion_groups = 0
    for g in by_flat.values():
        if len(g) < 2:
            continue
        smis = [x["canon_smi"] for x in g]
        if any("[O-]" in s for s in smis) and any(
                "C(=O)O" in s and "[O-]" not in s for s in smis):
            ion_groups += 1
    print(f"[4] Ionization-dup groups (anion+acid same scaffold): {ion_groups}")

    # 5. L/D variants
    ld = [a for a in analyzed if re.search(r"^(L|D)-", a["name"])]
    by_stem: dict[str, list[dict]] = defaultdict(list)
    for a in ld:
        stem = re.sub(r"^[LD]-", "", a["name"])
        by_stem[stem].append(a)
    mixed_stems = {s: g for s, g in by_stem.items()
                   if len({a["name"][0] for a in g}) > 1}
    print(f"[5] L/D stems with both variants present:             {len(mixed_stems)}")

    # 6. Undefined stereocenters
    undef = [a for a in analyzed if a["stereo_undefined"] > 0]
    print(f"[6] Molecules with ≥1 undefined stereocenter:         {len(undef)}")

    # Per-molecule audit TSV
    out_path = ensure_results_dir() / "stereo_audit.tsv"
    with out_path.open("w", newline="") as f:
        fieldnames = ["compound_name", "inchikey", "flat_inchikey", "formula",
                      "stereo_defined", "stereo_undefined", "canon_smiles"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for a in analyzed:
            writer.writerow({
                "compound_name": a["name"],
                "inchikey": a["inchikey"],
                "flat_inchikey": a["flat_inchikey"],
                "formula": a["formula"],
                "stereo_defined": a["stereo_defined"],
                "stereo_undefined": a["stereo_undefined"],
                "canon_smiles": a["canon_smi"],
            })
    print(f"\nWrote: {out_path.relative_to(out_path.parents[2])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
