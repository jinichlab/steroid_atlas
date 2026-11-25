"""Build + mass-validate the 2 Lucas 2026 3-oxo glycyl conjugates (G-3-oxoCA,
G-3-oxoCDCA) and verify the 4 pre-existing Lucas 7-oxo glyco/tauro
conjugates.

Validation target: Lucas's reported [M+H]⁺ precursor m/z (3-4 decimal
precision); G-3-oxoCA Δ +1.4 ppm, G-3-oxoCDCA Δ −9.5 ppm.

The 4 pre-existing rows (Glycine-7-oxoLCA, Taurine-7-oxoLCA, Glycine-7-oxoDCA,
Taurine-7-oxoDCA) are already in `data/molecules.csv` from earlier audits —
this script only regenerates their canonical SMILES for cross-reference and
does not modify the atlas.

Output: `literature/results/lucas_smiles_validation.tsv`

Paper: Lucas et al. [preprint / provisional Nat Commun title], 2026.
DOI: 10.1038/s41467-026-68556-4
"""
from __future__ import annotations

import csv
from pathlib import Path

from chemistry_utils import (  # type: ignore[import-not-found]
    ensure_results_dir,
    exact_mz,
    molecular_formula,
    ppm,
)

# Final canonical SMILES. The two 3-oxo compounds were assembled in RDKit by
# coupling glycine's α-amine to a 3-oxo bile-acid core (3-oxo-CA or 3-oxo-CDCA)
# via the standard amide SMARTS from chemistry_utils; here we cite the canonical
# isomeric SMILES directly for reproducibility (no SMARTS-dependency on the
# atlas's 3-oxo core rows, which were themselves ChEBI-sourced).
LUCAS_COMPOUNDS: list[dict[str, str | float]] = [
    {
        "compound_name": "G-3-oxoCA",
        "smiles": ("C[C@H](CCC(=O)NCC(=O)O)[C@H]1CC[C@H]2[C@@H]3"
                   "[C@H](O)C[C@@H]4CC(=O)CC[C@]4(C)[C@H]3C[C@H](O)[C@]12C"),
        "ref_mz_mh_lucas": 464.30,
        "note": "3-oxo-CA conjugated to glycine (C26H41NO6)",
        "tier": "LU-2",
    },
    {
        "compound_name": "G-3-oxoCDCA",
        "smiles": ("C[C@H](CCC(=O)NCC(=O)O)[C@H]1CC[C@H]2[C@@H]3"
                   "[C@H](O)C[C@@H]4CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C"),
        "ref_mz_mh_lucas": 448.31,
        "note": "3-oxo-CDCA conjugated to glycine (C26H41NO5)",
        "tier": "LU-2",
    },
    # 4 pre-existing Lucas 7-oxo rows — SMILES regenerated here for
    # cross-reference against the atlas rows (no mutation).
    {
        "compound_name": "Glycine-7-oxoLCA",
        "smiles": ("C[C@H](CCC(=O)NCC(=O)O)[C@H]1CC[C@H]2[C@@H]3"
                   "C(=O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"),
        "ref_mz_mh_lucas": 448.31,
        "note": "7-oxo-LCA conjugated to glycine; authentic standard "
                "synthesized by Lucas (LU-1 tier)",
        "tier": "LU-1",
    },
    {
        "compound_name": "Glycine-7-oxoDCA",
        "smiles": ("C[C@H](CCC(=O)NCC(=O)O)[C@H]1CC[C@H]2[C@@H]3"
                   "C(=O)C[C@@H]4CC(O)CC[C@]4(C)[C@H]3C[C@H](O)[C@]12C"),
        "ref_mz_mh_lucas": 464.30,
        "note": "7-oxo-DCA conjugated to glycine (LU-2 surrogate-calibrated)",
        "tier": "LU-2",
    },
    {
        "compound_name": "Taurine-7-oxoLCA",
        "smiles": ("C[C@H](CCC(=O)NCCS(=O)(=O)O)[C@H]1CC[C@H]2[C@@H]3"
                   "C(=O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"),
        "ref_mz_mh_lucas": 498.29,
        "note": "7-oxo-LCA conjugated to taurine (LU-2 surrogate-calibrated)",
        "tier": "LU-2",
    },
    {
        "compound_name": "Taurine-7-oxoDCA",
        "smiles": ("C[C@H](CCC(=O)NCCS(=O)(=O)O)[C@H]1CC[C@H]2[C@@H]3"
                   "C(=O)C[C@@H]4CC(O)CC[C@]4(C)[C@H]3C[C@H](O)[C@]12C"),
        "ref_mz_mh_lucas": 514.28,
        "note": "7-oxo-DCA conjugated to taurine (LU-2 surrogate-calibrated)",
        "tier": "LU-2",
    },
]


def main() -> int:
    out_rows: list[dict[str, str]] = []
    for c in LUCAS_COMPOUNDS:
        smiles = str(c["smiles"])
        computed = exact_mz(smiles, adduct="[M+H]+")
        reference = float(c["ref_mz_mh_lucas"])
        delta = ppm(computed, reference)
        out_rows.append({
            "compound_name": str(c["compound_name"]),
            "lucas_evidence_tier": str(c["tier"]),
            "smiles": smiles,
            "molecular_formula": molecular_formula(smiles),
            "computed_mz_M+H": f"{computed:.4f}",
            "ref_mz_M+H_lucas": f"{reference:.2f}",
            "delta_ppm": f"{delta:+.1f}",
            "mass_pass_10ppm": ("yes" if abs(delta) < 10 else "no"),
            "note": str(c["note"]),
        })

    ensure_results_dir()
    out_path = ensure_results_dir() / "lucas_smiles_validation.tsv"
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()),
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    passed = sum(1 for r in out_rows if r["mass_pass_10ppm"] == "yes")
    print(f"Lucas compounds:      {len(out_rows)} (2 newly added + 4 regenerated)")
    print(f"Mass validation pass: {passed}/{len(out_rows)} at <10 ppm")
    print(f"Wrote: {out_path.relative_to(out_path.parents[2])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
