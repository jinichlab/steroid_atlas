"""Build + mass-validate the 5 Lamichhane 2025 cysteinyl MCBAs.

Lamichhane 2025 quantified 110 MCBAs = 5 bile-acid cores (CA, CDCA, DCA, UDCA,
HDCA) × 22 amino acids. Cross-checking against the atlas (Gentry 2023's
8 × 21 = 168 MCBAs, which exclude Cys) identified only the 5 cysteinyl
conjugates as novel atlas additions:

  CA-Cys, CDCA-Cys, DCA-Cys, UDCA-Cys, HDCA-Cys

L-cysteine couples via its α-amine (not the thiol) to the BA C24 carboxyl.

Lamichhane does not publish a per-compound precursor-m/z table — they use
authentic synthetic standards with RT + MS/MS matching. The validation target
here is the theoretical [M+H]⁺ from the molecular formula; all 5 computed
values match within <10 ppm.

Output: `literature/results/lamichhane_smiles_validation.tsv`

Paper: Lamichhane et al. Microbially conjugated bile acids in early-life gut
microbiota and risk of islet autoimmunity. Nat Commun 16, 8932 (2025).
DOI: 10.1038/s41467-025-66619-6
"""
from __future__ import annotations

import csv
from pathlib import Path

from chemistry_utils import (  # type: ignore[import-not-found]
    amide_couple,
    ensure_results_dir,
    exact_mz,
    load_bile_acid_cores,
    molecular_formula,
    ppm,
)

L_CYSTEINE_SMILES = "N[C@@H](CS)C(=O)O"

CORES = ("CA", "CDCA", "DCA", "UDCA", "HDCA")

# Theoretical [M+H]+ values from molecular formulas:
#   CA + Cys  - H2O -> C27H45NO6S, [M+H]+ 512.3041
#   CDCA/DCA/UDCA/HDCA + Cys - H2O -> C27H45NO5S, [M+H]+ 496.3091
EXPECTED_MH_THEORETICAL: dict[str, float] = {
    "CA":   512.30,
    "CDCA": 496.31,
    "DCA":  496.31,
    "UDCA": 496.31,
    "HDCA": 496.31,
}


def main() -> int:
    ba_cores = load_bile_acid_cores()
    out_rows: list[dict[str, str]] = []

    for core in CORES:
        product = amide_couple(ba_cores[core], L_CYSTEINE_SMILES)
        assert product, f"coupling failed for {core}-Cys"
        computed = exact_mz(product, adduct="[M+H]+")
        reference = EXPECTED_MH_THEORETICAL[core]
        delta = ppm(computed, reference)
        out_rows.append({
            "bile_acid_core": core,
            "amino_acid": "L-cysteine",
            "compound_name": f"{core}-Cys",
            "smiles": product,
            "molecular_formula": molecular_formula(product),
            "computed_mz_M+H": f"{computed:.4f}",
            "theoretical_mz_M+H": f"{reference:.2f}",
            "delta_ppm": f"{delta:+.1f}",
            "mass_pass_10ppm": ("yes" if abs(delta) < 10 else "no"),
        })

    ensure_results_dir()
    out_path = ensure_results_dir() / "lamichhane_smiles_validation.tsv"
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()),
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    passed = sum(1 for r in out_rows if r["mass_pass_10ppm"] == "yes")
    print(f"Generated amidates:   {len(out_rows)}")
    print(f"Mass validation pass: {passed}/{len(out_rows)} at <10 ppm vs theoretical")
    print(f"Wrote: {out_path.relative_to(out_path.parents[2])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
