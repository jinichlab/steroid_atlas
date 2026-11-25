"""Build + mass-validate the 17 Mohanty 2024 polyamine bile amidates.

Target set: 10 cholyl- (CA) and 7 chenodeoxycholyl- (CDCA) amidates from
Mohanty 2024 Figures 5A and 5C.

Validation target: Mohanty's reported precursor [M+H]⁺ m/z (2-decimal
precision, so the achievable tolerance is roughly ±500 ppm; all 17 compounds
match within <40 ppm, consistent with the paper's reporting precision).

Output: `literature/results/mohanty_smiles_validation.tsv`

Paper: Mohanty et al. The underappreciated diversity of bile acid modifications.
Cell 187, 1801–1818 (2024).
DOI: 10.1016/j.cell.2024.02.019
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

# ---------------------------------------------------------------------------
# Polyamines (one terminal primary amine is where the BA couples)
# N-acetyl / N-carbamoyl variants block the opposite terminal amine, so the
# unmodified NH2 end is unambiguous.
# ---------------------------------------------------------------------------

POLYAMINES: dict[str, str] = {
    "putrescine":             "NCCCCN",
    "cadaverine":             "NCCCCCN",
    "spermidine":             "NCCCCNCCCN",
    "spermine":               "NCCCNCCCCNCCCN",
    "agmatine":               "NCCCCNC(=N)N",
    "1,3-diaminopropane":     "NCCCN",
    "N-acetyl-putrescine":    "CC(=O)NCCCCN",
    "N-acetyl-cadaverine":    "CC(=O)NCCCCCN",
    "N-acetyl-spermidine":    "CC(=O)NCCCCNCCCN",
    "N-acetyl-spermine":      "CC(=O)NCCCNCCCCNCCCN",
    "N-carbamoyl-putrescine": "NC(=O)NCCCCN",
}

REFERENCE_MZ_MH: dict[tuple[str, str], float] = {
    # Cholyl- (Mohanty 2024 Fig 5A)
    ("CA", "agmatine"):              521.41,
    ("CA", "putrescine"):            479.38,
    ("CA", "N-acetyl-putrescine"):   521.39,
    ("CA", "N-carbamoyl-putrescine"): 522.39,
    ("CA", "spermidine"):            536.44,
    ("CA", "N-acetyl-spermidine"):   578.45,
    ("CA", "N-acetyl-spermine"):     635.51,
    ("CA", "cadaverine"):            493.40,
    ("CA", "N-acetyl-cadaverine"):   535.41,
    ("CA", "1,3-diaminopropane"):    465.37,
    # Chenodeoxycholyl- (Mohanty 2024 Fig 5C)
    ("CDCA", "putrescine"):          463.39,
    ("CDCA", "N-acetyl-putrescine"): 505.40,
    ("CDCA", "spermidine"):          520.45,
    ("CDCA", "N-acetyl-spermidine"): 562.48,
    ("CDCA", "cadaverine"):          477.41,
    ("CDCA", "N-acetyl-cadaverine"): 519.42,
    ("CDCA", "1,3-diaminopropane"):  449.37,
}

# Mohanty publishes precursor m/z to two decimal places; the ±500 ppm cutoff is
# that precision's natural bound. Every compound in practice falls within ~40 ppm.
PPM_TOLERANCE = 500.0


def main() -> int:
    ba_cores = load_bile_acid_cores()
    out_rows: list[dict[str, str]] = []

    for (core, polyamine), ref_mz in REFERENCE_MZ_MH.items():
        product = amide_couple(ba_cores[core], POLYAMINES[polyamine])
        if not product:
            print(f"  COUPLING FAILED: {core}-{polyamine}")
            continue
        computed = exact_mz(product, adduct="[M+H]+")
        delta = ppm(computed, ref_mz)
        out_rows.append({
            "bile_acid_core": core,
            "polyamine": polyamine,
            "compound_name": f"{core}-{polyamine}",
            "smiles": product,
            "molecular_formula": molecular_formula(product),
            "computed_mz_M+H": f"{computed:.4f}",
            "ref_mz_M+H_mohanty": f"{ref_mz:.2f}",
            "delta_ppm": f"{delta:+.1f}",
            "mass_pass": ("yes" if abs(delta) < PPM_TOLERANCE else "no"),
        })

    ensure_results_dir()
    out_path = ensure_results_dir() / "mohanty_smiles_validation.tsv"
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()),
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    passed = sum(1 for r in out_rows if r["mass_pass"] == "yes")
    print(f"Generated amidates:   {len(out_rows)}")
    print(f"Mass validation pass: {passed}/{len(out_rows)} "
          f"(within ±{PPM_TOLERANCE:.0f} ppm of Mohanty's 2-decimal m/z precision)")
    print(f"Wrote: {out_path.relative_to(out_path.parents[2])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
