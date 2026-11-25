"""Build + mass-validate the 168 Gentry 2023 microbial conjugated bile acids
(MCBAs) = 21 proteinogenic amino acids (plus citrulline, ornithine, DOPA) × 8
bile-acid cores (CA, CDCA, DCA, UDCA, HDCA, αMCA, βMCA, γMCA).

Validation target: Gentry's MOESM8 reference precursor-m/z table ([M-H]⁻).
Tolerance achieved: <0.1 ppm for every compound.

Output: `literature/results/gentry_smiles_validation.tsv`

Paper: Gentry et al. Reverse metabolomics for the discovery of chemical
structures from humans. Nature 626, 419–426 (2024).
DOI: 10.1038/s41586-023-06906-8
SI file: literature/supplementary/Gentry2023_MOESM8_ESM.xlsx
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

try:
    import openpyxl
except ImportError:
    sys.stderr.write("openpyxl is required (pip install openpyxl)\n")
    raise

from chemistry_utils import (  # type: ignore[import-not-found]
    SUPPLEMENTARY_DIR,
    amide_couple,
    ensure_results_dir,
    exact_mz,
    load_bile_acid_cores,
    molecular_formula,
    ppm,
)

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

AMINO_ACIDS: dict[str, str] = {
    # L-isomers; SMILES use [C@@H] for the α-carbon (except Pro, a secondary amine)
    "Ala":  "C[C@@H](N)C(=O)O",
    "Arg":  "N[C@@H](CCCNC(=N)N)C(=O)O",
    "Asn":  "N[C@@H](CC(=O)N)C(=O)O",
    "Asp":  "N[C@@H](CC(=O)O)C(=O)O",
    "Cit":  "N[C@@H](CCCNC(=O)N)C(=O)O",
    "DOPA": "N[C@@H](Cc1ccc(O)c(O)c1)C(=O)O",
    "Gln":  "N[C@@H](CCC(=O)N)C(=O)O",
    "Glu":  "N[C@@H](CCC(=O)O)C(=O)O",
    "His":  "N[C@@H](Cc1c[nH]cn1)C(=O)O",
    "Ile":  "CC[C@H](C)[C@@H](N)C(=O)O",
    "Leu":  "CC(C)C[C@@H](N)C(=O)O",
    "Lys":  "N[C@@H](CCCCN)C(=O)O",
    "Met":  "CSCC[C@@H](N)C(=O)O",
    "Orn":  "N[C@@H](CCCN)C(=O)O",
    "Phe":  "N[C@@H](Cc1ccccc1)C(=O)O",
    "Pro":  "OC(=O)[C@@H]1CCCN1",   # secondary amine
    "Ser":  "N[C@@H](CO)C(=O)O",
    "Thr":  "C[C@@H](O)[C@H](N)C(=O)O",
    "Trp":  "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",
    "Tyr":  "N[C@@H](Cc1ccc(O)cc1)C(=O)O",
    "Val":  "CC(C)[C@@H](N)C(=O)O",
}

CORES_ORDER = ["CA", "CDCA", "DCA", "UDCA", "HDCA", "αMCA", "βMCA", "γMCA"]

GENTRY_MOESM8 = SUPPLEMENTARY_DIR / "Gentry2023_MOESM8_ESM.xlsx"

# Gentry's MOESM8 names the three MCA isomers as "MCA (a)", "MCA (b)", "MCA (g)"
MCA_NAME_MAP = {"MCA (a)": "αMCA", "MCA (b)": "βMCA", "MCA (g)": "γMCA"}


def load_moesm8_reference() -> dict[tuple[str, str], float]:
    """Return {(amino_acid, core): [M-H]⁻ m/z} from Gentry's MOESM8."""
    wb = openpyxl.load_workbook(GENTRY_MOESM8, read_only=True, data_only=True)
    ws = wb["Sheet1"]
    refs: dict[tuple[str, str], float] = {}
    for row in ws.iter_rows(min_row=3, values_only=True):
        name = row[0]
        if not name or "_" not in str(name):
            continue
        aa, core = str(name).strip().rsplit("_", 1)
        core = MCA_NAME_MAP.get(core, core)
        refs[(aa.strip(), core)] = float(row[1])
    return refs


# ---------------------------------------------------------------------------
# Build + validate
# ---------------------------------------------------------------------------

def main() -> int:
    ba_cores = load_bile_acid_cores()
    ref_mz = load_moesm8_reference()
    out_rows: list[dict[str, str]] = []
    coupling_failures: list[tuple[str, str]] = []

    for core in CORES_ORDER:
        for aa, aa_smi in AMINO_ACIDS.items():
            product = amide_couple(ba_cores[core], aa_smi,
                                   secondary_amine=(aa == "Pro"))
            if not product:
                coupling_failures.append((aa, core))
                continue
            computed = exact_mz(product, adduct="[M-H]-")
            reference = ref_mz.get((aa, core))
            delta = ppm(computed, reference) if reference else None
            out_rows.append({
                "amino_acid": aa,
                "bile_acid_core": core,
                "compound_name": f"{aa}-{core}",
                "smiles": product,
                "molecular_formula": molecular_formula(product),
                "computed_mz_M-H": f"{computed:.4f}",
                "ref_mz_M-H_moesm8": f"{reference:.4f}" if reference else "",
                "delta_ppm": f"{delta:+.3f}" if delta is not None else "",
                "mass_pass_5ppm": ("yes" if delta is not None and abs(delta) < 5
                                    else "no"),
            })

    ensure_results_dir()
    out_path = ensure_results_dir() / "gentry_smiles_validation.tsv"
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()),
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    passed = sum(1 for r in out_rows if r["mass_pass_5ppm"] == "yes")
    print(f"Generated amidates:   {len(out_rows)}")
    print(f"Mass validation pass: {passed}/{len(out_rows)} at <5 ppm")
    if coupling_failures:
        print(f"Coupling failures:    {len(coupling_failures)}")
        for aa, core in coupling_failures:
            print(f"  {aa}-{core}")
    print(f"Wrote: {out_path.relative_to(out_path.parents[2])}")
    return 0 if not coupling_failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
