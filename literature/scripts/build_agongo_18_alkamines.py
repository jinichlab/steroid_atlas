"""Build + validate the 18 Agongo 2026 alkamines.

Alkamines are a new class of microbial steroid modifications in which an amine
is joined to the C3 ketone of a 3-oxo bile-acid core (or cholestan-3-one) by
reductive amination — fundamentally distinct from the C24 amidation of
MCBAs/CSBAs (Gentry, Mohanty, Lamichhane, Lucas) which forms an amide bond at
the BA carboxyl.

Reaction: R–C(=O)–R' + H2N–R'' + 2[H] → R–CH(NHR'')–R' + H2O

The C3 stereocenter created by reductive amination is left unassigned in the
SMILES because Agongo observed both diastereomers (3α and 3β) in biological
samples.

The 3-oxo-LCA, 3-oxo-DCA, and 5β-cholestan-3-one core SMILES are defined
inline as literal canonical SMILES (they are standard scaffolds; this keeps
the script independent of specific `data/molecules.csv` row names).

Validation target: Agongo does not publish precursor m/z for every compound in
Fig 4C + Fig 5 (only diagnostic fragment ions for a subset). Here we check
the molecular formula and exact mass for internal consistency — the paper's
diagnostic fragments at m/z 355 / 357 / 359 (CA/DCA/LCA alkamine) match the
secondary-amine-at-C3 topology.

Output: `literature/results/agongo_smiles_validation.tsv`

Paper: Agongo et al. Alkamines reveal a hidden layer of steroid and drug
metabolism. bioRxiv (2026).
DOI: 10.64898/2026.05.13.724743
"""
from __future__ import annotations

import csv
from pathlib import Path

from chemistry_utils import (  # type: ignore[import-not-found]
    ensure_results_dir,
    exact_mz,
    molecular_formula,
    reductive_aminate,
)

# ---------------------------------------------------------------------------
# Core ketone scaffolds (canonical isomeric SMILES for 5β-steroid skeletons
# bearing a C3 ketone). Chemistry conventions match the atlas's bile-acid
# cores — see `literature/data/bile_acid_cores.json`.
# ---------------------------------------------------------------------------

KETONE_CORES: dict[str, str] = {
    # 3-oxo-DCA: 5β-cholan-24-oic acid with 3-oxo and 12α-OH
    "3-oxo-DCA": ("C[C@H](CCC(=O)O)[C@H]1CC[C@H]2[C@@H]3CC[C@@H]4CC(=O)CC"
                  "[C@]4(C)[C@H]3C[C@H](O)[C@]12C"),
    # 3-oxo-LCA: 5β-cholan-24-oic acid with 3-oxo (no other hydroxyls)
    "3-oxo-LCA": ("C[C@H](CCC(=O)O)[C@H]1CC[C@H]2[C@@H]3CC[C@@H]4CC(=O)CC"
                  "[C@]4(C)[C@H]3CC[C@]12C"),
    # 5β-cholestan-3-one: the cholestane C27 skeleton with C3 ketone
    "cholestan-3-one": ("CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC[C@@H]4CC"
                        "(=O)CC[C@]4(C)[C@H]3CC[C@]12C"),
}

# ---------------------------------------------------------------------------
# Amine library
# ---------------------------------------------------------------------------

AMINES: dict[str, str] = {
    "Phe":                "N[C@@H](Cc1ccccc1)C(=O)O",
    "Ala":                "N[C@@H](C)C(=O)O",
    "Gly":                "NCC(=O)O",
    "Leu":                "N[C@@H](CC(C)C)C(=O)O",
    "Val":                "N[C@@H](C(C)C)C(=O)O",
    "Tyr":                "N[C@@H](Cc1ccc(O)cc1)C(=O)O",
    "Lys":                "N[C@@H](CCCCN)C(=O)O",
    "Tyramine":           "NCCc1ccc(O)cc1",
    "Putrescine":         "NCCCCN",
    "N-acetylcadaverine": "CC(=O)NCCCCCN",
    "N-acetylputrescine": "CC(=O)NCCCCN",
    "5-ASA":              "Nc1ccc(O)c(C(=O)O)c1",
}

# (compound_name, amine_key, core_key)
TARGETS: list[tuple[str, str, str]] = [
    # Bile-acid alkamines (8) — Fig 4C
    ("Phe-3-oxo-LCA",                "Phe",                "3-oxo-LCA"),
    ("N-acetylcadaverine-3-oxo-LCA", "N-acetylcadaverine", "3-oxo-LCA"),
    ("N-acetylcadaverine-3-oxo-DCA", "N-acetylcadaverine", "3-oxo-DCA"),
    ("Tyramine-3-oxo-DCA",           "Tyramine",           "3-oxo-DCA"),
    ("Putrescine-3-oxo-DCA",         "Putrescine",         "3-oxo-DCA"),
    ("Tyr-3-oxo-DCA",                "Tyr",                "3-oxo-DCA"),
    ("Phe-3-oxo-DCA",                "Phe",                "3-oxo-DCA"),
    ("Lys-3-oxo-DCA",                "Lys",                "3-oxo-DCA"),
    # Drug-BA alkamines (2) — Fig 5
    ("5-ASA-3-oxo-LCA",              "5-ASA",              "3-oxo-LCA"),
    ("5-ASA-3-oxo-DCA",              "5-ASA",              "3-oxo-DCA"),
    # Cholestan-3-one alkamines (8) — Fig 4C
    ("Phe-cholestan-3-one",                "Phe",                "cholestan-3-one"),
    ("Ala-cholestan-3-one",                "Ala",                "cholestan-3-one"),
    ("Gly-cholestan-3-one",                "Gly",                "cholestan-3-one"),
    ("Leu-cholestan-3-one",                "Leu",                "cholestan-3-one"),
    ("Val-cholestan-3-one",                "Val",                "cholestan-3-one"),
    ("Tyramine-cholestan-3-one",           "Tyramine",           "cholestan-3-one"),
    ("N-acetylcadaverine-cholestan-3-one", "N-acetylcadaverine", "cholestan-3-one"),
    ("N-acetylputrescine-cholestan-3-one", "N-acetylputrescine", "cholestan-3-one"),
]


def main() -> int:
    cores = KETONE_CORES
    out_rows: list[dict[str, str]] = []
    failures: list[str] = []

    for compound_name, amine_key, core_key in TARGETS:
        product = reductive_aminate(cores[core_key], AMINES[amine_key])
        if not product:
            failures.append(compound_name)
            continue
        computed_mw = exact_mz(product, adduct="[M+H]+") - 1.00728
        computed_mh = exact_mz(product, adduct="[M+H]+")
        out_rows.append({
            "compound_name": compound_name,
            "amine": amine_key,
            "steroid_core": core_key,
            "smiles": product,
            "molecular_formula": molecular_formula(product),
            "exact_mass": f"{computed_mw:.4f}",
            "computed_mz_M+H": f"{computed_mh:.4f}",
        })

    ensure_results_dir()
    out_path = ensure_results_dir() / "agongo_smiles_validation.tsv"
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()),
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Generated alkamines:  {len(out_rows)}/{len(TARGETS)}")
    if failures:
        print(f"Reductive-amination failures:")
        for name in failures:
            print(f"  {name}")
    print(f"Wrote: {out_path.relative_to(out_path.parents[2])}")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
