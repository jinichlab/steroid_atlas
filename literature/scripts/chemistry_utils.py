"""Shared chemistry helpers for the literature-recruited SMILES build scripts.

Used by `build_gentry_168_mcbas.py`, `build_mohanty_17_polyamines.py`,
`build_lamichhane_5_cys_mcbas.py`, `build_lucas_g3oxo_mcbas.py`,
`build_agongo_18_alkamines.py`, and `stereo_audit_molecules.py`.

Each build script:
  1. Loads reference bile-acid cores from `literature/data/bile_acid_cores.json`
  2. Assembles target SMILES via reaction SMARTS
     (amide coupling or reductive amination)
  3. Validates the computed monoisotopic mass against the source paper's
     reported m/z (or theoretical value when the paper does not publish one)
  4. Writes a validation TSV to `literature/results/<paper>_smiles_validation.tsv`

The build scripts do NOT modify `data/molecules.csv`. Historical ingestion is
already recorded in `.bak_*` backups; these scripts exist so the per-compound
chemistry provenance can be regenerated from scratch.
"""
from __future__ import annotations

import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.Descriptors import ExactMolWt

REPO_ROOT = Path(__file__).resolve().parents[2]
BA_CORES_JSON = REPO_ROOT / "literature" / "data" / "bile_acid_cores.json"
RESULTS_DIR = REPO_ROOT / "literature" / "results"
SUPPLEMENTARY_DIR = REPO_ROOT / "literature" / "supplementary"
MOLECULES_CSV = REPO_ROOT / "data" / "molecules.csv"

PROTON_MASS = 1.00728


def load_bile_acid_cores() -> dict[str, str]:
    """Return the author-curated BA core SMILES dict (CA, CDCA, DCA, UDCA, HDCA,
    alpha/beta/gamma/omega-MCA). Keys include the Greek letters for the MCA
    isomers (αMCA, βMCA, γMCA, ωMCA)."""
    return json.loads(BA_CORES_JSON.read_text())


# ---------------------------------------------------------------------------
# Reaction SMARTS
# ---------------------------------------------------------------------------

AMIDE_COUPLE_RXN = AllChem.ReactionFromSmarts(
    "[CX3:1](=[O:2])[OX2H].[NX3H2:3][CX4:4] "
    ">> [CX3:1](=[O:2])[NX3H1:3][CX4:4]"
)
"""Amide bond formation between a carboxylic acid (-COOH) and a primary amine
on an sp3 carbon. Used for Gentry (21 proteinogenic amino acids × 8 BA cores),
Mohanty (polyamines), Lamichhane (cysteine), and Lucas (glycine 3-oxo
amidates)."""

PRO_AMIDE_COUPLE_RXN = AllChem.ReactionFromSmarts(
    "[CX3:1](=[O:2])[OX2H].[NX3H1:3]([CX4:4])[CX4:5] "
    ">> [CX3:1](=[O:2])[NX3:3]([CX4:4])[CX4:5]"
)
"""Amide bond formation when the amine is a secondary amine (proline)."""

REDUCTIVE_AMINATION_RXN = AllChem.ReactionFromSmarts(
    "[CX3:1](=[OX1])(-[CH2;R:2])(-[CH2;R:3])."
    "[NX3H2:4][*:5]"
    ">> [C:1](-[CH2;R:2])(-[CH2;R:3])-[NX3H1:4][*:5]"
)
"""Reductive amination of a ring ketone with a primary amine → secondary amine
+ H2O (net). Used for Agongo 2026 alkamines (C3 ketone of bile-acid 3-oxo
cores or cholestan-3-one). The ring-neighbor constraint prevents the SMARTS
from matching the BA C24 carboxyl."""


def amide_couple(acid_smiles: str, amine_smiles: str,
                 secondary_amine: bool = False) -> str:
    """Return the canonical isomeric SMILES for the amide formed between
    `acid_smiles` (any molecule with a -COOH) and `amine_smiles` (any molecule
    with a primary amine; pass `secondary_amine=True` for proline). Returns
    an empty string if the coupling SMARTS does not match."""
    rxn = PRO_AMIDE_COUPLE_RXN if secondary_amine else AMIDE_COUPLE_RXN
    acid = Chem.MolFromSmiles(acid_smiles)
    amine = Chem.MolFromSmiles(amine_smiles)
    products = rxn.RunReactants((acid, amine))
    if not products:
        return ""
    out = products[0][0]
    try:
        Chem.SanitizeMol(out)
        return Chem.MolToSmiles(out, isomericSmiles=True)
    except Exception:  # noqa: BLE001 — SMARTS can produce transient invalid mols
        return ""


def reductive_aminate(ketone_smiles: str, amine_smiles: str) -> str:
    """Return the canonical isomeric SMILES for the secondary amine formed by
    reductive amination of a ring ketone with a primary amine. The C3
    stereocenter is left unassigned (both 3α and 3β are biologically observed
    in Agongo 2026). Returns an empty string on no-match."""
    kmol = Chem.MolFromSmiles(ketone_smiles)
    amol = Chem.MolFromSmiles(amine_smiles)
    products = REDUCTIVE_AMINATION_RXN.RunReactants((kmol, amol))
    if not products:
        return ""
    for p_set in products:
        out = p_set[0]
        try:
            Chem.SanitizeMol(out)
            return Chem.MolToSmiles(out, isomericSmiles=True)
        except Exception:  # noqa: BLE001
            continue
    return ""


# ---------------------------------------------------------------------------
# Mass computation / validation
# ---------------------------------------------------------------------------

def exact_mz(smiles: str, adduct: str = "[M+H]+") -> float:
    """Return the monoisotopic m/z of `smiles` for the requested adduct.
    Supports `[M+H]+`, `[M-H]-`, and `[M+Na]+`."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles!r}")
    mass = ExactMolWt(mol)
    if adduct == "[M+H]+":
        return mass + PROTON_MASS
    if adduct == "[M-H]-":
        return mass - PROTON_MASS
    if adduct == "[M+Na]+":
        return mass + 22.98922
    raise ValueError(f"unsupported adduct: {adduct!r}")


def ppm(computed: float, reference: float) -> float:
    """Return the ppm deviation of a computed m/z from a reference m/z."""
    return (computed - reference) / reference * 1e6


def molecular_formula(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles!r}")
    return rdMolDescriptors.CalcMolFormula(mol)


def ensure_results_dir() -> Path:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    return RESULTS_DIR
