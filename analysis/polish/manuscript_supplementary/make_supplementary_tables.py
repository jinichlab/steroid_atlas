"""Generate manuscript Supplementary Tables S1a, S1b, and S1c (CSV only).

S1a — Literature-recruited proteins (13 rows)
     Columns: Accession, Paper, DOI, Species, Strain, Identifier type,
              Length (aa), Evidence level, Recruitment rationale, In atlas UMAP

S1b — Literature-recruited molecules (241 rows)
     Columns: Compound name, ChEBI ID, SMILES, Paper DOI(s), Evidence tier,
              Evidence note, Atlas inclusion rationale

S1c — Evidence-tier legend
     Defines every tier code used in S1b.

Sources:
  data/literature_recruited_proteins.csv
  data/molecules.csv  (is_literature_recruited == 1)
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]  # steroid_atlas/
OUT = Path(__file__).resolve().parent
PROTEINS_LR = ROOT / "data" / "literature_recruited_proteins.csv"
MOLECULES = ROOT / "data" / "molecules.csv"


# ---------------------------------------------------------------------------
# S1a — proteins
# ---------------------------------------------------------------------------

def build_s1a() -> pd.DataFrame:
    df = pd.read_csv(PROTEINS_LR)
    cols = [
        "Entry", "Paper_ShortRef", "Paper_DOI", "Species", "Strain",
        "Identifier_Type", "Sequence_Length", "Evidence_Level",
        "Recruitment_Rationale", "In_Atlas_UMAP",
    ]
    df = df[cols].copy()
    df = df.rename(columns={
        "Entry": "Accession",
        "Paper_ShortRef": "Paper",
        "Paper_DOI": "DOI",
        "Identifier_Type": "Identifier type",
        "Sequence_Length": "Length (aa)",
        "Evidence_Level": "Evidence level",
        "Recruitment_Rationale": "Recruitment rationale",
        "In_Atlas_UMAP": "In atlas UMAP",
    })
    return df.sort_values(["Paper", "Accession"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# S1b — molecules
# ---------------------------------------------------------------------------

def build_s1b() -> pd.DataFrame:
    m = pd.read_csv(MOLECULES)
    lit = m[m["is_literature_recruited"] == 1].copy()

    def fmt_chebi(x):
        if pd.isna(x) or str(x).strip() == "":
            return ""
        try:
            return f"CHEBI:{int(float(x))}"
        except (ValueError, TypeError):
            s = str(x).strip()
            if s.upper().startswith("CHEBI:"):
                return s.upper()
            return s

    lit["chebi_id"] = lit["chebi_id"].apply(fmt_chebi)

    def clean_urls(v):
        if pd.isna(v):
            return ""
        parts = [p.strip() for p in str(v).split("|") if p.strip()]
        return " | ".join(parts)

    lit["paper_url"] = lit["paper_url"].apply(clean_urls)

    out_cols = [
        "compound_name", "chebi_id", "smiles", "paper_url",
        "evidence_tier", "evidence_note",
    ]
    extras = [c for c in ("atlas_inclusion_rationale",) if c in lit.columns]
    out_cols += extras
    lit = lit[out_cols].copy()
    lit = lit.rename(columns={
        "compound_name": "Compound name",
        "chebi_id": "ChEBI ID",
        "smiles": "SMILES",
        "paper_url": "Paper DOI(s)",
        "evidence_tier": "Evidence tier",
        "evidence_note": "Evidence note",
        "atlas_inclusion_rationale": "Atlas inclusion rationale",
    })
    return lit.sort_values(["Paper DOI(s)", "Compound name"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# S1c — evidence-tier legend
# ---------------------------------------------------------------------------

TIER_LEGEND = [
    # Gentry 2023 — 5-level scale across 4 evidence streams
    ("1",     "Gentry 2023",     "Chemistry reference standard only (no biological detection)"),
    ("2",     "Gentry 2023",     "MASST match in public metabolomics data (isomer uncertain, MSI Level 2/3)"),
    ("3",     "Gentry 2023",     "Isomer-resolved detection in iHMP2 IBD cohort (RT + MS/MS)"),
    ("4-5",   "Gentry 2023",     "Cross-cohort detection (iHMP2 + Boland/PRISM IBD cohorts)"),
    ("5",     "Gentry 2023",     "Produced by HMP bacterial isolate(s) in culture (highest confidence)"),
    # Mohanty 2024 — polyamine bile amidates
    ("M-1",   "Mohanty 2024",    "Synthetic standard + RT/MS/MS match in BOTH mice and lion feces (strongest)"),
    ("M-2",   "Mohanty 2024",    "Synthetic standard + RT/MS/MS match in mice only"),
    ("M-3",   "Mohanty 2024",    "Synthetic standard + RT/MS/MS match in lion only"),
    ("M-4",   "Mohanty 2024",    "Same-mass peak detected but RT did not match standard (related isomer in vivo)"),
    # Lamichhane 2025 — DIABIMMUNE Cys-MCBAs
    ("L-1",   "Lamichhane 2025", "Targeted LC-MS/MS detection in 303-sample DIABIMMUNE cohort + statistically significant association with islet autoimmunity"),
    ("L-2",   "Lamichhane 2025", "Part of 110-MCBA targeted panel in DIABIMMUNE cohort (no significant clinical association)"),
    # Lucas 2026 — glyco-3-oxo bile acids
    ("LU-1",  "Lucas 2026",      "Authentic synthesized standard + detection in bacterial culture (highest Lucas confidence)"),
    ("LU-2",  "Lucas 2026",      "LC-MS/MS detection in bacterial culture; quantified via surrogate-standard calibration curve"),
    # Agongo 2026 — bile-acid alkamines (new chemistry class)
    ("AG-1",  "Agongo 2026",     "Level-1 identification (MS/MS + RT match to synthetic standard) + ion-mobility drift-time match"),
    ("AG-2",  "Agongo 2026",     "Level-1 identification (MS/MS + RT match to synthetic standard; no ion mobility)"),
    # Legacy / free-text
    ("biochem-substrate", "Various (biochem studies)", "Tool compound used as an enzyme substrate in a biochemical assay"),
    ("(empty)", "Rimal / Guzior / Arp / McCurry / etc.", "Pre-audit curation pass — no paper-specific tier assigned"),
]


def build_s1c() -> pd.DataFrame:
    return pd.DataFrame(TIER_LEGEND,
                        columns=["Evidence tier code", "Source paper", "Meaning"])


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    print("=== S1a — proteins ===")
    s1a = build_s1a()
    print(f"  rows: {len(s1a)}")
    s1a.to_csv(OUT / "S1a_literature_recruited_proteins.csv", index=False)
    print("  wrote: S1a_literature_recruited_proteins.csv")
    print()

    print("=== S1b — molecules ===")
    s1b = build_s1b()
    print(f"  rows: {len(s1b)}")
    s1b.to_csv(OUT / "S1b_literature_recruited_molecules.csv", index=False)
    print("  wrote: S1b_literature_recruited_molecules.csv")
    print()

    print("=== S1c — evidence-tier legend ===")
    s1c = build_s1c()
    print(f"  rows: {len(s1c)}")
    s1c.to_csv(OUT / "S1c_evidence_tier_legend.csv", index=False)
    print("  wrote: S1c_evidence_tier_legend.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
