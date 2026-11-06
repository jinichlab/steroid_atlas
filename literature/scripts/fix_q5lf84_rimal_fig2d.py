"""Apply Rimal 2024 Fig 2d corrections to the atlas for Q5LF84 (B. fragilis BSH).

Changes:
  1. Add `synonyms` column to molecules.csv (if absent).
  2. Merge each paper-shorthand MCBA row into the canonical ChEBI row:
       Ala-CA  -> L-alanocholate      (CHEBI:229710)
       Ser-CA  -> L-serocholate       (CHEBI:229711)
       Phe-CA  -> L-phenylalanocholate (CHEBI:229714)
     For each merged canonical row:
       * is_literature_recruited = 1
       * synonyms = paper-shorthand name (e.g. "Ala-CA")
       * paper_url: union of original canonical and shorthand paper_urls, plus
         Rimal 2024 (https://www.nature.com/articles/s41586-023-06990-w)
         when Fig 2d supports it.
       * interacting_protein_accessions: append Q5LF84 (preserves existing P0DXD2)
     Then delete the paper-shorthand row.
  3. Tyr-CA has no canonical ChEBI row in the atlas -> keep its row,
     append Q5LF84 to interacting_protein_accessions.
  4. Also append Q5LF84 to taurocholate (CHEBI:36257) reverse-map
     (Fig 2d substrate).
  5. Update proteins.csv for Q5LF84:
       * interacting_compounds: trim to 6 Fig-2d entities
         (cholate, taurocholate, L-alanocholate, L-phenylalanocholate,
          L-serocholate, Tyr-CA)
       * interacting_chebi_ids: 229710;229711;229714;36257
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/home/adsiordia/steroid_atlas")
MOL = ROOT / "data" / "molecules.csv"
PRO = ROOT / "data" / "proteins.csv"

RIMAL_URL = "https://www.nature.com/articles/s41586-023-06990-w"


def union_urls(*fields: str) -> str:
    """Union of ' | '-separated URL strings, order-preserving."""
    seen: list[str] = []
    for f in fields:
        if not isinstance(f, str):
            continue
        for u in f.split(" | "):
            u = u.strip()
            if u and u not in seen:
                seen.append(u)
    return " | ".join(seen)


def parse_prot_list(v) -> list[str]:
    """Interacting_protein_accessions is stored as a Python-repr list string."""
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


def main() -> int:
    m = pd.read_csv(MOL)
    p = pd.read_csv(PRO)

    if "synonyms" not in m.columns:
        m["synonyms"] = ""
        print("[molecules] added 'synonyms' column")

    pairs = [
        ("Ala-CA", "L-alanocholate",      229710.0),
        ("Ser-CA", "L-serocholate",       229711.0),
        ("Phe-CA", "L-phenylalanocholate", 229714.0),
    ]

    drop_idx: list[int] = []

    for short, long, chebi in pairs:
        short_mask = m["compound_name"] == short
        long_mask = (m["compound_name"].astype(str).str.strip() == long)
        if short_mask.sum() != 1 or long_mask.sum() != 1:
            print(f"ERROR {short}/{long}: counts short={short_mask.sum()} long={long_mask.sum()}",
                  file=sys.stderr)
            return 1
        s_idx = m.index[short_mask][0]
        l_idx = m.index[long_mask][0]

        # Union paper_url, always add Rimal (Fig 2d evidence)
        merged_url = union_urls(
            str(m.at[l_idx, "paper_url"]) if pd.notna(m.at[l_idx, "paper_url"]) else "",
            str(m.at[s_idx, "paper_url"]) if pd.notna(m.at[s_idx, "paper_url"]) else "",
            RIMAL_URL,
        )
        m.at[l_idx, "paper_url"] = merged_url

        m.at[l_idx, "is_literature_recruited"] = 1
        m.at[l_idx, "synonyms"] = short
        m.at[l_idx, "interacting_protein_accessions"] = append_protein(
            m.at[l_idx, "interacting_protein_accessions"], "Q5LF84"
        )

        drop_idx.append(s_idx)
        print(f"[merge] {short} -> {long} (CHEBI:{int(chebi)})  "
              f"synonyms='{short}'  paper_url+='Rimal 2024'  prots+='Q5LF84'")

    # Tyr-CA: no canonical, just annotate
    t_mask = m["compound_name"] == "Tyr-CA"
    if t_mask.sum() == 1:
        t_idx = m.index[t_mask][0]
        m.at[t_idx, "interacting_protein_accessions"] = append_protein(
            m.at[t_idx, "interacting_protein_accessions"], "Q5LF84"
        )
        m.at[t_idx, "synonyms"] = "L-tyrocholate (unofficial; not in ChEBI snapshot)"
        print("[tyr-ca] prots+='Q5LF84'  synonyms='L-tyrocholate (unofficial...)'")
    else:
        print(f"ERROR Tyr-CA rows: {t_mask.sum()}", file=sys.stderr)
        return 1

    # Taurocholate: substrate of Q5LF84 in Fig 2d
    tc_mask = m["chebi_id"] == 36257.0
    if tc_mask.sum() == 1:
        tc_idx = m.index[tc_mask][0]
        m.at[tc_idx, "interacting_protein_accessions"] = append_protein(
            m.at[tc_idx, "interacting_protein_accessions"], "Q5LF84"
        )
        print("[taurocholate CHEBI:36257] prots+='Q5LF84'")
    else:
        print(f"ERROR taurocholate rows: {tc_mask.sum()}", file=sys.stderr)
        return 1

    # Drop paper-shorthand rows
    before = len(m)
    m = m.drop(index=drop_idx).reset_index(drop=True)
    print(f"[molecules] dropped {before - len(m)} paper-shorthand rows "
          f"({before} -> {len(m)})")

    # -------- proteins.csv : Q5LF84 --------
    q_mask = p["accession"] == "Q5LF84"
    if q_mask.sum() != 1:
        print(f"ERROR Q5LF84 rows in proteins.csv: {q_mask.sum()}", file=sys.stderr)
        return 1
    q_idx = p.index[q_mask][0]

    # Trim interacting_compounds to the 6 Fig-2d-supported entities
    new_compounds = ";".join([
        "taurocholate",
        "L-alanocholate",
        "L-phenylalanocholate",
        "L-serocholate",
        "Tyr-CA",
    ])
    p.at[q_idx, "interacting_compounds"] = new_compounds
    # interacting_chebi_ids: products + substrate (taurocholate)
    p.at[q_idx, "interacting_chebi_ids"] = "36257;229710;229711;229714"
    print(f"[Q5LF84] interacting_compounds = '{new_compounds}'")
    print(f"[Q5LF84] interacting_chebi_ids = '36257;229710;229711;229714'")

    # -------- write out --------
    m.to_csv(MOL, index=False)
    p.to_csv(PRO, index=False)
    print(f"\nWrote {MOL} ({len(m)} rows) and {PRO}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
