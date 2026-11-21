"""
build_atlas_compound_docs.py — RAG documents for atlas compounds that the
ChEBI scraper can't cover.

small_molecule_scraper.py is keyed on ChEBI IDs, so compounds without one (most
literature-recruited microbial bile-acid conjugates, e.g. Ala-HDCA) — or with a
class-level ChEBI term that PubChem/ChEBI won't resolve — never get a corpus
document. This writes one document per such compound from what the atlas
itself knows: name, synonyms, SMILES, evidence tier + note, suggested
producers, source papers and interacting proteins.

Reads:   data/molecules.csv, data/proteins.csv,
         rag_data/rag_train/small_molecule/data_txt/   (which ChEBI docs exist)
Writes:  rag_data/rag_train/small_molecule/atlas_txt/*.txt   (regenerated each run)

    python rag/build_atlas_compound_docs.py
"""

import re
import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MOLECULES = ROOT / "data" / "molecules.csv"
PROTEINS = ROOT / "data" / "proteins.csv"
CHEBI_DIR = ROOT / "rag_data" / "rag_train" / "small_molecule" / "data_txt"
OUT_DIR = ROOT / "rag_data" / "rag_train" / "small_molecule" / "atlas_txt"

MAX_PROTEINS = 30

TIER_LABELS = {
    "1": "chemistry reference only — synthetic standard, no biological detection",
    "2": "MS/MS match in public data (isomer uncertain)",
    "3": "detected in human IBD stool (iHMP2, isomer-resolved)",
    "4-5": "cross-cohort validated in human IBD (iHMP2 + Boland)",
    "5": "produced by gut bacteria in culture",
}


def chebi_slug(raw: str) -> str | None:
    s = re.sub(r"(?i)^chebi[:_\s]*", "", str(raw).strip().strip("[]'\" "))
    try:
        return f"CHEBI_{int(float(s))}"
    except ValueError:
        return None


def file_slug(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")[:80] or "unnamed"


def split_list(s: str) -> list[str]:
    return [x.strip() for x in re.split(r"[;|,\n]", s.strip("[]")) if x.strip(" '\"")]


def main():
    mol = pd.read_csv(MOLECULES, dtype=str, keep_default_na=False)
    prot = pd.read_csv(PROTEINS, dtype=str, keep_default_na=False,
                       usecols=["accession", "protein_names", "organism"])
    prot_by_acc = prot.set_index("accession").to_dict("index")
    scraped = {p.stem for p in CHEBI_DIR.glob("*.txt")}

    if OUT_DIR.exists():
        shutil.rmtree(OUT_DIR)
    OUT_DIR.mkdir(parents=True)

    written, used = 0, set()
    for _, m in mol.iterrows():
        name = m["compound_name"].strip()
        slug = chebi_slug(m["chebi_id"])
        if slug in scraped or not name:
            continue  # covered by the ChEBI scraper (or nothing to call it)

        lines = [f"Compound: {name}"]
        if m["chebi_id"].strip():
            lines.append(f"ChEBI ID: CHEBI:{slug.split('_')[1] if slug else m['chebi_id']}")
        synonyms = [s.replace("aka ", "").strip() for s in m.get("synonyms", "").split("|")]
        synonyms = [s for s in synonyms if s and s.lower() != name.lower()]
        if synonyms:
            lines.append("Also known as: " + "; ".join(synonyms))
        if m["smiles"]:
            lines.append(f"SMILES: {m['smiles']}")
        lines.append("Source: Steroid Atlas compound table"
                     + (" (literature-recruited)" if m["is_literature_recruited"] == "1" else ""))

        tier = m.get("evidence_tier", "").strip()
        if tier:
            label = TIER_LABELS.get(tier, "")
            lines += ["", "== Evidence ==",
                      f"Evidence tier {tier}" + (f": {label}" if label else "")]
        if m.get("evidence_note", "").strip():
            lines += ["", m["evidence_note"].strip()]
        if m.get("gentry_producers", "").strip():
            lines += ["", "Suggested bacterial producers (HMP isolate culture screen, "
                          "enzyme not identified): " + m["gentry_producers"].replace(";", "; ")]
        papers = [u.strip() for u in m["paper_url"].split("|") if u.strip()]
        if papers:
            lines += ["", "Source paper(s): " + " ; ".join(papers)]

        accs = split_list(m["interacting_protein_accessions"])
        if accs:
            lines += ["", f"== Interacting proteins in the atlas ({len(accs)}) =="]
            for a in accs[:MAX_PROTEINS]:
                p = prot_by_acc.get(a, {})
                desc = " · ".join(x for x in (p.get("protein_names", "")[:100],
                                              p.get("organism", "")[:60]) if x)
                lines.append(f"{a}" + (f" — {desc}" if desc else ""))
            if len(accs) > MAX_PROTEINS:
                lines.append(f"… and {len(accs) - MAX_PROTEINS} more")
        else:
            lines += ["", "No protein in the atlas has a direct substrate attribution "
                          "to this compound."]

        fn = file_slug(name)
        while fn in used:
            fn += "_"
        used.add(fn)
        (OUT_DIR / f"{fn}.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
        written += 1

    print(f"wrote {written} atlas compound docs -> {OUT_DIR}")


if __name__ == "__main__":
    main()
