"""Fetch UniProt metadata for the 729 missing-candidate accessions surfaced
by script 20 (umbrella-GO coverage estimate). Categorize by likely
biological role so we can judge which are worth adding to the atlas vs
which are peripheral / regulatory / edge cases.

Uses the UniProt REST /uniprotkb/accessions batch endpoint (500-per-call).

Writes:
    analysis/umbrella_go_missing_details.tsv     one row per candidate
    analysis/umbrella_go_missing_summary.md      grouped summary + review notes
"""
from __future__ import annotations

import time
from pathlib import Path

import requests
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE

IN_TSV = HERE / "umbrella_go_missing.tsv"
BATCH = 500
FIELDS = (
    "accession,protein_name,gene_names,organism_name,length,ec,"
    "cc_function,go,keyword,cc_subcellular_location"
)
URL = "https://rest.uniprot.org/uniprotkb/accessions"


def fetch_batch(accs: list[str]) -> pd.DataFrame:
    params = {
        "accessions": ",".join(accs),
        "format": "tsv",
        "fields": FIELDS,
    }
    r = requests.get(URL, params=params, timeout=180)
    r.raise_for_status()
    from io import StringIO
    return pd.read_csv(StringIO(r.text), sep="\t")


accs = pd.read_csv(IN_TSV, sep="\t")["accession"].astype(str).tolist()
print(f"Fetching metadata for {len(accs):,} candidate accessions ...")
parts = []
for i in range(0, len(accs), BATCH):
    chunk = accs[i:i + BATCH]
    print(f"  batch {i // BATCH + 1}/{(len(accs) + BATCH - 1) // BATCH} "
          f"({len(chunk)} entries) ...", end=" ")
    df = fetch_batch(chunk)
    print(f"got {len(df)} rows")
    parts.append(df)
    time.sleep(0.4)

meta = pd.concat(parts, ignore_index=True)
meta = meta.fillna("")
print(f"\nFetched {len(meta):,} rows total")

# Categorize with simple heuristics so we can see the shape
def categorize(row) -> str:
    name = str(row.get("Protein names", "")).lower()
    kw = str(row.get("Keywords", "")).lower()
    fn = str(row.get("Function [CC]", "")).lower()
    ec = str(row.get("EC number", "")).strip()

    # 1) Peripheral / regulatory / cofactor donors
    if any(t in name for t in [
        "cytochrome p450 reductase", "cytochrome b5 reductase",
        "ferredoxin reductase", "adrenodoxin"
    ]):
        return "Peripheral · electron-donor"
    if any(t in name for t in [
        "steroidogenic acute regulatory",
        "star-related lipid transfer",
    ]):
        return "Peripheral · lipid transfer"
    if any(t in kw for t in [
        "transcription regulation", "dna-binding",
    ]) and "receptor" not in name:
        return "Regulatory · transcription"

    # 2) Genuine steroid enzyme (has an EC number & protein_name mentions a
    # steroid keyword)
    steroid_words = ("steroid", "sterol", "cholester", "bile", "androgen",
                     "estradiol", "estrogen", "cortisol", "progestin",
                     "progesterone", "aldos", "testoster", "hydroxysteroid",
                     "brassinost", "vitamin d", "calcitriol")
    if ec and any(w in name for w in steroid_words):
        return "Genuine · enzyme"
    if any(w in name for w in ["receptor"]) and \
            any(w in name for w in steroid_words):
        return "Genuine · receptor"
    if any(w in kw for w in ["transport"]) and \
            any(w in name for w in steroid_words):
        return "Genuine · transporter"
    if ec and not any(w in name for w in steroid_words):
        return "Enzyme · non-steroid name"

    # 3) Everything else — likely edge case
    return "Other · needs review"


meta["category"] = meta.apply(categorize, axis=1)
meta.to_csv(OUT_DIR / "umbrella_go_missing_details.tsv",
            sep="\t", index=False)
print(f"\nWrote analysis/umbrella_go_missing_details.tsv "
      f"({len(meta):,} rows)")

# Summary
print("\nBreakdown by category:")
cat_counts = meta["category"].value_counts()
for cat, n in cat_counts.items():
    print(f"  {cat:35s}  {n}")

print("\nTop 15 organisms:")
print(meta["Organism"].value_counts().head(15).to_string())

print("\nTop 15 EC-numbered proteins (by protein_name stem):")
ec_only = meta[meta["EC number"].str.strip().ne("")]
print(ec_only["Protein names"].str.slice(0, 60).value_counts().head(15).to_string())

# Markdown summary
lines = [
    "# Missing-candidate lookup — 729 accessions with steroid-related GO terms",
    "",
    "Fetched via UniProt `/uniprotkb/accessions` batch endpoint. All entries "
    "are Swiss-Prot annotation-score 5 and carry at least one of the four "
    "umbrella GO terms (steroid binding, steroid metabolic process, sterol "
    "metabolic process, bile acid transport), but were NOT recruited into "
    "the current 14,089-entry atlas via Rhea catalysis or literature audit.",
    "",
    "## Category breakdown",
    "",
    "| Category | Count | Notes |",
    "|---|---:|---|",
]
cat_notes = {
    "Genuine · enzyme":         "steroid enzyme with a steroid-named EC. High-priority recruit.",
    "Genuine · receptor":       "steroid receptor (nuclear or membrane). High-priority recruit.",
    "Genuine · transporter":    "steroid / bile-acid transporter. High-priority recruit.",
    "Enzyme · non-steroid name":"has EC but name doesn't mention a steroid. Needs case-by-case review — could be a general enzyme with a steroid-related side activity.",
    "Peripheral · electron-donor": "cytochrome P450 reductase, adrenodoxin, etc. — supply reducing equivalents to steroidogenic CYPs but do not bind steroids themselves. Debatable whether to include.",
    "Peripheral · lipid transfer": "StAR, StARD proteins — traffic cholesterol but not classical enzymes. Debatable.",
    "Regulatory · transcription":  "transcription factors regulating steroidogenic genes. Not steroid-interacting proteins by our definition.",
    "Other · needs review":       "does not match any of the above heuristics — likely a genuine edge case worth manual review.",
}
for cat, n in cat_counts.items():
    lines.append(f"| {cat} | {n} | {cat_notes.get(cat, '')} |")

lines += [
    "",
    "## Files",
    "- `analysis/umbrella_go_missing_details.tsv` — full 729-row metadata with a `category` column, ready for filtering / review before deciding what to recruit.",
]
(OUT_DIR / "umbrella_go_missing_summary.md").write_text(
    "\n".join(lines), encoding="utf-8"
)
print("\nWrote analysis/umbrella_go_missing_summary.md")
print("\nDone.")
