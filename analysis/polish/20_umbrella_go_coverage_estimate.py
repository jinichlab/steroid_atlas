"""Coverage estimate — how many extra steroid-related proteins does the
'umbrella GO term' pull surface versus the current 14,089-entry atlas?

Hits UniProt REST with:
  (go:0005496  OR  go:0008202  OR  go:0016125  OR  go:0015721)
  AND reviewed:true
  AND annotation_score:5
  AND (existence:1 OR existence:2)

UniProt propagates GO annotations along the ontology by default, so any
descendant term (e.g. cholesterol binding, aromatase activity) also
counts.

Reports:
    - total UniProt matches for each individual umbrella term
    - total UniProt matches for the union (all four combined)
    - how many are already in the atlas
    - how many are potentially MISSING
    - MISSING-accessions list  →  analysis/umbrella_go_missing.tsv

No writes to proteins.csv. Read-only recon.
"""
from __future__ import annotations

import time
from pathlib import Path

import requests
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUT_DIR = HERE

# Note: UniProt's `go:` search takes the NUMERIC id only (no `GO:` prefix)
UMBRELLA_TERMS = {
    "0005496": "steroid binding (molecular function)",
    "0008202": "steroid metabolic process",
    "0016125": "sterol metabolic process",
    "0015721": "bile acid and bile salt transport",
}

BASE_QUERY = "(reviewed:true) AND (annotation_score:5)"

SEARCH = "https://rest.uniprot.org/uniprotkb/search"
STREAM = "https://rest.uniprot.org/uniprotkb/stream"
HEADERS = {"Accept": "text/plain"}


def fetch_accessions(query: str) -> set[str]:
    """Stream every reviewed accession matching a UniProt query. Fast enough
    for tens of thousands of entries."""
    params = {"query": query, "format": "list"}
    r = requests.get(STREAM, params=params, headers=HEADERS, timeout=180)
    r.raise_for_status()
    accs = {line.strip() for line in r.text.splitlines() if line.strip()}
    return accs


def count_only(query: str) -> int:
    """Cheaper: HEAD the search endpoint to read the X-Total-Results
    header — one round-trip, no data downloaded."""
    params = {"query": query, "format": "list", "size": 1}
    r = requests.get(SEARCH, params=params, headers=HEADERS, timeout=60)
    r.raise_for_status()
    n = r.headers.get("X-Total-Results")
    return int(n) if n is not None else -1


# ── Per-umbrella-term counts (cheap) ─────────────────────────────────────
print("Per-umbrella-term UniProt counts (reviewed, annotation_score=5):")
per_term = {}
for go_id, label in UMBRELLA_TERMS.items():
    q = f"(go:{go_id}) AND {BASE_QUERY}"
    n = count_only(q)
    per_term[go_id] = n
    print(f"  GO:{go_id}  {n:>8,}   {label}")
    time.sleep(0.3)

union_q = (
    "(go:0005496 OR go:0008202 OR go:0016125 OR go:0015721) AND "
    f"{BASE_QUERY}"
)
n_union = count_only(union_q)
print(f"\nUnion of the 4 umbrella terms:  {n_union:,}")
print("Downloading the accession list ...")
union_accs = fetch_accessions(union_q)
print(f"  streamed {len(union_accs):,} accessions")

# ── Compare to atlas ─────────────────────────────────────────────────────
prot = pd.read_csv(ROOT / "data" / "proteins.csv", low_memory=False)
atlas = set(prot["accession"].astype(str))
print(f"\nAtlas size: {len(atlas):,} proteins")

already = union_accs & atlas
missing = union_accs - atlas
print(f"  In atlas already :  {len(already):,}")
print(f"  MISSING from atlas: {len(missing):,}")
print(f"  Atlas entries NOT covered by any umbrella term: "
      f"{len(atlas - union_accs):,}")

# Write the missing list for review
if missing:
    df = pd.DataFrame(sorted(missing), columns=["accession"])
    df.to_csv(OUT_DIR / "umbrella_go_missing.tsv", sep="\t", index=False)
    print(f"\nWrote analysis/umbrella_go_missing.tsv "
          f"({len(missing):,} rows)")

print("\nDone.")
