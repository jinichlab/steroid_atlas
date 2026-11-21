"""NCBI Taxonomy lookup for every organism in the protein atlas.

Produces authoritative taxonomic lineage (superkingdom → species) for every
protein by querying NCBI Taxonomy via Entrez. Reusable cache + per-protein
TSV let any downstream script (the Panels C/D classifier, the combined
taxomimic figure, future case studies) swap the hand-curated genus-set
classifier for NCBI-authoritative lineages.

Outputs (under analysis/polish/):
  ncbi_taxonomy_cache.json          binomial → {taxid, ranks…}
  ncbi_taxonomy_unresolved.tsv      organism strings that didn't resolve
  proteins_ncbi_taxonomy.tsv        per-protein lineage table for joining
                                     back to proteins.csv on accession
"""
from __future__ import annotations
import re, time, json, sys
from pathlib import Path
import pandas as pd
from Bio import Entrez
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUT  = HERE

Entrez.email = "research@example.com"   # generic — not the user's address
REQ_SLEEP   = 0.34                       # ≈3 req/sec without API key
BATCH_EFETCH = 100
CACHE_PATH   = OUT / "ncbi_taxonomy_cache.json"

RANKS = ["superkingdom", "kingdom", "phylum", "class",
         "order", "family", "genus", "species"]

# ── Helpers ────────────────────────────────────────────────────────
def binomial(org: str) -> str:
    """Extract 'Genus species' from a messy organism string."""
    if not org: return ""
    s = re.sub(r"\s*\([^)]*\)\s*", " ", org).strip()
    parts = s.split()
    return f"{parts[0]} {parts[1]}" if len(parts) >= 2 else s

def load_cache() -> dict:
    if not CACHE_PATH.exists():
        return {}
    try:
        raw = json.loads(CACHE_PATH.read_text())
    except Exception:
        return {}
    # Migrate legacy flat format ({binomial: taxid_string}) to new dict format
    migrated = {}
    for k, v in raw.items():
        if isinstance(v, dict):
            migrated[k] = v
        elif isinstance(v, str) or v is None:
            migrated[k] = {"taxid": v, "resolved_via": "binomial" if v else None}
        else:
            migrated[k] = {"taxid": None, "resolved_via": None}
    return migrated

def save_cache(cache: dict) -> None:
    CACHE_PATH.write_text(json.dumps(cache, indent=1))

def esearch_taxid(query: str) -> str | None:
    try:
        h = Entrez.esearch(db="taxonomy", term=query)
        r = Entrez.read(h); h.close()
        return r["IdList"][0] if r["IdList"] else None
    except Exception as e:
        print(f"    esearch error for {query!r}: {e}", file=sys.stderr)
        return None

def parse_lineage(taxon_elt) -> dict:
    """Return {rank: name} for a <Taxon> element from efetch."""
    out = {}
    # The queried taxon itself carries a rank (often "species" but may be
    # "genus" when the fallback hit the genus).
    rank = taxon_elt.findtext("Rank")
    name = taxon_elt.findtext("ScientificName")
    if rank in RANKS and name:
        out[rank] = name
    for lt in taxon_elt.findall("LineageEx/Taxon"):
        r = lt.findtext("Rank"); n = lt.findtext("ScientificName")
        if r in RANKS and n:
            out[r] = n
    return out

# ── Main ───────────────────────────────────────────────────────────
print("Loading proteins.csv …")
prot = pd.read_csv(ROOT / "data" / "proteins.csv", low_memory=False)
prot["organism"] = prot["organism"].fillna("").astype(str)
prot["binomial"] = prot["organism"].apply(binomial)

cache = load_cache()
unique_binoms = sorted({b for b in prot["binomial"].unique() if b})
print(f"Unique binomials: {len(unique_binoms)}   "
      f"(already cached: {sum(1 for b in unique_binoms if b in cache)})")

# ── Step 1: esearch each new binomial for a taxid ──────────────────
to_query = [b for b in unique_binoms
            if b not in cache or not cache[b].get("taxid")]
print(f"\nStep 1: esearch {len(to_query)} new binomials "
      f"(~{len(to_query)*REQ_SLEEP/60:.1f} min)")
for i, b in enumerate(to_query, 1):
    taxid = esearch_taxid(b)
    if taxid is None:
        # Fallback — try the genus only
        genus = b.split()[0]
        if genus and genus != b:
            time.sleep(REQ_SLEEP)
            taxid = esearch_taxid(genus)
            if taxid:
                cache[b] = {"taxid": taxid, "resolved_via": f"genus({genus})"}
    else:
        cache[b] = {"taxid": taxid, "resolved_via": "binomial"}
    if b not in cache:
        cache[b] = {"taxid": None, "resolved_via": None}
    time.sleep(REQ_SLEEP)
    if i % 50 == 0:
        print(f"  {i}/{len(to_query)}  (resolved so far: "
              f"{sum(1 for v in cache.values() if v.get('taxid'))})")
        save_cache(cache)
save_cache(cache)

# ── Step 2: efetch lineages in batches, update cache with ranks ────
need_lineage = [b for b, v in cache.items()
                if v.get("taxid") and not v.get("species") and not v.get("order")
                and not v.get("class")]
all_taxids = sorted({cache[b]["taxid"] for b in need_lineage})
if all_taxids:
    print(f"\nStep 2: efetch lineages for {len(all_taxids)} taxids")
    taxid_to_lineage: dict[str, dict] = {}
    for i in range(0, len(all_taxids), BATCH_EFETCH):
        batch = all_taxids[i:i+BATCH_EFETCH]
        try:
            h = Entrez.efetch(db="taxonomy", id=",".join(batch), retmode="xml")
            data = h.read(); h.close()
        except Exception as e:
            print(f"    efetch batch error: {e}", file=sys.stderr); continue
        root = ET.fromstring(data)
        for taxon in list(root):     # direct children of <TaxaSet>
            tid = taxon.findtext("TaxId")
            if tid:
                taxid_to_lineage[tid] = parse_lineage(taxon)
        time.sleep(REQ_SLEEP)
        print(f"  {min(i+BATCH_EFETCH, len(all_taxids))}/{len(all_taxids)}")
    # Merge lineages into the cache keyed by binomial
    for b, v in cache.items():
        tid = v.get("taxid")
        if tid and tid in taxid_to_lineage:
            v.update(taxid_to_lineage[tid])
    save_cache(cache)

# ── Reports ────────────────────────────────────────────────────────
resolved = {b: v for b, v in cache.items()
            if v.get("taxid") and (v.get("order") or v.get("species")
                                   or v.get("class"))}
unresolved = [b for b in unique_binoms if b not in resolved]
print(f"\nResolved {len(resolved)}/{len(unique_binoms)} "
      f"({len(resolved)/max(len(unique_binoms),1)*100:.1f}%)")
print(f"Unresolved: {len(unresolved)}")

if unresolved:
    pd.DataFrame({"binomial": unresolved}).to_csv(
        OUT / "ncbi_taxonomy_unresolved.tsv", sep="\t", index=False)
    print(f"  wrote {OUT / 'ncbi_taxonomy_unresolved.tsv'}")

# ── Build per-protein TSV ──────────────────────────────────────────
print("\nBuilding per-protein taxonomy TSV …")
rows = []
for _, r in prot.iterrows():
    b = r["binomial"]
    v = cache.get(b, {})
    row = {"accession": r["accession"],
           "organism":  r["organism"],
           "binomial":  b,
           "ncbi_taxid": v.get("taxid"),
           "resolved_via": v.get("resolved_via")}
    for rk in RANKS:
        row[f"ncbi_{rk}"] = v.get(rk)
    rows.append(row)
out_df = pd.DataFrame(rows)
tsv_path = OUT / "proteins_ncbi_taxonomy.tsv"
out_df.to_csv(tsv_path, sep="\t", index=False)
print(f"  wrote {tsv_path}  ({len(out_df)} rows)")

# ── Quick coverage summary per rank ────────────────────────────────
print("\nPer-rank coverage (fraction of proteins with that rank filled):")
for rk in RANKS:
    col = f"ncbi_{rk}"
    cov = out_df[col].notna().sum() / len(out_df) * 100
    print(f"  {rk:<13} {cov:5.1f}%")
