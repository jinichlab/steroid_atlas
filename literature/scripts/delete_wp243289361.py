"""Delete WP_243289361 (CpBSH/T) from the atlas.

The user reported that the Guzior 2024 paper attribution for this entry does
not hold up on re-audit — the protein isn't actually tied to that paper.
Removing it cleanly from every place it appears so a downstream re-projection
of the atlas can be run on the remaining 14 literature-recruited proteins.

Touched files:
  data/proteins.csv                                       (drop 1 row)
  data/literature_recruited_proteins.csv                  (drop 1 row)
  literature/sequences/literature_recruited.fasta         (drop 1 sequence)
  literature/sequences/literature_recruited_manifest.tsv  (drop 1 row)
  literature/embeddings/literature_recruited.h5           (drop 1 dataset)

Not touched (kept as historical provenance):
  literature/scripts/append_cpbsht_only.py
  literature/scripts/embed_cpbsht_and_project.py
  literature/scripts/fetch_guzior_sequences.py
  literature/scripts/parse_guzior_table1.py
  literature/scripts/annotate_guzior_evidence_tier.py
  literature/guzior2024_*.json
  literature/supplementary/Guzior2024_SupplTables.xlsx

molecules.csv is not touched — no molecule reverse-maps WP_243289361
(verified before running this script).
"""
from __future__ import annotations

import sys
from pathlib import Path

import h5py
import pandas as pd

ROOT = Path("/home/adsiordia/steroid_atlas")
TARGET = "WP_243289361"


def main() -> int:
    # 1. proteins.csv
    p_path = ROOT / "data" / "proteins.csv"
    p = pd.read_csv(p_path)
    before = len(p)
    p = p[p["accession"] != TARGET].reset_index(drop=True)
    print(f"[proteins.csv]  {before} -> {len(p)}  (removed {before - len(p)})")
    p.to_csv(p_path, index=False)

    # 2. literature_recruited_proteins.csv
    lrp_path = ROOT / "data" / "literature_recruited_proteins.csv"
    lrp = pd.read_csv(lrp_path)
    before = len(lrp)
    lrp = lrp[lrp["Entry"] != TARGET].reset_index(drop=True)
    print(f"[literature_recruited_proteins.csv]  {before} -> {len(lrp)}  "
          f"(removed {before - len(lrp)})")
    lrp.to_csv(lrp_path, index=False)

    # 3. FASTA
    fasta_path = ROOT / "literature" / "sequences" / "literature_recruited.fasta"
    text = fasta_path.read_text()
    records = text.split(">")[1:]  # first split entry is empty
    kept = []
    dropped = 0
    for rec in records:
        header_line = rec.splitlines()[0]
        if header_line.split("|")[0] == TARGET:
            dropped += 1
            continue
        kept.append(">" + rec)
    if dropped != 1:
        print(f"WARN: FASTA dropped={dropped} (expected 1)", file=sys.stderr)
    fasta_path.write_text("".join(kept))
    print(f"[FASTA]  {len(records)} -> {len(kept)}  (removed {dropped})")

    # 4. manifest.tsv
    mani_path = ROOT / "literature" / "sequences" / "literature_recruited_manifest.tsv"
    mani = pd.read_csv(mani_path, sep="\t")
    before = len(mani)
    # Figure out accession column
    acc_col = None
    for c in ("accession", "Entry", "entry", "id", "ID"):
        if c in mani.columns:
            acc_col = c
            break
    if acc_col is None:
        print(f"ERROR: no accession column in manifest "
              f"(columns: {list(mani.columns)})", file=sys.stderr)
        return 1
    mani = mani[mani[acc_col] != TARGET].reset_index(drop=True)
    print(f"[manifest.tsv]  {before} -> {len(mani)}  (removed {before - len(mani)})")
    mani.to_csv(mani_path, sep="\t", index=False)

    # 5. ProtT5 embeddings h5
    h5_path = ROOT / "literature" / "embeddings" / "literature_recruited.h5"
    with h5py.File(h5_path, "a") as f:
        before = len(f.keys())
        if TARGET in f:
            del f[TARGET]
        print(f"[H5]  {before} -> {len(f.keys())}  "
              f"(removed {TARGET!r})")

    print()
    print("Done. Downstream steps a human still needs to run:")
    print("  python web/scripts/build_atlas_data.py   # rebuild JSONs")
    print("  (later) re-project UMAP with 14 recruited proteins")
    return 0


if __name__ == "__main__":
    sys.exit(main())
