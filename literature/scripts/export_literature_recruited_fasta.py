"""Export the 15 literature-recruited protein sequences from data/proteins.csv
to a single consolidated FASTA — the canonical source for anyone wanting to
re-embed, BLAST, or cite these sequences.

Output:
    literature/sequences/literature_recruited.fasta
    literature/sequences/literature_recruited_manifest.tsv
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

ROOT = Path("/home/adsiordia/steroid_atlas")
PROTEINS_CSV = ROOT / "data" / "proteins.csv"
OUT_FASTA = ROOT / "literature" / "sequences" / "literature_recruited.fasta"
OUT_MANIFEST = ROOT / "literature" / "sequences" / "literature_recruited_manifest.tsv"


def main() -> None:
    prot = pd.read_csv(PROTEINS_CSV, low_memory=False)
    lit = prot[prot["is_literature_recruited"].fillna(False).astype(bool)].copy()
    lit = lit.sort_values("accession").reset_index(drop=True)
    print(f"literature-recruited rows: {len(lit)}")

    missing = []
    records = []
    for _, r in lit.iterrows():
        acc = str(r["accession"])
        seq = str(r.get("sequence", "")).strip()
        if not seq or seq.lower() == "nan":
            missing.append(acc)
            continue
        entry_name = str(r.get("entry_name") or "").strip()
        gene = str(r.get("gene_names") or "").strip().replace(" ", "_")
        organism = str(r.get("organism") or "").strip().replace(" ", "_")
        header = f">{acc}|{entry_name}|{gene}|{organism}|len={len(seq)}"
        records.append((header, seq))

    OUT_FASTA.parent.mkdir(parents=True, exist_ok=True)
    with OUT_FASTA.open("w") as f:
        for header, seq in records:
            f.write(header + "\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

    manifest = lit[[
        "accession", "entry_name", "protein_names", "gene_names", "organism",
        "length_aa", "umap_1", "umap_2", "cluster", "prott5_cluster",
        "identifier_type", "sequence_source", "paper_url",
    ]].copy()
    manifest.to_csv(OUT_MANIFEST, sep="\t", index=False)

    print(f"wrote {len(records)} records to {OUT_FASTA}  "
          f"({OUT_FASTA.stat().st_size/1024:.1f} KB)")
    print(f"wrote manifest to {OUT_MANIFEST}")
    if missing:
        print(f"MISSING sequences (not written): {missing}")


if __name__ == "__main__":
    main()
