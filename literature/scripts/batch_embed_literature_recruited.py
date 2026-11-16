"""Batch-embed all 15 literature-recruited proteins into a single H5.

Reads sequences directly from data/proteins.csv (rows where
is_literature_recruited == True). No external FASTA fetch needed — the atlas
build already pulled and curated these sequences into the main table.

Each entry is stored under its bare accession (e.g. "P0DXD2", "WP_243289361",
"MFU7515415") as a (1024,) float32 dataset, with attributes:
    entry_name, organism, paper, species, length_aa, source_url, umap_1,
    umap_2, cluster, prott5_cluster, identifier_type

Outputs
-------
literature/embeddings/literature_recruited.h5
literature/embeddings/literature_recruited_manifest.tsv
literature/embeddings/batch_embed_report.json

Run on a GPU node:
    LD_LIBRARY_PATH=/home/adsiordia/miniconda3/lib \\
      python literature/scripts/batch_embed_literature_recruited.py
"""
from __future__ import annotations
import json
import re
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
from transformers import T5EncoderModel, T5Tokenizer

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
PROTEINS_CSV = ROOT / "data" / "proteins.csv"
OUT_DIR = ROOT / "literature" / "embeddings"
OUT_DIR.mkdir(parents=True, exist_ok=True)
H5_PATH = OUT_DIR / "literature_recruited.h5"
TSV_PATH = OUT_DIR / "literature_recruited_manifest.tsv"
REPORT_PATH = OUT_DIR / "batch_embed_report.json"

MODEL_NAME = "Rostlab/prot_t5_xl_half_uniref50-enc"


def preprocess_sequence(seq: str) -> str:
    seq = seq.replace("*", "").upper()
    seq = re.sub(r"[UZOB]", "X", seq)
    return " ".join(list(seq))


def main() -> None:
    print(f"Reading {PROTEINS_CSV}")
    prot = pd.read_csv(PROTEINS_CSV, low_memory=False)
    lit = prot[prot["is_literature_recruited"].fillna(False).astype(bool)].copy()
    assert len(lit) > 0, "no literature-recruited rows found"
    print(f"literature-recruited rows: {len(lit)}")
    print(lit[["accession", "length_aa", "organism"]].to_string(index=False))

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"\nDevice: {device}")
    print(f"Loading {MODEL_NAME} ...")
    t0 = time.time()
    model = T5EncoderModel.from_pretrained(MODEL_NAME)
    if device.type == "cpu":
        model = model.to(torch.float32)
    model = model.to(device).eval()
    tokenizer = T5Tokenizer.from_pretrained(MODEL_NAME, do_lower_case=False)
    print(f"  model loaded in {time.time()-t0:.1f}s")

    report_rows = []
    with h5py.File(H5_PATH, "w") as h5:
        h5.attrs["model"] = MODEL_NAME
        h5.attrs["pooling"] = "mean over residue tokens (EOS trimmed)"
        h5.attrs["dim"] = 1024
        h5.attrs["dtype"] = "float32"
        h5.attrs["source_csv"] = str(PROTEINS_CSV.relative_to(ROOT))
        h5.attrs["created"] = time.strftime("%Y-%m-%d %H:%M:%S")

        for _, row in lit.iterrows():
            acc = str(row["accession"])
            seq = str(row.get("sequence", "")).strip()
            if not seq or seq.lower() == "nan":
                print(f"  SKIP {acc}: no sequence")
                report_rows.append({"accession": acc, "status": "skipped_no_seq"})
                continue
            seq_spaced = preprocess_sequence(seq)
            t0 = time.time()
            with torch.no_grad():
                tok = tokenizer(
                    [seq_spaced], add_special_tokens=True, padding=True,
                    return_tensors="pt",
                ).to(device)
                out = model(input_ids=tok["input_ids"],
                            attention_mask=tok["attention_mask"])
                per_res = out.last_hidden_state[0, :len(seq)].cpu().numpy().astype(np.float32)
                per_prot = per_res.mean(axis=0)  # (1024,)
            dt = time.time() - t0

            dset = h5.create_dataset(acc, data=per_prot)
            for key in ("entry_name", "protein_names", "gene_names", "organism",
                        "length_aa", "umap_1", "umap_2", "cluster",
                        "prott5_cluster", "identifier_type", "sequence_source",
                        "paper_url"):
                val = row.get(key, None)
                if val is None or (isinstance(val, float) and np.isnan(val)):
                    continue
                dset.attrs[key] = val if not isinstance(val, (np.integer, np.floating)) else val.item()

            report_rows.append({
                "accession": acc,
                "length_aa": int(row.get("length_aa", len(seq))),
                "elapsed_s": round(dt, 2),
                "status": "ok",
                "norm": float(np.linalg.norm(per_prot)),
            })
            print(f"  {acc:14s}  len={len(seq):4d}  {dt:5.2f}s  ||v||={np.linalg.norm(per_prot):.3f}")

    # Manifest TSV — easy human/reviewer inspection
    man = lit[[
        "accession", "entry_name", "protein_names", "organism", "length_aa",
        "umap_1", "umap_2", "cluster", "prott5_cluster", "identifier_type",
        "sequence_source", "paper_url",
    ]].copy()
    man.to_csv(TSV_PATH, sep="\t", index=False)

    REPORT_PATH.write_text(json.dumps({
        "model": MODEL_NAME,
        "n_entries": len(report_rows),
        "n_ok": sum(1 for r in report_rows if r["status"] == "ok"),
        "entries": report_rows,
        "h5_path": str(H5_PATH.relative_to(ROOT)),
    }, indent=2))

    print(f"\nWrote: {H5_PATH} ({H5_PATH.stat().st_size/1024:.1f} KB)")
    print(f"Wrote: {TSV_PATH}")
    print(f"Wrote: {REPORT_PATH}")


if __name__ == "__main__":
    main()
