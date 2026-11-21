"""
update_rag_store.py — bring the existing FAISS index in line with data/ without
re-embedding the whole corpus.

  · drops protein chunks whose accession is no longer in data/proteins.csv,
    and ChEBI compound chunks whose ChEBI ID is no longer in data/molecules.csv
  · embeds + appends ChEBI compound / protein docs not yet in the index
  · replaces all atlas compound docs (rag/build_atlas_compound_docs.py output;
    small, so they're simply re-embedded every run)

Vector i ↔ catalog line i is preserved: IndexFlatIP.remove_ids compacts in
order and the catalog is filtered with the same mask. The previous index and
catalog are kept as *.prev next to the new ones.

    python rag/update_rag_store.py            # shows the plan + cost, asks first
    python rag/update_rag_store.py --yes
"""

import argparse
import json
import os
import sys
from pathlib import Path

import faiss
import numpy as np
import pandas as pd
from openai import OpenAI

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_rag_store import (  # noqa: E402
    ATLAS_DIR, ATLAS_SECTION, DIM, MIN_CHUNK, MOLECULE_DIR, PROTEIN_DIR, ROOT, STORE,
    chunk, embed, load_names,
)


def chunks_for(path: Path, **meta) -> list[dict]:
    body = path.read_text(encoding="utf-8")
    if len(body) < MIN_CHUNK and meta.get("section") != ATLAS_SECTION:
        return []
    out = []
    for piece in chunk(body):
        # atlas docs are short by design — keep them even below MIN_CHUNK
        if len(piece) < MIN_CHUNK and meta.get("section") != ATLAS_SECTION:
            continue
        out.append({**meta, "page": None, "text": piece})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--yes", action="store_true", help="don't ask before embedding")
    args = ap.parse_args()

    key = os.getenv("OPENAI_API_KEY")
    if not key:
        sys.exit("No OPENAI_API_KEY in the environment.")

    mol = pd.read_csv(ROOT / "data" / "molecules.csv", dtype=str, keep_default_na=False)
    prot = pd.read_csv(ROOT / "data" / "proteins.csv", dtype=str,
                       keep_default_na=False, usecols=["accession"])
    want_prot = set(prot["accession"])
    want_chebi = set()
    for raw in mol["chebi_id"]:
        s = raw.strip().lower().replace("chebi:", "")
        try:
            want_chebi.add(f"CHEBI:{int(float(s))}")
        except ValueError:
            pass

    print("Loading index + catalog ...")
    index = faiss.read_index(str(STORE / "index.faiss"))
    with open(STORE / "catalog.jsonl", encoding="utf-8") as f:
        catalog = [json.loads(line) for line in f]
    assert index.ntotal == len(catalog), "index / catalog out of sync"

    # ── what to drop ──────────────────────────────────────────────────────
    def keep(r):
        if r["kind"] == "protein":
            return r["accession"] in want_prot
        if r["section"] == ATLAS_SECTION:
            return False  # always regenerated
        return r["chebi"] in want_chebi

    mask = np.array([keep(r) for r in catalog])
    dropped = [r for r, k in zip(catalog, mask) if not k]
    have_prot = {r["accession"] for r, k in zip(catalog, mask) if k and r["kind"] == "protein"}
    have_chebi = {r["chebi"] for r, k in zip(catalog, mask) if k and r["kind"] == "molecule"}

    # ── what to add ───────────────────────────────────────────────────────
    names = load_names()
    new = []
    for txt in sorted(MOLECULE_DIR.glob("*.txt")):
        chebi = txt.stem.replace("_", ":")
        if chebi in want_chebi and chebi not in have_chebi:
            new += chunks_for(txt, paper=names.get(txt.stem, txt.stem),
                              section="compound record", chebi=chebi,
                              accession=None, kind="molecule")
    for txt in sorted(PROTEIN_DIR.glob("*.txt")):
        if txt.stem in want_prot and txt.stem not in have_prot:
            new += chunks_for(txt, paper=txt.stem, section="protein record",
                              chebi=None, accession=txt.stem, kind="protein")
    for txt in sorted(ATLAS_DIR.glob("*.txt")):
        first = txt.read_text(encoding="utf-8").splitlines()[0]
        new += chunks_for(txt, paper=first.removeprefix("Compound: ").strip(),
                          section=ATLAS_SECTION, chebi=None, accession=None,
                          kind="molecule")

    def summary(rows):
        out = {}
        for r in rows:
            k = f"{r['kind']} / {r['section']}"
            out[k] = out.get(k, 0) + 1
        return out

    chars = sum(len(r["text"]) for r in new)
    print(f"\nindex: {index.ntotal:,} vectors")
    print(f"drop:  {len(dropped):,} chunks  {summary(dropped)}")
    for acc in sorted({r['accession'] for r in dropped if r['kind'] == 'protein'}):
        print(f"         protein {acc}")
    for c in sorted({r['chebi'] for r in dropped if r['kind'] == 'molecule' and r['chebi']}):
        print(f"         compound {c}")
    print(f"add:   {len(new):,} chunks  {summary(new)}")
    print(f"       {chars:,} chars · est. cost ${chars / 4 / 1e6 * 0.13:.3f}")
    if not dropped and not new:
        print("Nothing to do.")
        return
    if not args.yes and input("proceed? [y/N] ").strip().lower() != "y":
        sys.exit("aborted")

    vectors = embed([r["text"] for r in new], OpenAI(api_key=key)) if new else None
    if vectors is not None and vectors.shape[1] != DIM:
        sys.exit(f"Model returned {vectors.shape[1]} dims, expected {DIM}.")

    drop_ids = np.flatnonzero(~mask).astype("int64")
    if len(drop_ids):
        index.remove_ids(faiss.IDSelectorBatch(drop_ids))
    catalog = [r for r, k in zip(catalog, mask) if k]
    if vectors is not None:
        index.add(vectors)
        catalog += new
    assert index.ntotal == len(catalog)

    # write next to the old files, keep the old ones as *.prev
    tmp_idx, tmp_cat = STORE / "index.faiss.tmp", STORE / "catalog.jsonl.tmp"
    faiss.write_index(index, str(tmp_idx))
    with open(tmp_cat, "w", encoding="utf-8") as f:
        for r in catalog:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    for name, tmp in (("index.faiss", tmp_idx), ("catalog.jsonl", tmp_cat)):
        os.replace(STORE / name, STORE / f"{name}.prev")
        os.replace(tmp, STORE / name)

    print(f"\nwrote {index.ntotal:,} vectors -> {STORE}  (previous kept as *.prev)")
    print("Restart the RAG service (npm run dev:rag) to pick it up.")


if __name__ == "__main__":
    main()
