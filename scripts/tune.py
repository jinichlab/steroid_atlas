"""
tune_threshold.py — find a sensible MIN_SCORE for the chatbot's retriever.

Runs three sets of queries against your index:

  NAME LOOKUPS — compound names sampled straight from your own catalog.
                 These are the EASY case: the query is nearly a restatement
                 of the document title, so scores run high. Tuning on these
                 alone gives a threshold that's too strict.
  REALISTIC    — conceptual questions that name no compound. This is what
                 people actually type, and it's the group your threshold
                 has to accommodate.
  NEGATIVES    — questions with no possible answer in a steroid corpus.
                 They never score zero; how high they get is the noise floor.

MIN_SCORE goes between the negatives and the REALISTIC minimum — not the
name-lookup minimum.

Also reports near-ties: cases where the top two hits are different compounds
with almost identical scores, meaning retrieval is effectively a coin flip
and the chat may cite the wrong molecule.

    python scripts/tune_threshold.py
    python scripts/tune_threshold.py --n 12 --k 6
"""

import argparse
import json
import os
import random
import sys
from pathlib import Path

import faiss
import numpy as np
from dotenv import load_dotenv
from openai import OpenAI

ROOT = Path(__file__).resolve().parent.parent
STORE = ROOT / "data" / "rag_store"
MODEL = "text-embedding-3-large"

# Conceptual questions — no compound named. Edit these to match the kinds of
# things you and your lab actually ask; this group drives the threshold.
REALISTIC = [
    "why do gut bacteria deconjugate bile salts",
    "what makes a steroid an agonist versus an antagonist",
    "how does hydroxylation change a steroid's solubility",
    "which enzymes act on cholesterol first",
    "what is the difference between a primary and a secondary bile acid",
    "how are steroid hormones transported in blood",
    "what role do sterols play in membrane fluidity",
    "how does conjugation affect a bile acid's function",
]

NEGATIVES = [
    "how do I fix a flat bicycle tyre",
    "what is the offside rule in football",
    "best way to cook risotto",
    "explain the French Revolution",
    "how do I renew a passport",
    "what is the capital of Peru",
]

NEAR_TIE = 0.02      # top-2 within this margin = effectively indistinguishable

load_dotenv(ROOT / ".env")


def embed(texts, client):
    resp = client.embeddings.create(model=MODEL, input=texts)
    arr = np.array([d.embedding for d in resp.data], dtype="float32")
    arr /= np.linalg.norm(arr, axis=1, keepdims=True) + 1e-12
    return arr


def run_group(name, queries, index, catalog, client, k, show_spread=True):
    print(f"── {name} ──")
    vecs = embed(queries, client)
    scores, idxs = index.search(vecs, k)
    tops, ties = [], []
    for q, srow, irow in zip(queries, scores, idxs):
        top = float(srow[0])
        tops.append(top)
        hit = catalog[irow[0]]["paper"] if irow[0] != -1 else "—"
        print(f"  {top:.3f}  {q[:46]:48s} -> {hit[:30]}")
        if show_spread:
            print(f"         all-{k}: " + " ".join(f"{s:.3f}" for s in srow))
        # near-tie between two DIFFERENT source documents
        if len(srow) > 1 and irow[1] != -1:
            second = catalog[irow[1]]["paper"]
            if second != hit and (top - float(srow[1])) < NEAR_TIE:
                ties.append((q, hit, second, top, float(srow[1])))
    print()
    return tops, ties


def stats(label, xs):
    xs = sorted(xs)
    print(f"  {label:14s} min {xs[0]:.3f}  median {xs[len(xs)//2]:.3f}  max {xs[-1]:.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=10, help="name-lookup queries")
    ap.add_argument("--k", type=int, default=6, help="hits per query")
    args = ap.parse_args()

    key = os.getenv("OPENAI_API_KEY")
    if not key:
        sys.exit(f"No OPENAI_API_KEY in {ROOT / '.env'}")
    if not (STORE / "index.faiss").exists():
        sys.exit(f"No index at {STORE}. Run build_rag_index.py first.")

    index = faiss.read_index(str(STORE / "index.faiss"))
    catalog = [json.loads(l) for l in open(STORE / "catalog.jsonl", encoding="utf-8")]
    client = OpenAI(api_key=key)
    print(f"{index.ntotal:,} vectors · {len(catalog):,} catalog rows\n")

    labels = sorted({r["paper"] for r in catalog if not r["paper"].startswith("CHEBI")})
    if len(labels) < args.n:
        labels = sorted({r["paper"] for r in catalog})
    random.seed(0)
    lookups = [f"tell me about {n}" for n in random.sample(labels, min(args.n, len(labels)))]

    look_tops, look_ties = run_group("NAME LOOKUPS (easy case)", lookups,
                                     index, catalog, client, args.k)
    real_tops, real_ties = run_group("REALISTIC (drives the threshold)", REALISTIC,
                                     index, catalog, client, args.k)
    neg_tops, _ = run_group("NEGATIVES (noise floor)", NEGATIVES,
                            index, catalog, client, args.k, show_spread=False)

    print("── summary (top-1 score per query) ──")
    stats("name lookups", look_tops)
    stats("realistic", real_tops)
    stats("negatives", neg_tops)

    floor, ceiling = max(neg_tops), min(real_tops)
    print()
    if ceiling > floor:
        suggested = floor + (ceiling - floor) * 0.4      # bias toward permissive
        print(f"  usable window: {floor:.3f} (noise) .. {ceiling:.3f} (worst real question)")
        print(f"  -> MIN_SCORE = {suggested:.2f}")
        print("     Placed below the midpoint on purpose: a weak-but-real match is")
        print("     more useful than 'no relevant context', and the LLM is told to")
        print("     flag thin context anyway.")
    else:
        print(f"  NO GAP — worst realistic question {ceiling:.3f} sits at or below")
        print(f"  the noise floor {floor:.3f}. A fixed threshold can't separate them.")
        print("  Options:")
        print("    · re-index with CHUNK = 600 — long chunks dilute the signal")
        print("    · use the relative cutoff (keep hits within 75% of the top hit)")
        print("    · your corpus may not cover these concepts at all — check which")
        print("      documents the realistic queries actually matched above")

    all_ties = look_ties + real_ties
    print()
    if all_ties:
        print(f"── near-ties ({len(all_ties)}) — top two differ by < {NEAR_TIE} ──")
        for q, first, second, s1, s2 in all_ties:
            print(f"  {q[:42]:44s} {s1:.3f} {first[:24]}")
            print(f"  {'':44s} {s2:.3f} {second[:24]}")
        print("  Retrieval is a coin flip between these. Systematic IUPAC-style")
        print("  names look nearly identical to the embedding model. When the user")
        print("  has ticked a row, filter to that compound instead of trusting")
        print("  similarity — see the boost snippet in steroid_chat.")
    else:
        print("── no near-ties detected ──")


if __name__ == "__main__":
    main()