"""Sequence-similarity baseline for the atlas — Step 1 of the seq-vs-PLM story.

Two complementary characterizations of what pure sequence identity says
about the atlas, both computed on the current 35,117-protein set:

  A. PAIRWISE IDENTITY DISTRIBUTION
     - draw N random protein pairs (default 20 000)
     - local pairwise alignment (biopython PairwiseAligner, BLOSUM62)
     - report the distribution of percent identity across all sampled pairs
     - answers: "how sequence-diverse is the atlas?"

  B. REDUNDANCY VIA CD-HIT
     - cluster the atlas at 90 %, 70 %, 50 %, 40 % identity
     - report cluster count at each threshold
     - answers: "how much of the atlas would collapse if we dedup at X %?"
       — a canonical Table 1 entry that shows the atlas isn't just
       orthologous duplicates.

The figure produced is Fig X.A of the paper's sequence-vs-PLM section;
Fig X.B (the PLM comparison) is Step 2 of this narrative.

Reads:  ../data/proteins.csv
Writes: ../analysis/seq_identity_pairs.tsv          per-pair identity table
        ../analysis/seq_similarity_summary.txt      human-readable summary
        ../analysis/seq_similarity_figure.png       histogram + redundancy chart
        ../analysis/seq_similarity_figure.pdf       vector version
        ../analysis/tmp_cdhit/                      CD-HIT intermediate outputs
"""
from __future__ import annotations

import random
import subprocess
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
IN = ROOT / "data" / "proteins.csv"
OUT_TSV = HERE / "seq_identity_pairs.tsv"
OUT_TXT = HERE / "seq_similarity_summary.txt"
OUT_PNG = HERE / "seq_similarity_figure.png"
OUT_PDF = HERE / "seq_similarity_figure.pdf"
TMP = HERE / "tmp_cdhit"

N_PAIRS = 20000                 # random pairs to align
MAX_LEN = 2000                  # skip pairs where either seq is longer
CDHIT_THRESHOLDS = [0.90, 0.70, 0.50, 0.40]
RANDOM_SEED = 42


def _local_identity(aligner, s1: str, s2: str) -> float | None:
    """Local alignment identity = matches / alignment length. Returns None if either seq too long or empty."""
    if not s1 or not s2 or len(s1) > MAX_LEN or len(s2) > MAX_LEN:
        return None
    try:
        aln = aligner.align(s1, s2)[0]
    except Exception:
        return None
    counts = aln.counts()
    # counts.gaps + counts.mismatches + counts.identities = alignment length
    aln_len = counts.gaps + counts.mismatches + counts.identities
    if aln_len == 0:
        return None
    return counts.identities / aln_len


def compute_pairwise(seqs: list[str]) -> pd.DataFrame:
    print(f"\n[A] Sampling {N_PAIRS:,} random pairs and running local alignment...")
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1

    rng = random.Random(RANDOM_SEED)
    n = len(seqs)
    rows = []
    t0 = time.time()
    for k in range(N_PAIRS):
        i, j = rng.sample(range(n), 2)
        ident = _local_identity(aligner, seqs[i], seqs[j])
        if ident is not None:
            rows.append({"i": i, "j": j, "len_i": len(seqs[i]), "len_j": len(seqs[j]),
                         "identity": ident})
        if (k + 1) % 500 == 0:
            elapsed = time.time() - t0
            rate = (k + 1) / elapsed
            eta_s = (N_PAIRS - k - 1) / rate
            print(f"  {k+1:>6}/{N_PAIRS}  ({rate:.1f} pair/s, ETA {eta_s/60:.1f} min)", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"  → {len(df):,} successful alignments in {(time.time()-t0)/60:.1f} min")
    return df


def _write_fasta(accessions, seqs, path: Path):
    with path.open("w") as f:
        for acc, seq in zip(accessions, seqs):
            f.write(f">{acc}\n{seq}\n")


def _cdhit_word_size(thr: float) -> int:
    if thr >= 0.70: return 5
    if thr >= 0.60: return 4
    if thr >= 0.50: return 3
    return 2   # 0.40-0.50


def run_cdhit(accessions: list[str], seqs: list[str]) -> list[dict]:
    print(f"\n[B] Running CD-HIT at thresholds {CDHIT_THRESHOLDS}...")
    TMP.mkdir(exist_ok=True)
    fasta = TMP / "atlas.fasta"
    _write_fasta(accessions, seqs, fasta)
    print(f"  wrote FASTA: {fasta}  ({len(seqs):,} sequences)")

    results = []
    for thr in CDHIT_THRESHOLDS:
        out = TMP / f"cdhit_{int(thr*100)}"
        w = _cdhit_word_size(thr)
        cmd = ["cd-hit", "-i", str(fasta), "-o", str(out),
               "-c", str(thr), "-n", str(w), "-M", "8000", "-T", "8",
               "-d", "0"]
        print(f"  cd-hit -c {thr} -n {w} ...", flush=True)
        t0 = time.time()
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=1800)
        except subprocess.CalledProcessError as e:
            print(f"  ! cd-hit at {thr} failed: {e.stderr[:400]}")
            continue
        except subprocess.TimeoutExpired:
            print(f"  ! cd-hit at {thr} timed out (>30 min); skipping")
            continue
        # cd-hit outputs a .clstr file listing all clusters
        clstr = out.with_suffix(out.suffix + ".clstr")
        if not clstr.exists():
            clstr = Path(str(out) + ".clstr")
        n_clusters = 0
        with clstr.open() as f:
            for line in f:
                if line.startswith(">Cluster"):
                    n_clusters += 1
        dt = time.time() - t0
        print(f"    → {n_clusters:,} clusters at {int(thr*100)}% identity   ({dt:.1f}s)")
        results.append({
            "threshold_pct": int(thr * 100),
            "n_clusters": n_clusters,
            "reduction_pct": 100 * (1 - n_clusters / len(seqs)),
        })
    return results


def write_summary(pairs_df: pd.DataFrame, cdhit_res: list[dict], n_total: int) -> None:
    q = pairs_df["identity"].quantile([0.05, 0.25, 0.5, 0.75, 0.95])
    thresholds = [0.30, 0.50, 0.70, 0.90]
    frac_above = {t: (pairs_df["identity"] > t).mean() for t in thresholds}

    lines = []
    lines.append("=" * 74)
    lines.append(f"SEQUENCE-SIMILARITY BASELINE — {n_total:,} atlas proteins")
    lines.append("=" * 74)
    lines.append("")
    lines.append(f"[A] Pairwise identity distribution ({len(pairs_df):,} random pairs, "
                 "BLOSUM62 local alignment):")
    lines.append(f"    min      {pairs_df['identity'].min():.3f}")
    lines.append(f"    5%       {q[0.05]:.3f}")
    lines.append(f"    25%      {q[0.25]:.3f}")
    lines.append(f"    median   {q[0.50]:.3f}")
    lines.append(f"    75%      {q[0.75]:.3f}")
    lines.append(f"    95%      {q[0.95]:.3f}")
    lines.append(f"    max      {pairs_df['identity'].max():.3f}")
    lines.append(f"    mean     {pairs_df['identity'].mean():.3f}")
    lines.append("")
    lines.append(f"    Fraction of pairs above threshold:")
    for t in thresholds:
        lines.append(f"      > {int(t*100)}%   {frac_above[t]*100:5.2f}%")
    lines.append("")

    if cdhit_res:
        lines.append(f"[B] CD-HIT redundancy at multiple identity thresholds:")
        lines.append(f"    {'threshold':>10} {'clusters':>12} {'reduction':>12}")
        for r in cdhit_res:
            lines.append(f"      {r['threshold_pct']:>3}%    "
                         f"{r['n_clusters']:>12,}     {r['reduction_pct']:>6.1f}%")
        lines.append("")
        lines.append(f"    (Atlas has {n_total:,} sequences; the cluster count is the number of")
        lines.append("     non-redundant sequences that would remain after collapsing at each threshold.)")
    OUT_TXT.write_text("\n".join(lines) + "\n")
    print(f"\nWrote {OUT_TXT.name}")


def make_figure(pairs_df: pd.DataFrame, cdhit_res: list[dict], n_total: int) -> None:
    fig = plt.figure(figsize=(11.0, 4.5), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.4, 1])
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    BAR = "#1F4B99"
    REF = "#9CA3AF"

    # ── Panel A — pairwise identity distribution ──────────────────────────
    axA.hist(pairs_df["identity"] * 100, bins=np.arange(0, 105, 2.5),
             color=BAR, alpha=0.85, edgecolor="white", linewidth=0.5)
    axA.axvline(30, color=REF, linestyle="--", linewidth=0.9)
    axA.text(30, axA.get_ylim()[1] * 0.95, "  30% identity",
             color=REF, fontsize=9, va="top", ha="left")
    axA.set_xlabel("Pairwise sequence identity  (%, local alignment)", fontsize=10)
    axA.set_ylabel(f"Number of pairs (of {len(pairs_df):,})", fontsize=10)
    axA.set_title("A. Sequence-identity distribution across atlas",
                  fontsize=11, loc="left", pad=8)
    axA.grid(True, axis="y", alpha=0.25, linewidth=0.5)
    axA.spines["top"].set_visible(False)
    axA.spines["right"].set_visible(False)
    axA.tick_params(labelsize=9)

    # ── Panel B — CD-HIT redundancy ───────────────────────────────────────
    if cdhit_res:
        thrs = [r["threshold_pct"] for r in cdhit_res]
        counts = [r["n_clusters"] for r in cdhit_res]
        ypos = np.arange(len(thrs))[::-1]
        axB.barh(ypos, counts, color=BAR, height=0.6,
                 edgecolor="white", linewidth=0.4)
        for y, r in zip(ypos, cdhit_res):
            axB.text(r["n_clusters"] * 1.02, y,
                     f" {r['n_clusters']:,}  ({r['reduction_pct']:.1f}% reduction)",
                     va="center", ha="left", fontsize=9, color="#111827")
        axB.set_yticks(ypos)
        axB.set_yticklabels([f"CD-HIT {t}%" for t in thrs], fontsize=9)
        axB.axvline(n_total, color=REF, linestyle="--", linewidth=0.9)
        axB.text(n_total, len(thrs) - 0.4,
                 f"  {n_total:,} total atlas", color=REF,
                 fontsize=8.5, va="bottom", ha="left")
        axB.set_xlabel("Non-redundant sequence count", fontsize=10)
        axB.set_title("B. Redundancy at multiple identity thresholds",
                      fontsize=11, loc="left", pad=8)
        axB.set_xlim(0, n_total * 1.35)
        axB.tick_params(axis="x", labelsize=9)
        axB.tick_params(axis="y", length=0, pad=4)
        axB.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
        axB.grid(True, axis="x", alpha=0.25, linewidth=0.5)
        axB.spines["top"].set_visible(False)
        axB.spines["right"].set_visible(False)
        axB.spines["left"].set_visible(False)
    else:
        axB.text(0.5, 0.5, "CD-HIT step failed — see log",
                 ha="center", va="center", transform=axB.transAxes, fontsize=10)

    for f in (OUT_PNG, OUT_PDF):
        fig.savefig(f, dpi=300, bbox_inches="tight")
        print(f"Wrote {f.name}")


def main() -> int:
    print(f"Loading {IN.name}...")
    p = pd.read_csv(IN, low_memory=False)
    p = p[p["sequence"].notna() & (p["sequence"].astype(str).str.len() > 0)]
    p = p.reset_index(drop=True)
    accessions = p["accession"].astype(str).tolist()
    seqs = p["sequence"].astype(str).tolist()
    n_total = len(seqs)
    print(f"  {n_total:,} atlas sequences (length: median={int(p['length_aa'].median())}, "
          f"max={int(p['length_aa'].max())})")

    pairs_df = compute_pairwise(seqs)
    cdhit_res = run_cdhit(accessions, seqs)
    write_summary(pairs_df, cdhit_res, n_total)
    make_figure(pairs_df, cdhit_res, n_total)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
