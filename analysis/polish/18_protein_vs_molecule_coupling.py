"""Cross-view coupling analysis.

For every (protein cluster, molecule cluster) pair, compute the fraction of
proteins in the protein cluster that interact with at least one molecule
belonging to that molecule cluster. The result is an 82x8 matrix we call
the "chemistry fingerprint" of each protein cluster.

Outputs:
    analysis/protein_molecule_coupling_matrix.tsv
    analysis/protein_molecule_coupling_heatmap.{pdf,png}
    stdout — specialist/generalist rankings and molecule-cluster response breadth
"""
from __future__ import annotations

from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUT_DIR = ROOT / "analysis"

plt.rcParams.update({
    "font.family": ["DejaVu Sans", "Arial", "sans-serif"],
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 7,
    "figure.facecolor": "white",
    "axes.facecolor":   "white",
})

# ── Load ──────────────────────────────────────────────────────────────
prot = pd.read_csv(ROOT / "data" / "proteins.csv", low_memory=False)
mols = pd.read_csv(ROOT / "data" / "molecules.csv", low_memory=False)
prot["cluster"] = prot["cluster"].astype(int)
mols = mols.dropna(subset=["chebi_id", "cluster"])
mols["cluster"] = mols["cluster"].astype(int)
mols["chebi_str"] = mols["chebi_id"].astype(int).astype(str)

# ChEBI → molecule cluster
chebi_to_mcluster = dict(zip(mols["chebi_str"], mols["cluster"]))
mcluster_sizes = mols.groupby("cluster").size().to_dict()
mcluster_ids = sorted(mols["cluster"].unique())
pcluster_ids = sorted(prot["cluster"].unique())
print(f"Protein clusters: {len(pcluster_ids)}  "
      f"Molecule clusters: {len(mcluster_ids)}")
print(f"Molecule cluster sizes: {mcluster_sizes}")

# Map protein → set of molecule clusters it touches
prot["chebi_list"] = prot["interacting_chebi_ids"].fillna("").astype(str)
def mclusters_of(s: str) -> set[int]:
    out = set()
    for c in s.split(";"):
        c = c.strip()
        if c and c != "nan" and c in chebi_to_mcluster:
            out.add(chebi_to_mcluster[c])
    return out

prot["mclusters"] = prot["chebi_list"].apply(mclusters_of)

# ── 82×8 coupling matrix: fraction of P's proteins that touch ≥1 ChEBI in M ──
coupling = pd.DataFrame(0.0, index=pcluster_ids, columns=mcluster_ids)
for pc, sub in prot.groupby("cluster"):
    n = len(sub)
    for mc in mcluster_ids:
        hits = sub["mclusters"].apply(lambda s: mc in s).sum()
        coupling.loc[pc, mc] = hits / n

coupling.to_csv(OUT_DIR / "protein_molecule_coupling_matrix.tsv", sep="\t")
print(f"\nCoupling matrix saved: {OUT_DIR / 'protein_molecule_coupling_matrix.tsv'}")

# ── Rankings ──────────────────────────────────────────────────────────
# A. Specialist / generalist score = number of molecule clusters a protein cluster
#    touches at ≥50% of its proteins
breadth = (coupling >= 0.5).sum(axis=1)
print("\n=== Protein-cluster chemical breadth (number of molecule clusters at ≥50% coverage) ===")
print("Distribution:")
print(breadth.value_counts().sort_index().to_string())
print("\nMost specialist (breadth = 0 or 1):")
specialists = breadth[breadth <= 1].index.tolist()
print(f"  {len(specialists)} clusters")
print("\nMost generalist (breadth ≥ 3):")
generalists = breadth[breadth >= 3].sort_values(ascending=False)
for pc in generalists.index[:15]:
    row = coupling.loc[pc]
    sz = (prot["cluster"] == pc).sum()
    touched = [f"M{int(m)}" for m in mcluster_ids if row[m] >= 0.5]
    print(f"  C{pc+1} (n={sz}): touches {', '.join(touched)}")

# B. Molecule-cluster response breadth = number of protein clusters that
#    have ≥50% of members touching this M
response = (coupling >= 0.5).sum(axis=0)
print("\n=== Molecule-cluster response breadth (number of protein clusters at ≥50% coverage) ===")
for mc in mcluster_ids:
    sz = mcluster_sizes[mc]
    n_p = int(response[mc])
    print(f"  M{mc} (n={sz} molecules): engaged by {n_p} protein clusters")

# ── Figure: heatmap ────────────────────────────────────────────────────
# Sort protein clusters by their total coupling score for a readable heatmap
total_coupling = coupling.sum(axis=1)
sort_idx = total_coupling.sort_values(ascending=False).index.tolist()
coupling_sorted = coupling.loc[sort_idx]

# Translate raw IDs to display IDs (1-based)
display_rows = [f"C{int(c)+1}" for c in coupling_sorted.index]
display_cols = [f"M{int(c)}" for c in coupling_sorted.columns]

fig, ax = plt.subplots(figsize=(8, 14))
# pcolormesh gives explicit edges between cells so the heatmap reads as a grid
mat = coupling_sorted.to_numpy()
nrows, ncols = mat.shape
im = ax.pcolormesh(
    np.arange(ncols + 1), np.arange(nrows + 1), mat,
    cmap="Blues", vmin=0, vmax=1,
    edgecolors="white", linewidth=0.7,
)
ax.set_xticks(np.arange(ncols) + 0.5, labels=display_cols, rotation=35,
              ha="right", rotation_mode="anchor")
ax.set_yticks(np.arange(nrows) + 0.5, labels=display_rows)
ax.set_xlim(0, ncols)
ax.set_ylim(nrows, 0)  # top-to-bottom
ax.set_xlabel("Molecule cluster  (LIPID MAPS–inspired labels; M-ID in parentheses)")
ax.set_ylabel("Protein cluster  (ranked by total molecule-cluster coverage)")
ax.set_title("Protein × molecule coupling\n(fraction of P's proteins touching ≥1 ChEBI in M)", pad=10)
cbar = plt.colorbar(im, ax=ax, shrink=0.4)
cbar.set_label("Fraction of protein-cluster members\ninteracting with molecule cluster", fontsize=8)
# Reduce y-tick density (show every other label)
for i, lbl in enumerate(ax.get_yticklabels()):
    if i % 2 == 1:
        lbl.set_visible(False)

fig.tight_layout()
pdf = OUT_DIR / "protein_molecule_coupling_heatmap.pdf"
png = OUT_DIR / "protein_molecule_coupling_heatmap.png"
fig.savefig(pdf, bbox_inches="tight")
fig.savefig(png, dpi=300, bbox_inches="tight")
print(f"\nWrote: {pdf}")
print(f"Wrote: {png}")
