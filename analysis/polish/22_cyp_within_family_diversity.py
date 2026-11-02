"""Within-family diversity across the 12 majority-P450 clusters.

For each of the 12 CYP clusters we compute:
  - # proteins, # unique substrates
  - median protein length + 10-90 percentile range
  - top 3 organisms (taxonomic composition)
  - domain-of-life mix (mammal / fish / reptile / plant / bacteria / etc.)
  - UMAP centroid
  - mean pairwise substrate Tanimoto (chemistry)
  - substrate Jaccard vs its same-family siblings

Then we build a comparison figure GROUPED BY CYP FAMILY so the reader can
see: for the families that split into multiple clusters (CYP46A × 4,
CYP39A × 3), what makes each cluster distinct? Is it taxonomy, length,
or substrate scope?

Outputs (under analysis/):
    cyp_within_family_diversity.tsv    per-cluster comparison table (12 rows)
    cyp_within_family_jaccard.tsv      pairwise substrate Jaccard for
                                         same-family cluster pairs
    figure_cyp_within_family.{pdf,png}  multi-panel comparison figure
    cyp_within_family_report.md        interpretation notes
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
except ImportError as e:
    raise SystemExit("rdkit not installed. Run under the conda python.") from e

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUT_DIR = HERE

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "#ffffff",
})

# The 12 majority-P450 clusters (display_id → internal cluster id + label)
CYP_CLUSTERS = {
    6:  ("CYP46A", "Cholesterol 24-OH"),
    9:  ("CYP39A", "24-OH-Chol 7α-OH"),
    18: ("CYP46A", "Cholesterol 24-OH"),
    20: ("CYP1A",  "Cytochrome P450 1A"),
    26: ("CYP46A", "Cholesterol 24-OH"),
    29: ("CYP51A", "Lanosterol 14α-DM"),
    34: ("CYP39A", "24-OH-Chol 7α-OH"),
    46: ("CYP17A", "Steroid 17α-OH"),
    47: ("CYP90B", "Cytochrome P450 85A1"),
    48: ("CYP46A", "Cholesterol 24-OH"),
    57: ("CYP3A",  "Cytochrome P450 3A"),
    71: ("CYP39A", "(3R)-3-OH-acyl-CoA DH (borderline)"),
}
internal = {d - 1: (d, sub, name) for d, (sub, name) in CYP_CLUSTERS.items()}

# ── Load substrate fingerprints ───────────────────────────────────────
print("Loading substrate fingerprints ...")
mol = pd.read_csv(ROOT / "data" / "molecules.csv", low_memory=False)
mol = mol[["chebi_id", "smiles"]].copy()
mol["chebi_id"] = pd.to_numeric(mol["chebi_id"], errors="coerce").astype("Int64")
mol = mol[mol["smiles"].astype(str).str.strip().ne("")].dropna(subset=["chebi_id"])
fp_by_chebi = {}
for _, r in mol.iterrows():
    m = Chem.MolFromSmiles(str(r["smiles"]))
    if m is not None:
        fp_by_chebi[int(r["chebi_id"])] = AllChem.GetMorganFingerprintAsBitVect(
            m, 2, nBits=2048
        )

# ── Load atlas + restrict to the 12 CYP clusters ─────────────────────
prot = pd.read_csv(ROOT / "data" / "proteins.csv", low_memory=False)
prot["cluster"] = prot["cluster"].astype(int)
prot["organism"] = prot["organism"].fillna("").astype(str)
prot["interacting_chebi_ids"] = prot["interacting_chebi_ids"].fillna("").astype(str)
prot["length_aa"] = pd.to_numeric(prot["length_aa"], errors="coerce")

sub = prot[prot["cluster"].isin(internal.keys())].copy()
_NUM_RE = re.compile(r"\d+")
sub["chebi_set"] = sub["interacting_chebi_ids"].map(
    lambda s: {int(x) for x in _NUM_RE.findall(s)}
)


# Rough kingdom / domain-of-life guess from organism string
def kingdom(org: str) -> str:
    o = org.lower()
    if any(t in o for t in ("homo sapiens", "sapiens", "mus musculus", "rattus",
                             "bos taurus", "sus scrofa", "canis lupus",
                             "macaca", "callithrix", "cricetulus", "papio")):
        return "Mammal"
    if any(t in o for t in ("gallus", "coturnix", "meleagris", "chick")):
        return "Bird"
    if any(t in o for t in ("thamnophis", "python", "gecko", "anolis",
                             "chelonia", "trachemys")):
        return "Reptile"
    if any(t in o for t in ("xenopus", "silurana", "rhinella", "bufo")):
        return "Amphibian"
    if any(t in o for t in ("danio rerio", "oryzias", "poecilia",
                             "salmo", "oncorhynchus", "scleropages",
                             "esox", "cyprinus", "zebrafish", "trout",
                             "salmon")):
        return "Fish"
    if any(t in o for t in ("saccharomyces", "candida", "aspergillus",
                             "penicillium", "neurospora", "schizosac",
                             "kluyveromyces")):
        return "Fungi"
    if any(t in o for t in ("arabidopsis", "oryza", "zea mays",
                             "solanum", "nicotiana", "vitis", "medicago")):
        return "Plant"
    if any(t in o for t in ("escherichia", "staphyloco", "mycobact",
                             "bacillus", "streptoco", "clostrid",
                             "salmonella", "pseudomonas", "corynebacterium",
                             "listeria", "helicobact")):
        return "Bacteria"
    if any(t in o for t in ("drosophila", "aedes", "anopheles",
                             "caenorhabditis", "ciona")):
        return "Invertebrate"
    return "Other"


sub["kingdom"] = sub["organism"].map(kingdom)


def pairwise_tanimoto(chebis: list[int]) -> np.ndarray:
    n = len(chebis)
    if n < 2:
        return np.array([], dtype=np.float32)
    fps = [fp_by_chebi[c] for c in chebis]
    sims = []
    for i in range(n - 1):
        sims.extend(DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:]))
    return np.asarray(sims, dtype=np.float32)


# ── Per-cluster comparison table ─────────────────────────────────────
rows = []
cluster_chebis: dict[int, set[int]] = {}
for cid, s in sub.groupby("cluster"):
    disp, subfam, label = internal[cid]
    all_c: set[int] = set()
    for st in s["chebi_set"]:
        all_c.update(st)
    known = sorted(c for c in all_c if c in fp_by_chebi)
    cluster_chebis[cid] = set(known)
    sims = pairwise_tanimoto(known)
    top_orgs = s["organism"].value_counts().head(3)
    kingdom_frac = s["kingdom"].value_counts(normalize=True).round(3)
    top_king = kingdom_frac.index[0] if not kingdom_frac.empty else ""
    length = s["length_aa"].dropna()
    rows.append(dict(
        display_id=disp,
        subfamily=subfam,
        family_label=label,
        n_proteins=len(s),
        n_substrates=len(known),
        median_length=int(length.median()) if len(length) else np.nan,
        length_p10=int(np.percentile(length, 10)) if len(length) else np.nan,
        length_p90=int(np.percentile(length, 90)) if len(length) else np.nan,
        umap_cx=round(s["umap_1"].mean(), 2),
        umap_cy=round(s["umap_2"].mean(), 2),
        top_kingdom=top_king,
        top_kingdom_pct=round(100 * kingdom_frac.iloc[0], 1)
                        if not kingdom_frac.empty else 0,
        top_org1=f"{top_orgs.index[0][:40]} ({int(top_orgs.iloc[0])})"
                 if len(top_orgs) >= 1 else "",
        top_org2=f"{top_orgs.index[1][:40]} ({int(top_orgs.iloc[1])})"
                 if len(top_orgs) >= 2 else "",
        top_org3=f"{top_orgs.index[2][:40]} ({int(top_orgs.iloc[2])})"
                 if len(top_orgs) >= 3 else "",
        mean_tanimoto=round(float(sims.mean()), 3) if sims.size else np.nan,
        median_tanimoto=round(float(np.median(sims)), 3) if sims.size else np.nan,
    ))

comp = pd.DataFrame(rows).sort_values(["subfamily", "display_id"])
comp.to_csv(OUT_DIR / "cyp_within_family_diversity.tsv",
            sep="\t", index=False)
print(f"\nWrote cyp_within_family_diversity.tsv ({len(comp)} rows)")
print(comp[["display_id", "subfamily", "n_proteins", "n_substrates",
            "median_length", "top_kingdom", "top_kingdom_pct",
            "mean_tanimoto"]].to_string(index=False))

# ── Pairwise Jaccard for same-family cluster pairs ───────────────────
jac_rows = []
for subfam in comp["subfamily"].unique():
    cids = comp[comp["subfamily"] == subfam]["display_id"].to_list()
    if len(cids) < 2:
        continue
    ints = [c - 1 for c in cids]
    for i, ci in enumerate(ints):
        for cj in ints[i + 1:]:
            A = cluster_chebis[ci]
            B = cluster_chebis[cj]
            u = len(A | B)
            jac_rows.append(dict(
                subfamily=subfam,
                cluster_a=f"C{ci + 1}",
                cluster_b=f"C{cj + 1}",
                n_shared=len(A & B),
                jaccard=round(len(A & B) / u, 3) if u else 0.0,
            ))
jac = pd.DataFrame(jac_rows)
jac.to_csv(OUT_DIR / "cyp_within_family_jaccard.tsv",
           sep="\t", index=False)
print(f"\nSubstrate Jaccard for same-family pairs:")
print(jac.to_string(index=False))

# ── Figure: multi-panel comparison, grouped by CYP family ────────────
print("\nRendering multi-panel comparison figure ...")

# One row per CYP subfamily (7 rows), each showing its clusters side-by-side
subfam_order = ["CYP46A", "CYP39A", "CYP1A", "CYP17A", "CYP51A", "CYP3A", "CYP90B"]
subfam_present = [s for s in subfam_order if s in comp["subfamily"].values]

kingdom_colors = {
    "Mammal":    "#c026d3",
    "Fish":      "#0369a1",
    "Reptile":   "#059669",
    "Amphibian": "#0891b2",
    "Bird":      "#eab308",
    "Fungi":     "#a16207",
    "Plant":     "#16a34a",
    "Bacteria":  "#dc2626",
    "Invertebrate": "#7c3aed",
    "Other":     "#94a3b8",
}

fig = plt.figure(figsize=(15, 3 * len(subfam_present)))
gs = GridSpec(len(subfam_present), 3, figure=fig,
              width_ratios=[1.0, 1.2, 1.6], hspace=0.55, wspace=0.35)

for row_idx, subfam in enumerate(subfam_present):
    sub_rows = comp[comp["subfamily"] == subfam].sort_values("display_id")
    n_cl = len(sub_rows)

    # Panel 1: length distribution (KDE-style histogram)
    ax = fig.add_subplot(gs[row_idx, 0])
    for j, r in enumerate(sub_rows.itertuples()):
        cid_int = r.display_id - 1
        L = sub[sub["cluster"] == cid_int]["length_aa"].dropna()
        ax.hist(L, bins=25, alpha=0.5,
                label=f"C{r.display_id} (n={r.n_proteins})",
                edgecolor="white", linewidth=0.4)
    ax.set_xlabel("Protein length (aa)")
    ax.set_ylabel("count")
    ax.set_title(f"{subfam}  ·  length", fontsize=10, pad=4)
    ax.legend(fontsize=7.5, loc="upper right")
    ax.grid(axis="y", alpha=0.3)

    # Panel 2: kingdom composition (stacked bars, one per cluster)
    ax = fig.add_subplot(gs[row_idx, 1])
    kg = (sub[sub["cluster"].isin(sub_rows["display_id"] - 1)]
             .groupby(["cluster", "kingdom"]).size().unstack(fill_value=0))
    kg = kg.reindex(sub_rows["display_id"] - 1)
    pcts = kg.div(kg.sum(axis=1), axis=0) * 100
    ord_kingdoms = pcts.sum().sort_values(ascending=False).index.to_list()
    pcts = pcts[ord_kingdoms]
    xpos = np.arange(n_cl)
    bottoms = np.zeros(n_cl)
    for king in ord_kingdoms:
        vals = pcts[king].values
        ax.bar(xpos, vals, bottom=bottoms,
               color=kingdom_colors.get(king, "#94a3b8"),
               edgecolor="black", linewidth=0.4, label=king)
        bottoms += vals
    ax.set_xticks(xpos)
    ax.set_xticklabels([f"C{r.display_id}" for r in sub_rows.itertuples()])
    ax.set_ylabel("% of proteins")
    ax.set_title(f"{subfam}  ·  kingdom mix", fontsize=10, pad=4)
    ax.set_ylim(0, 100)
    if row_idx == 0:
        ax.legend(fontsize=6.5, loc="upper right", framealpha=0.95,
                  ncol=2, handlelength=1.2, handletextpad=0.4)
    ax.grid(axis="y", alpha=0.3)

    # Panel 3: substrate scope summary (bars for n_substrates + mean Tanimoto)
    ax = fig.add_subplot(gs[row_idx, 2])
    xpos = np.arange(n_cl)
    ax.bar(xpos - 0.2, sub_rows["n_substrates"], width=0.4,
           color="#0369a1", edgecolor="black", linewidth=0.4,
           label="# substrates")
    ax.set_xticks(xpos)
    ax.set_xticklabels([f"C{r.display_id}" for r in sub_rows.itertuples()])
    ax.set_ylabel("# unique substrates", color="#0369a1")
    ax.tick_params(axis="y", labelcolor="#0369a1")

    ax2 = ax.twinx()
    ax2.bar(xpos + 0.2, sub_rows["mean_tanimoto"], width=0.4,
            color="#c026d3", edgecolor="black", linewidth=0.4,
            label="mean Tanimoto")
    ax2.set_ylabel("mean pairwise Tanimoto", color="#c026d3")
    ax2.tick_params(axis="y", labelcolor="#c026d3")
    ax2.set_ylim(0, 1)
    ax.set_title(f"{subfam}  ·  substrate chemistry", fontsize=10, pad=4)
    ax.grid(axis="y", alpha=0.3)

fig.suptitle(
    "Within-CYP-family diversity across the 12 majority-P450 clusters",
    fontsize=13, fontweight="bold", y=0.995,
)
fig.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig(OUT_DIR / "figure_cyp_within_family.pdf")
fig.savefig(OUT_DIR / "figure_cyp_within_family.png",
            dpi=300, bbox_inches="tight")
print(f"Wrote figure_cyp_within_family.{{pdf,png}}")

# ── Interpretation notes ─────────────────────────────────────────────
lines = [
    "# Within-CYP-family diversity — the 12 majority-P450 clusters",
    "",
    "## What we asked",
    "",
    "For the CYP families whose members split into multiple atlas clusters "
    "(**CYP46A × 4**, **CYP39A × 3**), what actually differentiates each "
    "cluster from its same-family siblings? Is it:",
    "- **Taxonomy** (different organisms → different sequence lineage)?",
    "- **Length / domain architecture** (different protein length distributions)?",
    "- **Substrate scope** (each cluster acts on a distinct chemical set)?",
    "",
    "## Per-cluster comparison table",
    "",
    "| C# | Subfamily | n_prot | n_sub | median len | top kingdom | mean Tanimoto |",
    "|---|---|---:|---:|---:|---|---:|",
]
for _, r in comp.iterrows():
    lines.append(
        f"| C{int(r['display_id'])} | {r['subfamily']} | "
        f"{int(r['n_proteins'])} | {int(r['n_substrates'])} | "
        f"{int(r['median_length'])} ({int(r['length_p10'])}–{int(r['length_p90'])}) | "
        f"{r['top_kingdom']} ({int(r['top_kingdom_pct'])}%) | "
        f"{r['mean_tanimoto']} |"
    )

lines += [
    "",
    "## Substrate Jaccard between same-family cluster pairs",
    "",
    "| Subfamily | Cluster A | Cluster B | # shared substrates | Jaccard |",
    "|---|---|---|---:|---:|",
]
for _, r in jac.iterrows():
    lines.append(
        f"| {r['subfamily']} | {r['cluster_a']} | {r['cluster_b']} | "
        f"{int(r['n_shared'])} | {r['jaccard']} |"
    )

lines += [
    "",
    "## Reading the figure",
    "",
    "The multi-panel figure lays out one row per CYP subfamily. Each row has:",
    "- **Length distribution** — histogram of protein length per cluster.",
    "- **Kingdom mix** — stacked bar of the taxonomic composition per cluster.",
    "- **Substrate chemistry** — number of unique substrates (blue) and "
    "mean pairwise Tanimoto (magenta) per cluster.",
    "",
    "Where the length + kingdom bars differ but the substrate bars are "
    "similar → the clusters are **taxonomic variants of the same enzyme**.",
    "Where the substrate bars differ → the clusters actually act on "
    "distinct chemical sets.",
    "",
    "## Files",
    "- `cyp_within_family_diversity.tsv` — the 12-row comparison table",
    "- `cyp_within_family_jaccard.tsv` — pairwise Jaccard for same-family pairs",
    "- `figure_cyp_within_family.{pdf,png}` — the multi-panel figure",
]
(OUT_DIR / "cyp_within_family_report.md").write_text("\n".join(lines),
                                                       encoding="utf-8")
print("\nWrote cyp_within_family_report.md")
