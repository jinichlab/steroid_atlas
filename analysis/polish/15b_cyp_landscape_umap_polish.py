"""Polished CYP-landscape UMAP figure — subfamily-level, app-style labels.

Improvements over analysis/15_cyp_superfamily_landscape.py:

- Dominant CYP designator is at SUBFAMILY level (CYP46A vs CYP46B), not
  just family. Biochemically distinct enzymes stop being lumped.
- Each highlighted cluster carries a pillow-style label with BOTH the
  cluster id AND the family name (`C6 · Cholesterol 24-hydroxylase`),
  matching the label style the web app draws when you zoom in.
- P450 dots are enlarged so they pop against the grey background.
- Legend lists subfamilies present in the highlighted clusters (with a
  count of clusters per subfamily in parentheses).
- Palette taken from the app's spatially-aware distinct-color generator
  so the figure reads like an app screenshot.

Output:
    analysis/cyp_landscape_umap_polished.{pdf,png}
"""
from __future__ import annotations

import colorsys
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUT_DIR = HERE

MAJORITY_THRESHOLD = 0.50

plt.rcParams.update({
    "font.family": ["-apple-system", "BlinkMacSystemFont", "Arial",
                     "sans-serif"],
    "font.size": 10,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "#f8fafc",
    "axes.facecolor":   "#ffffff",
})

# ── Load atlas + cluster fingerprints ───────────────────────────────────
prot = pd.read_csv(ROOT / "data" / "proteins.csv", low_memory=False)
prot["gene_names"] = prot["gene_names"].fillna("").astype(str)
prot["protein_names"] = prot["protein_names"].fillna("").astype(str)
prot["cluster"] = prot["cluster"].astype(int)

fp = pd.read_csv(HERE / "cluster_fingerprints.tsv", sep="\t")
name_by_cluster: dict[int, str] = {
    int(r["cluster"]): str(r["dominant_stem"])[:44]
    for _, r in fp.iterrows()
}

# ── CYP subfamily extraction (family digit + subfamily letter) ─────────
CYP_RE = re.compile(r"CYP(\d+)([A-Z])?(\d*)")
CYP_NAME_RE = re.compile(r"cytochrome\s+p[-\s]?450\s+(\d+)([A-Z])?")


def cyp_hit(row) -> tuple[bool, str]:
    """Return (is_p450, subfamily) — subfamily like `CYP46A` or `CYP3A`,
    falls back to family (`CYP46`) when subfamily letter is missing."""
    m = CYP_RE.search(row["gene_names"].upper())
    if m:
        fam = f"CYP{m.group(1)}"
        return True, f"{fam}{m.group(2)}" if m.group(2) else fam
    if "cytochrome p450" in row["protein_names"].lower() or \
            "cytochrome p-450" in row["protein_names"].lower():
        m2 = CYP_NAME_RE.search(row["protein_names"].lower())
        if m2:
            fam = f"CYP{m2.group(1)}"
            return True, f"{fam}{m2.group(2).upper()}" if m2.group(2) else fam
        return True, ""
    return False, ""


hits = prot.apply(cyp_hit, axis=1)
prot["is_cyp"] = [h[0] for h in hits]
prot["cyp_sub"] = [h[1] for h in hits]

# Rescue unlabeled entries: if a family has one dominant subfamily letter
# (>= 70% of labeled entries), assume unlabeled fall under it
letter_dominant: dict[str, str] = {}
labeled = prot[prot["cyp_sub"].str.contains(r"^CYP\d+[A-Z]$", regex=True,
                                            na=False)]
for _, sub in labeled.groupby(labeled["cyp_sub"].str.extract(r"^(CYP\d+)")[0]):
    counts = sub["cyp_sub"].str.extract(r"^CYP\d+([A-Z])$")[0].value_counts()
    if not counts.empty and counts.iloc[0] / len(sub) >= 0.7:
        letter_dominant[sub["cyp_sub"].iloc[0].rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")] = counts.index[0]


def resolve_sub(s: str) -> str:
    if not s:
        return ""
    if re.fullmatch(r"CYP\d+", s):
        letter = letter_dominant.get(s, "")
        return f"{s}{letter}" if letter else s
    return s


prot["cyp_sub"] = prot["cyp_sub"].apply(resolve_sub)

# ── Per-cluster: is this a majority-P450 cluster + which subfamily? ──
maj: list[dict] = []
for cid, sub in prot.groupby("cluster"):
    n = len(sub)
    n_cyp = int(sub["is_cyp"].sum())
    if n_cyp / n < MAJORITY_THRESHOLD:
        continue
    # Dominant subfamily = most frequent non-empty subfamily label
    sf_counts = (sub[sub["cyp_sub"] != ""]["cyp_sub"]
                 .value_counts())
    dom_sub = sf_counts.index[0] if not sf_counts.empty else ""
    maj.append(dict(
        cluster=int(cid),
        display_id=int(cid) + 1,
        family_name=name_by_cluster.get(int(cid), ""),
        n_proteins=n,
        n_cyps=n_cyp,
        pct_cyp=round(100 * n_cyp / n, 1),
        dominant_subfamily=dom_sub,
        cx=float(sub["umap_1"].mean()),
        cy=float(sub["umap_2"].mean()),
    ))
maj_df = pd.DataFrame(maj).sort_values("pct_cyp", ascending=False)
print(f"Highlighted clusters: {len(maj_df)}")
print(maj_df[["display_id", "family_name", "dominant_subfamily",
              "n_proteins", "pct_cyp"]].to_string(index=False))

# ── Palette — spatially-aware assignment so neighbors get contrasting hues
def hsl_rgb(h, s, l):
    return colorsys.hls_to_rgb(h, l, s)


def sub_color_for(sub: str, i: int) -> tuple[float, float, float]:
    """Golden-angle palette with alternating lightness/saturation."""
    h = ((i * 137.508) / 360) % 1
    l = 0.55 if i % 2 == 0 else 0.42
    s = [0.9, 0.72, 0.6][i % 3]
    return hsl_rgb(h, s, l)


# Assign a color per unique subfamily in the highlighted set
subs_present = sorted(maj_df["dominant_subfamily"].unique())
color_by_sub = {sub: sub_color_for(sub, i) for i, sub in enumerate(subs_present)}
sub_counts = maj_df["dominant_subfamily"].value_counts().to_dict()

# ── Figure ─────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(12, 9))

# Grey background = every atlas protein
ax.scatter(prot["umap_1"], prot["umap_2"], s=2, c="#e2e8f0",
           alpha=0.55, linewidth=0, zorder=1)

# Foreground: each highlighted cluster's P450s, colored by subfamily
for _, row in maj_df.iterrows():
    cid = int(row["cluster"])
    sub = prot[prot["cluster"] == cid]
    color = color_by_sub[row["dominant_subfamily"]]
    ax.scatter(sub["umap_1"], sub["umap_2"], s=14, c=[color],
               alpha=0.85, edgecolor="black", linewidth=0.25, zorder=2)

# Compact labels at each centroid — small font + abbreviated family name
# so the underlying colored dots stay visible.
FAMILY_ABBR = {
    "Cholesterol 24-hydroxylase":            "Chol 24-OH",
    "24-hydroxycholesterol 7-alpha-hydroxylase": "24-OH-Chol 7α-OH",
    "Lanosterol 14-alpha demethylase":       "Lano 14α-demethylase",
    "Steroid 17-alpha-hydroxylase":          "Steroid 17α-OH",
    "Cytochrome P450 85A1":                  "P450 85A1",
    "Cytochrome P450 1A":                    "P450 1A",
    "Cytochrome P450 3A":                    "P450 3A",
    "(3R)-3-hydroxyacyl-CoA dehydrogenase":  "(3R)-3-OH-acyl-CoA DH",
}


def short(name: str) -> str:
    return FAMILY_ABBR.get(name.strip(), name[:20])


for _, row in maj_df.iterrows():
    text = f"C{int(row['display_id'])} · {short(row['family_name'])}"
    ax.annotate(
        text,
        xy=(row["cx"], row["cy"]),
        xytext=(0, 12),   # small offset so the badge sits just above the dots
        textcoords="offset points",
        fontsize=6.5,     # much smaller than before
        fontweight="600",
        ha="center", va="center",
        color="#0f172a",
        zorder=4,
        bbox=dict(
            boxstyle="round,pad=0.22",
            facecolor="white",
            edgecolor="#475569",
            linewidth=0.4,
            alpha=0.92,
        ),
        arrowprops=dict(
            arrowstyle="-",
            color="#475569",
            lw=0.35,
            shrinkA=0.5,
            shrinkB=2,
        ),
    )

# Legend — one row per subfamily present, with a count of clusters
from matplotlib.patches import Patch
handles = []
for sub in sorted(subs_present, key=lambda s: (-sub_counts[s], s)):
    n_cl = sub_counts[sub]
    handles.append(Patch(
        facecolor=color_by_sub[sub], edgecolor="black",
        label=f"{sub}  ({n_cl} cluster{'s' if n_cl != 1 else ''})",
    ))
handles.append(Patch(facecolor="#e2e8f0", edgecolor="none",
                     label="Other proteins (non-P450 majority)"))
ax.legend(handles=handles, title="Dominant CYP subfamily",
          title_fontsize=10, fontsize=9, loc="upper left",
          framealpha=0.95, edgecolor="#cbd5e1")

ax.set_xlabel("UMAP-1")
ax.set_ylabel("UMAP-2")
ax.set_title(
    f"Cytochrome P450 clusters across the Steroid Atlas UMAP  ·  "
    f"{len(maj_df)} majority-P450 clusters of the {prot['cluster'].nunique()}",
    pad=12,
)
ax.grid(alpha=0.22, color="#94a3b8", linewidth=0.4)

fig.tight_layout()
fig.savefig(OUT_DIR / "cyp_landscape_umap_polished.pdf")
fig.savefig(OUT_DIR / "cyp_landscape_umap_polished.png",
            dpi=300, bbox_inches="tight")
print(f"\nWrote analysis/cyp_landscape_umap_polished.{{pdf,png}}")
