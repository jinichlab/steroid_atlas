"""Panel A for the CYP case study — color-by-subfamily version.

Each atlas cluster that is majority-P450 gets a cluster-ID badge on the
plot. Clusters are coloured by their dominant CYP SUBFAMILY so that
multiple clusters sharing a subfamily (e.g. the four CYP46A clusters)
show up in the same colour, visually revealing within-family fragmentation.
Legend is placed in the upper-left of the UMAP and uses the human-readable
function name as the color key.

Output:
    analysis/figure2_panelA_cyp_landscape.{pdf,png}

Convention: cluster IDs are 1-based (display_id = cluster + 1) to match the
web app's labels.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUT_DIR = ROOT / "analysis"

MAJORITY_THRESHOLD = 0.50
USE_ONE_BASED = True

plt.rcParams.update({
    "font.family": ["DejaVu Sans", "Arial", "sans-serif"],
    "font.size": 10,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "axes.facecolor":   "white",
})

# ── CYP subfamily → function name + color (hand-picked, visually distinct) ──
SUBFAMILY_LABEL = {
    "CYP46A": "Cholesterol 24-hydroxylase",
    "CYP51A": "Sterol 14α-demethylase",
    "CYP17A": "Steroid 17α-hydroxylase",
    "CYP1A":  "Estrogen / xenobiotic hydroxylase",
    "CYP3A":  "Xenobiotic hydroxylase",
    "CYP39A": "Oxysterol 7α-hydroxylase",
    "CYP90B": "Brassinosteroid oxidase",
    "CYP85A": "Brassinosteroid oxidase",
    "CYP11B": "11β-steroid hydroxylase",
    "CYP11A": "Cholesterol side-chain cleavage",
    "CYP27A": "Sterol 27-hydroxylase",
    "CYP7A":  "Cholesterol 7α-hydroxylase",
    "CYP7B":  "Oxysterol 7α-hydroxylase",
    "CYP19A": "Aromatase",
    "CYP21A": "Steroid 21-hydroxylase",
}

SUBFAMILY_COLOR = {
    # Legend order (by cluster count then alphabetical) → requested palette
    "CYP46A": "#2563eb",   # 1. blue       — 4 clusters
    "CYP39A": "#f97316",   # 2. orange     — 3 clusters
    "CYP17A": "#9333ea",   # 3. purple     — 1 cluster
    "CYP1A":  "#ca8a04",   # 4. mustard    — 1 cluster  (second "orange" interpretation)
    "CYP3A":  "#16a34a",   # 5. green      — 1 cluster
    "CYP51A": "#dc2626",   # 6. crimson    — 1 cluster
    "CYP90B": "#0d9488",   # 7. teal       — 1 cluster  (distinct from crimson)
    # Fallbacks for subfamilies that could enter the top-12 at different k
    "CYP85A": "#be185d",
    "CYP11B": "#14b8a6",
    "CYP11A": "#0d9488",
    "CYP27A": "#4338ca",
}


def family_description(sub: str) -> str:
    if sub in SUBFAMILY_LABEL:
        return SUBFAMILY_LABEL[sub]
    fam = re.match(r"(CYP\d+)", sub or "")
    if fam:
        key = fam.group(1)
        if key in SUBFAMILY_LABEL:
            return SUBFAMILY_LABEL[key]
    return "Cytochrome P450"


# ── Load atlas ──────────────────────────────────────────────────────────
prot = pd.read_csv(ROOT / "data" / "proteins.csv", low_memory=False)
prot["gene_names"] = prot["gene_names"].fillna("").astype(str)
prot["protein_names"] = prot["protein_names"].fillna("").astype(str)
prot["cluster"] = prot["cluster"].astype(int)

CYP_PATTERNS = [
    re.compile(r"CYP(\d+)([A-Z])?\d*"),
    re.compile(r"CYTOCHROME\s+P-?450\s+(\d+)([A-Z])?"),
]


def cyp_info(gene: str, pname: str) -> tuple[bool, str]:
    for src in (gene.upper(), pname.upper()):
        for pat in CYP_PATTERNS:
            m = pat.search(src)
            if m:
                sub = f"CYP{m.group(1)}{m.group(2) or ''}"
                return True, sub
    if "cytochrome p450" in pname.lower():
        return True, ""
    return False, ""


hits = prot.apply(lambda r: cyp_info(r["gene_names"], r["protein_names"]), axis=1)
prot["is_cyp"] = [h[0] for h in hits]
prot["cyp_sub"] = [h[1] for h in hits]

labeled = prot[prot["cyp_sub"].str.match(r"^CYP\d+[A-Z]$", na=False)]
letter_dominant: dict[str, str] = {}
for fam, sub in labeled.groupby(labeled["cyp_sub"].str.extract(r"^(CYP\d+)")[0]):
    counts = sub["cyp_sub"].str.extract(r"^CYP\d+([A-Z])$")[0].value_counts()
    if not counts.empty and counts.iloc[0] / len(sub) >= 0.7:
        letter_dominant[fam] = counts.index[0]


def resolve_sub(s: str) -> str:
    if not s:
        return ""
    if re.fullmatch(r"CYP\d+", s):
        letter = letter_dominant.get(s, "")
        return f"{s}{letter}" if letter else s
    return s


prot["cyp_sub"] = prot["cyp_sub"].apply(resolve_sub)

# ── Identify majority-P450 clusters ─────────────────────────────────────
maj: list[dict] = []
for cid, sub in prot.groupby("cluster"):
    n = len(sub)
    n_cyp = int(sub["is_cyp"].sum())
    if n_cyp / n < MAJORITY_THRESHOLD:
        continue
    sf_counts = sub[sub["cyp_sub"] != ""]["cyp_sub"].value_counts()
    dom_sub = sf_counts.index[0] if not sf_counts.empty else ""
    maj.append(dict(
        cluster=int(cid),
        display_id=int(cid) + 1 if USE_ONE_BASED else int(cid),
        n_proteins=n,
        n_cyps=n_cyp,
        pct_cyp=round(100 * n_cyp / n, 1),
        dominant_subfamily=dom_sub,
        cx=float(sub["umap_1"].mean()),
        cy=float(sub["umap_2"].mean()),
    ))

maj_df = (pd.DataFrame(maj)
          .sort_values(["dominant_subfamily", "n_proteins"],
                       ascending=[True, False])
          .reset_index(drop=True))

print(f"\nHighlighted clusters: {len(maj_df)}")
print(maj_df[["display_id", "cluster", "n_proteins", "pct_cyp",
              "dominant_subfamily"]].to_string(index=False))

# ── Figure ─────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(11, 9))

# Grey background
other = prot[~prot["cluster"].isin(maj_df["cluster"].to_list())]
ax.scatter(other["umap_1"], other["umap_2"], s=1.6, c="#e2e8f0",
           alpha=0.55, linewidth=0, zorder=1)

# Foreground: each highlighted cluster, coloured by subfamily
for _, row in maj_df.iterrows():
    cid = int(row["cluster"])
    sub_dom = row["dominant_subfamily"]
    color = SUBFAMILY_COLOR.get(sub_dom, "#64748b")
    sub = prot[prot["cluster"] == cid]
    ax.scatter(sub["umap_1"], sub["umap_2"], s=16, c=[color],
               alpha=0.88, edgecolor="black", linewidth=0.25, zorder=2)

# On-plot labels: just cluster ID
for _, row in maj_df.iterrows():
    ax.annotate(
        f"C{int(row['display_id'])}",
        xy=(row["cx"], row["cy"]),
        xytext=(0, 11),
        textcoords="offset points",
        fontsize=8.5,
        fontweight="700",
        ha="center", va="center",
        color="#0f172a",
        zorder=4,
        bbox=dict(
            boxstyle="round,pad=0.22",
            facecolor="white",
            edgecolor="#334155",
            linewidth=0.5,
            alpha=0.95,
        ),
        arrowprops=dict(arrowstyle="-", color="#334155", lw=0.4,
                        shrinkA=0.5, shrinkB=2),
    )

# Legend: one row per subfamily, showing function name + cluster-count
sub_cluster_counts = maj_df["dominant_subfamily"].value_counts().to_dict()
# Preserve order by cluster count (largest group first), then alphabetical
sub_order = sorted(sub_cluster_counts.keys(),
                   key=lambda s: (-sub_cluster_counts[s], s))

handles = []
for sub_dom in sub_order:
    color = SUBFAMILY_COLOR.get(sub_dom, "#64748b")
    label_name = family_description(sub_dom)
    ncl = sub_cluster_counts[sub_dom]
    handles.append(Patch(
        facecolor=color, edgecolor="black", linewidth=0.4,
        label=f"{sub_dom}",
    ))
handles.append(Patch(facecolor="#e2e8f0", edgecolor="none",
                     label="Non-P450 proteins"))

ax.legend(handles=handles,
          title="Dominant CYP subfamily",
          title_fontsize=10,
          fontsize=9,
          loc="lower left",
          framealpha=0.95,
          edgecolor="#cbd5e1",
          borderaxespad=0.8)

ax.set_xlabel("UMAP-1")
ax.set_ylabel("UMAP-2")
# title removed by user request
ax.grid(alpha=0.22, color="#94a3b8", linewidth=0.4)

fig.tight_layout()
pdf_path = OUT_DIR / "figure2_panelA_cyp_landscape.pdf"
png_path = OUT_DIR / "figure2_panelA_cyp_landscape.png"
fig.savefig(pdf_path, bbox_inches="tight")
fig.savefig(png_path, dpi=300, bbox_inches="tight")

print(f"\nWrote: {pdf_path}")
print(f"Wrote: {png_path}")
