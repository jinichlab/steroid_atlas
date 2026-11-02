"""Panel B for the CYP case study — zoomed UMAP of the CYP46 region,
dots coloured by their specific CYP family.

The four CYP46 clusters (C6, C18, C26, C48) are badged on the plot.
Each dot is coloured by the CYP family it carries — so the predominantly
CYP46A clusters (C6, C18, C26) appear mostly blue, and the mixed C48
cluster appears visually fragmented into CYP46A + CYP11 + CYP27 + CYP1A
+ Halloween-gene subpopulations. No formal legend; family colour coding
is explained in the figure caption.

Output:
    analysis/figure2_panelB_cyp46_zoom.{pdf,png}
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUT_DIR = ROOT / "analysis"

plt.rcParams.update({
    "font.family": ["DejaVu Sans", "Arial", "sans-serif"],
    "font.size": 10,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "axes.facecolor":   "white",
})

CYP46_CLUSTERS = {
    5:  "C6",
    17: "C18",
    25: "C26",
    47: "C48",
}

# ── CYP family detection ──────────────────────────────────────────────
CYP_PATTERNS = [
    re.compile(r"CYP(\d+)([A-Z])?\d*"),
    re.compile(r"CYTOCHROME\s+P-?450\s+(\d+)([A-Z])?"),
]

def cyp_info(gene: str, pname: str) -> tuple[bool, str]:
    for src in (gene.upper(), pname.upper()):
        for pat in CYP_PATTERNS:
            m = pat.search(src)
            if m:
                return True, f"CYP{m.group(1)}{m.group(2) or ''}"
    if "cytochrome p450" in pname.lower():
        return True, ""
    return False, ""


# Specific family → colour. Any CYP46 family becomes blue (dominant).
# Palette matches Panel A for the seven "major" subfamilies;
# minor CYP families (CYP11, CYP27, Halloween genes) get complementary distinct colours.
FAMILY_COLOR = {
    "CYP46":  "#2563eb", "CYP46A": "#2563eb",   # Panel A blue
    "CYP39":  "#f97316", "CYP39A": "#f97316",   # Panel A orange
    "CYP17":  "#9333ea", "CYP17A": "#9333ea",   # Panel A purple
    "CYP1":   "#ca8a04", "CYP1A":  "#ca8a04", "CYP1B":  "#ca8a04", "CYP1D": "#ca8a04",  # Panel A mustard
    "CYP3":   "#16a34a", "CYP3A":  "#16a34a",   # Panel A green
    "CYP51":  "#dc2626", "CYP51A": "#dc2626",   # Panel A crimson
    "CYP90":  "#64748b", "CYP90B": "#64748b", "CYP85A": "#64748b",   # all treated as "other" (not focus of Panel B)
    # Minor families (new colours, distinct from Panel A palette)
    "CYP11":  "#64748b", "CYP11A": "#475569", "CYP11B": "#64748b",   # slate
    "CYP27":  "#be185d", "CYP27A": "#be185d",   # dark magenta
    # Insect Halloween genes
    "CYP302": "#8b5cf6", "CYP302A": "#8b5cf6",
    "CYP314": "#8b5cf6", "CYP314A": "#8b5cf6",
    "CYP315": "#8b5cf6", "CYP315A": "#8b5cf6",
    "CYP306": "#8b5cf6",
}


def color_for_family(fam: str) -> str:
    if not fam:
        return "#94a3b8"  # unlabeled P450 — muted slate
    if fam in FAMILY_COLOR:
        return FAMILY_COLOR[fam]
    # Fall back to family (drop subfamily letter)
    m = re.match(r"CYP\d+", fam)
    if m and m.group(0) in FAMILY_COLOR:
        return FAMILY_COLOR[m.group(0)]
    return "#64748b"  # other CYP (unknown family in palette)


# ── Load atlas ─────────────────────────────────────────────────────────
prot = pd.read_csv(ROOT / "data" / "proteins.csv", low_memory=False)
prot["gene_names"] = prot["gene_names"].fillna("").astype(str)
prot["protein_names"] = prot["protein_names"].fillna("").astype(str)
prot["cluster"] = prot["cluster"].astype(int)

hits = prot.apply(lambda r: cyp_info(r["gene_names"], r["protein_names"]), axis=1)
prot["is_cyp"] = [h[0] for h in hits]
prot["cyp_fam"] = [h[1] for h in hits]

# Zoom bounds based on the four CYP46 clusters
target = prot[prot["cluster"].isin(CYP46_CLUSTERS.keys())]
PAD = 0.3
x_lo = target["umap_1"].min() - PAD
x_hi = target["umap_1"].max() + PAD
y_lo = target["umap_2"].min() - PAD
y_hi = target["umap_2"].max() + PAD
print(f"Zoom window: UMAP-1 [{x_lo:.2f}, {x_hi:.2f}], UMAP-2 [{y_lo:.2f}, {y_hi:.2f}]")

in_window = prot[
    prot["umap_1"].between(x_lo, x_hi) & prot["umap_2"].between(y_lo, y_hi)
].copy()
in_window["point_color"] = in_window["cyp_fam"].apply(color_for_family)
# Non-P450 proteins get light grey
in_window.loc[~in_window["is_cyp"], "point_color"] = "#e2e8f0"
print(f"Atlas proteins in zoom window: {len(in_window)}  "
      f"(P450: {int(in_window['is_cyp'].sum())}, "
      f"non-P450: {int((~in_window['is_cyp']).sum())})")

# Breakdown of CYP families in the window
fam_counts = in_window.loc[in_window["is_cyp"], "cyp_fam"].value_counts()
print("CYP families in zoom window:")
print(fam_counts.head(20).to_string())

# ── Figure ─────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 8))

# Plot non-P450 first (bottom layer), then P450s on top so they're not hidden
non_p = in_window[~in_window["is_cyp"]]
ax.scatter(non_p["umap_1"], non_p["umap_2"], s=4,
           c="#e2e8f0", alpha=0.55, linewidth=0, zorder=1)

p = in_window[in_window["is_cyp"]]
ax.scatter(p["umap_1"], p["umap_2"], s=18,
           c=p["point_color"].to_list(),
           alpha=0.9, edgecolor="black", linewidth=0.25, zorder=2)

# Cluster-ID badges at each centroid of the 4 CYP46 clusters
for cid, label in CYP46_CLUSTERS.items():
    sub = prot[prot["cluster"] == cid]
    cx = float(sub["umap_1"].mean())
    cy = float(sub["umap_2"].mean())
    ax.annotate(
        label,
        xy=(cx, cy),
        xytext=(0, 14),
        textcoords="offset points",
        fontsize=10,
        fontweight="700",
        ha="center", va="center",
        color="#0f172a",
        zorder=4,
        bbox=dict(
            boxstyle="round,pad=0.25",
            facecolor="white",
            edgecolor="#334155",
            linewidth=0.55,
            alpha=0.97,
        ),
        arrowprops=dict(arrowstyle="-", color="#334155", lw=0.5,
                        shrinkA=0.5, shrinkB=3),
    )


# ── Legend ─────────────────────────────────────────────────────────────
from matplotlib.patches import Patch
legend_entries = [
    ("CYP46A",            "#2563eb"),
    ("CYP39A",            "#f97316"),
    ("CYP1A",             "#ca8a04"),
    ("CYP11 (A/B)",       "#475569"),
    ("CYP27A",            "#be185d"),
    ("CYP302A/314A/315A (ecdysone P450s)", "#8b5cf6"),
    ("Other / unlabeled P450", "#64748b"),
    ("Non-P450 proteins", "#e2e8f0"),
]
handles = [Patch(facecolor=c, edgecolor="black", linewidth=0.4, label=n)
           for n, c in legend_entries]
ax.legend(handles=handles,
          title="CYP family",
          title_fontsize=10,
          fontsize=8.5,
          loc="upper left",
          framealpha=0.95,
          edgecolor="#cbd5e1",
          borderaxespad=0.8)

ax.set_xlim(x_lo, x_hi)
ax.set_ylim(y_lo, y_hi)
ax.set_xlabel("UMAP-1")
ax.set_ylabel("UMAP-2")
ax.grid(alpha=0.22, color="#94a3b8", linewidth=0.4)

fig.tight_layout()
pdf_path = OUT_DIR / "figure2_panelB_cyp46_zoom.pdf"
png_path = OUT_DIR / "figure2_panelB_cyp46_zoom.png"
fig.savefig(pdf_path, bbox_inches="tight")
fig.savefig(png_path, dpi=300, bbox_inches="tight")
print(f"\nWrote: {pdf_path}")
print(f"Wrote: {png_path}")
