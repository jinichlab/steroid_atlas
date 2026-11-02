"""Publication-ready CYP-landscape UMAP figure.

Fixes from the previous iteration:
  - Labels on the plot are tiny numeric-only badges (`C6`, `C18` …) so the
    underlying colored dots stay fully visible.
  - Legend is horizontal and sits BELOW the plot (along the x-axis).
  - A two-column key table underneath maps `C# → family name` so the
    reader knows what each cluster is without labels blocking the dots.
  - Distinct hand-picked palette (no two greens, no two blues).
  - Larger colored dots so P450 clusters pop against the grey background.

Output:
    analysis/cyp_landscape_publication.{pdf,png}
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUT_DIR = HERE

MAJORITY_THRESHOLD = 0.50

plt.rcParams.update({
    "font.family": ["DejaVu Sans", "sans-serif"],
    "font.size": 10,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "#ffffff",
    "axes.facecolor":   "#ffffff",
})

# Hand-picked categorical palette — max contrast, no accidental twins
SUB_PALETTE = {
    "CYP46A": "#0369a1",   # deep blue
    "CYP39A": "#7c3aed",   # violet
    "CYP1A":  "#059669",   # emerald
    "CYP51A": "#b91c1c",   # crimson
    "CYP17A": "#ea580c",   # orange
    "CYP3A":  "#facc15",   # gold
    "CYP90B": "#db2777",   # pink
    # fallbacks
    "_other1": "#0891b2",  # cyan
    "_other2": "#6b7280",  # slate
    "_other3": "#a16207",  # amber-dark
}

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

CYP_RE = re.compile(r"CYP(\d+)([A-Z])?(\d*)")
CYP_NAME_RE = re.compile(r"cytochrome\s+p[-\s]?450\s+(\d+)([A-Z])?")


def cyp_hit(row) -> tuple[bool, str]:
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

# Rescue unlabeled entries where a family has one dominant subfamily
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

# ── Majority-P450 clusters ─────────────────────────────────────────────
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
        display_id=int(cid) + 1,
        family_name=name_by_cluster.get(int(cid), ""),
        n_proteins=n,
        n_cyps=n_cyp,
        pct_cyp=round(100 * n_cyp / n, 1),
        dominant_subfamily=dom_sub,
        cx=float(sub["umap_1"].mean()),
        cy=float(sub["umap_2"].mean()),
    ))
maj_df = pd.DataFrame(maj).sort_values("display_id")
print(f"Highlighted clusters: {len(maj_df)}")

# ── Color per subfamily ─────────────────────────────────────────────────
subs_present = sorted(maj_df["dominant_subfamily"].unique(),
                      key=lambda s: -maj_df["dominant_subfamily"].value_counts()[s])
color_by_sub = {}
fallback_idx = 1
for sub in subs_present:
    if sub in SUB_PALETTE:
        color_by_sub[sub] = SUB_PALETTE[sub]
    else:
        color_by_sub[sub] = SUB_PALETTE[f"_other{fallback_idx}"]
        fallback_idx += 1
sub_counts = maj_df["dominant_subfamily"].value_counts().to_dict()

# ── Figure with GridSpec: plot on top, legend + key table below ─────────
fig = plt.figure(figsize=(12, 11))
gs = GridSpec(3, 1, figure=fig, height_ratios=[10, 0.6, 3.4],
              hspace=0.15)

ax = fig.add_subplot(gs[0])
ax_legend = fig.add_subplot(gs[1])
ax_key = fig.add_subplot(gs[2])

# ── Main scatter ────────────────────────────────────────────────────────
ax.scatter(prot["umap_1"], prot["umap_2"], s=2, c="#e5e7eb",
           alpha=0.55, linewidth=0, zorder=1)
for _, row in maj_df.iterrows():
    cid = int(row["cluster"])
    sub = prot[prot["cluster"] == cid]
    color = color_by_sub[row["dominant_subfamily"]]
    ax.scatter(sub["umap_1"], sub["umap_2"], s=22, c=[color],
               alpha=0.9, edgecolor="black", linewidth=0.3, zorder=2)

# Small numeric-only badges — offset slightly so dots stay visible
for _, row in maj_df.iterrows():
    ax.annotate(
        f"C{int(row['display_id'])}",
        xy=(row["cx"], row["cy"]),
        xytext=(9, 9), textcoords="offset points",
        fontsize=9, fontweight="bold",
        color="#0f172a",
        ha="center", va="center",
        bbox=dict(
            boxstyle="circle,pad=0.28",
            facecolor="white",
            edgecolor="#0f172a",
            linewidth=0.9,
            alpha=0.97,
        ),
        arrowprops=dict(
            arrowstyle="-",
            color="#0f172a",
            lw=0.6,
            shrinkA=1,
            shrinkB=6,
        ),
        zorder=5,
    )

ax.set_xlabel("UMAP-1")
ax.set_ylabel("UMAP-2")
ax.set_title(
    f"Cytochrome P450 clusters across the Steroid Atlas UMAP\n"
    f"{len(maj_df)} majority-P450 clusters of the "
    f"{prot['cluster'].nunique()} · colored by dominant CYP subfamily",
    pad=10,
)
ax.grid(alpha=0.22, color="#94a3b8", linewidth=0.4)

# ── Legend row (horizontal, below plot) ─────────────────────────────────
ax_legend.axis("off")
handles = []
for sub in sorted(subs_present, key=lambda s: (-sub_counts[s], s)):
    n_cl = sub_counts[sub]
    handles.append(Patch(
        facecolor=color_by_sub[sub], edgecolor="black",
        label=f"{sub}  ({n_cl})",
    ))
handles.append(Patch(facecolor="#e5e7eb", edgecolor="none",
                     label="non-P450 majority"))
ax_legend.legend(
    handles=handles, loc="center",
    title="Dominant CYP subfamily  (# clusters in atlas)",
    title_fontsize=10, fontsize=9,
    ncol=len(handles), frameon=False,
    columnspacing=1.0, handletextpad=0.5,
)

# ── Key table below legend — 2-column mapping of C# → family name ───────
ax_key.axis("off")
maj_sorted = maj_df.sort_values("display_id")
rows = []
for _, r in maj_sorted.iterrows():
    rows.append([
        f"C{int(r['display_id'])}",
        r["family_name"],
        r["dominant_subfamily"],
        f"n={int(r['n_proteins'])}",
    ])
# Split into 2 columns for compactness
half = (len(rows) + 1) // 2
col_A = rows[:half]
col_B = rows[half:]
n_rows = max(len(col_A), len(col_B))
# Pad
while len(col_A) < n_rows: col_A.append(["", "", "", ""])
while len(col_B) < n_rows: col_B.append(["", "", "", ""])

table_data = []
for a, b in zip(col_A, col_B):
    table_data.append(a + [""] + b)

table = ax_key.table(
    cellText=table_data,
    colLabels=(
        ["Cluster", "Family name (dominant stem)", "Subfamily", "Size", "",
         "Cluster", "Family name (dominant stem)", "Subfamily", "Size"]
    ),
    cellLoc="left",
    colLoc="left",
    loc="center",
    edges="horizontal",
)
table.auto_set_font_size(False)
table.set_fontsize(8.5)
table.scale(1, 1.25)

# Style header row
for i in range(9):
    c = table[(0, i)]
    c.set_text_props(fontweight="bold", color="#334155")
    c.set_facecolor("#f1f5f9")
    if i == 4:  # spacer column
        c.set_facecolor("#ffffff")

# Colour the Cluster cell background to match its subfamily
for row_i, r_data in enumerate(table_data, start=1):
    for col_offset in (0, 5):
        cid_str = r_data[col_offset]
        if not cid_str:
            continue
        # find the subfamily for this cid
        cid_num = int(cid_str[1:])
        row = maj_sorted[maj_sorted["display_id"] == cid_num]
        if row.empty: continue
        sub = row["dominant_subfamily"].iloc[0]
        col = color_by_sub[sub]
        table[(row_i, col_offset)].set_facecolor(col)
        table[(row_i, col_offset)].set_text_props(color="white",
                                                  fontweight="bold")

# Column widths — narrow for cluster/subfamily/size, wide for name
col_widths = [0.05, 0.32, 0.07, 0.05, 0.02, 0.05, 0.32, 0.07, 0.05]
for i, w in enumerate(col_widths):
    for row_i in range(len(table_data) + 1):
        table[(row_i, i)].set_width(w)

fig.tight_layout()
fig.savefig(OUT_DIR / "cyp_landscape_publication.pdf")
fig.savefig(OUT_DIR / "cyp_landscape_publication.png",
            dpi=300, bbox_inches="tight")
print(f"Wrote analysis/cyp_landscape_publication.{{pdf,png}}")
