"""P450 subfamily-level substrate diversity — notable trends only.

Same idea as script 18 but at SUBFAMILY granularity (CYP7A vs CYP7B,
CYP11A vs CYP11B, CYP19A vs CYP17A, etc.) so biochemically distinct
enzymes don't get merged.

Filter: only subfamilies with ≥5 proteins AND ≥5 substrates (drops single-
protein or single-substrate rows that are just noise). Ranks by median
substrate Tanimoto; colored by biological role.

Outputs (under analysis/):
    cyp_subfamily_substrate_tanimoto.tsv
    figure_cyp_subfamily_ranking.{pdf,png}
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
except ImportError as e:
    raise SystemExit("rdkit not installed. Run under the conda python.") from e

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUT_DIR = HERE

# Filters — chosen so every remaining row has genuine biological weight
MIN_PROTEINS = 5
MIN_SUBSTRATES = 5

# Biology-role coloring (curated from CYP nomenclature literature)
SUBFAM_ROLE = {
    # Steroidogenic
    "CYP11A": "steroidogenic",  "CYP11B": "steroidogenic",
    "CYP17A": "steroidogenic",  "CYP19A": "steroidogenic",
    "CYP21":  "steroidogenic",  "CYP21A": "steroidogenic",
    # Bile acid / oxysterol
    "CYP7A":  "bile_acid",      "CYP7B":  "bile_acid",
    "CYP8B":  "bile_acid",      "CYP27A": "bile_acid",
    "CYP39A": "bile_acid",      "CYP46A": "bile_acid",
    # Sterol biosynthesis
    "CYP51":  "sterol_bio",     "CYP51A": "sterol_bio",
    # Xenobiotic
    "CYP1A":  "xenobiotic",     "CYP1B":  "xenobiotic",
    "CYP2A":  "xenobiotic",     "CYP2B":  "xenobiotic",
    "CYP2C":  "xenobiotic",     "CYP2D":  "xenobiotic",
    "CYP2E":  "xenobiotic",     "CYP3A":  "xenobiotic",
    "CYP4A":  "xenobiotic",     "CYP4F":  "xenobiotic",
    # Plant brassinosteroid
    "CYP85A": "plant_br",       "CYP90A": "plant_br",
    "CYP90B": "plant_br",       "CYP90C": "plant_br",
    "CYP90D": "plant_br",
}
ROLE_COLOR = {
    "steroidogenic": "#dc2626",
    "bile_acid":     "#0369a1",
    "sterol_bio":    "#16a34a",
    "xenobiotic":    "#c026d3",
    "plant_br":      "#eab308",
    "other":         "#94a3b8",
}
ROLE_LABEL = {
    "steroidogenic": "Steroidogenic (cortisol · androgen · estrogen · aldo)",
    "bile_acid":     "Bile acid / oxysterol",
    "sterol_bio":    "Sterol biosynthesis",
    "xenobiotic":    "Xenobiotic metabolism",
    "plant_br":      "Plant brassinosteroid",
    "other":         "Other / mixed",
}

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 12.5,
    "axes.labelsize": 11,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "legend.fontsize": 8.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
})
NULL_COLOR = "#7f1d1d"
RANDOM_MEDIAN = 0.20

# ── Fingerprints ────────────────────────────────────────────────────────
print("Loading fingerprints ...")
mol = pd.read_csv(ROOT / "data" / "molecules.csv", low_memory=False)
mol = mol[["chebi_id", "smiles"]].copy()
mol = mol[mol["smiles"].astype(str).str.strip().ne("")]
mol["chebi_id"] = pd.to_numeric(mol["chebi_id"], errors="coerce").astype("Int64")
mol = mol.dropna(subset=["chebi_id"])
fp_by_chebi = {}
for _, r in mol.iterrows():
    m = Chem.MolFromSmiles(str(r["smiles"]))
    if m is not None:
        fp_by_chebi[int(r["chebi_id"])] = AllChem.GetMorganFingerprintAsBitVect(
            m, 2, nBits=2048
        )
print(f"  {len(fp_by_chebi):,} substrate fingerprints")

# ── Proteins + CYP subfamily assignment ─────────────────────────────────
print("Assigning CYP subfamilies ...")
prot = pd.read_csv(ROOT / "data" / "proteins.csv", low_memory=False)
prot["gene_names"] = prot["gene_names"].fillna("").astype(str)
prot["interacting_chebi_ids"] = prot["interacting_chebi_ids"].fillna("").astype(str)

CYP_RE = re.compile(r"CYP(\d+)([A-Z])?(\d*)")


def cyp_family_and_letter(row) -> tuple[str, str]:
    """Return (family, subfamily_letter). letter='' when the gene name
    lacks the subfamily letter (e.g. old-style `CYP51` instead of
    `CYP51A1`)."""
    m = CYP_RE.search(row["gene_names"].upper())
    if not m:
        return "", ""
    return f"CYP{m.group(1)}", (m.group(2) or "")


fam_letter = prot.apply(cyp_family_and_letter, axis=1)
prot["cyp_family_only"] = [t[0] for t in fam_letter]
prot["cyp_subfam_letter"] = [t[1] for t in fam_letter]

# For each CYP family, if one subfamily accounts for >70% of resolved
# entries (the ones WITH a letter), assume unlabeled proteins belong to
# that dominant subfamily — a defensible way to rescue older annotations
# like plain `CYP51` that pre-date subfamily assignment.
letter_dominant: dict[str, str] = {}
for fam, sub in prot[prot["cyp_family_only"] != ""].groupby("cyp_family_only"):
    labeled = sub[sub["cyp_subfam_letter"] != ""]
    if labeled.empty:
        continue
    counts = labeled["cyp_subfam_letter"].value_counts()
    top_letter, top_n = counts.index[0], counts.iloc[0]
    if top_n / len(labeled) >= 0.70:
        letter_dominant[fam] = top_letter


def resolve_subfamily(row) -> str:
    fam = row["cyp_family_only"]
    if not fam:
        return ""
    letter = row["cyp_subfam_letter"]
    if not letter and fam in letter_dominant:
        letter = letter_dominant[fam]
    return f"{fam}{letter}" if letter else fam


prot["cyp_subfamily"] = prot.apply(resolve_subfamily, axis=1)
p450 = prot[prot["cyp_subfamily"] != ""].copy()
print(f"  {len(p450):,} P450s with a resolvable subfamily "
      f"({p450['cyp_subfamily'].nunique()} distinct subfamilies)")
if letter_dominant:
    print(f"  merged unlabeled entries into dominant subfamily where "
          f">70%: {letter_dominant}")

_NUM_RE = re.compile(r"\d+")
p450["chebi_set"] = p450["interacting_chebi_ids"].map(
    lambda s: {int(x) for x in _NUM_RE.findall(s)}
)


def pairwise_tanimoto(chebis: list[int]) -> np.ndarray:
    n = len(chebis)
    if n < 2:
        return np.array([], dtype=np.float32)
    fps = [fp_by_chebi[c] for c in chebis]
    sims = []
    for i in range(n - 1):
        sims.extend(DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:]))
    return np.asarray(sims, dtype=np.float32)


# ── Per-subfamily aggregation ───────────────────────────────────────────
print("Computing per-subfamily substrate Tanimoto ...")
records = []
for subfam, sub in p450.groupby("cyp_subfamily"):
    all_c = set().union(*sub["chebi_set"].to_list())
    known = sorted(c for c in all_c if c in fp_by_chebi)
    n_prot = len(sub)
    n_sub = len(known)
    if n_prot < MIN_PROTEINS or n_sub < MIN_SUBSTRATES:
        continue
    sims = pairwise_tanimoto(known)
    records.append(dict(
        subfamily=subfam,
        family=f"CYP{re.match(r'CYP(\d+)', subfam).group(1)}",
        n_proteins=n_prot,
        n_clusters=int(sub["cluster"].nunique()),
        n_substrates=n_sub,
        mean=round(float(sims.mean()), 4),
        median=round(float(np.median(sims)), 4),
        std=round(float(sims.std()), 4),
        sims=sims,
        role=SUBFAM_ROLE.get(subfam, "other"),
    ))
records.sort(key=lambda r: r["median"], reverse=True)
print(f"  {len(records)} subfamilies pass filters "
      f"(n_proteins ≥ {MIN_PROTEINS} AND n_substrates ≥ {MIN_SUBSTRATES})")

summary = pd.DataFrame([{k: v for k, v in r.items() if k != "sims"}
                         for r in records])
summary.to_csv(OUT_DIR / "cyp_subfamily_substrate_tanimoto.tsv",
               sep="\t", index=False)
print("  wrote cyp_subfamily_substrate_tanimoto.tsv")
print()
print(summary.to_string(index=False))

# ── Figure ──────────────────────────────────────────────────────────────
print("\nRendering subfamily ranking violin ...")
fig, ax = plt.subplots(figsize=(11, max(5, 0.55 * len(records))))
pos = np.arange(len(records))
data = [r["sims"] for r in records]
parts = ax.violinplot(data, positions=pos, vert=False, widths=0.85,
                     showmedians=True, showextrema=False)
for i, body in enumerate(parts["bodies"]):
    body.set_facecolor(ROLE_COLOR[records[i]["role"]])
    body.set_alpha(0.78)
    body.set_edgecolor("black")
    body.set_linewidth(0.5)
if "cmedians" in parts:
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.8)

labels = [
    f"{r['subfamily']}   ·   "
    f"n={r['n_proteins']}, subs={r['n_substrates']}, clusters={r['n_clusters']}"
    for r in records
]
ax.set_yticks(pos)
ax.set_yticklabels(labels, fontsize=9)
ax.invert_yaxis()

ax.axvline(RANDOM_MEDIAN, color=NULL_COLOR, linestyle="--", linewidth=1.3)
ax.set_xlabel("Pairwise substrate Tanimoto (ECFP4)")
ax.set_xlim(0, 1)
ax.set_title(
    f"Notable P450 subfamily trends in the Steroid Atlas\n"
    f"({len(records)} subfamilies ranked by substrate coherence, "
    f"filtered to n_proteins≥{MIN_PROTEINS} & n_substrates≥{MIN_SUBSTRATES})",
    pad=10,
)
ax.grid(axis="x", alpha=0.3)

# Annotate extremes — placed OUTSIDE the plot at the top-right / bottom-right
# so they don't overlap violins
ax.text(1.01, 0.98, "◄ SPECIALIST", fontsize=10, color="#0369a1",
        fontweight="bold", va="top", ha="left",
        transform=ax.transAxes)
ax.text(1.01, 0.02, "◄ BROAD-\n   SPECIFICITY", fontsize=10, color="#c026d3",
        fontweight="bold", va="bottom", ha="left",
        transform=ax.transAxes)

# Legend — placed BELOW the plot so it doesn't cover any violins
from matplotlib.patches import Patch
roles_present = list({r["role"] for r in records})
role_order = ["steroidogenic", "bile_acid", "sterol_bio",
              "xenobiotic", "plant_br", "other"]
handles = [Patch(facecolor=ROLE_COLOR[r], edgecolor="black",
                 label=ROLE_LABEL[r])
           for r in role_order if r in roles_present]
handles.append(Patch(facecolor="none", edgecolor=NULL_COLOR,
                     label=f"Atlas random-pair Tanimoto = {RANDOM_MEDIAN:.2f}"))
ax.legend(handles=handles, loc="upper center",
          bbox_to_anchor=(0.5, -0.09), ncol=3, framealpha=0.95,
          title="P450 biological role", title_fontsize=10)

fig.tight_layout(rect=[0, 0.02, 1, 1])
fig.savefig(OUT_DIR / "figure_cyp_subfamily_ranking.pdf")
fig.savefig(OUT_DIR / "figure_cyp_subfamily_ranking.png",
            dpi=300, bbox_inches="tight")
print("  wrote figure_cyp_subfamily_ranking.{pdf,png}")
print("\nDone.")
