"""Per-CYP-family substrate diversity ranking.

For each cytochrome P450 family (CYP1, CYP3, CYP7, CYP11, CYP17, CYP19,
CYP21, CYP39, CYP46, CYP51, CYP85, CYP90, ...), aggregate every unique
substrate its member proteins bind, compute pairwise ECFP4 Tanimoto over
the aggregated pool, and rank families by their substrate coherence.

Panel structure (horizontal violin, ranked by median):
  * SPECIALIST families at the top  — narrow substrate range, high Tanimoto
  * BROAD-SPECIFICITY families at the bottom — many substrates, low Tanimoto

Also writes a per-family table with n_proteins / n_substrates / mean /
median Tanimoto, plus a size-adjusted mean (bootstrap subsample to N=10).

Outputs (under analysis/):
    cyp_family_substrate_tanimoto.tsv
    figure_cyp_family_ranking.{pdf,png}
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

# Style
plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
})
NULL_COLOR = "#7f1d1d"
RANDOM_MEDIAN = 0.20

# ── Load substrate fingerprints ────────────────────────────────────────
print("Loading substrate fingerprints ...")
mol = pd.read_csv(ROOT / "data" / "molecules.csv", low_memory=False)
mol = mol[["compound_name", "chebi_id", "smiles"]].copy()
mol = mol[mol["smiles"].astype(str).str.strip().ne("")]
mol["chebi_id"] = pd.to_numeric(mol["chebi_id"], errors="coerce").astype("Int64")
mol = mol.dropna(subset=["chebi_id"])
fp_by_chebi: dict[int, object] = {}
for _, r in mol.iterrows():
    m = Chem.MolFromSmiles(str(r["smiles"]))
    if m is None:
        continue
    fp_by_chebi[int(r["chebi_id"])] = AllChem.GetMorganFingerprintAsBitVect(
        m, 2, nBits=2048
    )
print(f"  {len(fp_by_chebi):,} substrate fingerprints")

# ── Load proteins + classify CYP family ─────────────────────────────────
print("Loading proteins + classifying P450 family ...")
prot = pd.read_csv(ROOT / "data" / "proteins.csv", low_memory=False)
prot["gene_names"] = prot["gene_names"].fillna("").astype(str)
prot["protein_names"] = prot["protein_names"].fillna("").astype(str)
prot["interacting_chebi_ids"] = prot["interacting_chebi_ids"].fillna("").astype(str)
prot["cluster"] = prot["cluster"].astype(int)

_CYP_GENE_RE = re.compile(r"CYP(\d+)([A-Z])?(\d*)")
_CYP_NAME_RE = re.compile(r"cytochrome\s+p[-\s]?450\s+(\d+)")


def cyp_family(row) -> str:
    m = _CYP_GENE_RE.search(row["gene_names"].upper())
    if m:
        return f"CYP{m.group(1)}"
    m2 = _CYP_NAME_RE.search(row["protein_names"].lower())
    if m2:
        return f"CYP{m2.group(1)}"
    return ""


prot["cyp_family"] = prot.apply(cyp_family, axis=1)
p450 = prot[prot["cyp_family"] != ""].copy()
print(f"  {len(p450):,} P450 proteins across {p450['cyp_family'].nunique()} "
      f"resolvable CYP families")

_NUM_RE = re.compile(r"\d+")
p450["chebi_set"] = p450["interacting_chebi_ids"].map(
    lambda s: {int(x) for x in _NUM_RE.findall(s)}
)


def pairwise_tanimoto(chebis: list[int]) -> np.ndarray:
    """Upper-triangle pairwise Tanimoto over the given ChEBI list."""
    n = len(chebis)
    if n < 2:
        return np.array([], dtype=np.float32)
    fps = [fp_by_chebi[c] for c in chebis]
    sims = []
    for i in range(n - 1):
        sims.extend(DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:]))
    return np.asarray(sims, dtype=np.float32)


# ── Per-family aggregation ──────────────────────────────────────────────
print("Aggregating substrates per CYP family ...")
family_data = {}
for fam, sub in p450.groupby("cyp_family"):
    all_c = set().union(*sub["chebi_set"].to_list())
    known = sorted(c for c in all_c if c in fp_by_chebi)
    n_prot = len(sub)
    n_clusters = int(sub["cluster"].nunique())
    n_sub = len(known)
    sims = pairwise_tanimoto(known)
    family_data[fam] = dict(
        family=fam, n_proteins=n_prot, n_clusters=n_clusters,
        n_substrates=n_sub, sims=sims,
        mean=float(sims.mean()) if sims.size else np.nan,
        median=float(np.median(sims)) if sims.size else np.nan,
        std=float(sims.std()) if sims.size else np.nan,
    )

# Size-adjusted mean (bootstrap resample to N=10 substrates)
SUBSAMPLE_N = 10
N_BOOTSTRAP = 300
rng = np.random.default_rng(20260824)
for fam, d in family_data.items():
    if d["n_substrates"] < SUBSAMPLE_N:
        d["mean_adj"] = np.nan
        continue
    all_chebis = sorted({c for _, s in p450[p450["cyp_family"] == fam][
        "chebi_set"].items() for c in s if c in fp_by_chebi})
    means = []
    for _ in range(N_BOOTSTRAP):
        pick = rng.choice(len(all_chebis), size=SUBSAMPLE_N, replace=False)
        chebis = [all_chebis[i] for i in pick]
        means.append(pairwise_tanimoto(chebis).mean())
    d["mean_adj"] = float(np.mean(means))

# ── Save summary ────────────────────────────────────────────────────────
summary = pd.DataFrame([
    {k: v for k, v in d.items() if k != "sims"}
    for d in family_data.values()
])
summary = summary.sort_values("median", ascending=False, na_position="last")
for c in ("mean", "median", "std", "mean_adj"):
    summary[c] = summary[c].astype(float).round(4)
summary.to_csv(OUT_DIR / "cyp_family_substrate_tanimoto.tsv",
               sep="\t", index=False)
print(f"\n  wrote cyp_family_substrate_tanimoto.tsv "
      f"({len(summary)} families)")
print(summary.to_string(index=False))

# ── Figure: horizontal violin ranked by median ──────────────────────────
print("\nRendering per-CYP-family violin ranking ...")

# Filter to families with enough substrates for a meaningful distribution
MIN_SUB = 4
rows = [d for d in family_data.values() if d["n_substrates"] >= MIN_SUB]
rows.sort(key=lambda d: d["median"], reverse=True)

# Color families by their "specialty function" (rough, hand-curated groups)
CYP_ROLE_COLOR = {
    "CYP1": "#c026d3",   # magenta — xenobiotic
    "CYP3": "#c026d3",   # magenta — xenobiotic
    "CYP4": "#c026d3",   # magenta — xenobiotic
    "CYP7": "#0369a1",   # blue — bile acid biosynthesis
    "CYP8": "#0369a1",   # blue — bile acid
    "CYP27": "#0369a1",  # blue — bile acid
    "CYP39": "#0369a1",  # blue — oxysterol
    "CYP46": "#0369a1",  # blue — brain oxysterol
    "CYP11": "#dc2626",  # red — steroidogenic (cortisol/aldo)
    "CYP17": "#dc2626",  # red — steroidogenic (androgen)
    "CYP19": "#dc2626",  # red — steroidogenic (aromatase)
    "CYP21": "#dc2626",  # red — steroidogenic (cortisol precursor)
    "CYP51": "#16a34a",  # green — sterol biosynthesis
    "CYP85": "#eab308",  # yellow — plant BR
    "CYP90": "#eab308",  # yellow — plant BR
    "CYP2":  "#94a3b8",  # slate — mixed / other
}
LEGEND_ROLES = [
    ("#dc2626", "Steroidogenic (cortisol / androgen / estrogen)"),
    ("#0369a1", "Bile acid / oxysterol"),
    ("#16a34a", "Sterol biosynthesis"),
    ("#eab308", "Plant brassinosteroid"),
    ("#c026d3", "Xenobiotic metabolism"),
    ("#94a3b8", "Other / mixed"),
]

fig, ax = plt.subplots(figsize=(11, max(4.5, 0.55 * len(rows))))
pos = np.arange(len(rows))
data = [r["sims"] for r in rows]
parts = ax.violinplot(data, positions=pos, vert=False, widths=0.85,
                     showmedians=True, showextrema=False)
for i, body in enumerate(parts["bodies"]):
    body.set_facecolor(CYP_ROLE_COLOR.get(rows[i]["family"], "#94a3b8"))
    body.set_alpha(0.78)
    body.set_edgecolor("black")
    body.set_linewidth(0.5)
if "cmedians" in parts:
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.4)

# Labels + annotation
labels = [
    f"{r['family']}\n"
    f"n_prot={r['n_proteins']}, "
    f"n_sub={r['n_substrates']}, "
    f"clusters={r['n_clusters']}"
    for r in rows
]
ax.set_yticks(pos)
ax.set_yticklabels(labels, fontsize=8.5)
ax.invert_yaxis()

ax.axvline(RANDOM_MEDIAN, color=NULL_COLOR, linestyle="--", linewidth=1.2,
           label=f"Atlas random-pair median = {RANDOM_MEDIAN:.2f}")
ax.set_xlabel("Pairwise substrate Tanimoto (ECFP4)")
ax.set_xlim(0, 1)
ax.grid(axis="x", alpha=0.3)
ax.set_title(
    f"P450 substrate diversity per family  ·  "
    f"{len(rows)} CYP families with ≥{MIN_SUB} substrates, ranked by median",
    pad=10,
)

# Legend for the role colors
from matplotlib.patches import Patch
legend_handles = [Patch(facecolor=c, edgecolor="black", label=lbl)
                  for c, lbl in LEGEND_ROLES]
legend_handles.append(
    Patch(facecolor="none", edgecolor=NULL_COLOR,
          label=f"Random baseline = {RANDOM_MEDIAN:.2f}")
)
ax.legend(handles=legend_handles, loc="lower right", framealpha=0.95,
          fontsize=8, title="P450 role", title_fontsize=9)

fig.tight_layout()
fig.savefig(OUT_DIR / "figure_cyp_family_ranking.pdf")
fig.savefig(OUT_DIR / "figure_cyp_family_ranking.png",
            dpi=300, bbox_inches="tight")
print("  wrote figure_cyp_family_ranking.{pdf,png}")
print("\nDone.")
