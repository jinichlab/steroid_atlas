"""Panels C and D for the CYP case study.

Panel C — Taxonomic composition of each CYP46 cluster (horizontal stacked bar).
Panel D — CYP family composition of each CYP46 cluster (horizontal stacked bar).

Produces one combined 2-panel figure, saved to:
    analysis/figure2_panelsCD_cyp46_composition.{pdf,png}
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUT_DIR = ROOT / "analysis"

plt.rcParams.update({
    "font.family": ["DejaVu Sans", "Arial", "sans-serif"],
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10.5,
    "xtick.labelsize": 9,
    "ytick.labelsize": 10,
    "legend.fontsize": 8.5,
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
# Order top-to-bottom (reverse so the pure clusters appear at top of the stack)
# Panels use different cluster orderings:
#   Panel C — original left-to-right (C6, C18, C26, C48)
#   Panel D — by decreasing CYP46A fraction (C26 100%, C6 99.6%, C18 69%, C48 60%)
CLUSTER_ORDER_C = [25, 5, 17, 47]  # same as Panel D
CLUSTER_ORDER_D = [25, 5, 17, 47]
CLUSTER_ORDER = CLUSTER_ORDER_C  # retained for backward compat where used downstream

# ── CYP detection ──────────────────────────────────────────────────────
CYP_PATTERNS = [re.compile(r"CYP(\d+)([A-Z])?\d*"),
                re.compile(r"CYTOCHROME\s+P-?450\s+(\d+)([A-Z])?")]

def cyp_info(gene: str, pname: str) -> tuple[bool, str]:
    for src in (gene.upper(), pname.upper()):
        for pat in CYP_PATTERNS:
            m = pat.search(src)
            if m: return True, f"CYP{m.group(1)}{m.group(2) or ''}"
    if "cytochrome p450" in pname.lower(): return True, ""
    return False, ""

# ── Kingdom classifier (expanded) ───────────────────────────────────────
def kingdom(org: str) -> str:
    o = (org or "").lower()
    if any(k in o for k in ['salmon','trout','danio','zebrafish','oryzias','killi',
            'poecilia','fundulus','xiphophor','tilapia','teleost','cichlid','guppy',
            'medaka','gasterost','pike','catfish','perch','cod','tuna','carp',
            'minnow','stickleb','platyfish','sturgeon','eel','rockfish','flounder',
            'halibut','sole','notothenioid','cyprinod']): return 'Fish'
    if any(k in o for k in ['alligator','crocod','turtle','lizard','snake','gecko',
            'anolis','python','tortoise','iguana','skink','monitor','viper','testudin',
            'squamat','sauria','pogona']): return 'Reptile'
    if any(k in o for k in ['chicken','duck','pigeon','sparrow','eagle','owl','hawk',
            'parrot','falcon','columba','gallus','corvus','passer','taenio','ficedul',
            'geospiza','melopsitta','junco','meleag','anser','anas','aquila','pelecan',
            'phoeni','struthio','pavo','cuculus','pico','trochil','psittac','phasian',
            'picidae','charadri','laridae','sulidae','ciconiidae','gruidae','accipitridae',
            'anatidae','opistho','threski','pterodroma','apteryx','emu','ostrich',
            'penguin','aptenodytes','pygoscelis','spheniscus','balearica','egretta',
            'corvinella','erpornis','ciconia','zapornia','corythaixoides','circaetus',
            'piaya','eudyptes','cardinalis','bucero']): return 'Bird'
    if any(k in o for k in ['xenopus','frog','salamander','newt','axolot','rana ',
            'hyla','bufo','pipa','ambystoma','andrias','anura',
            'leptobrachium']): return 'Amphibian'
    if any(k in o for k in ['mus musculus','mouse','rat','rattus','human','homo',
            'primate','monkey','macaca','pan ','gorilla','pig','sus scrofa','cattle',
            'bos ','sheep','ovis','horse','equus','dog','canis','cat ','felis','rabbit',
            'orycto','bear','ursus','whale','dolph','bat ','rodent','hamster',
            'mesocric','pero','squirrel','otter','lemur','mole','elephant','marsupial',
            'dasyp','armad','tupaia','erinaceus','sorex','cavia','guinea','mustelid',
            'vicugna','camelus','capra','mustela','proca','nomascus','rhinop','mandrill',
            'colobus','callith','saimiri','aotus','trichechus','odocoileus','cervus',
            'rangif','capreolus','tragulus','balaenop','orcinus','tursiops','lipotes',
            'pteropus','myotis','rhinolophus','rousettus','monodelphis','sarcophilus',
            'phascolarctos','trichosurus','vombatus','macropus','wallab','potorous',
            'phalanger','cebus','sapajus','papio','erythrocebus','chlorocebus',
            'trachypithecus','ceratotherium','sotalia','phocoena','mustel','ailuropoda',
            'manis','octodon','castor','fukomys','jaculus','spalax','dipodomys',
            'spermophilus','marmota','mesocricetus','suncus','mandrillus','rhinopithecus',
            'cnephaeus','molossus','balaenoptera','sousa','monodon','lynx','panthera',
            'carlito','pipistrellus','prionailurus']): return 'Mammal'
    if any(k in o for k in ['arabidopsis','oryza','sativa','medicago','glycine','zea',
            'populus','vitis','solanum','brachypodium']): return 'Plant'
    if any(k in o for k in ['bacillus','escherichia','streptomyces','mycobac','pseudo',
            'clostr','staph','strept','rhodo','archae','sulfolob','halo','methano',
            'eggerthella','bacteroides','bifidobacterium']): return 'Microbe'
    if any(k in o for k in ['drosophila','caenorhabditis','elegans','insect','arthropod',
            'bombyx','apis','anopheles','aedes','manduca','tribolium','acyrthosiphon',
            'nasonia']): return 'Invertebrate'
    return 'Other'


# ── Load atlas ─────────────────────────────────────────────────────────
prot = pd.read_csv(ROOT / "data" / "proteins.csv", low_memory=False)
prot["gene_names"] = prot["gene_names"].fillna("").astype(str)
prot["protein_names"] = prot["protein_names"].fillna("").astype(str)
prot["cluster"] = prot["cluster"].astype(int)

hits = prot.apply(lambda r: cyp_info(r["gene_names"], r["protein_names"]), axis=1)
prot["is_cyp"] = [h[0] for h in hits]
prot["cyp_fam"] = [h[1] for h in hits]
prot["kingdom"] = prot["organism"].fillna("").apply(kingdom)

# ── Taxonomic composition per cluster ──────────────────────────────────
TAX_ORDER = ['Mammal', 'Bird', 'Reptile', 'Amphibian', 'Fish',
             'Invertebrate', 'Plant', 'Microbe', 'Other']
# Taxonomy palette — Set2-inspired, distinct from Panel B family colours
TAX_COLOR = {
    'Mammal':       '#8c6d31',   # warm brown
    'Bird':         '#fc8d62',   # salmon
    'Reptile':      '#66c2a5',   # light teal
    'Amphibian':    '#e78ac3',   # pink
    'Fish':         '#8da0cb',   # periwinkle
    'Invertebrate': '#a6d854',   # lime
    'Plant':        '#17b36b',   # saturated green
    'Microbe':      '#b3b3b3',   # neutral grey
    'Other':        '#f1f5f9',   # near-white slate
}

tax_table = pd.DataFrame(0, index=CLUSTER_ORDER, columns=TAX_ORDER, dtype=int)
for cid in CLUSTER_ORDER:
    counts = prot[prot.cluster == cid]['kingdom'].value_counts()
    for k, v in counts.items():
        if k in TAX_ORDER:
            tax_table.loc[cid, k] = v
# Normalise to fractions
tax_frac = tax_table.div(tax_table.sum(axis=1), axis=0)

# ── Family composition per cluster ─────────────────────────────────────
FAM_BUCKETS = {
    'CYP46A':      lambda f: f.startswith('CYP46'),
    'CYP1A':       lambda f: f.startswith('CYP1A') or f == 'CYP1',
    'CYP11A/B':    lambda f: f.startswith('CYP11'),
    'CYP27A':      lambda f: f.startswith('CYP27'),
    'CYP17A':      lambda f: f.startswith('CYP17'),
    'CYP302/314/315A': lambda f: any(f.startswith(x) for x in ['CYP302', 'CYP314', 'CYP315']),
    'Other P450':  None,  # catch-all for remaining P450s
    'Non-P450':    None,  # non-P450 proteins
}
# Matches Panel B family palette exactly
FAM_COLOR = {
    'CYP46A':            '#2563eb',  # blue           — Panel B
    'CYP1A':             '#ca8a04',  # mustard        — Panel B
    'CYP11A/B':          '#475569',  # slate          — Panel B
    'CYP27A':            '#be185d',  # dark magenta   — Panel B
    'CYP17A':            '#9333ea',  # purple         (Panel A, used if present)
    'CYP302/314/315A':   '#8b5cf6',  # violet         — Panel B
    'Other P450':        '#64748b',  # slate-grey     — Panel B
    'Non-P450':          '#e2e8f0',  # light grey     — Panel B
}

def family_bucket(row) -> str:
    if not row['is_cyp']:
        return 'Non-P450'
    for name, test in FAM_BUCKETS.items():
        if test is not None and test(row['cyp_fam']):
            return name
    return 'Other P450'


prot['fam_bucket'] = prot.apply(family_bucket, axis=1)
FAM_ORDER = list(FAM_BUCKETS.keys())

fam_table = pd.DataFrame(0, index=CLUSTER_ORDER, columns=FAM_ORDER, dtype=int)
for cid in CLUSTER_ORDER:
    counts = prot[prot.cluster == cid]['fam_bucket'].value_counts()
    for k, v in counts.items():
        if k in FAM_ORDER:
            fam_table.loc[cid, k] = v
fam_frac = fam_table.div(fam_table.sum(axis=1), axis=0)

# ── Print summary ───────────────────────────────────────────────────────
print("Taxonomic composition (counts):")
print(tax_table.to_string())
print("\nFamily composition (counts):")
print(fam_table.to_string())

# ── Figure ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5))

# Panel C — Taxonomic composition
axC = axes[0]
y_C = np.arange(len(CLUSTER_ORDER_C))
labels_C = [CYP46_CLUSTERS[c] for c in CLUSTER_ORDER_C]
tax_frac_C = tax_frac.loc[CLUSTER_ORDER_C]
left = np.zeros(len(CLUSTER_ORDER_C))
for tax in TAX_ORDER:
    vals = tax_frac_C[tax].to_numpy()
    axC.barh(y_C, vals, left=left, color=TAX_COLOR[tax], edgecolor='white',
             linewidth=0.3, label=tax if vals.sum() > 0 else None)
    left += vals
axC.set_yticks(y_C, labels=labels_C)
axC.invert_yaxis()
axC.set_xlim(0, 1.0)
axC.set_xlabel("Fraction of cluster members")
axC.set_title("C  Taxonomic composition", loc='left', fontweight='bold')
axC.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
axC.set_xticklabels(['0', '25%', '50%', '75%', '100%'])
# Only show legend for taxa actually present
present = [t for t in TAX_ORDER if tax_frac_C[t].sum() > 0]
handles = [plt.Rectangle((0, 0), 1, 1, color=TAX_COLOR[t]) for t in present]
axC.legend(handles, present, loc='center left', bbox_to_anchor=(1.0, 0.5),
           frameon=False, title='Kingdom', title_fontsize=9)

# Panel D — Family composition
axD = axes[1]
y_D = np.arange(len(CLUSTER_ORDER_D))
labels_D = [CYP46_CLUSTERS[c] for c in CLUSTER_ORDER_D]
fam_frac_D = fam_frac.loc[CLUSTER_ORDER_D]
left = np.zeros(len(CLUSTER_ORDER_D))
for fam in FAM_ORDER:
    vals = fam_frac_D[fam].to_numpy()
    axD.barh(y_D, vals, left=left, color=FAM_COLOR[fam], edgecolor='white',
             linewidth=0.3, label=fam if vals.sum() > 0 else None)
    left += vals
axD.set_yticks(y_D, labels=labels_D)
axD.invert_yaxis()
axD.set_xlim(0, 1.0)
axD.set_xlabel("Fraction of cluster members")
axD.set_title("D  Protein-family composition", loc='left', fontweight='bold')
axD.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
axD.set_xticklabels(['0', '25%', '50%', '75%', '100%'])
present_f = [f for f in FAM_ORDER if fam_frac_D[f].sum() > 0]
handles_f = [plt.Rectangle((0, 0), 1, 1, color=FAM_COLOR[f]) for f in present_f]
axD.legend(handles_f, present_f, loc='center left', bbox_to_anchor=(1.0, 0.5),
           frameon=False, title='Family', title_fontsize=9)

fig.tight_layout()
pdf = OUT_DIR / "figure2_panelsCD_cyp46_composition.pdf"
png = OUT_DIR / "figure2_panelsCD_cyp46_composition.png"
fig.savefig(pdf, bbox_inches='tight')
fig.savefig(png, dpi=300, bbox_inches='tight')
print(f"\nWrote: {pdf}")
print(f"Wrote: {png}")
