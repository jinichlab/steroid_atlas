"""Panels C and D for the CYP case study.

Panel C — CYP family composition of each CYP46 cluster (horizontal stacked bar).
Panel D — Taxonomic composition of each CYP46 cluster (horizontal stacked bar).

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

# ── Kingdom classifier (genus-based, comprehensive) ────────────────────
# Genus-level lookup built from the full list of unique organism strings
# present in the atlas across the 12 majority-P450 clusters. More accurate
# than substring matching, which misses ~140 non-passerine bird genera and
# ~95 songbird (Passerida) genera in CYP46 clusters C26 and C48.

_BIRD_GENERA = set("""
Accipiter Acanthisitta Acrocephalus Aegithalos Aegotheles Agelaius Alaudala
Alca Alectura Aleadryas Amazona Anas Anhinga Anser Anseranas Anthoscopus
Antrostomus Aphelocoma Apteryx Aptenodytes Aquila Ardea Ardeotis Arenaria
Asarcornis Astur Atlantisia Atrichornis Aythya Balaeniceps Balearica Baryphthengus
Bombycilla Brachypodius Brachypteracias Branta Bucco Buceros Bubo Buphagus
Buteo Callaeas Caloenas Calonectris Calypte Calcarius Calidris Campylorhamphus
Cardinalis Casuarius Cathartes Catharus Centropus Cepphus Certhia Cettia Ceyx
Chaetops Chaetorhynchus Chaetura Charadrius Chionis Chlamydotis Chloroceryle
Chloropsis Chordeiles Chroicocephalus Chrysolophus Chunga Ciccaba Cinclus
Circaetus Ciconia Cisticola Climacteris Cnemophilus Cochlearius Columba Columbina
Contopus Corvus Corythaeola Corythaixoides Coturnix Crypturellus Crotophaga
Cuculus Cyanistes Cyanocitta Cyanocorax Cyanoderma Dasyornis Delichon
Dendrocolaptes Dendroica Dicaeum Dicrurus Donacobius Drymodes Dromaius Dromas
Dryobates Dryoscopus Dumetella Edolisoma Egretta Elachura Emberiza Empidonax
Eolophus Erithacus Erpornis Erythrocercus Eudromia Eudyptes Eudyptula Eulacestoma
Eurystomus Falco Falcunculus Ficedula Formicarius Fregata Fregetta Fringilla
Fulmarus Furnarius Galbula Gallus Garrulus Gavia Geococcyx Geospiza Geothlypis
Glaucidium Glareola Grallaria Grantiella Grus Gymnorhina Haliaeetus Halcyon
Helopsaltes Hemiprocne Herpetotheres Himantopus Hippolais Hirundo Horornis Hylia
Hylocichla Hypocryptadius Ibidorhyncha Icterus Ifrita Illadopsis Indicator Irena
Jacana Junco Lamprotornis Lanius Larus Leiothrix Leptosomus Lepidocolaptes
Lepidothrix Leptocoma Leucopsar Limosa Lonchura Lophotis Lophophanes Loxia
Luscinia Machaerirhynchus Malurus Manacus Melopsittacus Melospiza Meleagris
Melanocharis Menura Merops Mesembrinibis Mesitornis Mimus Mionectes Mniotilta
Mohoua Molothrus Motacilla Mycerobas Myiagra Myiarchus Mystacornis Nesospiza
Nestor Neodrepanis Neopipo Nicator Nipponia Notiomystis Nothoprocta Nothocercus
Nucifraga Nyctibius Nyctiprogne Oceanites Oceanodroma Odontophorus Oenanthe
Onychorhynchus Opisthocomus Oreocharis Oreotrochilus Oriolus Origma Orthonyx Otus
Oxylabes Oxyruncus Pachyramphus Pampusana Pandion Panurus Pardalotus Parus
Passer Passerina Patagioenas Pedionomus Pelecanus Pelecanoides Penelope Peucedramus
Phaethon Phaetusa Phainopepla Phalacrocorax Phaps Phasianus Pheucticus
Phoenicopterus Phoenicurus Pica Picathartes Picoides Pipra Pipilo Piprites Pitangus
Pitta Piranga Platysteira Ploceus Pluvianellus Podargus Podiceps Podilymbus Poecile
Polioptila Pomatorhinus Pomatostomus Probosciger Progne Promerops Prunella
Psilopogon Psittacus Psophia Pterocles Pteruthius Ptilonorhynchus Ptilorrhoa
Pycnonotus Pygoscelis Pyrocephalus Pyrrhocorax Pyrrhula Quiscalus Ramphastos
Regulus Rhabdornis Rhagologus Rhadina Rhea Rhipidura Rhodinocichla Rhynochetos
Riparia Rissa Rostratula Rynchops Sagittarius Sakesphorus Sapayoa Saxicola Scopus
Sclerurus Scytalopus Semnornis Serilophus Serinus Setophaga Sialia Sitta Smithornis
Smutsornis Spheniscus Spiza Spinus Spizaetus Spizella Steatornis Sterrhoptilus
Streptopelia Strigops Struthidea Struthio Sturnus Sinosuthora Sula Sylvia Sylvietta
Syrrhaptes Tachuris Taeniopygia Terpsiphone Thalassarche Thinocorus Thryothorus
Tichodroma Toxostoma Trogon Tricholaema Troglodytes Turdus Turnix Tyrannus
Tychaedon Tyto Upupa Uria Urocolius Urocynchramus Vermivora Vidua Vireo Zapornia
Zenaida Zonotrichia Zosterops
""".split())

_MAMMAL_GENERA = set("""
Homo Pan Pongo Gorilla Macaca Papio Cercopithecus Chlorocebus Mandrillus
Trachypithecus Rhinopithecus Colobus Callithrix Cebus Sapajus Saimiri Aotus
Prolemur Lemur Microcebus Tarsius Carlito Nomascus Hylobates Erythrocebus
Mus Rattus Peromyscus Microtus Cricetulus Mesocricetus Marmota Spermophilus
Urocitellus Ictidomys Dipodomys Nannospalax Acomys Jaculus Octodon Chinchilla
Cavia Hydrochoerus Castor Fukomys Heterocephalus Ochotona Lepus Oryctolagus
Ornithorhynchus Tachyglossus Zaglossus Monodelphis Sarcophilus Trichosurus
Phascolarctos Vombatus Macropus Potorous Phalanger Dasypus Chrysochloris Manis
Loxodonta Trichechus Procavia Elephantulus Bos Bubalus Bison Ovis Capra Capreolus
Odocoileus Cervus Rangifer Dama Alces Antilocapra Giraffa Camelus Vicugna Lama
Sus Phacochoerus Equus Ceratotherium Diceros Tapirus Rhinolophus Myotis
Pteropus Rousettus Pipistrellus Nyctalus Carollia Artibeus Molossus Noctilio
Eptesicus Cnephaeus Balaenoptera Physeter Mesoplodon Delphinapterus Neophocaena
Globicephala Phocoena Tursiops Orcinus Lipotes Sousa Sotalia Monodon Lagenorhynchus
Eschrichtius Megaptera Ailuropoda Enhydra Lutra Mustela Mellivora Gulo Ursus Panthera
Lynx Puma Felis Prionailurus Suricata Hyaena Crocuta Vulpes Canis Nyctereutes
Odobenus Phoca Neomonachus Leptonychotes Zalophus Eumetopias Mirounga
""".split())

_REPTILE_GENERA = set("""
Crocodylus Alligator Gavialis Caiman Chelonia Chrysemys Pelodiscus Trachemys
Terrapene Dermochelys Caretta Gopherus Anolis Pogona Sceloporus Zootoca Podarcis
Ophiophagus Naja Python Varanus Aspidoscelis Gekko Lacerta Elgaria Eublepharis
Paroedura Protobothrops Crotalus Thamnophis Pseudonaja Sphenodon Chelydra
Mauremys Pantherophis Platysternon Phrynocephalus Salvator
""".split())

_AMPHIBIAN_GENERA = set("""
Xenopus Rana Nanorana Microcaecilia Rhinella Bufo Hyla Lithobates Hymenochirus
Pelophylax Pleurodeles Ambystoma Caecilia Geotrypetes Gegeneophis Eleutherodactylus
Engystomops Leptobrachium Pipa Andrias Pyxicephalus Pelobates
""".split())

_FISH_GENERA = set("""
Danio Danionella Oryzias Takifugu Oreochromis Salmo Esox Gadus Scyliorhinus
Latimeria Callorhinchus Petromyzon Astyanax Poecilia Hippocampus Lepisosteus
Protopterus Monopterus Oncorhynchus Clupea Ictalurus Scleropages Nothobranchius
Gasterosteus Fundulus Xiphophorus Kryptolebias Dicentrarchus Larimichthys
Pundamilia Maylandia Astatotilapia Pygocentrus Hippoglossus Paralichthys
Cynoglossus Solea Amphiprion Stegastes Acanthochromis Neolamprologus Haplochromis
Simochromis Boleophthalmus Periophthalmus Anabas Mastacembelus Sparus Seriola
Tetraodon Fugu Mola Myripristis Anguilla Hucho Salvelinus Thymallus Coregonus
Osmerus Coilia Engraulis Sardina Dorosoma Chanos Megalops Elops Albula Hemibagrus
Pangasianodon Pelteobagrus Clarias Silurus Rhincodon Carcharodon Stegostoma
Chiloscyllium Erpetoichthys Polypterus Neogobius Sphaeramia Cirrhinus
Synaphobranchus Hemitripterus Sinocyclocheilus Electrophorus Amphilophus Salarias
Cyprinus Cyprinodon Cyclopterus Lates Mugilogobius Carassius Betta Gambusia
Acipenser Triplophysa Knipowitschia Umbra Scomber Alosa Crenichthys Scophthalmus
Paramormyrops Labrus Labeo Collichthys Perca Notothenia Cottoperca Parambassis
Gymnodraco Atractosteus Pleuronectes Stegastes Eleginops Phoxinus Xyrichtys
Zoarces Culter Pagothenia Dissostichus Onychostoma Gouania Austrofundulus
Channa Denticeps Sander Echeneis Pygocentrus Pocillopora Austrobatrachus
Champsocephalus
""".split())

_INVERTEBRATE_GENERA = set("""
Drosophila Caenorhabditis Nematostella Branchiostoma Ciona Apis Bombyx Tribolium
Daphnia Strongylocentrotus Octopus Crassostrea Hydra Aplysia Lottia Haliotis
Pinctada Mizuhopecten Euprymna Styela Stylophora Pocillopora Acropora Orbicella
Acanthaster Actinia Patiria Holothuria Stichopus Desmophyllum
""".split())

_PLANT_GENERA = set("""
Arabidopsis Oryza Zea Solanum Nicotiana Glycine Vitis Populus Cucumis Brassica
Hordeum Triticum Sorghum Medicago Phaseolus Pisum Brachypodium
""".split())

_MICROBE_GENERA = set("""
Escherichia Bacillus Streptomyces Mycobacterium Pseudomonas Clostridium
Staphylococcus Streptococcus Rhodobacter Sulfolobus Haloarcula Methanothermobacter
Eggerthella Bacteroides Bifidobacterium Saccharomyces Candida Aspergillus
Neurospora Cryptococcus Schizosaccharomyces Ustilago Pichia
""".split())

def kingdom(org: str) -> str:
    """Return the vertebrate / invertebrate / plant / microbe group of the organism,
    based on the first word (genus) of the organism string."""
    if not org: return 'Unknown'
    g = org.split()[0].strip()
    if g in _BIRD_GENERA: return 'Bird'
    if g in _MAMMAL_GENERA: return 'Mammal'
    if g in _REPTILE_GENERA: return 'Reptile'
    if g in _AMPHIBIAN_GENERA: return 'Amphibian'
    if g in _FISH_GENERA: return 'Fish'
    if g in _INVERTEBRATE_GENERA: return 'Invertebrate'
    if g in _PLANT_GENERA: return 'Plant'
    if g in _MICROBE_GENERA: return 'Microbe'
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

# Panel C — Protein-family composition (left)
axC = axes[0]
y_C = np.arange(len(CLUSTER_ORDER_D))
labels_C = [CYP46_CLUSTERS[c] for c in CLUSTER_ORDER_D]
fam_frac_C = fam_frac.loc[CLUSTER_ORDER_D]
left = np.zeros(len(CLUSTER_ORDER_D))
for fam in FAM_ORDER:
    vals = fam_frac_C[fam].to_numpy()
    axC.barh(y_C, vals, left=left, color=FAM_COLOR[fam], edgecolor='white',
             linewidth=0.3, label=fam if vals.sum() > 0 else None)
    left += vals
axC.set_yticks(y_C, labels=labels_C)
axC.invert_yaxis()
axC.set_xlim(0, 1.0)
axC.set_xlabel("Fraction of cluster members")
axC.set_title("C  Protein-family composition", loc='left', fontweight='bold')
axC.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
axC.set_xticklabels(['0', '25%', '50%', '75%', '100%'])
present_f = [f for f in FAM_ORDER if fam_frac_C[f].sum() > 0]
handles_f = [plt.Rectangle((0, 0), 1, 1, color=FAM_COLOR[f]) for f in present_f]
axC.legend(handles_f, present_f, loc='center left', bbox_to_anchor=(1.0, 0.5),
           frameon=False, title='Family', title_fontsize=9)

# Panel D — Taxonomic composition (right)
axD = axes[1]
y_D = np.arange(len(CLUSTER_ORDER_C))
labels_D = [CYP46_CLUSTERS[c] for c in CLUSTER_ORDER_C]
tax_frac_D = tax_frac.loc[CLUSTER_ORDER_C]
left = np.zeros(len(CLUSTER_ORDER_C))
for tax in TAX_ORDER:
    vals = tax_frac_D[tax].to_numpy()
    axD.barh(y_D, vals, left=left, color=TAX_COLOR[tax], edgecolor='white',
             linewidth=0.3, label=tax if vals.sum() > 0 else None)
    left += vals
axD.set_yticks(y_D, labels=labels_D)
axD.invert_yaxis()
axD.set_xlim(0, 1.0)
axD.set_xlabel("Fraction of cluster members")
axD.set_title("D  Taxonomic composition", loc='left', fontweight='bold')
axD.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
axD.set_xticklabels(['0', '25%', '50%', '75%', '100%'])
present = [t for t in TAX_ORDER if tax_frac_D[t].sum() > 0]
handles = [plt.Rectangle((0, 0), 1, 1, color=TAX_COLOR[t]) for t in present]
axD.legend(handles, present, loc='center left', bbox_to_anchor=(1.0, 0.5),
           frameon=False, title='Kingdom', title_fontsize=9)

fig.tight_layout()
pdf = OUT_DIR / "figure2_panelsCD_cyp46_composition.pdf"
png = OUT_DIR / "figure2_panelsCD_cyp46_composition.png"
fig.savefig(pdf, bbox_inches='tight')
fig.savefig(png, dpi=300, bbox_inches='tight')
print(f"\nWrote: {pdf}")
print(f"Wrote: {png}")
