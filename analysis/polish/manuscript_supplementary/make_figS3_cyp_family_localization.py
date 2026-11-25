"""Taxonomic composition of selected CYP-enriched clusters (manuscript Supp. Fig. S3).

Promoted from analysis/polish/_archive/23_cyp_singlecluster_taxonomic_signature.py
with manuscript label edits; data numbers are unchanged (verified at k=82).

Five CYP subfamilies land in exactly one majority-P450 cluster each and
contain the vast majority of that subfamily's members in the atlas:

    CYP1A   → C20    (Cytochrome P450 1A)
    CYP51A  → C29    (Lanosterol 14α-demethylase)
    CYP17A  → C46    (Steroid 17α-hydroxylase)
    CYP90B  → C47    (Brassinosteroid 6-oxidase, plant)
    CYP3A   → C57    (Cytochrome P450 3A)

Unlike CYP46A (4 clusters) and CYP39A (3 clusters), these families don't
split across the atlas. This script asks whether that "cleanness" is
reflected taxonomically: do the five clusters draw from the same pool
of organisms, or does each pick up a distinct taxonomic order?

Outputs:
  figS3a_cyp_family_localization.{pdf,png}
  figS3a_cyp_family_localization_data.tsv
"""
from __future__ import annotations
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]  # steroid_atlas/
OUT  = HERE

# ── CYP family detection (regex, matches the Fig 2 panel-A/B script) ──
CYP_PAT = [re.compile(r"CYP(\d+)([A-Z])?\d*"),
           re.compile(r"CYTOCHROME\s+P-?450\s+(\d+)([A-Z])?")]
def cyp_subfamily(gene: str, pname: str) -> str:
    for src in (str(gene).upper(), str(pname).upper()):
        for pat in CYP_PAT:
            m = pat.search(src)
            if m:
                return f"CYP{m.group(1)}{m.group(2) or ''}"
    return ""

# ── Single-cluster CYP families ──
# key = proteins.csv.cluster (0-indexed); display = cluster + 1
# NOTE on C47: only 4 of 11 members are CYP90B (the steroid C-22
# hydroxylases, EC 1.14.14.178). The other 7 are sister subfamilies in
# the SAME plant brassinosteroid-biosynthesis pathway (CYP85A, CYP90A/D/G,
# CYP724, CYP72). So this row shows the whole BR-pathway cluster — not
# only the CYP90B subset — and is relabelled accordingly.
SINGLE_CLUSTER = [
    # cid, display, subfamily_filter (None = all cluster members), label, desc
    (19, 20, "CYP1A",  "CYP1A",   "Cytochrome P450 1A"),
    (28, 29, "CYP51A", "CYP51A",  "Lanosterol 14α-DM"),
    (45, 46, "CYP17A", "CYP17A",  "Steroid 17α-OH"),
    (46, 47, None,     "BR-path", "Plant brassinosteroid P450s"),
    (56, 57, "CYP3A",  "CYP3A",   "Cytochrome P450 3A"),
]

# ── Taxonomic classifier (regex on genus names) ──
# High-level kingdom/class groups only — Mammals (all orders, placental
# + marsupials + monotremes), Birds, Reptiles, Amphibians, Fish, Other
# vertebrates (jawless fish etc.), Invertebrates, Plants, Fungi, Bacteria.
TAXON_PATTERNS = [
    ("Plants", re.compile(r"\b(Arabidopsis|Oryza|Paris|Musa|Trigonella|"
        r"Nicotiana|Zea|Solanum|Spinacia|Vitis|Gossypium|Glycine|"
        r"Brachypodium|Hordeum|Triticum|Sorghum|Setaria|Medicago|"
        r"Cucumis|Fragaria|Malus|Prunus|Populus|Eucalyptus|Pinus|Picea|"
        r"Physcomitrella|Marchantia|Selaginella|Chara|Chlamydomonas|"
        r"Ostreococcus|Micromonas|Volvox|Salvia|Catharanthus|Mentha)\b")),
    ("Mammals", re.compile(r"\b("
        # Primates
        r"Homo|Pan|Gorilla|Pongo|Macaca|Papio|Saimiri|Aotus|Callithrix|"
        r"Saguinus|Rhinopithecus|Prolemur|Lemur|Nomascus|Hylobates|"
        r"Cercopithecus|Cercocebus|Chlorocebus|Eulemur|Trachypithecus|"
        r"Mandrillus|Theropithecus|Hapalemur|Microcebus|Tarsius|"
        r"Daubentonia|Carlito|Nycticebus|Otolemur|Loris|Propithecus|"
        r"Pithecia|Cebus|Erythrocebus|Alouatta|Ateles|Nasalis|Colobus|"
        r"Piliocolobus|Presbytis|Semnopithecus|Oedipomidas|"
        # Rodentia
        r"Mus|Rattus|Spermophilus|Cavia|Nannospalax|Spalax|Microtus|"
        r"Cricetulus|Peromyscus|Mesocricetus|Dipodomys|Chinchilla|"
        r"Castor|Marmota|Octodon|Hystrix|Sciurus|Jaculus|Fukomys|"
        r"Heterocephalus|Ictidomys|Urocitellus|Perognathus|Myodes|"
        r"Sigmodon|Onychomys|Apodemus|Neotoma|Myocastor|Cricetomys|"
        r"Meriones|Gerbillus|Dasyprocta|Hydrochoerus|Ondatra|"
        r"Lasiopodomys|Tupaia|"  # vole; Scandentia (tree shrews)

        # Carnivora
        r"Canis|Felis|Panthera|Ursus|Neovison|Mustela|Vulpes|Ailuropoda|"
        r"Procyon|Odobenus|Ailurus|Lontra|Enhydra|Meles|Zalophus|"
        r"Eumetopias|Pusa|Phoca|Leopardus|Lynx|Hyaena|Suricata|Nasua|"
        r"Puma|Acinonyx|Gulo|Taxidea|Nyctereutes|Thalarctos|Neomonachus|"
        r"Monachus|Leptonychotes|Mirounga|Halichoerus|Callorhinus|"
        r"Crocuta|Spilogale|Mephitis|Conepatus|Bassariscus|Potos|"
        r"Arctocephalus|Paradoxurus|Herpestes|Genetta|Civettictis|"
        # Cetacea
        r"Balaenoptera|Delphinus|Tursiops|Megaptera|Orcinus|Globicephala|"
        r"Lagenorhynchus|Stenella|Physeter|Monodon|Phocoena|Eubalaena|"
        r"Berardius|Balaena|Neophocaena|Mesoplodon|Kogia|Inia|"
        r"Delphinapterus|Sousa|Eschrichtius|Lagenodelphis|Peponocephala|"
        r"Hyperoodon|Platanista|Pontoporia|Ziphius|Steno|Lipotes|"
        # Artiodactyla
        r"Bos|Capra|Ovis|Sus|Camelus|Vicugna|Bubalus|Cervus|Odocoileus|"
        r"Rangifer|Giraffa|Hippopotamus|Taurotragus|Tragelaphus|Alces|"
        r"Dama|Antilope|Boselaphus|Syncerus|Tayassu|Lama|Antilocapra|"
        r"Ovibos|Addax|Kobus|Hemitragus|Pecari|Pantholops|Moschus|"
        r"Catagonus|Muntiacus|Hyemoschus|Tragulus|Potamochoerus|Bison|"
        r"Phacochoerus|Hylochoerus|Babyrousa|Okapia|Gazella|Saiga|"
        r"Alcelaphus|Hippotragus|Oryx|Damaliscus|"
        # Perissodactyla
        r"Equus|Tapirus|Ceratotherium|Diceros|Rhinoceros|"
        # Chiroptera
        r"Rhinolophus|Rousettus|Myotis|Pteropus|Hipposideros|Miniopterus|"
        r"Pipistrellus|Eptesicus|Desmodus|Phyllostomus|Artibeus|"
        r"Carollia|Molossus|Tadarida|Vespertilio|Eonycteris|Macroglossus|"
        r"Nyctalus|Lasiurus|Noctilio|Sturnira|Vampyrum|Pteronotus|"
        # Lagomorpha
        r"Oryctolagus|Lepus|Ochotona|Sylvilagus|"
        # Sirenia, Dermoptera, Afrosoricida, Eulipotyphla, Tubulidentata,
        # Xenarthra/Cingulata, Proboscidea, Hyracoidea
        r"Trichechus|Dugong|Galeopterus|Cynocephalus|Echinops|"
        r"Chrysochloris|Microgale|Tenrec|Galemys|Erinaceus|Sorex|"
        r"Solenodon|Condylura|Scalopus|Talpa|Suncus|Orycteropus|"
        r"Choloepus|Bradypus|Dasypus|Tolypeutes|Myrmecophaga|Tamandua|"
        r"Cyclopes|Loxodonta|Elephas|Mammuthus|Procavia|Heterohyrax|"
        r"Dendrohyrax|"
        # Marsupials
        r"Sarcophilus|Dasyurus|Antechinus|Thylacinus|Phascolarctos|"
        r"Vombatus|Trichosurus|Macropus|Wallabia|Dendrolagus|Setonix|"
        r"Petaurus|Petauroides|Pseudocheirus|Tarsipes|Cercartetus|"
        r"Acrobates|Monodelphis|Didelphis|Thylamys|Marmosa|Caluromys|"
        r"Perameles|Isoodon|Macrotis|"
        # Monotremes
        r"Ornithorhynchus|Tachyglossus|Zaglossus"
        r")\b")),
    ("Birds", re.compile(r"\b("
        # Base list
        r"Gallus|Taeniopygia|Zonotrichia|Cyanistes|Columba|Anas|Falco|"
        r"Haliaeetus|Melopsittacus|Serinus|Corvus|Pavo|Struthio|Parus|"
        r"Passer|Ficedula|Catharus|Meleagris|Coturnix|Anser|Branta|"
        r"Nothoprocta|Tyto|Cuculus|Eudyptula|Poephila|Aquila|Buteo|"
        r"Dromaius|Apteryx|Casuarius|"
        # Added from CYP46A1 clusters
        r"Zosterops|Eudyptes|Amazona|Geospiza|Chrysolophus|Loxia|"
        r"Cyanoderma|Nothocercus|Sylvia|Grus|Chloropsis|Anhinga|Prunella|"
        r"Certhia|Malurus|Nyctibius|Centropus|Cairina|Vidua|Pterocles|"
        r"Aptenodytes|Lepidothrix|Pygoscelis|Chaetops|Buphagus|Oenanthe|"
        r"Zapornia|Jacana|Pluvianellus|Indicator|Serilophus|Syrrhaptes|"
        r"Pelecanoides|Sitta|Pardalotus|Thalassarche|Drymodes|Turnix|"
        r"Horornis|Nycticryphes|Oreotrochilus|Xiphorhynchus|Cettia|"
        r"Cepphus|Illadopsis|Gymnorhina|Uria|Molothrus|Erpornis|"
        r"Hypocryptadius|Brachypteracias|Balaeniceps|Calonectris|Aramus|"
        r"Thinocorus|Helopsaltes|Cisticola|Oriolus|Bombycilla|Himantopus|"
        r"Scytalopus|Leiothrix|Heliornis|Rynchops|Alaudala|Anthoscopus|"
        r"Quiscalus|Sagittarius|Semnornis|Cinclus|Smutsornis|Hippolais|"
        r"Rhadina|Pycnonotus|Phainopepla|Oxylabes|Neodrepanis|"
        r"Mystacornis|Pomatostomus|Atlantisia|Passerina|Fregetta|"
        r"Spheniscus|Larus|Nicator|Tricholaema|Spizella|Urocolius|"
        r"Pteruthius|Atrichornis|Accipiter|Bubo|Calidris|Junco|Otus|"
        r"Ploceus|Peucedramus|Ramphastos|Circaetus|Piaya|Phaetusa|"
        r"Nyctiprogne|Eurystomus|Pampusana|Podargus|Scopus|Pomatorhinus|"
        r"Ceuthmochares|Callaeas|Glareola|Ceyx|Lamprotornis|Fregata|"
        r"Aegithalos|Calcarius|Tichodroma|Sula|Picathartes|Dryoscopus|"
        r"Pitta|Donacobius|Corythaeola|Campylorhamphus|Penelope|"
        r"Tychaedon|Elachura|Eolophus|Halcyon|Oxyruncus|Edolisoma|"
        r"Tyrannus|Regulus|Crotophaga|Chunga|Probosciger|Polioptila|"
        r"Pachyramphus|Motacilla|Furnarius|Tachuris|Mionectes|Emberiza|"
        r"Sinosuthora|Columbina|Neopipo|Melospiza|Crypturellus|Geococcyx|"
        r"Toxostoma|Vireo|Alectura|Upupa|Falcunculus|Oreocharis|"
        r"Machaerirhynchus|Dasyornis|Grantiella|Origma|Ptilonorhynchus|"
        r"Onychorhynchus|Cardinalis|Dromas|Hylia|Platysteira|"
        r"Urocynchramus|Lanius|Erythrocercus|Chroicocephalus|Phaethon|"
        r"Podiceps|Manacus|Fulmarus|Dryobates|Gavia|Opisthocomus|"
        r"Nipponia|Mesitornis|Merops|Acanthisitta|Cathartes|Chlamydotis|"
        r"Egretta|Calypte|Buceros|Antrostomus|Charadrius|Lonchura|Pipra|"
        r"Strigops|Phasianus|Chloebia|Limosa|Panurus|Pedionomus|Ifrita|"
        r"Rhabdornis|Rhipidura|Chloroceryle|Galbula|Sterrhoptilus|Irena|"
        r"Balearica|Oceanodroma|Dicaeum|Sclerurus|Myiagra|Bucco|"
        r"Baryphthengus|Hemiprocne|Psilopogon|Grallaria|Sakesphorus|"
        r"Psophia|Odontophorus|Dicrurus|Ciconia|Leptocoma|Chordeiles|"
        r"Podilymbus|Glaucidium|Setophaga|Mesembrinibis|Formicarius|"
        r"Sylvietta|Chauna|Piprites|Arenaria|Corythaixoides|Trogon|"
        r"Rostratula|Spizaetus|Eubucco|Ciccaba|Climacteris|Chionis|"
        r"Sapayoa|Nesospiza|Acrocephalus|Brachypodius|Asarcornis|"
        r"Agelaius|Erithacus|Pheucticus|Oceanites|Alca|Promerops|"
        r"Steatornis|Notiomystis|Aegotheles|Caloenas|Rhynochetos|"
        r"Ibidorhyncha|Eudromia|Mohoua|Thryothorus|Rhodinocichla|"
        r"Smithornis|Cochlearius|Ptilorrhoa|Ardeotis|Lophotis|"
        r"Chaetorhynchus|Aleadryas|Orthonyx|Struthidea|Leucopsar|"
        r"Eulacestoma|Rhagologus|Cnemophilus|Melanocharis|Menura"
        r")\b")),
    ("Reptiles", re.compile(r"\b("
        r"Thamnophis|Python|Boa|Crotalus|Alligator|Crocodylus|Testudo|"
        r"Trachemys|Varanus|Gekko|Anolis|Pogona|Ophiophagus|Naja|"
        r"Lacerta|Pelodiscus|Gallotia|Elgaria|Heloderma|Eublepharis|"
        r"Chrysemys|Chelonia|Caretta|Dermochelys|Sphenodon|Iguana|"
        # Added from CYP46A1 clusters
        r"Pantherophis|Podarcis|Gopherus|Terrapene|Phrynocephalus|"
        r"Salvator|Mauremys|Pseudonaja|Platysternon"
        r")\b")),
    ("Amphibians", re.compile(r"\b("
        r"Xenopus|Rana|Pleurodeles|Engystomops|Bufo|Ambystoma|Hyla|Pipa|"
        r"Cynops|Ceratophrys|Microhyla|Odorrana|Eleutherodactylus|"
        r"Litoria|Lissotriton|Salamandra|Taricha|Triturus|Nanorana|"
        r"Dendrobates|Mantella|"
        # Added from CYP46A1 clusters
        r"Leptobrachium|Pyxicephalus|Hymenochirus|Pelobates|Geotrypetes|"
        r"Microcaecilia"
        r")\b")),
    ("Fish", re.compile(r"\b("
        r"Oncorhynchus|Mastacembelus|Danio|Takifugu|Oryzias|Gasterosteus|"
        r"Esox|Thunnus|Salmo|Clupea|Gadus|Cyprinus|Astyanax|Latimeria|"
        r"Lepisosteus|Xiphophorus|Poecilia|Tetraodon|Scomber|Sparus|"
        r"Dicentrarchus|Lates|Morone|Perca|Carassius|Ictalurus|Acipenser|"
        r"Callorhinchus|Rhincodon|Scyliorhinus|Protopterus|Neoceratodus|"
        # Added from CYP46A1 clusters
        r"Sinocyclocheilus|Neogobius|Oreochromis|Myripristis|Sphaeramia|"
        r"Scleropages|Hucho|Electrophorus|Maylandia|Amphiprion|Salarias|"
        r"Cyclopterus|Amphilophus|Astatotilapia|Mugilogobius|Anabas|"
        r"Nothobranchius|Kryptolebias|Echeneis|Erpetoichthys|Pygocentrus|"
        r"Austrofundulus|Denticeps|Cirrhinus|Danionella|Channa|"
        r"Neolamprologus|Gouania|Umbra|Triplophysa|Knipowitschia|Chanos|"
        r"Haplochromis|Fundulus|Pundamilia|Periophthalmus|Onychostoma|"
        r"Dissostichus|Betta|Cyprinodon|Crenichthys|Gambusia|"
        r"Acanthochromis|Seriola|Alosa|Monopterus|Albula|Anguilla|"
        r"Coregonus|Champsocephalus|Collichthys|Labrus|Hippocampus|"
        r"Cynoglossus|Paramormyrops|Scophthalmus|Pleuronectes|"
        r"Synaphobranchus|Stegastes|Eleginops|Phoxinus|Xyrichtys|Solea|"
        r"Zoarces|Culter|Gymnodraco|Parambassis|Pagothenia|Notothenia|"
        r"Cottoperca"
        r")\b")),
    # Jawless fish (lampreys, hagfish) and other paraphyletic vertebrates
    ("Other vertebrates", re.compile(
        r"\b(Petromyzon|Lampetra|Eptatretus|Myxine|Branchiostoma|Ciona)\b")),
    ("Invertebrates", re.compile(r"\b("
        r"Drosophila|Anopheles|Aedes|Culex|Apis|Bombus|Bombyx|Tribolium|"
        r"Caenorhabditis|Octopus|Aplysia|Crassostrea|Lymnaea|"
        r"Biomphalaria|Strongylocentrotus|Nematostella|Hydra|Trichoplax|"
        r"Capitella|Daphnia|Ixodes|Rhipicephalus|Schistosoma|Brugia|Loa|"
        r"Onchocerca|Trichinella|Anisakis|Haemonchus|"
        # Added: corals, sea stars, cucumbers, anemones from CYP46A1
        r"Pocillopora|Actinia|Acanthaster|Holothuria|Stylophora|Patiria|"
        r"Stichopus|Desmophyllum|Acropora"
        r")\b")),
    ("Fungi", re.compile(r"\b(Saccharomyces|Candida|Aspergillus|Neurospora|"
        r"Schizosaccharomyces|Cryptococcus|Ustilago|Yarrowia|Debaryomyces|"
        r"Kluyveromyces|Komagataella|Penicillium|Fusarium|Magnaporthe|"
        r"Trichoderma|Trichophyton|Malassezia|Rhizopus|Mucor)\b")),
    ("Bacteria", re.compile(r"\b(Escherichia|Mycobacterium|Streptomyces|"
        r"Bacillus|Clostridium|Pseudomonas|Staphylococcus|Salmonella|"
        r"Bacteroides|Lactobacillus|Enterococcus|Bifidobacterium|"
        r"Streptococcus|Vibrio|Helicobacter|Rhodopseudomonas|Rhodobacter|"
        r"Synechocystis|Nostoc|Anabaena|Thermus|Deinococcus|"
        r"Acinetobacter|Klebsiella|Yersinia|Francisella|Burkholderia|"
        r"Chlamydia|Rickettsia|Treponema|Borrelia|Leptospira|"
        r"Methanococcus|Methanobacterium|Sulfolobus|Pyrococcus|"
        r"Halobacterium)\b")),
]

def classify(org: str) -> str:
    for name, pat in TAXON_PATTERNS:
        if pat.search(org):
            return name
    return "Other"

TAXON_ORDER = [
    "Mammals", "Birds", "Reptiles", "Amphibians", "Fish",
    "Other vertebrates", "Invertebrates",
    "Plants", "Fungi", "Bacteria", "Other",
]
TAXON_COLOR = {
    "Mammals":           "#dc2626",
    "Birds":             "#f97316",
    "Reptiles":          "#eab308",
    "Amphibians":        "#65a30d",
    "Fish":              "#2563eb",
    "Other vertebrates": "#64748b",
    "Invertebrates":     "#a855f7",
    "Plants":            "#16a34a",
    "Fungi":             "#be185d",
    "Bacteria":          "#0d9488",
    "Other":             "#cbd5e1",
}
LIGHT_TAXA = {"Reptiles"}  # dark text on these

# ── Load ──
print("Loading proteins.csv ...")
prot = pd.read_csv(ROOT / "data" / "proteins.csv", low_memory=False)
prot["organism"]      = prot["organism"].fillna("").astype(str)
prot["gene_names"]    = prot["gene_names"].fillna("").astype(str)
prot["protein_names"] = prot["protein_names"].fillna("").astype(str)
prot["cyp_sub"] = prot.apply(
    lambda r: cyp_subfamily(r["gene_names"], r["protein_names"]), axis=1)
prot["taxon"]   = prot["organism"].apply(classify)

# ── Build the per-cluster tally ──
# For CYP1A/3A/17A/51A clusters we filter to the dominant subfamily so
# the bar shows that family's taxonomy. For C47 (BR-pathway), we take
# every cluster member because the cluster spans multiple sister
# subfamilies in the same plant pathway.
rows = []
for cid, disp, sub_fam, label, desc in SINGLE_CLUSTER:
    if sub_fam is None:
        sub = prot[prot["cluster"] == cid]
        total_fam = len(sub)         # no "elsewhere in atlas" comparator
        confined  = 100.0
    else:
        sub = prot[(prot["cluster"] == cid) & (prot["cyp_sub"] == sub_fam)]
        total_fam = int((prot["cyp_sub"] == sub_fam).sum())
        confined  = len(sub) / total_fam * 100 if total_fam else 0
    row = {"cid": cid, "disp": disp, "subfamily": label, "desc": desc,
           "n": len(sub), "n_subfamily_total": total_fam,
           "confinement_pct": round(confined, 1)}
    for t in TAXON_ORDER:
        row[t] = int((sub["taxon"] == t).sum())
    rows.append(row)
df = pd.DataFrame(rows)
df.to_csv(OUT / "figS3a_cyp_family_localization_data.tsv", sep="\t", index=False)
print(df[["disp", "subfamily", "n", "n_subfamily_total",
          "confinement_pct", "Other"]].to_string(index=False))

# ── Figure: horizontal stacked bars ──
plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 12, "axes.labelsize": 11,
    "xtick.labelsize": 9, "ytick.labelsize": 10, "legend.fontsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.facecolor": "white",
})
fig, ax = plt.subplots(figsize=(11.5, 4.6))
# Top row is the first entry of SINGLE_CLUSTER → reverse y so C20 is at top.
y_positions = np.arange(len(df))[::-1]

labels = [f"C{r.disp} · {r.subfamily}\n({r.desc}, n={r.n})"
          if r.subfamily != "BR-path"
          else f"C{r.disp} · Plant BR P450s\n(CYP85A/90A/B/D/G/724/72, n={r.n})"
          for r in df.itertuples()]

left = np.zeros(len(df))
for t in TAXON_ORDER:
    vals = df[t].values / df["n"].values * 100
    ax.barh(y_positions, vals, left=left, color=TAXON_COLOR[t],
            edgecolor="white", linewidth=0.6, label=t, height=0.78)
    text_color = "#1e293b" if t in LIGHT_TAXA else "white"
    for yi, v, lf in zip(y_positions, vals, left):
        if v >= 7:
            ax.text(lf + v/2, yi, f"{v:.0f}%",
                    ha="center", va="center", fontsize=8.5,
                    color=text_color, fontweight="700")
    left = left + vals

# Right-edge "x% of the subfamily lives here" note (skip for BR-path
# since it's not a single CYP subfamily to compare against)
for yi, r in zip(y_positions, df.itertuples()):
    if r.subfamily == "BR-path":
        note = "multi-subfamily\nplant BR cluster"
    else:
        note = f"{r.confinement_pct:.0f}% of atlas\n{r.subfamily} proteins"
    ax.text(101.5, yi, note, ha="left", va="center",
            fontsize=8.5, color="#334155")

ax.set_yticks(y_positions); ax.set_yticklabels(labels)
ax.set_xlim(0, 100)
ax.set_xlabel("Composition (%)")
ax.set_title(
    "Taxonomic composition of selected CYP-enriched clusters",
    loc="left", pad=10)
ax.grid(axis="x", alpha=0.25, linewidth=0.4)

# Legend lists only taxa that actually appear
present = [t for t in TAXON_ORDER if df[t].sum() > 0]
handles = [Patch(facecolor=TAXON_COLOR[t], edgecolor="black",
                 linewidth=0.3, label=t) for t in present]
ax.legend(handles=handles, loc="center left",
          bbox_to_anchor=(1.18, 0.5), frameon=False,
          fontsize=9, title="Taxonomic group", title_fontsize=10)

fig.tight_layout()
fig.savefig(OUT / "figS3a_cyp_family_localization.pdf",
            bbox_inches="tight")
fig.savefig(OUT / "figS3a_cyp_family_localization.png",
            dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  wrote {OUT / 'figS3a_cyp_family_localization.{pdf,png}'}")
