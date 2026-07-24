# Methods — Steroid-interacting protein corpus construction

Reconstructed from the dated scripts under `/home/adsiordia/steroid_core_classifier/` and `/home/adsiordia/marimo_visualizer/MarimoSteroidVisualizer/`. Every step below is reproducible from the referenced script.

## Timeline

| Date | Step | Script | Input → Output |
|---|---|---|---|
| **2025-11-20** | Original seed pipeline (slow, scraping-based) | `20251120_steroid_uniprot.py` | `all_steroid_uniprot_unique_annot4plus.xlsx` → `all_steroid_uniprot_with_steroid_rhea_chebi_smiles.xlsx` |
| **2025-11-20** | Second-pass filter to rescue omitted entries | `20231120_omitted_secondpass.py` | ↑ + `omitted_nonsteroid_chebi_rhea_entries.xlsx` → `rescued_steroid_entries_from_second_pass.xlsx`, `..._FINAL.xlsx` |
| **2026-05-11** | Comprehensive rewrite (bulk Rhea TSV) | `20260511_comprehensive_steroid_uniprot.py` | Rhea + UniProt bulk downloads → `all_steroid_uniprot_comprehensive.xlsx` (34,131 rows) |
| **2026-05-11** | Union with old + binder entries | `20260511_add_binders_and_old_only.py` | Comprehensive + old_only + KW-0754/GO:0005496 binders → `all_steroid_uniprot_comprehensive_v2.xlsx` (**37,391 rows**) |
| **2026-05-11** | FASTA extraction | `20260511_make_fasta.py` | v2 xlsx → `all_steroid_uniprot_comprehensive_v2.fasta` |
| **2026-05-11** | ProtT5 embedding (SLURM) | `script_prot5.sh` calling `prott5_embedder.py` | v2 FASTA → `all_steroid_uniprot_comprehensive_v2.h5` (37,391 × 1024) |
| **May–Jun 2026** | UMAP fit + HDBSCAN clustering + CSV assembly | (in-notebook, `refit_full_umap.py`) | h5 → `protein_sequence_embedding.csv` (37,381 rows with 2D coords + clusters) |
| **~Jul 2026** | Exact-sequence deduplication | `demo/_dedupe_sequences.py` | 37,381 → **35,835 unique sequences** (1,546 duplicates removed, 4.1%) |
| **Jul 2026** | Literature-recruited additions + audit | `literature/scripts/*.py` | +15 papers-recruited entries; **35,834 final** |

## Detailed methodology

### Step 1 — Steroid ChEBI identification (`20260511_comprehensive_steroid_uniprot.py`, section B)

We downloaded the Rhea bulk chemistry file (`rhea-chebi-smiles.tsv`, from ftp.expasy.org/databases/rhea/tsv/), which lists every ChEBI participant of any Rhea reaction with its SMILES. Each of the **14,173 unique ChEBI SMILES** was tested against a sterane substructure query:

- **Sterane core query:** SMILES `C1CCC2CCC3C4CCCC4CCC3C2C1` (the 4-fused-ring steroid backbone, all single bonds, no heteroatoms).
- **Molecule pre-processing** (function `modify_molecule`): before matching, every ChEBI SMILES was RDKit-deprotonated + de-aromatized, and all N, O, S atoms were replaced with C. This lets the sterane match tolerate steroid derivatives with heteroatom substitutions (bile acids with -OH, steroid sulfates, N-containing analogs, etc.).
- **Match:** RDKit `HasSubstructMatch(sterane_core, useChirality=False)`.

**Result:** 875 ChEBI IDs classified as steroid-like → saved as `steroid_chebis.tsv`.

### Step 2 — Steroid Rhea reaction identification (same script, section C)

We downloaded `rhea-reaction-smiles.tsv` (36,014 Rhea reactions with reaction SMILES) and `rhea-directions.tsv` (master ID + LR + RL + BI direction IDs per reaction). For each reaction:

1. Split the reaction SMILES on `>>` to get reactants and products
2. Split each side on `.` to get individual molecule SMILES
3. Canonicalize each with RDKit
4. Match against the canonicalized steroid ChEBI SMILES set

**Result:** 2,210 Rhea reactions have ≥1 steroid participant → expanded across all 4 direction IDs → **4,420 steroid Rhea IDs** → saved as `steroid_rheas.tsv`.

### Step 3 — Protein corpus (same script, section D)

We downloaded `rhea2uniprot_sprot.tsv` (SwissProt) and `rhea2uniprot_trembl.tsv.gz` (TrEMBL) from Rhea. Filtering these to Rhea IDs that matched our steroid set:

- SwissProt: 8,791 rows
- TrEMBL: 586,934 rows
- **Union of unique accessions: 95,322**

### Step 4 — UniProt metadata fetching (same script, section E)

Metadata retrieved via the UniProt REST bulk-accessions endpoint (`https://rest.uniprot.org/uniprotkb/accessions`) in batches of 300, with automatic retry on 429/5xx. Fields requested: `accession, id, protein_name, gene_names, organism_name, length, sequence, annotation_score, rhea`.

**Filter:** `annotation_score ≥ 4` (UniProt-curated quality threshold — 4 or 5 stars). This retained **34,131 rows** → `all_steroid_uniprot_comprehensive.xlsx`.

### Step 5 — Union with additional sources (`20260511_add_binders_and_old_only.py`)

Two additional sources were unioned into the corpus:

**(a) Old-only rescue** (514 entries): entries present in the 2025-11-20 `all_steroid_uniprot_with_steroid_rhea_chebi_smiles_FINAL.xlsx` but *not* rediscovered by the 2026 comprehensive pass. These come from earlier sources (`source` column values):

| source tag | count |
|---|---:|
| `ChEBI; old_only_rescue` | 401 |
| `PDB_Uniprot; old_only_rescue` | 65 |
| `Manual; old_only_rescue` | 42 |
| misc combined tags | 6 |

**(b) Steroid binders** (2,746 entries): proteins annotated as steroid-binding in UniProt without necessarily catalyzing a steroid reaction. Query:

```
(keyword:KW-0754 OR go:0005496) AND (annotation_score:4 OR annotation_score:5)
```

- `KW-0754` = UniProt keyword "Steroid-binding"
- `GO:0005496` = Gene Ontology term "steroid binding"

**Union deduplicated by UniProt accession:** 37,391 rows → `all_steroid_uniprot_comprehensive_v2.xlsx`.

### Step 6 — FASTA generation (`20260511_make_fasta.py`)

Sequences extracted from the v2 xlsx, whitespace stripped, deduplicated by `Entry`, non-standard amino-acid characters flagged (but retained). Output: `all_steroid_uniprot_comprehensive_v2.fasta` (24.8 MB, 37,391 sequences).

### Step 7 — ProtT5 embedding (`script_prot5.sh`)

**Model:** `Rostlab/prot_t5_xl_half_uniref50-enc` (Elnaggar et al., IEEE TPAMI 2022, doi:10.1109/TPAMI.2021.3095381), encoder-only, half-precision.

**Embedding script:** `prott5_embedder.py` by Michael Heinzinger (bundled with the Rostlab ProtTrans repository).

**Compute:** SLURM job, 1 GPU, 20 CPUs, 50 GB RAM. Per-protein (mean-pooled) embeddings, 1024 dimensions.

**Output:** `all_steroid_uniprot_comprehensive_v2.h5` (37,391 × 1024 float32, ~150 MB).

### Step 8 — UMAP + HDBSCAN

2D UMAP embedding fit on the 37,391 × 1024 ProtT5 matrix (n_neighbors=15, min_dist=0.1). HDBSCAN clustering on the 2D output. Coordinates + cluster labels merged with the v2 metadata to produce `protein_sequence_embedding.csv` (37,381 rows; a few entries dropped due to missing sequences before embedding).

### Step 9 — Exact-sequence deduplication (`demo/_dedupe_sequences.py`)

Priority-based dedup preserving:

1. `is_new=1` rows (literature-recruited entries never dropped in favor of a UniProt sibling)
2. Row with the most non-null metadata columns (Rhea ID, ChEBI ID, EC, etc.)
3. Alphabetical `Entry` for tiebreaks

**Result:** 37,381 → **35,835 unique sequences** (1,546 exact-sequence duplicates collapsed, 4.1% reduction) → `protein_sequence_embedding.DEDUP.csv`.

Cluster-level dedup impact (top 5):

| Cluster | Before | After | Removed |
|---:|---:|---:|---:|
| 68 (BSH family) | 307 | 234 | −73 |
| 39 | 140 | 97 | −43 |
| 49 | 126 | 84 | −42 |
| 44 | 426 | 386 | −40 |
| 67 | 315 | 277 | −38 |

### Step 10 — Literature-recruited additions

+15 protein entries manually recruited from 5 audited papers with paper-specific evidence tiers (see `data/literature_recruited_proteins.csv` and `literature/scripts/`). Full audit corrected one wrong attribution (`A0A8F5DVT9` removed).

**Final atlas: 35,834 protein sequences.**

## Data-source citations

- **Rhea** — Bansal et al., *Nucleic Acids Res*. 50(D1):D622–D631 (2022). https://www.rhea-db.org
- **UniProt** — The UniProt Consortium, *Nucleic Acids Res*. 51(D1):D523–D531 (2023). https://www.uniprot.org
- **ChEBI** — Hastings et al., *Nucleic Acids Res*. 44(D1):D1214–D1219 (2016). https://www.ebi.ac.uk/chebi
- **ProtT5** — Elnaggar et al., *IEEE TPAMI* (2022). doi:10.1109/TPAMI.2021.3095381
- **RDKit** — https://www.rdkit.org
- **UMAP** — McInnes et al., *arXiv* 1802.03426 (2018)
- **HDBSCAN** — Campello et al., *ACM TKDD* 10(1):5 (2015)

## Reproducibility — exact commands

```bash
cd /home/adsiordia/steroid_core_classifier

# Step 1-4 — corpus + metadata
python 20260511_comprehensive_steroid_uniprot.py

# Step 5 — union with old + binders
python 20260511_add_binders_and_old_only.py

# Step 6 — extract FASTA
python 20260511_make_fasta.py

# Step 7 — embed (SLURM)
sbatch script_prot5.sh

# Step 8-9 — done in the marimo project directory:
cd /home/adsiordia/marimo_visualizer/MarimoSteroidVisualizer
# (UMAP + HDBSCAN in refit_full_umap.py / notebook)
python demo/_dedupe_sequences.py
```
