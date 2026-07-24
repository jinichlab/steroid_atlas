# Methods — Steroid-interacting protein corpus construction

The atlas draws steroid-interacting proteins from **two complementary evidence sources** — Rhea-catalogued steroid-catalyzing reactions and UniProt-annotated steroid-binding proteins — augmented with a small set of literature-recruited entries from a targeted paper audit. The union is deduplicated by UniProt accession, ProtT5-embedded, 2D-UMAP projected, and deduplicated a second time by exact sequence identity.

## 1. Sterane substructure classifier — chemistry filter

Every candidate steroid compound is identified by an **RDKit substructure match against a sterane backbone**. The classifier is intentionally permissive so it catches modified steroids (bile acids, sulfates, glucuronides, N-containing steroid alkaloids, oxidized/reduced ring variants) rather than only exact carbon skeletons.

**Query molecule** — the fully-saturated sterane core, all single bonds, no heteroatoms:

```
SMILES: C1CCC2CCC3C4CCCC4CCC3C2C1
```

This is the 4-fused-ring backbone (6-6-6-5 cyclopentanoperhydrophenanthrene) shared by every steroid.

**Pre-processing applied to each candidate compound before matching.** For each ChEBI SMILES we build an RDKit `Mol` object and then normalize it in three steps:

1. **De-aromatize every atom and bond** — steroid rings in some ChEBI records are drawn aromatic (e.g., estrone's A ring). The query is fully saturated, so aromaticity would prevent a match.
2. **Replace every N, O, S atom with C** — bile acids carry hydroxyls, steroid sulfates carry sulfate esters, steroid alkaloids carry ring nitrogens. Substituting heteroatoms with carbon lets the query recognize these as *"same skeleton, different substituents"* rather than rejecting them for atomic-identity mismatch.
3. **Reduce every non-single bond to a single bond** — double bonds in specific ring positions (e.g., Δ4-3-ketosteroids) would prevent an exact match to the fully-saturated query.

**Match call**: `mol_modified.HasSubstructMatch(sterane_core, useChirality=False)`. Ignoring chirality is intentional — the classifier is a "does this molecule contain a steroid skeleton?" filter, not a stereochemistry check.

**Applied to Rhea's chemistry catalog**: the `rhea-chebi-smiles.tsv` bulk file lists every ChEBI compound appearing as a participant in any Rhea reaction, together with its canonical SMILES.

| Input | Count |
|---|---:|
| Total ChEBI participants in Rhea | 14,173 |
| After sterane substructure filter | **875 steroid ChEBIs** |

The 875 steroid ChEBIs are exported to `steroid_chebis.tsv` for downstream reaction matching.

## 2. Protein evidence sources

### Source A — Steroid-catalyzing enzymes via Rhea

Rhea is an expert-curated database of biochemical reactions where every reaction is written as a reaction SMILES built from ChEBI components. The path from "steroid ChEBI" → "steroid-catalyzing enzyme" has four sub-steps:

**Step A1 — Identify steroid reactions.** The `rhea-reaction-smiles.tsv` bulk file lists every Rhea reaction ID together with its reaction SMILES in `reactants>>products` form. Each reaction SMILES is:

1. Split on `>>` to get the reactant side and product side
2. Each side split on `.` to get individual molecule SMILES
3. Each molecule canonicalized with RDKit (`Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)`)
4. Canonical SMILES compared against the canonicalized steroid ChEBI SMILES from §1

A reaction is retained if at least one of its molecules (reactant or product) matches a steroid ChEBI.

| Input | Count |
|---|---:|
| Total Rhea reactions | 36,014 |
| Reactions with ≥1 steroid participant | **2,210** |

**Step A2 — Expand across reaction directions.** Every Rhea reaction has up to four IDs registered in `rhea-directions.tsv`: `MASTER` (undirected), `LR` (left-to-right), `RL` (right-to-left), and `BI` (bidirectional). Downstream Rhea→UniProt mappings can reference any of these IDs, so we expand the 2,210 reactions into the full set of direction-IDs any of them can carry.

| Input | Count |
|---|---:|
| Master IDs (unique reactions) | 2,210 |
| Expanded to all 4 directions | **4,420 IDs** |

**Step A3 — Map Rhea → UniProt.** Rhea ships two bulk protein mapping files:

- `rhea2uniprot_sprot.tsv` — Swiss-Prot reviewed entries (higher quality)
- `rhea2uniprot_trembl.tsv.gz` — TrEMBL unreviewed entries (larger coverage)

Both are filtered to rows whose Rhea ID is in the 4,420-ID steroid set from step A2:

| Source | Rows matching steroid Rhea | Notes |
|---|---:|---|
| Swiss-Prot (`rhea2uniprot_sprot`) | 8,791 | curated |
| TrEMBL (`rhea2uniprot_trembl`) | 586,934 | auto-annotated |
| Union of unique UniProt accessions | **95,322** | one row per accession |

**Step A4 — Fetch metadata + quality filter.** For each of the 95,322 candidate accessions we retrieve full metadata from the UniProt REST bulk endpoint:

- **Endpoint**: `https://rest.uniprot.org/uniprotkb/accessions`
- **Batch size**: 300 accessions per request (URL-length constrained)
- **Fields**: `accession, id, protein_name, gene_names, organism_name, length, sequence, annotation_score, rhea`
- **Format**: TSV
- **Retry policy**: on HTTP 429, 500, 502, 503, 504 → exponential backoff (2, 4, 6, 8 s), max 4 attempts

Metadata rows are accumulated across all batches (incremental save to `steroid_uniprot_metadata.tsv` every 50 batches for crash recovery).

The final filter is **`annotation_score ≥ 4`** — UniProt's curated-quality score, on a 1-5 scale. Levels 4 and 5 correspond to "high" or "highest" evidence (multiple experimental references, curator review, or complete Swiss-Prot annotation). This filter is the main reason the TrEMBL half of the raw union (586,934 rows → mostly auto-annotated homology hypotheses) collapses down: only a small fraction of TrEMBL entries score ≥ 4.

Each retained protein is then augmented with **per-protein steroid annotations**:

- `Steroid_Rhea_numeric` — the specific Rhea IDs (of its full reaction set) that involve a steroid
- `Steroid_ChEBI_grouped` — the ChEBI IDs of the steroid participants in those reactions, grouped by reaction
- `Steroid_SMILES_grouped` — the corresponding SMILES

So every row carries not just "this enzyme touches a steroid" but "this enzyme catalyzes reactions R1, R2, R3, in which steroid ChEBIs C1, C2 appear."

| Result | Count |
|---|---:|
| Candidate accessions (before quality filter) | 95,322 |
| **After `annotation_score ≥ 4` filter** | **34,131** |

Source A contributes 91.3% of the atlas.

### Source B — Annotated steroid-binding proteins

UniProt REST query:

```
(keyword:KW-0754 OR go:0005496) AND (annotation_score:4 OR annotation_score:5)
```

- `KW-0754` = UniProt keyword *Steroid-binding*
- `GO:0005496` = Gene Ontology term *steroid binding*

Captures proteins with annotated steroid-binding function that may not catalyze a Rhea-registered steroid reaction (receptors, transporters, carrier proteins).

→ **2,746 binding entries.**

### Union

The two automated sources are unioned by UniProt accession, with duplicates collapsed.

→ **36,877 unique UniProt accessions.**

## 3. Sequence embedding

Sequences are extracted to FASTA, whitespace stripped, deduplicated by accession, and embedded with **ProtT5** (`Rostlab/prot_t5_xl_half_uniref50-enc`, encoder-only, half-precision). Per-residue embeddings are mean-pooled to a **1,024-dimensional per-protein vector**.

Compute: single GPU, 20 CPU threads, 50 GB RAM.

## 4. Dimensionality reduction + clustering

- **UMAP** (n_neighbors=15, min_dist=0.1, 2D output) fit on the 36,877 × 1,024 ProtT5 matrix.
- **HDBSCAN** clustering on the 2D UMAP output produces cluster labels.

## 5. Exact-sequence deduplication

The corpus is deduplicated by exact amino-acid sequence identity. When multiple entries share a sequence, the surviving row is selected by:

1. Literature-recruited entries (see §7) are always preserved,
2. Otherwise the row with the most non-null metadata columns wins,
3. Alphabetical UniProt accession as tiebreaker.

→ **35,835 unique sequences** (1,546 exact duplicates removed, 4.1%).

## 6. Literature-recruited additions

15 protein entries were manually recruited from an audit of five recent papers characterizing new steroid-metabolizing activities (Rimal 2024, Guzior 2024, McCurry 2024, Jacoby 2025, Arp 2025). Each is annotated with an evidence tier reflecting the specificity of the paper's demonstration for that individual protein — biochemical characterization, in vivo genetic evidence, or bioinformatic assignment. Sequences are fetched from UniProt, RefSeq, NCBI Genome, or the papers' supplementary tables; ProtT5-embedded; and projected into the existing UMAP via cosine-nearest-neighbor in ProtT5 space (the UMAP model is inherited rather than refit to preserve the reference layout).

Full provenance is in `data/literature_recruited_proteins.csv`.

## Final atlas

**35,349 protein sequences** (final after literature audit removed one incorrectly attributed entry).

## Corpus composition by evidence source

| Source | Entries | % of union |
|---|---:|---:|
| A. Steroid-catalyzing (Rhea) | 34,131 | 92.6% |
| B. Steroid-binding (UniProt KW-0754 / GO:0005496) | 2,746 | 7.4% |
| **Union (dedup by accession)** | **36,877** | 100% |
| After exact-sequence dedup | 35,349 | |
| After literature audit (+15 STARs) | **35,349** | |

## Data sources cited

- **Rhea** — Bansal et al., *Nucleic Acids Res*. 50(D1):D622–D631 (2022). https://www.rhea-db.org
- **UniProt** — The UniProt Consortium, *Nucleic Acids Res*. 51(D1):D523–D531 (2023). https://www.uniprot.org
- **ChEBI** — Hastings et al., *Nucleic Acids Res*. 44(D1):D1214–D1219 (2016). https://www.ebi.ac.uk/chebi
- **PDB** — Berman et al., *Nucleic Acids Res*. 28(1):235–242 (2000). https://www.rcsb.org
- **ProtT5** — Elnaggar et al., *IEEE TPAMI* (2022). doi:10.1109/TPAMI.2021.3095381
- **RDKit** — https://www.rdkit.org
- **UMAP** — McInnes et al., *arXiv* 1802.03426 (2018)
- **HDBSCAN** — Campello et al., *ACM TKDD* 10(1):5 (2015)
