# Methods — Steroid-interacting protein corpus construction

The atlas draws steroid-interacting proteins from **five complementary evidence sources**, each independently defined. The union is deduplicated by UniProt accession, ProtT5-embedded, 2D-UMAP projected, and deduplicated a second time by exact sequence identity.

## 1. Sterane substructure classifier

The chemistry backbone. Every ChEBI compound and every candidate Rhea reaction is filtered through the same RDKit-based **sterane substructure query**:

- **Sterane query:** SMILES `C1CCC2CCC3C4CCCC4CCC3C2C1` — the 4-fused-ring steroid backbone, all single bonds, no heteroatoms.
- **Molecule pre-processing before matching**: RDKit is used to (i) de-aromatize every atom and bond, (ii) replace all N, O, S atoms with C, and (iii) reduce every bond to a single bond. This is deliberately permissive — it lets the sterane query match steroid derivatives with heteroatom substitutions (bile acids with hydroxyls, steroid sulfates, N-containing steroid alkaloids, and other modified steroid analogs) that would otherwise fail an exact match on the pure-carbon skeleton.
- **Match:** RDKit `HasSubstructMatch(sterane_core, useChirality=False)`.

Applied to every ChEBI participant in the Rhea chemistry catalog (14,173 unique ChEBI SMILES) this classifier identified **875 steroid ChEBIs**.

## 2. Protein evidence sources

### Source A — Steroid-catalyzing enzymes via Rhea

Every Rhea reaction (36,014 total) is parsed at the reactant + product SMILES level and its constituent molecules are canonicalized with RDKit. Reactions with ≥1 molecule matching the steroid ChEBI catalog from §1 are retained.

- **2,210 steroid Rhea reactions** (expanded to **4,420 IDs** across the four reaction directions: master, LR, RL, BI)
- Intersected with the full `rhea2uniprot_sprot` + `rhea2uniprot_trembl` mappings from Rhea's bulk TSVs
- **95,322 candidate UniProt accessions** whose annotated Rhea reactions include at least one steroid participant
- UniProt metadata fetched via the REST bulk-accessions endpoint (batches of 300, retry on 429/5xx)
- Filtered to `annotation_score ≥ 4` (UniProt curated-quality threshold — 4 or 5 stars)

→ **34,131 catalytic entries.**

### Source B — Annotated steroid-binding proteins

UniProt REST query:

```
(keyword:KW-0754 OR go:0005496) AND (annotation_score:4 OR annotation_score:5)
```

- `KW-0754` = UniProt keyword *Steroid-binding*
- `GO:0005496` = Gene Ontology term *steroid binding*

Captures proteins with annotated steroid-binding function that may not catalyze a Rhea-registered steroid reaction (receptors, transporters, carrier proteins).

→ **2,746 binding entries.**

### Source C — ChEBI cross-referenced proteins without a Rhea reaction

UniProt entries whose CC line cross-references a steroid ChEBI from §1 but which are not linked to any Rhea reaction. These are proteins with documented steroid interactions in UniProt (from published biochemistry or structural evidence) whose reactions are not represented in the Rhea catalog.

→ **407 entries.**

### Source D — PDB structural evidence

UniProt entries with a PDB cross-reference where the co-crystallized ligand is a steroid ChEBI from §1. Captures proteins whose steroid interaction is established by a solved structure but not by a Rhea reaction annotation.

→ **65 entries.**

### Source E — Expert manual curation

Hand-selected entries added by expert review of the literature — steroid-related proteins where automated cross-references miss the connection (fusion proteins, orphan enzymes, recently reannotated entries).

→ **42 entries.**

### Union

The five sources are unioned by UniProt accession, with duplicates collapsed. Any entry appearing in Source A takes precedence for source labeling, followed by B → C → D → E for provenance clarity.

→ **37,391 unique UniProt accessions.**

## 3. Sequence embedding

Sequences are extracted to FASTA, whitespace stripped, deduplicated by accession, and embedded with **ProtT5** (`Rostlab/prot_t5_xl_half_uniref50-enc`, encoder-only, half-precision). Per-residue embeddings are mean-pooled to a **1,024-dimensional per-protein vector**.

Compute: single GPU, 20 CPU threads, 50 GB RAM.

## 4. Dimensionality reduction + clustering

- **UMAP** (n_neighbors=15, min_dist=0.1, 2D output) fit on the 37,391 × 1,024 ProtT5 matrix.
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

**35,834 protein sequences** (final after literature audit removed one incorrectly attributed entry).

## Corpus composition by evidence source

| Source | Entries | % of union |
|---|---:|---:|
| A. Steroid-catalyzing (Rhea) | 34,131 | 91.3% |
| B. Steroid-binding (UniProt KW/GO) | 2,746 | 7.3% |
| C. ChEBI cross-reference (no Rhea) | 407 | 1.1% |
| D. PDB structural evidence | 65 | 0.2% |
| E. Expert manual curation | 42 | 0.1% |
| **Union (dedup by accession)** | **37,391** | 100% |
| After exact-sequence dedup | 35,835 | |
| After literature audit | **35,834** | |

## Data sources cited

- **Rhea** — Bansal et al., *Nucleic Acids Res*. 50(D1):D622–D631 (2022). https://www.rhea-db.org
- **UniProt** — The UniProt Consortium, *Nucleic Acids Res*. 51(D1):D523–D531 (2023). https://www.uniprot.org
- **ChEBI** — Hastings et al., *Nucleic Acids Res*. 44(D1):D1214–D1219 (2016). https://www.ebi.ac.uk/chebi
- **PDB** — Berman et al., *Nucleic Acids Res*. 28(1):235–242 (2000). https://www.rcsb.org
- **ProtT5** — Elnaggar et al., *IEEE TPAMI* (2022). doi:10.1109/TPAMI.2021.3095381
- **RDKit** — https://www.rdkit.org
- **UMAP** — McInnes et al., *arXiv* 1802.03426 (2018)
- **HDBSCAN** — Campello et al., *ACM TKDD* 10(1):5 (2015)
