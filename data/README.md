# Data files

Three CSVs form the atlas. Every row is human-readable and traceable back to a paper or public database.

## `proteins.csv` — 35,349 rows × 21 columns

One row per protein sequence. Deduplicated (no identical sequences repeated).

| Column | Description |
|---|---|
| `accession` | Primary identifier (UniProt accession, RefSeq WP_ ID, or locus tag) |
| `entry_name` | UniProt-style entry name (e.g. `CBH_BIFL2`) |
| `protein_names` | Full protein name(s) with EC numbers as UniProt records them |
| `gene_names` | Gene name + locus tag(s) |
| `organism` | Species + strain designation |
| `length_aa` | Sequence length (amino acids) |
| `sequence` | Full protein sequence |
| `ec_numbers` | EC numbers **only** for reactions experimentally confirmed by the cited paper for this specific protein. Blank if no per-protein EC evidence. |
| `rhea_reactions` | Rhea reaction database IDs |
| `interacting_chebi_ids` | `;`-separated ChEBI IDs of steroid substrates/products |
| `interacting_compounds` | `;`-separated human-readable names matching the ChEBI IDs |
| `reaction_descriptions` | Free-text description of the reactions this protein catalyzes |
| `umap_1`, `umap_2` | 2D UMAP embedding of the ProtT5 sequence embedding |
| `cluster` | Cluster ID from HDBSCAN clustering of the UMAP |
| `prott5_cluster` | Original ProtT5-space cluster ID |
| `is_literature_recruited` | 1 if this row was added from a specific paper's audit; 0 if from the original Rhea/UniProt corpus |
| `paper_url` | DOI URL of the paper this entry was recruited from (only if `is_literature_recruited=1`) |
| `annotation` | Human-readable annotation citing the paper |
| `sequence_source` | Where the sequence was fetched from + evidence trail |
| `identifier_type` | UniProt / RefSeq_WP / Locus_Tag / NCBI_Genome |
| `atlas_inclusion_rationale` | **One-line summary of WHY this protein is in the Steroid Atlas.** Three patterns: (1) `Literature-recruited {role} → {DOI}` for literature-recruited rows; (2) `Rhea/UniProt enzyme-substrate match (Rhea: XXXX) | GO: …` for base-corpus enzymes with Rhea substrate links; (3) `GO: steroid binding; cholesterol binding | KW: Steroid-binding` for receptors, transporters, binders (no enzymatic substrate link, but inclusion is justified by an explicit steroid-family GO term or UniProt keyword). Coverage: 14,085/14,087 (99.98%). The 2 Prominin-1 rows (CD133) are flagged as `inclusion trigger unclear` — cholesterol-rich-membrane markers without an explicit steroid GO/KW term; manual review candidate. |

## `molecules.csv` — 864 rows × 14 columns

One row per steroid/bile-acid small molecule.

| Column | Description |
|---|---|
| `compound_name` | Human-readable name (common or IUPAC). The atlas convention is to use the **canonical ChEBI name** when a ChEBI entry exists (e.g. `L-alanocholate`, not `Ala-CA`). |
| `chebi_id` | ChEBI database ID (numeric, no `CHEBI:` prefix). Blank if no ChEBI entry known — never fabricated. |
| `smiles` | SMILES string (canonical or as-published) |
| `umap_1`, `umap_2` | 2D UMAP embedding of the RDKit Morgan fingerprint. For literature-recruited additions made after the original corpus (Gentry, Mohanty, …), coordinates come from Tanimoto-nearest-neighbor projection + small Gaussian jitter (σ=0.15). To be recomputed from scratch in the end-of-audit cleanup pass. |
| `cluster` | Cluster ID from HDBSCAN clustering of the UMAP (inherited from nearest neighbor for post-hoc additions) |
| `is_literature_recruited` | 1 if added from a specific paper; 0 if from original corpus |
| `paper_url` | Paper DOI (only if `is_literature_recruited=1`). May be a pipe-separated list when multiple papers describe the same compound. |
| `interacting_protein_accessions` | List of protein accessions known to act on this compound |
| `chebi_url` | Direct link to the ChEBI page (blank if no ChEBI ID) |
| `synonyms` | Alternative name(s) used in the literature for the same molecule — primarily the paper-side shorthand when the canonical ChEBI name is used in `compound_name`. Example: for `L-alanocholate` the `synonyms` field holds `Ala-CA` (the name used in Rimal 2024 / Guzior 2024 / Lucas 2023). Blank if no alternative is recorded. |
| `evidence_tier` | Paper-specific confidence code for the detection/validation level in the source paper. See [`literature/literature_recruited_molecules.md`](../literature/literature_recruited_molecules.md) for the full per-paper schema. Shorthand: `1`/`2`/`3`/`4-5`/`5` for Gentry 2023; `M-1`/`M-2`/`M-3`/`M-4` for Mohanty 2024. Blank for molecules without a per-paper evidence audit. |
| `evidence_note` | Free-text caveat combining the tier label with sample-specific details (e.g., `"Mohanty: synthesized standard + RT+MS/MS match in BOTH mice AND lion feces"`). Rendered in the UI as a color-coded banner. |
| `gentry_producers` | For Gentry 2023 compounds only — semicolon-separated list of HMP bacterial isolates that produce the compound in vitro (from Gentry 2023 MOESM6). Blank otherwise. |

**Molecule-side provenance:** every literature-recruited molecule row is documented in [`literature/literature_recruited_molecules.md`](../literature/literature_recruited_molecules.md) with per-paper SMILES derivation, mass-validation evidence, and evidence-tier assignments. Ingestion scripts live under the scratchpad for each session.

## `literature_recruited_proteins.csv` — 15 rows × 12 columns

Provenance table: every row is a protein that was added to `proteins.csv` from a specific paper (i.e., `is_literature_recruited=1`).

| Column | Description |
|---|---|
| `Entry` | Same accession as in `proteins.csv` |
| `Paper_ShortRef` | Short citation (e.g. "Rimal 2024") |
| `Paper_Citation` | Full citation string |
| `Paper_DOI` | DOI |
| `Species` | Organism the enzyme is from |
| `Strain` | Strain designation |
| `Identifier_Type` | UniProt / RefSeq_WP / Locus_Tag / NCBI_Genome |
| `Sequence_Length` | aa |
| `Recruitment_Rationale` | Human-readable description of why this specific enzyme was added — what activity did the paper report |
| `Evidence_Level` | Biochemical / genetic / bioinformatic — nature of the evidence in the cited paper |
| `Sequence_Fetch_Source` | Which database the sequence came from + any special notes |
| `In_Atlas_UMAP` | "yes" if the entry has UMAP coordinates in `proteins.csv` |

A human-readable Markdown version is also provided at `../literature/literature_recruited_proteins.md`.

## Curation principles

1. **No fabricated identifiers.** If a compound has no known ChEBI ID, its `chebi_id` field is blank. Same for proteins that lack a UniProt accession — they use their locus tag instead of a made-up UniProt-like ID, and the `identifier_type` field tells you which.
2. **EC numbers reflect what the paper actually tested.** `ec_numbers` in `proteins.csv` lists only ECs directly demonstrated by the cited paper for that specific protein. Sequence-similarity inferences are excluded.
3. **Every literature-recruited row is auditable.** The `annotation` and `sequence_source` fields, plus the provenance table, form a complete chain from paper → gene → sequence → row.

## Optional: RAG AI Agent
Include your OpenAI API in the landing page of the atlas tool to use an AI agent that can answer questions about the data. The agent is trained on all available data on the small molecules visualized. It can provide answers for proteins based on general information available. 
