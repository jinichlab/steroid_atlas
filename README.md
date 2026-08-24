# Steroid Atlas

An interactive UMAP atlas of steroid- and bile-acid-metabolizing enzymes and their small-molecule substrates, curated from public databases (UniProt, Rhea, ChEBI, NCBI RefSeq) and augmented with hand-audited literature recruitments from recent papers.

- **14,089 protein sequences** — Swiss-Prot annotation-score 5, active in UniProtKB, steroid-associated by either Rhea catalysis or GO / ChEBI / keyword annotation. Each entry carries ProtT5 embeddings, 2D UMAP coordinates, k-means clusters (k = 82, chosen by silhouette + composite score in [50, 95]), and PubMed provenance.
- **677 steroid small molecules** with SMILES, verified ChEBI IDs, and RDKit-drawn structures
- **2,889 natural + synthetic steroid entries** for chemical-space comparison
- **15 literature-recruited proteins** from 5 audited papers, each with full evidence-based provenance
- **Two visualizers under one roof:**
  - **`marimo/`** — Python + marimo + Altair + RDKit notebook, the original interactive tool: clickable steroid catalogue, per-tile inline protein panels, rich detail cards with GO / keyword chips, PubMed refs, EC / Rhea reactions, UniProt + AlphaFold links.
  - **`web/`** — Next.js 14 + deck.gl web app, a redesigned frontend consuming the same data: WebGL scatter, zoom-conditional cluster labels, marimo-style hover tooltips, multi-column search (name / gene / GO / EC / ChEBI / …), pan-or-select toggle with persistent lasso, structure-thumbnail catalogue tiles, multi-select detail cards.
- **AI chatbot** with retrieval over PubChem, Europe PMC, and ChEBI records for compound / protein questions  
  *(requires a personal OpenAI API key)*

## Quick start

Clone once, run either visualizer:
```bash
git clone https://github.com/jinichlab/steroid_atlas.git
cd steroid_atlas
```

**Marimo (Python):**
```bash
pip install -r requirements.txt
./marimo/run.sh
# open http://localhost:2730
```

**Web (Next.js):**
```bash
cd web
npm install
npm run dev -- --port 3000
# open http://localhost:3000
```

If on a remote server, tunnel first — e.g. `ssh -N -L 2730:localhost:2730 -L 3000:localhost:3000 <you>@<server>`.

## Repo layout

```
steroid-atlas/
├── README.md                              # you are here
├── LICENSE                                # MIT
├── CITATION.cff                           # how to cite this atlas
├── requirements.txt                       # Python deps
├── .gitignore
│
├── data/                                  # canonical data — all CSVs
│   ├── proteins.csv                       # 14,089 rows × 33 cols
│   ├── molecules.csv                      # 677 rows × 10 cols
│   ├── natural_synthetic_steroids.csv     # 2,889 rows × 7 cols
│   ├── literature_recruited_proteins.csv  # 15 rows × 12 cols — provenance
│   ├── literature_recruited_proteins.md   # human-readable version
│   └── README.md                          # data dictionary
│
├── marimo/                                # original marimo visualizer (Python)
│   ├── visualizer.py                      # marimo notebook
│   ├── run.sh                             # launcher script
│   └── README.md
│
├── web/                                   # Next.js 14 + deck.gl web app
│   ├── src/                               # App Router pages + client components
│   ├── public/atlas/                      # pre-built JSON + structure PNGs
│   ├── scripts/                           # build_atlas_data.py, build_structures.py
│   ├── package.json
│   └── README.md
│
├── literature/                            # provenance materials
│   ├── scripts/                           # audit + fetch scripts (reproducible)
│   ├── supplementary/                     # SI Excel files from cited papers
│   ├── sequences/                         # fetched FASTAs organized per paper
│   ├── embeddings/                        # ProtT5 embeddings + projection reports
│   └── README.md
│
└── docs/
    ├── methodology.md                     # how the atlas was built
    └── changelog.md
```

## Data — what's in it

Three CSVs in `data/`. Every row is human-readable and traceable back to a paper or public database. See `data/README.md` for the full data dictionary.

**`proteins.csv`** — one row per protein. Columns include `accession`, `organism`, `sequence`, `ec_numbers` (only those experimentally confirmed by the cited paper), `interacting_chebi_ids`, `umap_1`/`umap_2`, and `is_literature_recruited`.

**`molecules.csv`** — one row per steroid/bile-acid. Columns include `compound_name`, `chebi_id` (numeric, no prefix; blank if unknown — **never fabricated**), `smiles`, `umap_1`/`umap_2`.

**`literature_recruited_proteins.csv`** — provenance table for the 15 entries added from specific papers. Every row names the paper, the specific evidence in that paper, and where the sequence was fetched from.

## Curation principles

1. **No fabricated identifiers.** ChEBI IDs and UniProt accessions are only listed when verified against the source database. Compounds without a known ChEBI ID have a blank field.
2. **EC numbers reflect only what the paper tested.** For literature-recruited entries, `ec_numbers` in `proteins.csv` lists only ECs directly demonstrated by the cited paper for that specific protein. Sequence-similarity inferences are excluded.
3. **Every literature recruitment is auditable.** The `annotation`, `sequence_source`, and provenance-table entries form a chain from paper → gene → sequence → row.

## Papers audited and recruited from

| Paper | Enzymes added | Focus |
|---|---|---|
| [Rimal 2024 (Nature)](https://doi.org/10.1038/s41586-023-06990-w) | 2 | Bile salt hydrolase amine N-acyltransferase activity forms BBAAs |
| [Guzior 2024 (Nature)](https://doi.org/10.1038/s41586-024-07017-8) | 1 | CpBSH/T bifunctional MCBA acyltransferase kinetics |
| [McCurry 2024 (Cell)](https://doi.org/10.1016/j.cell.2024.05.005) | 4 | *Eggerthella lenta* 21-dehydroxylation cluster |
| [Jacoby 2025 (Cell Host Microbe)](https://doi.org/10.1016/j.chom.2025.09.014) | 3 | OsrABC oxidative steroid reduction pathway |
| [Arp 2025 (Nat Commun)](https://doi.org/10.1038/s41467-025-61425-6) | 4 | Bile-acid ring reductases (5β, 3β-HSDH/Δ5-4, Δ6) |
| Bacteroides theta background | 1 | B. theta BSH (Q8A6H3) |

Full provenance in `data/literature_recruited_proteins.csv` and human-readable `.md`.

## Reproducibility

All fetch + audit scripts under `literature/scripts/` are runnable. See `literature/README.md` for the end-to-end pipeline (download SI files → parse → fetch sequences → verify → append to atlas → embed → project into UMAP).

## Citation

If you use this atlas, please cite the underlying source papers (see `CITATION.cff`) and this repository.

## License

MIT (see `LICENSE`). Data is aggregated from public databases and open-access publications.

## Contact

Adriana Siordia · <adsiordia@ucsd.edu> · Jinich Lab, UC San Diego
Akshay Uppal · <akuppal@ucsd.edu> · Jinich Lab, UC San Diego
