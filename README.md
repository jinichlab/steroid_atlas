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
│   ├── scripts/                           # audit + fetch + embedding-export scripts
│   │   └── export_literature_recruited_fasta.py   # regenerates FASTA + manifest from proteins.csv
│   ├── supplementary/                     # SI Excel files from cited papers
│   ├── sequences/
│   │   ├── literature_recruited.fasta             # canonical 15-record FASTA
│   │   ├── literature_recruited_manifest.tsv      # per-protein metadata (paper, UMAP, cluster)
│   │   └── guzior/                                # per-paper candidate FASTAs
│   ├── embeddings/
│   │   └── literature_recruited.h5                # ProtT5 1024-d vectors (15 proteins, 72 KB)
│   └── README.md
│
├── analysis/
│   └── polish/                            # finalized paper-figure + CYP-superfamily scripts
│       ├── 15d_cyp_landscape_cluster_ids.py       # Figure 2 Panel A
│       ├── 16_panelB_cyp46_zoom.py                # Figure 2 Panel B
│       ├── 17_panelsCD_cyp46_composition.py       # Figure 2 Panels C + D
│       ├── 18_protein_vs_molecule_coupling.py     # protein × molecule coupling heatmap
│       └── (18–22) CYP family rankings, within-family diversity, umbrella GO coverage
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

**Literature-recruitment pipeline.** All fetch + audit scripts under `literature/scripts/` are runnable. See `literature/README.md` for the end-to-end pipeline (download SI files → parse → fetch sequences → verify → append to atlas → embed → project into UMAP).

The 15 recruited proteins' ProtT5 1024-d per-protein vectors are tracked at `literature/embeddings/literature_recruited.h5` (72 KB) — no need to recompute to use the atlas. To regenerate from scratch:

```bash
# 1. canonical FASTA (reads from data/proteins.csv)
python literature/scripts/export_literature_recruited_fasta.py

# 2. ProtT5 embeddings (requires GPU; prot_t5_xl_half_uniref50-enc)
python prott5_embedder.py \
    -i literature/sequences/literature_recruited.fasta \
    -o literature/embeddings/literature_recruited.h5 \
    --per_protein 1
```

**Paper figures (Figure 2 and supporting analyses).** Finalized scripts live under `analysis/polish/`:

| Script | Produces |
|---|---|
| `15d_cyp_landscape_cluster_ids.py` | Figure 2 Panel A — full atlas, CYP cluster palette |
| `16_panelB_cyp46_zoom.py` | Figure 2 Panel B — CYP46 region zoom |
| `17_panelsCD_cyp46_composition.py` | Figure 2 Panels C + D — taxonomy + CYP family breakdown |
| `18_protein_vs_molecule_coupling.py` | 82×8 protein × molecule coupling heatmap |
| `18_cyp_family_substrate_ranking.py` | CYP family substrate-Tanimoto ranking |
| `19_cyp_subfamily_notable_trends.py` | CYP subfamily trend analysis |
| `20_umbrella_go_coverage_estimate.py` | Umbrella-GO coverage estimate |
| `21_umbrella_go_missing_lookup.py` | Missing umbrella-GO lookup |
| `22_cyp_within_family_diversity.py` | Within-family Tanimoto + Jaccard diversity |

Each script reads `data/proteins.csv` + `data/molecules.csv` and writes PDFs / PNGs / TSVs into `analysis/` on run (outputs are regenerable and not tracked).

## Citation

If you use this atlas, please cite the underlying source papers (see `CITATION.cff`) and this repository.

## License

MIT (see `LICENSE`). Data is aggregated from public databases and open-access publications.

## Contact

Adriana Siordia · <adsiordia@ucsd.edu> · Jinich Lab, UC San Diego
Akshay Uppal · <akuppal@ucsd.edu> · Jinich Lab, UC San Diego
