# Visualizer app

An interactive marimo notebook that renders the atlas as three UMAP views:

- **Protein-centric** — 35,834 protein sequences colored by cluster, with search + detail cards
- **Small-molecule-centric** — 681 steroid/bile-acid molecules with 2D structures + interacting-enzyme lists
- **Natural + Synthetic** — 2,889 natural and synthetic steroids from Rhea

## Run

From the repo root:

```bash
./app/run.sh
```

The app binds to `0.0.0.0:2730` by default. Override with env vars:

```bash
PORT=3000 HOST=127.0.0.1 ./app/run.sh
```

If you're on a remote server, tunnel from your laptop:

```bash
ssh -L 2730:localhost:2730 <you>@<server>
```

Then browse to <http://localhost:2730>.

## Data source

The app reads `../data/*.csv`. Column names on disk are snake_case for downstream users; the app internally re-maps them to its historical column names via a small shim in each load cell.

## What's inside `visualizer.py`

A marimo reactive notebook (~900 lines). Every cell is decorated with `@app.cell` and returns a tuple of values consumed by subsequent cells. The main sections:

1. Imports + RDKit availability check
2-4. Data loading (molecules, natural+synthetic, proteins) with column rename shim
5. Structure cache — precompute 2D depictions of every SMILES via RDKit
6-8. UI controls (view selector, search, method radio)
9-11. Altair UMAP chart with STAR markers for literature-recruited entries
12-14. Selection table
15-19. Detail panel — protein or molecule card with structures + interacting entities

To modify: edit `visualizer.py` in an editor (marimo notebooks are plain Python), or open it in `marimo edit` mode for the reactive-development UI.
