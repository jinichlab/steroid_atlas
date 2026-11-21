# Intra-cluster Tanimoto — steroid molecule clusters

- Total molecules with valid SMILES: **826**
- Random-pair Tanimoto: mean **0.257**, median **0.214** (n = 19,965 sampled pairs)
- All intra-cluster pairs: mean **0.368**, median **0.314** (n = 44,468 pairs)

## Per-cluster

| Cluster | n | mean | median | std | min | max |
|---|---|---|---|---|---|---|
| C1 | 50 | 0.453 | 0.437 | 0.150 | 0.105 | 0.935 |
| C2 | 66 | 0.647 | 0.716 | 0.212 | 0.141 | 1.000 |
| C3 | 52 | 0.448 | 0.276 | 0.293 | 0.084 | 1.000 |
| C4 | 109 | 0.288 | 0.234 | 0.158 | 0.082 | 1.000 |
| C5 | 181 | 0.279 | 0.246 | 0.122 | 0.090 | 1.000 |
| C6 | 88 | 0.392 | 0.373 | 0.152 | 0.109 | 1.000 |
| C7 | 36 | 0.663 | 0.671 | 0.128 | 0.455 | 1.000 |
| C8 | 46 | 0.515 | 0.470 | 0.171 | 0.247 | 1.000 |
| C9 | 148 | 0.393 | 0.373 | 0.150 | 0.102 | 1.000 |
| C10 | 50 | 0.708 | 0.714 | 0.112 | 0.398 | 1.000 |

**Reading:** higher mean = the cluster is chemically tight (its members share substructure). Big gap between mean and min = the cluster contains some outliers RDKit sees as chemically distant despite sharing the UMAP neighborhood.

## Files
- `intra_cluster_tanimoto.tsv` — per-cluster summary
- `per_molecule_tanimoto.tsv` — per-molecule centrality inside its cluster
- `intra_cluster_tanimoto_distribution.{pdf,png}` — box plot
- `intra_vs_random_tanimoto.{pdf,png}` — KDE overlay vs null