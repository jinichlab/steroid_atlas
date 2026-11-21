# High-dimensional k-means vs. 2D UMAP k-means

Same k as the published clusters; only the space differs. ARI = adjusted Rand index (1 = identical partitions, 0 = chance).

## Steroids (826 molecules, k = 8)

- Space: Morgan fingerprints, radius 2, 2048 bits
- ARI vs 2D UMAP clusters: **0.394**

| High-dim cluster | Size | Main 2D cluster (share) | 2D clusters spanned |
|---:|---:|---|---:|
| 5 | 175 | 9 (42%) | 6 |
| 7 | 155 | 5 (55%) | 7 |
| 4 | 140 | 6 (55%) | 6 |
| 1 | 98 | 2 (62%) | 5 |
| 8 | 91 | 5 (89%) | 3 |
| 6 | 84 | 4 (89%) | 3 |
| 2 | 51 | 1 (96%) | 2 |
| 3 | 32 | 3 (100%) | 1 |

## Proteins (14,087 of 14,087 embedded, k = 82)

- Space: ProtT5 (`prott5.h5`, 1024-d, L2-normalised)
- ARI vs 2D UMAP clusters: **0.647**
- Assigned to nearest centroid (vectors from `literature_recruited.h5`, not used to fit): 12 — C8WL28, C8WL29, C8WL30, C8WL31, MFU7515415, MFU7516964, MFU7517346, dw0526, cp1309, mf2052, mf0519, Q5LF84
- No vector (cluster_hd = -1): 0

| High-dim cluster | Size | Main 2D cluster (share) | 2D clusters spanned |
|---:|---:|---|---:|
| 3 | 838 | 3 (33%) | 4 |
| 47 | 524 | 6 (98%) | 2 |
| 25 | 517 | 20 (85%) | 6 |
| 58 | 496 | 33 (96%) | 3 |
| 32 | 412 | 14 (43%) | 5 |
| 66 | 392 | 28 (68%) | 3 |
| 30 | 378 | 22 (87%) | 4 |
| 48 | 334 | 2 (83%) | 3 |
| 9 | 334 | 26 (60%) | 3 |
| 24 | 320 | 9 (78%) | 6 |

## Settings

- `KMeans(random_state=42, n_init=10)`
- scikit-learn 1.7.1, rdkit 2026.3.4, numpy 2.2.6
