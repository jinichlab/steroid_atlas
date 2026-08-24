# Cluster substrate-Tanimoto significance

## Method

For each cluster we take the observed mean pairwise Tanimoto over its unique substrate set, then draw **2,000 random substrate samples of the same size** from the full molecule pool (n = 592) and compute the same mean each time to build a null distribution. The empirical p-value is `(#null ≥ observed + 1) / (N_perm + 1)` (specialist direction) or `(#null ≤ observed + 1) / (N_perm + 1)` (promiscuous direction). We correct each direction separately by Benjamini-Hochberg FDR.

**Z-score sign convention:** positive = observed > null (specialist), negative = observed < null (promiscuous). Larger magnitudes mean the cluster is farther from what random sampling would produce.

## Summary

- Testable clusters (≥2 substrates): **77 of 82**
- Significantly SPECIALIST (q_specialist < 0.05, z > 0): **72**
- Significantly PROMISCUOUS (q_promiscuous < 0.05, z < 0): **0**

---

## Top 10 significantly SPECIALIST protein families
*(substrates are chemically MORE similar than random pools of the same size)*

| Protein family | # substrates | Observed | Null mean | Z | q |
|---|---|---|---|---|---|
| **Sterol O-acyltransferase 1** | 61 | 0.401 | 0.231 | +19.02 | 7.7e-04 |
| **Acyl-coenzyme A diphosphatase NUDT19** | 17 | 0.574 | 0.231 | +17.70 | 7.7e-04 |
| **Cytochrome P450 85A1** | 29 | 0.453 | 0.231 | +15.66 | 7.7e-04 |
| **Bile salt export pump** | 10 | 0.581 | 0.231 | +12.70 | 7.7e-04 |
| **Cholesterol 24-hydroxylase** | 20 | 0.434 | 0.231 | +11.59 | 7.7e-04 |
| **Sterol carrier protein 2** | 50 | 0.348 | 0.231 | +11.31 | 7.7e-04 |
| **Phosphatidylcholine-sterol acyltransferase** | 38 | 0.367 | 0.231 | +11.23 | 7.7e-04 |
| **Delta(24)-sterol reductase** | 58 | 0.327 | 0.232 | +10.22 | 7.7e-04 |
| **Cholesterol 24-hydroxylase** | 26 | 0.378 | 0.231 | +9.85 | 7.7e-04 |
| **Glucosylceramidase** | 4 | 0.749 | 0.232 | +9.83 | 7.7e-04 |

## Top 10 significantly PROMISCUOUS protein families
*(substrates are chemically MORE diverse than random pools of the same size)*

| Protein family | # substrates | Observed | Null mean | Z | q |
|---|---|---|---|---|---|

---

## Files
- `cluster_tanimoto_significance.tsv` — full 82-row table with z-scores + q-values in both directions.
- `cluster_tanimoto_significance_volcano.{pdf,png}` — Z vs −log₁₀(q).