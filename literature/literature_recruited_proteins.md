# Literature-recruited proteins in the Steroid Atlas

This document lists every protein added to `data/proteins.csv` beyond the initial Rhea-based UniProt corpus (rows where `is_literature_recruited=1`). Each entry records:

- the **paper** it came from (DOI),
- the **specific figure or table** that provides the enzyme-activity evidence,
- the **SI file in this repo** (if we have it locally),
- the **protein-ID attribution strategy** — how the paper's description was matched to the identifier we use in the atlas (direct mention vs. strain + annotation inference).

**Total:** 13 entries across 4 papers (Rimal 2024, McCurry 2024, Jacoby 2025, Arp 2025). Guzior 2024 and Lucas 2026 are **molecule-side only** — no protein rows are recruited from them.

> **Previously** this document listed 15 entries across 6 papers. `WP_243289361` (CpBSH/T, *Clostridium perfringens* ATCC 13124), previously attributed to Guzior 2024, was removed on re-audit after the paper attribution did not hold up. The Guzior 2024 paper now contributes **molecules only** (its 9 MCBA products) to the atlas. The Guzior-side scripts and provenance JSONs under `literature/scripts/` and `literature/guzior2024_*.json` are kept as historical audit trail.

Evidence we have walked through end-to-end is marked **✓ walked through**. The others are marked **pending — not yet walked through**, meaning the rows in the atlas were accepted based on the paper's overall claim but the exact figure/panel attribution hasn't been audited yet.

---

## Rimal 2024 ✓ walked through

**Rimal, B. et al. Bile salt hydrolase catalyses formation of amine-conjugated bile acids. *Nature* 626, 859–863 (2024).**
DOI: [10.1038/s41586-023-06990-w](https://doi.org/10.1038/s41586-023-06990-w)
SI file in repo: `literature/supplementary/Rimal_2024_bsh_bfagilis.xlsx`

### Naming convention — paper shorthand ↔ atlas canonical name

Rimal 2024 (and the downstream MCBA literature) names amino-acid-conjugated cholates with a short code like `Ala-CA`. The atlas stores these under their **canonical ChEBI name** in `data/molecules.csv:compound_name` and carries the paper shorthand in the `synonyms` column (see `data/README.md` for the schema). Mapping for this paper's MCBAs:

| Paper shorthand | Atlas `compound_name` | `chebi_id` | `synonyms` |
|---|---|---|---|
| `Ala-CA` | `L-alanocholate` | 229710 | `Ala-CA` |
| `Ser-CA` | `L-serocholate` | 229711 | `Ser-CA` |
| `Phe-CA` | `L-phenylalanocholate` | 229714 | `Phe-CA` |
| `Tyr-CA` | `Tyr-CA` *(no canonical ChEBI entry; see note on row)* | — | — |
| `His-CA` | `L-histidocholate` | 229712 | — *(not a product of Q5LF84 by Fig 2d; tracked separately)* |
| `Thr-CA` | — *(not currently in the atlas as a dedicated row)* | — | — |

A browser search on `Ala-CA` / `Ser-CA` / `Phe-CA` is intended to resolve to the canonical row via the `synonyms` column.

### `P0DXD2` — *Bifidobacterium longum* subsp. *longum* (NCTC 11818)

- **Identifier type:** UniProt (Swiss-Prot curated)
- **Length:** 317 aa
- **Evidence type:** Biochemical (purified enzyme + kinetics + mutagenesis)
- **Paper evidence location:**
  - **Fig 2b** — in-vitro BlBSH assay → 6 MCBA products: Ala-CA, His-CA, Phe-CA, Ser-CA, Thr-CA, Tyr-CA  *(SI tab "Fig 2b")*
  - Fig 2a — BlBSH inhibitor screen  *(SI tab "Fig 2a")*
  - Fig 2c — BlBSH in E. coli (heterologous expression)  *(SI tab "Fig 2c")*
  - Structure-guided mutagenesis (C1A, R17A, N172A, R225A, T2A, W21F, T171A, K264A, D267A)
- **Protein-ID attribution strategy:** Directly named as *B. longum* BSH; Swiss-Prot entry `P0DXD2` (`CBH_BIFL2`). **High confidence — direct match.**
- **In atlas UMAP:** yes

### `Q5LF84` — *Bacteroides fragilis* (NCTC 9343)

- **Identifier type:** UniProt
- **Length:** 331 aa
- **Evidence type:** Genetic (KO / complementation, in vivo)
- **Paper evidence location:**
  - **Fig 2d** — *B. fragilis* WT vs Δbsh-KO vs KI complementation → 4 MCBA products: Ala-CA, Phe-CA, Ser-CA, Tyr-CA from CA and TCA substrates. KO abolishes all four; KI restores all four.  *(SI tab "Fig 2d")*
- **Protein-ID attribution strategy:** The paper uses *B. fragilis* NCTC 9343 with a Δbsh knockout but does **not** name a specific UniProt accession. `Q5LF84` was **tentatively** matched to the paper's BSH by:
  1. strain match (*B. fragilis* NCTC 9343),
  2. locus tag `BF9343_1433` → UniProt `Q5LF84`,
  3. BSH / choloylglycine hydrolase annotation (`EC 2.3.1.-`, GO:0045302).

  **Confidence: strain-level + annotation-level association, not a direct paper mention.** Flagged for independent biochemical confirmation if the paper is re-audited.
- **In atlas UMAP:** yes

---

## Guzior 2024 — molecule-side only

**Guzior, D.V. et al. Bile salt hydrolase acyltransferase activity expands bile acid diversity. *Nature* 626, 852–858 (2024).**
DOI: [10.1038/s41586-024-07017-8](https://doi.org/10.1038/s41586-024-07017-8)
SI file in repo: `literature/supplementary/Guzior2024_SupplTables.xlsx`

This paper contributes 9 MCBA products to `data/molecules.csv` (Ala-, Ser-, Glu-, His-, Ile-, Val-, Lys-, Ala-, Glu-CDCA — see the per-molecule `paper_url` filter). **No protein rows** are recruited from Guzior 2024. A previous version of this doc listed `WP_243289361` (CpBSH/T) under Guzior 2024; that row was removed on re-audit because the paper attribution did not hold up.

---

## McCurry 2024 ✓ walked through

**McCurry, M.D. et al. Gut bacteria convert glucocorticoids into progestins in the presence of hydrogen gas. *Cell* 187, 2952–2968 (2024).**
DOI: [10.1016/j.cell.2024.05.005](https://doi.org/10.1016/j.cell.2024.05.005)
SI file in repo: *not yet downloaded* (main-text PDF at `literature/papers/McCurry_2024_gluccocorticoids.pdf`)

### Cluster-level attribution — convention established by this paper

McCurry 2024 identifies a **four-gene cluster** (`Elen_2451–Elen_2454`) in *Eggerthella lenta* and closely related *Gordonibacter pamelaeae* that is required for 21-dehydroxylation of 3α/3β,5α-THDOC to 3α/3β,5α-THP (allopregnanolone / isoallopregnanolone — the former is brexanolone, the FDA-approved postpartum-depression drug). The paper does **not** biochemically isolate a single catalytic enzyme; the catalytic unit is left to future work (Discussion: *"Future studies under strict anaerobic conditions using the purified enzyme complex are needed to determine the exact mechanism"*).

Accordingly, each of the four cluster members is represented in the atlas with:

- the **same substrate / product list** (the cluster-level reaction),
- the `annotation` field restricted to **what the paper says** — cluster membership + the UniProt / pfam domain annotation for that member. No inference of catalytic role, electron carriage, or dispensability beyond what the paper explicitly demonstrates,
- a `putative_role` banner (new column, rendered at the top of the protein-centric card in the web app) reading **"Putative steroid 21-dehydroxylase (gene-cluster-level functional evidence)"** — to make clear to the reader that individual catalytic activity was not biochemically confirmed,
- `cluster_level_evidence = yes` (new boolean column) for filtering.

This convention is intended as the template for future gene-cluster papers where the catalytic unit is not isolated. Note that **not every pathway paper fits this template** — see Jacoby 2025 / OsrABC below, where each of OsrA/B/C was individually purified and assayed, so the cluster-level banner does **not** apply there.

### Was this the only cluster the paper considered?

Within the paper's scope, **yes — the cluster identity is settled, only the mechanism is deferred.** The authors:

1. Screened a **26-strain *Eggerthellaceae* library** and found Elen_2451–2454 segregates perfectly with 21-dehydroxylation activity: 15/15 producers carry the cluster, 13/13 non-producers lack it, and the one near-exception (*E. lenta* DSM 15644) has a 2-bp frame-shifting truncation in Elen_2453 and does not produce THPs (Fig 4A/B).
2. Heterologously tested **two alternative molybdopterin-oxidoreductase candidate clusters** in *E. lenta* — Elen_3048–3051 and Elen_2917–2921 — and found **neither** supports 21-dehydroxylation (Fig 4G).
3. Showed that other 21-dehydroxylation-inducing bacteria (*E. coli* Nissle, *C. scindens* VPI 12780, *C. perfringens* ATCC 13124, *C. rodentium* ATCC 8090, *C. freundii* ICC 168) act as **H₂-gas inducers of *E. lenta* activity**, not as alternative catalysts (Fig 3M).

No competing cluster or standalone enzyme is left unresolved in this paper.

### Substrates and products (shared across all 4 cluster members)

| Role | Atlas `compound_name` | ChEBI | Paper source |
|---|---|---|---|
| Substrate | `3α,5α-tetrahydrodeoxycorticosterone (THDOC)` | 34461 | Fig 1B/C/D, Fig 2E/F |
| Substrate | `3β,5α-tetrahydrodeoxycorticosterone` | *(ChEBI lookup needed — may use `synonyms`-only convention like Tyr-CA)* | Fig 1B/C, Fig 4F |
| Product | `allopregnanolone (3α,5α-THP, brexanolone)` | 50169 | Fig 1C/D/G, Fig 2E/F |
| Product | `isoallopregnanolone (3β,5α-THP)` | 11909 | Fig 1C, Fig 4F |

> `17650` (cortisol) was previously present in `interacting_chebi_ids` for all 4 cluster members. It was **removed on walk-through** because the paper never tested cortisol — it tested the 3α- and 3β,5α-reduced THDOC derivatives specifically. The paper mentions cortisol only as part of the broader glucocorticoid class.

### `C8WL28` — *Eggerthella lenta* (DSM 2243), Elen_2451

- **Identifier type:** UniProt (from 2009 *E. lenta* genome deposit, DOI 10.4056/sigs.33592)
- **Length:** 281 aa
- **Evidence type:** Genetic (cluster identification + heterologous activity)
- **Paper evidence location:**
  - Fig 4A/B — cluster conservation across the producer strain library
  - Fig 4G — cumate-inducible homologous expression of the full Elen_2451–2454 cluster in *E. lenta* DSM 2243 restores 21-dehydroxylation of 3α5αTHDOC → 3α5α/3β5α-THP
  - Fig 4I — heterologous expression of the Elen_2451–2453 subset in *Gordonibacter urolithinfaciens* (a non-producer) restores activity
- **Annotation (paper-only):** *Member of the Elen_2451–2454 cluster responsible for 21-dehydroxylation of 3α/3β,5α-THDOC to 3α/3β,5α-THP (McCurry et al., Cell 2024). UniProt annotation: FdhD-family protein.*
- **Putative-role banner:** *Putative steroid 21-dehydroxylase (gene-cluster-level functional evidence)*
- **Protein-ID attribution strategy:** Direct — the locus tag `Elen_2451` in the paper corresponds one-to-one with UniProt `C8WL28` from the 2009 *E. lenta* DSM 2243 genome deposit. **High confidence.**
- **In atlas UMAP:** yes

### `C8WL29` — *Eggerthella lenta* (DSM 2243), Elen_2452

- **Identifier type:** UniProt
- **Length:** 203 aa
- **Evidence type:** Genetic (cluster identification + heterologous activity)
- **Paper evidence location:** same as `C8WL28` (Fig 4A/B, 4G, 4I)
- **Annotation (paper-only):** *Member of the Elen_2451–2454 cluster responsible for 21-dehydroxylation of 3α/3β,5α-THDOC to 3α/3β,5α-THP (McCurry et al., Cell 2024). UniProt annotation: 4Fe-4S ferredoxin iron-sulfur binding domain protein.*
- **Putative-role banner:** *Putative steroid 21-dehydroxylase (gene-cluster-level functional evidence)*
- **Protein-ID attribution strategy:** Direct — locus tag `Elen_2452` ↔ UniProt `C8WL29`. **High confidence.**
- **In atlas UMAP:** yes

### `C8WL30` — *Eggerthella lenta* (DSM 2243), Elen_2453

- **Identifier type:** UniProt
- **Length:** 909 aa
- **Evidence type:** Genetic (cluster identification + heterologous activity)
- **Paper evidence location:** same as `C8WL28` (Fig 4A/B, 4G, 4I)
- **Annotation (paper-only):** *Member of the Elen_2451–2454 cluster responsible for 21-dehydroxylation of 3α/3β,5α-THDOC to 3α/3β,5α-THP (McCurry et al., Cell 2024). UniProt annotation: molybdopterin oxidoreductase.*
- **Putative-role banner:** *Putative steroid 21-dehydroxylase (gene-cluster-level functional evidence)*
- **Protein-ID attribution strategy:** Direct — locus tag `Elen_2453` ↔ UniProt `C8WL30`. **High confidence.**
- **In atlas UMAP:** yes

### `C8WL31` — *Eggerthella lenta* (DSM 2243), Elen_2454

- **Identifier type:** UniProt
- **Length:** 307 aa
- **Evidence type:** Genetic (cluster identification; shown by the paper to be absent from the Elen_2451–2453 heterologous construct that reconstitutes activity in *G. urolithinfaciens*, Fig 4I)
- **Paper evidence location:**
  - Fig 4A/B — cluster conservation across the producer strain library
  - Fig 4G — cumate-inducible homologous expression of the full Elen_2451–2454 cluster in *E. lenta* DSM 2243 restores activity
  - Fig 4I — the Elen_2451–2453 subset alone (i.e. **without Elen_2454**) restores activity when expressed in *G. urolithinfaciens*
- **Annotation (paper-only):** *Member of the Elen_2451–2454 cluster responsible for 21-dehydroxylation of 3α/3β,5α-THDOC to 3α/3β,5α-THP (McCurry et al., Cell 2024). UniProt annotation: band 7 / SPFH family protein.*
- **Putative-role banner:** *Putative steroid 21-dehydroxylase (gene-cluster-level functional evidence)*
- **Protein-ID attribution strategy:** Direct — locus tag `Elen_2454` ↔ UniProt `C8WL31`. **High confidence.**
- **In atlas UMAP:** yes

---

## Jacoby 2025 ✓ walked through

**Jacoby, C. et al. Gut bacteria metabolize natural and synthetic steroid hormones via the reductive OsrABC pathway. *Cell Host & Microbe* 33, 1873–1885 (Nov 12, 2025).**
DOI: [10.1016/j.chom.2025.09.014](https://doi.org/10.1016/j.chom.2025.09.014)
Paper PDF in repo: `literature/papers/Jaccoby2024_OsrABC.pdf` *(filename predates the published year; paper is Jacoby 2025)*

### Identifier convention — GenBank protein, not UniProt

The three Osr enzymes are NOT UniProt-curated entries. The paper deposits each as a **GenBank protein accession** (CDS translation, with explicit `.1` version) on the newly-sequenced *C. steroidoreducens* HCS.1 genome (GenBank chromosome `CP170704.1`, plasmids `CP170705.1`/`CP170706.1`):

| Atlas accession | Paper's cited identifier | Locus tag (where given) | Length |
|---|---|---|---|
| `MFU7515415` | **GenBank: MFU7515415.1** (OsrC) | `ACFYH6_05340` | 254 aa |
| `MFU7516964` | **GenBank: MFU7516964.1** (OsrA) | — | 620 aa |
| `MFU7517346` | **GenBank: MFU7517346.1** (OsrB) | — | 645 aa |

Because these are not UniProt entries:

- `identifier_type` in `data/proteins.csv` and `data/literature_recruited_proteins.csv` is set to `GenBank_Protein` (was `NCBI_Genome`).
- The web app's protein-centric view branches on `identifier_type` and renders an **NCBI Protein** link (`https://www.ncbi.nlm.nih.gov/protein/{accession}.1`) instead of the UniProt + AlphaFold buttons — AlphaFold has no entry for GenBank-only accessions.
- No AlphaFold structure is listed; the paper does not report a solved structure for OsrA/OsrC and uses a model built for OsrB (Fig 3B, FadH 2,4-dienoyl-CoA reductase PDB 1PS9) only as a ligand-binding comparison.

### These are individually characterized enzymes — NOT a McCurry-style cluster

Unlike the McCurry 21-dehydroxylation cluster (above), each of OsrA, OsrB, OsrC was heterologously expressed **individually** in *E. coli* and purified. The paper reports **k_cat for each enzyme separately** on multiple substrates (Fig 3E for OsrB, Fig 4B for OsrC, Fig 5D for OsrA). As a result, each atlas row carries its **own** paper-grounded substrate/product list — not a shared pathway list — and no `putative_role` banner is used.

### Pathway (Fig 5E)

```
prednisolone  ──(OsrA, Δ¹-reductase)──▶  cortisol  ──(OsrB, Δ⁴-reductase)──▶  5β-DHC  ──(OsrC, 3-oxo-5β-steroid 3β-reductase)──▶  3β,5β-THC
```

The three enzymes act sequentially. In vivo (p. 1881), a 3α-hydroxysteroid dehydrogenase (not part of OsrABC) can further isomerize 3β,5β-THC to 3α,5β-THC; this is a downstream activity and is NOT attributed to OsrC on its row.

### Was this the only reductive pathway the paper considered?

Within the paper's scope, **yes** — OsrABC is proposed as the single pathway responsible for the observed 3β,5β-tetrahydrosteroid phenotype in *C. steroidoreducens* HCS.1. OsrA and OsrB were identified via RNA-seq (Data S1): both are induced ≥2-fold when cortisol is present (Fig 3A), and heterologous expression of either reduced cortisol (OsrB) or prednisolone (OsrA). For OsrC, the paper screened **all 12 SDR-domain genes in the HCS.1 genome** (Fig 4A) and found that only ACFYH6_05340 (= MFU7515415.1) produced 3β,5β-THC in E. coli; the others produced different stereoisomers or no product. No competing pathway or standalone enzyme is left unassigned.

### `MFU7516964` — *Clostridium steroidoreducens* HCS.1, OsrA

- **Identifier type:** `GenBank_Protein` — GenBank accession `MFU7516964.1`, translated from the HCS.1 chromosome deposit (`CP170704.1`)
- **Length:** 620 aa
- **Protein family:** Fe-S flavoenzyme, FadH-like (Pfam PF00724 + PF07992 — 2,4-dienoyl-CoA reductase family)
- **Evidence type:** Biochemical (purified enzyme + k_cat on multiple substrates)
- **Paper evidence location:**
  - Fig 3A — RNA-seq induction by cortisol (log2 fold change, ‑log10 p-value both significant)
  - Fig 5B — heterologous expression in *C. steroidoreducens* / *E. coli* reduces synthetic glucocorticoids (prednisolone, methylprednisolone, prednisone, dexamethasone, budesonide) to Δ⁴ and/or 3β,5β-tetrahydro products
  - Fig 5C — co-expression of *osrA* + *osrB* reduces prednisolone fully; *osrA* alone converts prednisolone → cortisol (then OsrB takes over)
  - Fig 5D — purified OsrA k_cat on prednisolone (highest), prednisone, methylprednisolone, dexamethasone, budesonide. **No detectable activity on cortisol** (confirms strict specificity for Δ¹)
- **Substrates:** prednisolone (CHEBI:8378), prednisone (CHEBI:8382), methylprednisolone (CHEBI:6888), dexamethasone (CHEBI:41879), budesonide (CHEBI:3207), boldenone (CHEBI:3027)
- **Products:** cortisol (CHEBI:17650), cortisone (CHEBI:16962), testosterone (CHEBI:17347 — from boldenone)
- **Protein-ID attribution strategy:** Direct — the paper cites `OsrA (GenBank: MFU7516964.1)` by name (p. 1876). **High confidence — direct match.**
- **Atlas entry link:** [NCBI Protein — MFU7516964.1](https://www.ncbi.nlm.nih.gov/protein/MFU7516964.1)
- **In atlas UMAP:** yes

### `MFU7517346` — *Clostridium steroidoreducens* HCS.1, OsrB

- **Identifier type:** `GenBank_Protein` — GenBank accession `MFU7517346.1`
- **Length:** 645 aa
- **Protein family:** Fe-S flavoenzyme, FadH-like (Pfam PF00724 + PF07992)
- **Evidence type:** Biochemical (purified enzyme + k_cat + structural homology model)
- **Paper evidence location:**
  - Fig 3A — RNA-seq induction by cortisol (most significant hit; `osrB` is the top-induced gene)
  - Fig 3B — ligand-binding model built on *E. coli* FadH 2,4-dienoyl-CoA reductase scaffold (PDB 1PS9); shows FMN + Fe-S + FAD + NADP⁺ cofactors
  - Fig 3C — percent sequence identity of OsrA/OsrB to characterized Fe-S flavoenzymes (FadH, BaiH, BaiCD) — in the 27–33% range
  - Fig 3D — heterologous expression of *osrB* (but not *osrA*) in *E. coli* converts cortisol → 5β-DHC
  - Fig 3E — purified OsrB k_cat on cortisol, corticosterone, cortisone, testosterone, progesterone. Poor on 4-cholesten-3-one and 4-cholenic-3-one
- **Substrates:** cortisol (CHEBI:17650), cortisone (CHEBI:16962), corticosterone (CHEBI:16827), testosterone (CHEBI:17347), progesterone (CHEBI:17026), 4-cholesten-3-one (CHEBI:16175, poor)
- **Products:** 5β-dihydrocortisol / 5β-DHC (CHEBI:732), 5β-dihydroprogesterone / 5β-DHP (CHEBI:30154), 5β-dihydrotestosterone (CHEBI:2150; atlas label matches)
- **Protein-ID attribution strategy:** Direct — the paper cites `OsrB (GenBank: MFU7517346.1)` by name (p. 1876). **High confidence — direct match.**
- **Atlas entry link:** [NCBI Protein — MFU7517346.1](https://www.ncbi.nlm.nih.gov/protein/MFU7517346.1)
- **In atlas UMAP:** yes

### `MFU7515415` — *Clostridium steroidoreducens* HCS.1, OsrC

- **Identifier type:** `GenBank_Protein` — GenBank accession `MFU7515415.1`; also carries the locus tag `ACFYH6_05340` (used in Fig 4A)
- **Length:** 254 aa
- **Protein family:** Short-chain dehydrogenase/reductase (SDR) superfamily
- **Evidence type:** Biochemical (purified enzyme + k_cat on multiple substrates; whole-genome SDR screen excluding alternatives)
- **Paper evidence location:**
  - Fig 4A — phylogenetic analysis of all 12 SDR-domain genes in HCS.1; heterologous expression in *E. coli*. Only `ACFYH6_12520` and `ACFYH6_05340` (= OsrC) produced 3β,5β-THC; `ACFYH6_12160` and `ACFYH6_06930` produced 3α,5β-THC instead; none of the others showed activity
  - Fig 4B — purified OsrC k_cat on DHT, 5α-DHP, 5α-DHC, DLCA (dehydrolithocholic acid), 5β-DHP, 5β-DHC. Prefers 5β-dihydrosteroids over the bile acid analog DLCA
- **Substrates:** 5β-dihydrocortisol / 5β-DHC (CHEBI:732), 5β-dihydroprogesterone / 5β-DHP (CHEBI:30154), 5α-dihydroprogesterone / 5α-DHP (CHEBI:28952), 5β-dihydrotestosterone / 5β-DHT (CHEBI:2150)
- **Products:** 3β,5β-tetrahydroprogesterone (epipregnanolone, CHEBI:16229); 3β,5β-THC (no canonical atlas ChEBI — the paper identifies the stereoisomer by reference standards, p. 1878)
- **NOT a product in vitro:** 3α,5β-THC / 3α,5β-THP (CHEBI:1712, pregnanolone) — produced by **other** SDRs (ACFYH6_12160, ACFYH6_06930) per Fig 4A. In vivo, pregnanolone-type 3α,5β-THC does appear (p. 1881) but is attributed to a separate downstream 3α-HSDH, not OsrC.
- **Protein-ID attribution strategy:** Direct — the paper renames `ACFYH6_05340` to OsrC and cites `GenBank: MFU7515415.1` (p. 1878). **High confidence — direct match.**
- **Atlas entry link:** [NCBI Protein — MFU7515415.1](https://www.ncbi.nlm.nih.gov/protein/MFU7515415.1)
- **In atlas UMAP:** yes

---

## Arp 2025 ✓ walked through

**Arp, G. et al. Identification of gut bacteria reductases that biotransform steroid hormones. *Nature Communications* 16, 6285 (2025).**
DOI: [10.1038/s41467-025-61425-6](https://doi.org/10.1038/s41467-025-61425-6)
Paper PDF in repo: `literature/papers/Arp2025_reductases.pdf`
SI files in repo: `literature/supplementary/Arp2025_MOESM{1..8}_ESM.{pdf,xlsx}` — Supp Data 1 (5β-RDs, 318 entries), Supp Data 2 (Δ6-RDs, 94 entries), Supp Data 4 (strain list), Supp Data 5 (cloning primers)

### Identifier convention — Arp SI Protein ID, not UniProt or NCBI

The four Arp enzymes are stored in the atlas with the paper's short names as the primary `accession` (e.g. `cp1309`, `mf2052`). These short names are **Arp-internal identifiers derived from the paper's own re-annotation** (Prokka-like gene numbering) and do **NOT** correspond to NCBI RefSeq locus tags.

The authoritative public anchor for each is the **"Protein ID" from Supplementary Data 1 / 2**, which has the format `{GCA}_{numeric suffix}`. All four atlas sequences have been **byte-for-byte verified** against the SI. Mapping:

| Atlas accession | Arp SI Protein ID | Length | SI note | Real NCBI locus (where known) |
|---|---|---|---|---|
| `cp1309` | `GCA_900447045.1_01309` (Supp Data 1) | 645 aa | — | `NCTC11833_01490` / `SQB90530.1` (NADH:flavin oxidoreductase) |
| `dw0526` | `GCA_005121165.3_00526` (Supp Data 1) | 652 aa | — | `EIO64_02605` / `QCI58255.1` (NAD(P)/FAD-dependent oxidoreductase) |
| `mf2052` | `GCA_001312505.1_02052` (Supp Data 1) | 1155 aa | **"Fused with 3β-HSD"** | not uniquely resolvable by length in GCA_001312505.1 |
| `mf0519` | `GCA_000153925.1_00519` (Supp Data 2) | 641 aa | — | 3 same-length candidates in GCA (needs BLAST to disambiguate) |

Important: the numeric suffix in the SI Protein ID (e.g. `01309`) is **not the same as the NCBI locus tag number**. For `cp1309`, the NCBI locus tag is `NCTC11833_01490`, not `NCTC11833_01309`. The SI number is from Arp's internal re-annotation. This is why the atlas's `uniprot_url` (`.../uniprotkb/cp1309/entry`) returns 404 — these short names are not UniProt accessions. The protein-centric view renders the `Locus_Tag` chip ("no external DB — see paper") instead.

### These are a mix of biochemically tested and bioinformatic entries — NOT a cluster, NOT a 3-step pathway

Arp 2025 identifies three enzyme **families** (Δ⁴-5β-reductase, 3β-HSDH/Δ⁵⁻⁴ isomerase, Δ⁶-3-ketosteroid reductase) by phylogenetic reconstruction + mutagenesis + biochemical characterization of representative members. The four atlas rows are a mix:

| Row | Evidence type in Arp 2025 |
|---|---|
| `mf2052` | **Biochemical** — purified + expressed in E. coli (Fig 3b/c/e); pregnenolone → epipregnanolone end-to-end |
| `dw0526` | **Biochemical** (Fig 3e) — but 5β-RD domain alone shows **no activity on pregnenolone**, consistent with needing the upstream 3β-HSDH/I activity (encoded separately by dw0523) to make its real progesterone substrate |
| `cp1309` | **Bioinformatic / phylogenetic** — homolog of biochemically characterized `ci2350`; not individually tested. Whole-culture *C. paraputrificum* 5β-reductase activity previously reported (Bokkenheuser 1983) |
| `mf0519` | **Bioinformatic / phylogenetic** — one of 94 Δ⁶-RD family homologs in Supp Data 2; not individually tested; the Δ⁶-RD prototype is `ci2349` |

### Prior-literature status (fork-checked Nov 2025)

All four genes are **newly identified at the gene/sequence level in Arp 2025**. One reaction (whole-culture *C. paraputrificum* progesterone 5β-reduction) was known since Bokkenheuser et al. 1983 (*J Clin Microbiol* 18:500–504), but the specific gene was not cloned. *Dysosmobacter welbionis* itself was only described in 2019–2020. The 3β-HSDH/Δ⁵⁻⁴ isomerase on steroid hormones is explicitly claimed by Arp as the first bacterial enzyme of this class acting on steroid hormones (vs. bile acids).

### What's conspicuously missing from the atlas

The paper's primary biochemically characterized **prototype** enzymes are **not in the atlas**, because Arp never published their standalone sequences (they appear only in figures as shorthand `ci2350`, `ci2349`, `dw0523`, `fp1651`):

| Missing prototype | Role in paper | Strongest evidence |
|---|---|---|
| `ci2350` (C. innocuum 6_1_30) | **primary 5β-reductase** — ITC Kd on progesterone / cortisol / 6-DHP / hydrocortisone; substrate scope on 4 natural + 8 synthetic progestins | Fig 1, Fig 2 |
| `ci2349` (C. innocuum 6_1_30) | **primary Δ⁶-reductase** — full conversion of 6-DHP → progesterone | Fig 5c/d |
| `dw0523` (D. welbionis J115) | **3β-HSDH/Δ⁵⁻⁴ isomerase** (separated companion to dw0526) | Fig 3e |
| `fp1651` (F. porci LKV-178-WT-2G) | **second fused enzyme** — pregnenolone → epipregnanolone | Fig 3e |

The paper's shorthand names (`ci*`, `dw*`, `fp*`) do NOT resolve to NCBI locus tags one-to-one (e.g. `HMPREF9406_2350` in GCA_000183585.2 is a 169-aa signal peptidase, not the ci2350 5β-reductase). Adding these prototypes would require either (i) direct communication with the Arp authors to get the sequences, or (ii) BLAST + phylogenetic inference against SI candidates — flagged as a future follow-up.

### Substrate / product summary — per-row

| Row | Substrates | Products |
|---|---|---|
| `mf2052` | pregnenolone (CHEBI:16581) | progesterone (17026) → 5β-DHP (30154) → epipregnanolone (16229) |
| `dw0526` | progesterone (17026) *when paired with dw0523* | 5β-DHP (30154) |
| `cp1309` | progesterone (17026), cortisol (17650), cortisone (16962), corticosterone (16827), testosterone (17347) *(inferred by homology; not tested)* | 5β-DHC (732), 5β-DHP (30154) |
| `mf0519` | 6-dehydroprogesterone *(not currently in atlas)*; progesterone (17026) as product class | progesterone (17026) *(inferred; not tested)* |

### `mf2052` — *Mediterraneibacter faecis* JCM 15917 (fused enzyme)

- **Arp SI Protein ID:** `GCA_001312505.1_02052` (Supp Data 1; SI note: *"Fused with 3β-HSD"*)
- **Length:** 1155 aa (two-domain fusion protein)
- **Identifier type:** `Locus_Tag` — UI renders "no external DB — see paper" chip
- **Evidence type:** **Biochemical** — purified + E. coli expression (Fig 3b/c/e)
- **Paper evidence location:**
  - Fig 3b — gene architecture (two-domain MerR-family protein; N-term 5β-RD / C-term 3β-HSDH)
  - Fig 3c — AlphaFold-predicted structure showing both catalytic domains
  - Fig 3e — pregnenolone substrate → 3 sequential products (progesterone, 5β-DHP, epipregnanolone)
- **Attribution strategy:** Direct — SI Protein ID byte-for-byte matches atlas sequence (MD5 confirmed). **High confidence.**
- **In atlas UMAP:** yes

### `dw0526` — *Dysosmobacter welbionis* J115 (5β-RD alone)

- **Arp SI Protein ID:** `GCA_005121165.3_00526` (Supp Data 1)
- **Length:** 652 aa
- **Identifier type:** `Locus_Tag`
- **Evidence type:** **Biochemical** (Fig 3e) — nuance: 5β-RD domain alone shows **no activity on pregnenolone** because it requires progesterone as substrate, which is generated by the companion 3β-HSDH/I enzyme (`dw0523`) encoded separately in the same *D. welbionis* genome
- **Paper evidence location:**
  - Fig 3b — gene architecture (*D. welbionis* has the fusion split into two genes, dw0523 + dw0526, vs. the single fused mf2052/fp1651)
  - Fig 3e — pregnenolone + dw0526 alone → no products (negative result consistent with needing upstream dw0523)
- **Attribution strategy:** Direct — SI Protein ID byte-for-byte matches atlas sequence. **High confidence.**
- **Companion gene not in atlas:** `dw0523` (3β-HSDH/Δ⁵⁻⁴ isomerase) — flagged for future addition
- **In atlas UMAP:** yes

### `cp1309` — *Clostridium paraputrificum* NCTC 11833 (5β-RD homolog)

- **Arp SI Protein ID:** `GCA_900447045.1_01309` (Supp Data 1)
- **Length:** 645 aa
- **Identifier type:** `Locus_Tag`
- **Evidence type:** **Bioinformatic** (phylogenetic member of Δ⁴-5β-RD clade; biochemically characterized prototype is ci2350, not cp1309 directly)
- **Paper evidence location:**
  - Fig 1b — whole-culture *C. paraputrificum* reduces progesterone → 5β-DHP
  - Fig 1c — phylogenetic tree placing cp1309 in same OYE clade as ci2350
  - Fig 3b — gene architecture (AcrR-family-adjacent)
- **Attribution strategy:** Direct — SI Protein ID byte-for-byte matches atlas sequence. **High confidence on sequence identity; medium confidence on individual activity (never tested).**
- **Prior literature:** Whole-culture *C. paraputrificum* 5β-reductase reported by Bokkenheuser et al. 1983 (*J Clin Microbiol* 18:500–504, doi:10.1128/jcm.18.3.500-504.1983). The SPECIFIC gene (cp1309) was first identified in Arp 2025.
- **In atlas UMAP:** yes

### `mf0519` — *Mediterraneibacter torques* ATCC 27756 (Δ⁶-RD homolog) — SPECIES CORRECTED

- **Arp SI Protein ID:** `GCA_000153925.1_00519` (Supp Data 2)
- **Length:** 641 aa
- **Identifier type:** `Locus_Tag`
- **Evidence type:** **Bioinformatic** (phylogenetic member of Δ⁶-RD steroid-hormone clade c2; biochemically characterized prototype is ci2349)
- **Paper evidence location:**
  - Fig 5a — phylogenetic tree placing mf0519 in Δ⁶-RD(c2) clade (steroid hormone, baiH-related)
  - Supp Data 2 — 94 Δ⁶-RD candidates; mf0519 is one of them
- **Attribution strategy:** Direct — SI Protein ID byte-for-byte matches atlas sequence. **High confidence on sequence identity; LOW confidence on individual activity (never tested).**
- **⚠ SPECIES CORRECTION:** The atlas previously labeled this as *Mediterraneibacter faecis* ATCC BAA-2716 — this was **wrong**. The assembly `GCA_000153925.1` is `[Ruminococcus] torques ATCC 27756`, reclassified as `Mediterraneibacter torques` by Togo et al. 2018. Arp SI GTDB taxonomy confirms `s__Mediterraneibacter torques`. The atlas's `mf` prefix was previously mis-expanded as *M. faecis*; both species are *Mediterraneibacter* so share the "mf" two-letter code in Arp's shorthand, but they are different species.
- **In atlas UMAP:** yes

---

## Lucas 2026 — molecule-side only

**Lucas, L. N. et al. Investigation of bile salt hydrolase activity in human gut bacteria reveals production of conjugated secondary bile acids. *Nature Communications* 17, 3077 (2026).**
DOI: [10.1038/s41467-026-68556-4](https://doi.org/10.1038/s41467-026-68556-4)
Paper PDF in repo: `literature/papers/Lucas2026_BSH.pdf`
SI files in repo: `literature/supplementary/Lucas2026_MOESM{1,3,4}_ESM.pdf`

**No protein rows are recruited from Lucas 2026.** The paper is a whole-strain survey (77 human gut bacterial strains × 5 conjugated bile acids via LC-MS/MS) with one genetic KO (B. thetaiotaomicron Δhsd, Fig 6/7). BSH identities across the 77 strains are bioinformatic HMM calls; activity is measured at the whole-culture level without individual enzyme purification, cloning, or kinetics. This evidence depth is below the atlas's bar for literature-recruited protein rows.

The single genetic experiment targets **B. theta HsdA / BT_2086** (a canonical UniProt Swiss-Prot 7α-HSD) and shows the *hsdA* gene is required for in vivo G-7-oxoLCA production from GCDCA (Δhsd knockout experiment, Fig 6/7). Lucas's claim that the enzyme acts *directly* on conjugated primary BAs is an inference from the genetic knockout — **no purified-enzyme assay was performed**. The paper's enzyme-level claim is a proposed refinement of previously known activity (Sherrod & Hylemon 1977, ref 47; McMillan 2023, ref 52), not a new-sequence discovery.

Lucas 2026's contribution to the atlas is **molecule-side only** — the novel conjugated secondary bile acids (CSBAs) and related -oxo compounds detected in the 77-strain survey and the mouse colonization experiment. See `data/molecules.csv` for the specific entries (filter by `paper_url` containing `10.1038/s41467-026-68556-4`).

> **Q8A6H3 history**: Previously this document listed Q8A6H3 as a literature-recruited protein under an "Attribution flagged for re-audit" section; in a pre-session re-audit it was moved back to base Rhea/UniProt corpus (`is_literature_recruited=0`). In a later pass during this session, Q8A6H3 was **removed from the atlas entirely** because its `protein_names` field had been annotated with an unsupported extension (*"— directly transforms conjugated bile acids"*) based on Lucas's Δhsd knockout — a claim that fails the atlas's purified-enzyme bar. Total atlas protein count dropped by 1 (14088 → 14087). The Glycine-7-oxoLCA `evidence_note` in `data/molecules.csv` keeps the generic descriptor "B. theta HsdA / 7α-HSD" without naming any specific UniProt accession.

---

## Atlas inclusion rationale — every protein row's "why"

The atlas carries 14,087 protein rows total. Every row has a one-line explanation of why it is in the Steroid Atlas, stored in the `atlas_inclusion_rationale` column of `data/proteins.csv` (and surfaced to the web UI). Three patterns cover the entire set:

| # | Pattern | Count | Example |
|---|---|---:|---|
| 1 | **Literature-recruited** — this document's 13 rows, each with paper DOI + mechanism-level annotation | 13 | `Literature-recruited OsrA (GenBank: MFU7516964.1) — oxosteroid reductase A, Fe-S flavoenzyme → 10.1016/j.chom.2025.09.014` |
| 2 | **Rhea/UniProt enzyme–substrate match** — the base corpus; a Rhea reaction annotates the protein with a specific steroid/sterol/bile-acid substrate | 12,812 | `Rhea/UniProt enzyme-substrate match (Rhea: 11064, 11065) | GO: cholesterol binding; steroid 11-alpha-monooxygenase activity` |
| 3 | **Receptor / binder / transporter** — no enzymatic substrate link, but the row has at least one steroid-family GO term or UniProt keyword that justifies inclusion | 1,260 | `GO: nuclear glucocorticoid receptor activity; steroid hormone binding | KW: Steroid-binding` |

**Coverage: 14,085/14,087 (99.98%).** The only 2 rows flagged as `inclusion trigger unclear` are both **Prominin-1 (CD133)** — tetraspanin membrane markers known to localize to cholesterol-rich membrane microdomains but lacking an explicit `steroid`, `sterol`, or `cholester…` term in their current GO/KW annotation. Flagged for manual review (not deleted — their atlas placement predates this audit and may reflect a legitimate but under-annotated inclusion trigger such as `cholesterol-binding membrane domain` or a specific PubMed-citation link).

### Why this rationale column exists

Before this audit, 1,262 "receptor/binder" rows had empty `interacting_compounds` and empty `paper_url` — which looked like orphan data at first glance. The audit confirmed they were classified as steroid-relevant via UniProt's own GO/KW annotation at ingestion time, but the atlas didn't surface that classification trail. The `atlas_inclusion_rationale` column makes the trail explicit:
- UI users clicking a protein card can immediately see the inclusion basis
- Automated tools can distinguish enzymatic rows (Rhea-based) from receptor/binder rows (GO/KW-based) from literature-recruited rows (paper-attributed)
- Future audits can re-run the rationale computation to catch any new rows lacking a clear trigger

### Expanded steroid-family vocabulary

The rationale builder recognizes these term stems (case-insensitive substring match on `go_labels` and `keyword_labels`):

- Core: `steroid`, `sterol`, `cholester`, `oxysterol`
- Mammalian hormones: `androgen`, `estrogen`, `progester`, `testost`, `dihydrotest`, `corticoster`, `glucocorticoid`, `mineralocorticoid`
- Bile acids: `bile acid`, `bile salt`, `cholate`, `deoxycholate`, `lithocholate`, `chenodeoxycholate`, `ursodeoxycholate`, `hyodeoxycholate`, `taurocholate`, `glycocholate`, `muricholic`, `hyocholic`
- Vitamin D family: `vitamin d`, `calcitri`, `secosteroid`, `1,25-dihydroxyvitamin`
- Invertebrate: `ecdys`, `ecdysteroid`, `dafach`
- Plant/fungal: `brassinosteroid`, `ergost`, `phytosterol`
- Pregnane family: `pregnan`, `pregnen`, `pregnolone`
- Specific enzyme classes: `hydroxysteroid dehydrogenase`, `3-ketosteroid`, `cyp7`, `cyp27`, `cyp11`, `cyp19`, `cyp17`, `cyp21`, `cyp51`, `cyp46`, `aromatase`

When additional steroid-family vocabulary surfaces, extend this list and re-run `scratchpad/add_inclusion_rationale.py`.

## One-file sequence bundle

All 13 sequences are in one place for re-embedding or BLASTing:

- FASTA: `literature/sequences/literature_recruited.fasta` (13 entries after Q8A6H3 un-recruitment and subsequent removal from the atlas; header format `>accession|entry_name|gene|organism|len=N`)
- Per-entry manifest: `literature/sequences/literature_recruited_manifest.tsv`
- Pre-computed ProtT5 embeddings (1024-dim per protein): `literature/embeddings/literature_recruited.h5`

To regenerate the FASTA from `data/proteins.csv`:

```bash
python literature/scripts/export_literature_recruited_fasta.py
```
