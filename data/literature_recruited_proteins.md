# Literature-recruited proteins in the Steroid Atlas

This document lists every protein added to `data/proteins.csv` beyond the initial Rhea-based UniProt corpus (rows where `is_literature_recruited=1`). Each entry records:

- the **paper** it came from (DOI),
- the **specific figure or table** that provides the enzyme-activity evidence,
- the **SI file in this repo** (if we have it locally),
- the **protein-ID attribution strategy** — how the paper's description was matched to the identifier we use in the atlas (direct mention vs. strain + annotation inference).

**Total:** 15 entries across 6 papers.

Evidence we have walked through end-to-end is marked **✓ walked through**. The others are marked **pending — not yet walked through**, meaning the rows in the atlas were accepted based on the paper's overall claim but the exact figure/panel attribution hasn't been audited yet.

---

## Rimal 2024 ✓ walked through

**Rimal, B. et al. Bile salt hydrolase catalyses formation of amine-conjugated bile acids. *Nature* 626, 859–863 (2024).**
DOI: [10.1038/s41586-023-06990-w](https://doi.org/10.1038/s41586-023-06990-w)
SI file in repo: `literature/supplementary/Rimal_2024_bsh_bfagilis.xlsx`

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

## Guzior 2024 — pending re-audit of SI tables

**Guzior, D.V. et al. Bile salt hydrolase acyltransferase activity expands bile acid diversity. *Nature* 626, 852–858 (2024).**
DOI: [10.1038/s41586-024-07017-8](https://doi.org/10.1038/s41586-024-07017-8)
SI file in repo: `literature/supplementary/Guzior2024_SupplTables.xlsx`

### `WP_243289361` — *Clostridium perfringens* (ATCC 13124)

- **Identifier type:** RefSeq_WP
- **Length:** 329 aa
- **Evidence type:** Biochemical (kinetics + mutagenesis + structure)
- **Paper evidence location:** Main text + Supp Table 1 (CpBSH/T). Structure: PDB 2BJG. *(exact SI panel/tab pending — not yet walked through)*
- **Why recruited:** CpBSH/T — the only Guzior 2024 protein with direct biochemical evidence. Bifunctional BSH (EC 3.5.1.24, pH 3–7) AND novel amine N-acyltransferase (EC 2.3.1.-, pH optimum 5.3). Peak acyl transfer = 7% of hydrolysis rate. Substrate scope: 16/20 amino acids from TCA, 12/20 from CA; Pro and Asp never conjugated. Mutants: C2A abolishes both activities; N82Y preserves activity but shifts specificity. Other Guzior 2024 Table 1 entries were deliberately NOT added — no protein-level biochemical evidence.
- **Protein-ID attribution strategy:** pending — not yet walked through
- **In atlas UMAP:** yes

---

## McCurry 2024 — pending re-audit

**McCurry, M.D. et al. Gut bacteria convert glucocorticoids into progestins in the presence of hydrogen gas. *Cell* 187, 2949–2963.e19 (2024).**
DOI: [10.1016/j.cell.2024.05.005](https://doi.org/10.1016/j.cell.2024.05.005)
SI file in repo: *not yet downloaded*

### `C8WL28` — *Eggerthella lenta* (DSM 2243)

- **Identifier type:** UniProt (2009 *E. lenta* genome deposit, DOI 10.4056/sigs.33592)
- **Length:** 281 aa
- **Evidence type:** Genetic (cluster identification + heterologous activity)
- **Why recruited:** Elen_2451 — one of four proteins in the *elen_2451–2454* cluster responsible for 21-dehydroxylation of glucocorticoids. Cluster identified via ATc-inducible expression + LC-MS activity assay in *E. lenta*.
- **Protein-ID attribution strategy:** pending — not yet walked through
- **In atlas UMAP:** yes

### `C8WL29` — *Eggerthella lenta* (DSM 2243)

- **Identifier type:** UniProt
- **Length:** 203 aa
- **Evidence type:** Genetic (cluster identification)
- **Why recruited:** Elen_2452 — ferredoxin (4Fe-4S) component of the *elen_2451–2454* cluster.
- **Protein-ID attribution strategy:** pending — not yet walked through
- **In atlas UMAP:** yes

### `C8WL30` — *Eggerthella lenta* (DSM 2243)

- **Identifier type:** UniProt
- **Length:** 909 aa
- **Evidence type:** Genetic (cluster identification)
- **Why recruited:** Elen_2453 — molybdopterin oxidoreductase/dehydrogenase in the *elen_2451–2454* cluster.
- **Protein-ID attribution strategy:** pending — not yet walked through
- **In atlas UMAP:** yes

### `C8WL31` — *Eggerthella lenta* (DSM 2243)

- **Identifier type:** UniProt
- **Length:** 307 aa
- **Evidence type:** Genetic (cluster identification)
- **Why recruited:** Elen_2454 — SPFH / band-7 / stomatin-prohibitin family membrane protein in the *elen_2451–2454* cluster.
- **Protein-ID attribution strategy:** pending — not yet walked through
- **In atlas UMAP:** yes

---

## Jacoby 2025 — pending re-audit

**Jacoby, C. et al. *Cell Host & Microbe* (2025).**
DOI: [10.1016/j.chom.2025.09.014](https://doi.org/10.1016/j.chom.2025.09.014)
SI file in repo: *not yet downloaded*

### `MFU7515415` — *Clostridium steroidoreducens* (HCS.1)

- **Identifier type:** NCBI_Genome
- **Length:** 254 aa
- **Evidence type:** Genetic (pathway identification)
- **Why recruited:** OsrA — one of three genes (OsrABC) in the oxidative steroid reduction pathway identified in *C. steroidoreducens* HCS.1.
- **Protein-ID attribution strategy:** pending — not yet walked through
- **In atlas UMAP:** yes

### `MFU7516964` — *Clostridium steroidoreducens* (HCS.1)

- **Identifier type:** NCBI_Genome
- **Length:** 620 aa
- **Evidence type:** Genetic (pathway identification)
- **Why recruited:** OsrB — component of the OsrABC pathway from *C. steroidoreducens* HCS.1.
- **Protein-ID attribution strategy:** pending — not yet walked through
- **In atlas UMAP:** yes

### `MFU7517346` — *Clostridium steroidoreducens* (HCS.1)

- **Identifier type:** NCBI_Genome
- **Length:** 645 aa
- **Evidence type:** Genetic (pathway identification)
- **Why recruited:** OsrC — component of the OsrABC pathway from *C. steroidoreducens* HCS.1.
- **Protein-ID attribution strategy:** pending — not yet walked through
- **In atlas UMAP:** yes

---

## Arp 2025 — pending re-audit

**Arp, G. et al. *Nature Communications* (2025).**
DOI: [10.1038/s41467-025-61425-6](https://doi.org/10.1038/s41467-025-61425-6)
SI file in repo: *not yet downloaded (Supp Data 1 and Supp Data 2 referenced)*

### `dw0526` — *Dysosmobacter welbionis* (J115)

- **Identifier type:** Locus_Tag
- **Length:** 652 aa
- **Evidence type:** Biochemical (Arp 2025 Supp Data 1)
- **Sequence source:** GCA_005121165.3 genome
- **Why recruited:** Δ4-3-ketosteroid 5β-reductase. Bile-acid-ring reduction activity characterized in Supp Data 1.
- **Protein-ID attribution strategy:** pending — not yet walked through
- **In atlas UMAP:** yes

### `cp1309` — *Clostridium paraputrificum* (NCTC 11833)

- **Identifier type:** Locus_Tag
- **Length:** 645 aa
- **Evidence type:** Biochemical (Arp 2025 Supp Data 1)
- **Sequence source:** GCA_900447045.1 genome
- **Why recruited:** Δ4-3-ketosteroid 5β-reductase.
- **Protein-ID attribution strategy:** pending — not yet walked through
- **In atlas UMAP:** yes

### `mf2052` — *Mediterraneibacter faecis* (JMC 15917)

- **Identifier type:** Locus_Tag
- **Length:** 1155 aa
- **Evidence type:** Biochemical (Arp 2025 Supp Data 2)
- **Sequence source:** GCA_001312505.1 genome
- **Why recruited:** Fused 3β-HSDH / Δ5-4 isomerase (1155 aa fused enzyme).
- **Protein-ID attribution strategy:** pending — not yet walked through
- **In atlas UMAP:** yes

### `mf0519` — *Mediterraneibacter faecis* (ATCC BAA-2716)

- **Identifier type:** Locus_Tag
- **Length:** 641 aa
- **Evidence type:** Biochemical (Arp 2025 Supp Data 2)
- **Sequence source:** GCA_000153925.1 genome
- **Why recruited:** Δ6-3-ketosteroid reductase.
- **Protein-ID attribution strategy:** pending — not yet walked through
- **In atlas UMAP:** yes

---

## Attribution flagged for re-audit

**Paper attribution uncertain — re-audit before relying on this row.**
DOI (currently in `paper_url`): [10.1038/s41467-026-68556-4](https://doi.org/10.1038/s41467-026-68556-4)
SI file in repo: *not yet downloaded*

### `Q8A6H3` — *Bacteroides thetaiotaomicron* (VPI-5482)

- **Identifier type:** UniProt
- **Length:** 259 aa
- **Evidence type:** Bioinformatic (needs re-audit)
- **Why recruited:** *B. theta* BSH — canonical gut Bacteroides bile salt hydrolase.
- **Protein-ID attribution strategy:** pending — **paper attribution itself flagged**. The DOI currently in `paper_url` is `10.1038/s41467-026-68556-4` (Nat Commun 2026) but the entry label in prior docs read "Yao 2018 / see note" — this mismatch must be resolved before this row is cited.
- **In atlas UMAP:** yes

---

## One-file sequence bundle

All 15 sequences are in one place for re-embedding or BLASTing:

- FASTA: `literature/sequences/literature_recruited.fasta` (15 entries, header format `>accession|entry_name|gene|organism|len=N`)
- Per-entry manifest: `literature/sequences/literature_recruited_manifest.tsv`
- Pre-computed ProtT5 embeddings (1024-dim per protein): `literature/embeddings/literature_recruited.h5`

To regenerate the FASTA from `data/proteins.csv`:

```bash
python literature/scripts/export_literature_recruited_fasta.py
```
