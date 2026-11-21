# Literature-recruited molecules in the Steroid Atlas

This document lists every molecule added to `data/molecules.csv` beyond the initial Rhea/UniProt/HMDB-based molecule corpus. Each entry records:

- the **paper** it came from (DOI),
- the **specific figure or table** that provides the structural/detection evidence,
- the **SI file in this repo** (if we have it locally),
- the **chemistry provenance** — how the SMILES was derived and validated against the paper's reported m/z,
- the **evidence tier** — explicit confidence level based on the source paper's own hierarchy (chemistry only → cross-cohort validated).

This doc is the molecule-side counterpart to `literature_recruited_proteins.md`. Both use the same "✓ walked through" marker when the per-paper audit has been completed.

**Current totals (data/molecules.csv):** 864 molecule rows (was 675 before literature audits). 268 rows are flagged `is_literature_recruited=1`.

---

## Paper inventory (molecule side)

| Paper | Mol. count | Status | Notes |
|---|---:|---|---|
| Rimal 2024 (Nature) | 10 | audit via protein-side walk-through | MCBAs (bile amidates on CA core) |
| Guzior 2024 (Nature) | 9 | audit via protein-side walk-through | MCBAs on CDCA core |
| Arp 2025 (Nat Commun) | 15 | audit via protein-side walk-through | Steroid hormones + synthetic progestins |
| **Lucas 2026 (Nat Commun)** | **4 + 2 added = 6** | **✓ walked through** | CSBAs (conjugated secondary bile acids); added G-3-oxoCA + G-3-oxoCDCA from C. intestinalis cultures |
| McCurry 2024 (Cell) | — | — | 21-dehydroxylation substrates/products |
| **Gentry 2023 (Nature)** | **168** | **✓ walked through** | Bile amidates — 21 AAs × 8 BA cores combinatorial library |
| **Mohanty 2024 (Cell)** | **17** | **✓ walked through** | Polyamine bile amidates (Cholyl- + Chenodeoxycholyl-) |
| **Lamichhane 2025 (Nat Commun)** | **27 + 5 added** | **✓ walked through** | Early-life + islet autoimmunity MCBAs; added 5 Cys-conjugates missing from Gentry's library |
| Lamichhane preprint (bioRxiv 2025) | 9 | pending | Preprint of the above |
| **Xue 2025 (Gut Microbes)** | **22 enriched (0 new)** | **✓ walked through** | OSA / atherosclerosis mouse study — annotation-only pass; all compounds already present |
| **Aatsinki 2025 preprint (bioRxiv)** | **34 enriched (0 new)** | **✓ walked through** | Early-life stool microbiome — annotation-only pass; all MCBAs named in Results overlap Gentry/Lamichhane |
| **Agongo 2026 preprint** | **18 added** | **✓ walked through** | **Alkamines** — NEW reductive-amination chemistry at C3 (not C24 amidate); 8 BA-alkamines + 2 drug-BA (5-ASA) + 8 cholestan-3-one alkamines |
| **Reilly 2026 (Cell Systems)** | **45 enriched (0 new)** | **✓ walked through** | BA rhythmicity across 14 mouse organs — annotation-only pass |
| Kvitne 2025 · Kostenko 2024 · Mohanty 2025 preprint | 1 each | pending | Single-compound cites |

---

## Evidence-tier convention (schema shared across all paper walk-throughs)

Every literature-recruited molecule row carries three fields:

| Field | Values | Meaning |
|---|---|---|
| `evidence_tier` | paper-specific code (see each walk-through) | A short code for the detection/validation level |
| `evidence_note` | free text | Longer caveat — the tier label plus sample-specific details |
| `gentry_producers` | ;-separated species list, or blank | Only populated for Gentry 2023 Tier 5; names bacterial strains that produce the compound in culture |

UI treatment: each `evidence_tier` renders as a color-coded banner (red=chemistry only, amber=MASST only, sky blue=single cohort, emerald=cross-cohort or culture-validated, deep blue=longitudinal clinical cohort).

---

## Naming convention for literature-recruited rows

The atlas inherits mixed naming conventions from its sources — ChEBI uses one style (`L-alanocholate`), source papers use another (`Ala-CA`), and the atlas currently stores one `compound_name` per row. To keep display options open, **every literature-recruited row has its `synonyms` field structured as three ordered parts separated by ` | `**:

```
<full chemical name> | aka <paper-notation name> | <chemistry metadata>
```

| Part | Purpose | Example |
|---|---|---|
| **Full chemical name** | Descriptive, human-readable, unambiguous as to chemistry class | `phenylalanine-3-oxolithocholic acid alkamine` |
| **Paper-notation name** | The exact spelling the source paper uses (underscore or hyphen notation) | `aka Phenylalanine_3-oxo-LCA` |
| **Chemistry metadata** | Formula, mass, caveats on stereochem, synthesis route, Gentry producers, etc. | `C33H49NO4, [M+H]+ 524.37; C3 stereochem unspecified (both R and S exist in biology)` |

The **`compound_name` field keeps the short acronym form** (as it appears in paper figures and the current atlas UI) because:
1. Compact rendering in the UMAP scatter plot + search typeahead
2. Direct match to paper figure labels (helpful when cross-referencing)
3. The full name is preserved in `synonyms[0]`, so the UI can swap the two if desired later

### Per-class suffix conventions

| Chemistry class | Acronym style (compound_name) | Full-name suffix (synonyms part 1) |
|---|---|---|
| **MCBAs / bile amidates** (Gentry, Lamichhane, Mohanty — C24 carboxyl → amide) | `{AA}-{BA}` (e.g., `Ala-CA`, `CA-Cys`) | `{amino acid}yl-{BA full name} amidate` (e.g., `L-cysteinyl-cholic acid amidate`) |
| **CSBAs** (Lucas — glycyl / tauryl + 3-oxo BA) | `G-3-oxo{BA}` or `T-3-oxo{BA}` | `glyco-3-oxo{BA full name} amidate` |
| **Alkamines** (Agongo — C3 ketone → secondary amine) | `{amine}-{ketone-core}` (e.g., `Phe-3-oxo-LCA`) | `{amine full name}-{ketone-core full name} alkamine` (e.g., `phenylalanine-3-oxolithocholic acid alkamine`) |

### Known inconsistencies flagged for the deferred cleanup pass

- Mixed naming for the same compound (e.g., `Phe-CA` AND `L-phenylalanocholate` both exist in the atlas as separate rows — Rimal/Guzior era vs Gentry era)
- Ionization-form duplicates (e.g., `cholate` vs `cholic acid`, `taurocholate` vs `taurocholic acid`)
- Stereochem-form duplicates (e.g., L-D amino-acid variants from the Gentry re-run)

The naming convention above applies to all **newly-added literature-recruited rows (Lamichhane, Lucas completion, Agongo)** as of this audit. Pre-existing rows from earlier ingestion passes will be normalized in the end-of-audit cleanup.

---

## Gentry 2023 ✓ walked through

**Gentry, E.C. et al. Reverse metabolomics for the discovery of chemical structures from humans. *Nature* 626, 419–428 (2024).**
DOI: [10.1038/s41586-023-06906-8](https://doi.org/10.1038/s41586-023-06906-8)
Paper PDF in repo: `literature/papers/Gentry2023_reversemetabolomics_steroidmolecules.pdf`
SI files in repo: `literature/supplementary/Gentry2023_MOESM{2..10}_ESM.{pdf,xlsx,csv}`

### What Gentry 2023 contributes

A combinatorial reverse-metabolomics library of **168 bile amidates** = 21 amino acids × 8 bile acid cores (CA, CDCA, DCA, UDCA, HDCA, α-/β-/γ-muricholic acid). Each compound was synthesized as a reference standard, MS/MS spectrum acquired, then searched against public metabolomics data (MASST) and two independent human IBD cohorts (iHMP2 + Boland/PRISM) for isomer-resolved detection.

### Chemistry provenance (how we got the SMILES)

- **8 bile acid cores:** 5 from PubChem (CA CID 221493, CDCA 10133, DCA 222528, UDCA 31401, HDCA 5283820), 3 from user-approved ChEBI reference (α-MCA, β-MCA, γ-MCA — PubChem had stereochem inconsistencies for the muricholic acids; the ChEBI-consensus SMILES were used instead)
- **21 L-amino acids:** standard SMILES; L-proline handled specially for its secondary α-nitrogen
- **Amide coupling:** RDKit reaction SMARTS `[CX3:1](=O)[OH].[NX3H2:3][CX4] >> [CX3:1](=O)[NH:3][CX4]` between BA C24 carboxyl and amino acid α-amine
- **Mass validation:** all 168 compounds computed [M-H]⁻ monoisotopic mass matched Gentry's MOESM8 reported precursor m/z within ±0.1 ppm
- **Ingestion script:** `scratchpad/build_smiles_table.py` → final table at `scratchpad/gentry_amidate_final.csv`

### Evidence tiers for Gentry 2023 (the 4 streams in the paper)

| Tier | Count | Meaning | UI banner |
|---|---:|---|---|
| **1** | 83 | Chemistry reference standard only — no biological detection reported in Gentry 2023 | red "⚠ Chemistry reference only" |
| **2** | 44 | MASST public-data match (isomer uncertain, MSI Level 2/3) | amber "⚠ MS/MS match in public data" |
| **3** | 15 | Isomer-resolved detection in iHMP2 IBD cohort (RT+MS/MS matching, Fig 3a) | sky-blue "✓ Detected in human IBD stool" |
| **4-5** | 13 | Cross-cohort replicated in iHMP2 + Boland IBD cohorts (strongest human evidence) | emerald "✓ Cross-cohort validated" |
| **5** | 13 | Produced in bacterial culture by ≥1 HMP isolate (Fig 4, Supp Data 4; MOESM6 producer grid) | emerald "✓ Produced by gut bacteria in culture" |

Any row with `evidence_tier ∈ {3, 4-5, 5}` is considered biologically-validated per Gentry 2023. The 83 Tier-1 rows (including all 63 MCA conjugates) are present in the atlas as chemistry reference only — the loud red banner warns users that no biological detection was reported.

### Where to find each evidence stream in the SI

| Stream | File | Content |
|---|---|---|
| Stream 1 (chemistry) | `Gentry2023_MOESM8_ESM.xlsx` | 168-row table: compound, precursor m/z, RT, CCS, mass error, MS2 fragment |
| Stream 2 (MASST) | `Gentry2023_MOESM4_ESM.csv` | 2,792 MS/MS matches across all compound classes with USI links |
| Stream 3 (iHMP2) | `Gentry2023_MOESM5_ESM.pdf` | iHMP2 cohort detection list — 28 bile amidates |
| Stream 3 (Boland) | `Gentry2023_MOESM9_ESM.pdf` | Boland IBD cohort detection list — 12 bile amidates |
| Stream 4 (bacterial culture) | `Gentry2023_MOESM6_ESM.csv` | Bile acid × HMP species grid (28 compounds × 72 species) |

### Audit script trail

Build + validate SMILES + attach evidence tiers: `scratchpad/add_gentry_molecules.py`

---

## Mohanty 2024 ✓ walked through

**Mohanty, I. et al. The underappreciated diversity of bile acid modifications. *Cell* 187, 1801–1818 (2024).**
DOI: [10.1016/j.cell.2024.02.019](https://doi.org/10.1016/j.cell.2024.02.019)
Paper PDF in repo: `literature/papers/Mohanty2024_bileacidmoddiversity.pdf`
SI files in repo: `literature/supplementary/Mohanty2024_mmc{1..5}.{pdf,xlsx}`

### What Mohanty 2024 contributes

**17 polyamine bile amidates** synthesized as reference standards and validated by MS/MS + retention-time matching against detections in public metabolomics data (mice, lion, and other carnivore feces):

- **10 Cholyl- amidates** (CA core, Figure 5A)
- **7 Chenodeoxycholyl- amidates** (CDCA core, Figure 5C)

### Chemistry provenance

- **BA cores:** CA and CDCA from the same user-approved authoritative source used for Gentry (see above)
- **11 polyamine building blocks:** putrescine, cadaverine, spermidine, spermine, agmatine, 1,3-diaminopropane, N-acetyl-putrescine, N-acetyl-cadaverine, N-acetyl-spermidine, N-acetyl-spermine, N-carbamoyl-putrescine — SMILES constructed with the terminal primary amine explicit for amide coupling
- **Amide coupling:** same RDKit SMARTS as Gentry, targeting the polyamine's α-primary-amine only (asymmetric polyamines like spermidine are coupled at the propyl terminus; N-acetyl variants have the acetyl group on the opposite terminal amine)
- **Mass validation:** all 17 compounds' computed [M+H]⁺ matched Mohanty's Figure 5A/5C reported m/z within ±40 ppm (within the ±500 ppm tolerance of Mohanty's 2-decimal m/z precision)
- **Build script:** `scratchpad/mohanty_amidates.py` → final table at `scratchpad/mohanty_amidate_final.csv`
- **Ingestion script:** `scratchpad/add_mohanty_amidates.py` (updates 4 existing + appends 13 new)

### Evidence tiers for Mohanty 2024

| Tier | Count | Meaning | Evidence source |
|---|---:|---|---|
| **M-1** | 6 | Synthesized standard + RT+MS/MS match in **both mice AND lion feces** (strongest) | Fig 5B/D |
| **M-2** | 6 | Synthesized standard + RT+MS/MS match in **mice only** | Fig S5A |
| **M-3** | 2 | Synthesized standard + RT+MS/MS match in **lion only** | Fig S5B |
| **M-4** | 3 | Synthesized standard + same-mass peak detected in lion samples but RT did NOT match — suggests a related BA isomer present in vivo rather than the exact CDCA standard | Fig 5D, S5B, S6C |

3 of the M-1 compounds also carry **ICH-validated absolute quantification** (mmc5): Cholyl-putrescine, Cholyl-cadaverine, Cholyl-N-acetyl-putrescine — the paper's highest evidence tier.

### Per-compound assignments

**Cholyl- series (10):**
- M-1: Cholyl-putrescine, Cholyl-spermidine, Cholyl-cadaverine, Cholyl-N-acetyl-cadaverine
- M-2: Cholyl-agmatine, Cholyl-N-acetyl-putrescine, Cholyl-N-carbamoyl-putrescine, Cholyl-N-acetyl-spermidine, Cholyl-N-acetyl-spermine, Cholyl-1,3-diaminopropane

**Chenodeoxycholyl- series (7):**
- M-1: Chenodeoxycholyl-N-acetyl-putrescine, Chenodeoxycholyl-N-acetyl-spermidine
- M-3: Chenodeoxycholyl-putrescine, Chenodeoxycholyl-N-acetyl-cadaverine
- M-4: Chenodeoxycholyl-cadaverine, Chenodeoxycholyl-spermidine, Chenodeoxycholyl-1,3-diaminopropane

### Where to find each evidence stream in the SI

| Stream | File | Content |
|---|---|---|
| Chemistry (references) | `Mohanty2024_mmc2.xlsx` sheet `Inclusion_list` | 21,550 candidate compounds — explicitly-named polyamines are a subset |
| ICH-validated quantification | `Mohanty2024_mmc5.xlsx` sheet `Method_validation_ICH_guideline` | 5 compounds with linearity, LOD, LOQ, precision, accuracy per ICH guidelines |
| Polyamines in big cats | `Mohanty2024_mmc4.xlsx` sheet `Polyamines_in_felines` | Cholyl-putrescine / -cadaverine / -N-acetyl-putrescine quantified in lion, serval, cheetah, jaguar, jaguarundi, tiger |
| Mirror plots + spectral links | `Mohanty2024_mmc4.xlsx` sheet `Mirror_plot_USIs` | GNPS-dashboard USI links for mice/lion detections per compound |

### Known issues (to fix in next cleanup pass)

- The row `Spermidine-CA (3-trihydroxy-spermidine bile amidate)` and the row `Cholyl-spermidine` now have the same canonical SMILES — duplicate representing the same molecule under two legacy name conventions
- Two pre-existing atlas rows attributed to Mohanty (`3-acetyl-CA`, `3-succinyl-CA`) are NOT in Fig 5A/5C synthesized set — their Mohanty attribution should be re-sourced (they likely came from a different section of the paper or earlier Dorrestein-lab work)

---

## Lamichhane 2025 ✓ walked through

**Lamichhane, S. et al. Microbiome-derived bile acid signatures in early life and their association with islet autoimmunity. *Nature Communications* (2026, in press).**
DOI: [10.1038/s41467-025-66619-6](https://doi.org/10.1038/s41467-025-66619-6)
Paper PDF in repo: `literature/papers/Lamichhane2025_bileacid.pdf`
SI files in repo: `literature/supplementary/Lamichhane2025_MOESM{1,2}_ESM.pdf` (MOESM1 = supplementary figures + Tables S1–S7; MOESM2 = peer review file)

### What Lamichhane 2025 contributes

Longitudinal fecal metabolomics of the DIABIMMUNE cohort — 303 stool samples from 74 Finnish children at 3, 6, 12, 18, 24, and 36 months, split into CTR (seroconverted-negative), P1Ab (single islet autoantibody), and P2Ab (multiple autoantibodies → T1D-progression-prone) groups. Targeted LC-MS/MS of **110 microbially conjugated bile acids (MCBAs)** = 5 bile-acid cores (CA, CDCA, DCA, UDCA, HDCA) × 22 amino acids (Ala, Arg, Asn, Asp, Cit, **Cys**, Gln, Glu, His, Ile, Leu, Lys, L-DOPA, Met, Orn, Phe, Pro, Ser, Thr, Trp, Tyr, Val). Gentry is a co-author; the targeted panel builds directly on their synthesized reference library — except for cysteine, which Gentry did not include. 78/110 MCBAs were detectable in the cohort.

### Cross-check vs existing atlas

| Lamichhane's 110-MCBA universe | In atlas via | Overlap |
|---|---|---:|
| 5 cores (CA/CDCA/DCA/UDCA/HDCA) × 21 shared AAs | Gentry 2023 (8×21 library) | 105 compounds |
| 5 cores × **Cys** (not in Gentry's 21-AA panel) | **NEW** | **5 compounds** |

So Lamichhane's only genuinely atlas-novel compounds are the **5 Cys-conjugate MCBAs**: **CA-Cys**, **CDCA-Cys**, **DCA-Cys**, **UDCA-Cys**, **HDCA-Cys**. The other 105 overlap with Gentry's library already present.

### Chemistry provenance

- **5 BA cores:** reused from Gentry walk-through (`scratchpad/final_ba_cores.json`)
- **L-Cysteine SMILES:** `N[C@@H](CS)C(=O)O` — α-amine couples via amide bond to BA C24 carboxyl; thiol (-SH) stays free (consistent with biological conjugation chemistry; Gentry/Dorrestein-lab coupling methodology)
- **Amide coupling:** RDKit reaction SMARTS targeting primary α-amine on chiral center; thiol unperturbed
- **Mass validation:** all 5 computed [M+H]⁺ masses matched theoretical values within <10 ppm
  - CA-Cys: C₂₇H₄₅NO₆S, [M+H]⁺ = 512.3040 (expected 512.30)
  - CDCA/DCA/UDCA/HDCA-Cys: C₂₇H₄₅NO₅S, [M+H]⁺ = 496.3091 (expected 496.31)
- Lamichhane does not publish a per-compound m/z table in MOESM1 — they use authentic synthetic standards with RT+MS/MS matching (same Gentry-lab targeted methodology). Theoretical masses are the rigorous validation target.
- **Ingestion scripts:** `scratchpad/lamichhane_cys_mcbas.py` → `scratchpad/lamichhane_cys_mcbas.csv`; `scratchpad/add_lamichhane_cys_mcbas.py` appended to `data/molecules.csv`
- **Structure rendering:** `scratchpad/render_lamichhane_cys_structures.py` generated 5 PNGs + updated `structures_index.json`

### Evidence tiers for Lamichhane 2025

| Tier | Count | Meaning | UI banner |
|---|---:|---|---|
| **L-1** | 1 | Targeted LC-MS/MS detection in 303-sample DIABIMMUNE cohort **+ statistically significant association with islet autoimmunity** (linear mixed-effect model, p<0.05) | deep-blue "✓ Associated with islet autoimmunity" |
| **L-2** | 4 | Part of 110-MCBA targeted panel in 303-sample DIABIMMUNE cohort (no significant clinical association reported) | sky-blue "✓ Measured in DIABIMMUNE cohort" |

Per-compound assignments (based on MOESM1 Supplementary Tables S1–S4 linear mixed-effect model outputs):

| Compound | Tier | Finding |
|---|---|---|
| CA-Cys | **L-1** | Elevated in two-antibody (TwoAb) case group vs controls (coef = +0.095, p = 0.045, q = 0.291; MOESM1 Table S2). Also referenced in breastfeeding-length analyses. |
| CDCA-Cys | L-2 | Panel-detected, no significant association in Tables S1–S4 |
| DCA-Cys | L-2 | Panel-detected, no significant association in Tables S1–S4 |
| UDCA-Cys | L-2 | Panel-detected, no significant association in Tables S1–S4 |
| HDCA-Cys | L-2 | Panel-detected, no significant association in Tables S1–S4; HDCA-conjugate trajectories shown in Fig S1 (figure focuses on amino-acid conjugates generally, not Cys specifically) |

### UMAP placement

All 5 Cys-MCBAs placed by Tanimoto-nearest-neighbor projection + Gaussian jitter (σ=0.15). NN in every case is the corresponding Ser-BA (Gentry library), which is chemically exact: cysteine and serine differ only in the terminal thiol vs hydroxyl. Tanimoto similarity ≥ 0.85 for every pairing. All inherit cluster 0 (the main amidate cluster).

### What Lamichhane 2025 does NOT add (caveat)

Of the 27 pre-existing rows in `molecules.csv` that carry the Lamichhane DOI as `paper_url`, roughly half look miscategorized and should be revisited in the end-of-audit cleanup pass:

- **Legitimate MCBAs** (14 rows): Tyr-CDCA, Met-DCA, Phe-DCA, Thr-CA, His-DCA, Asn-UDCA, Trp-UDCA, Asp-UDCA, Glu-UDCA, Ser-CDCA, Met-CDCA, Tyr-HDCA, Tyr-UDCA, Ser-UDCA — these are part of Lamichhane's panel and the attribution is correct
- **Secondary BAs likely from context, not primary focus** (10 rows): 3-oxoDCA, 3-oxo-CA, 3-oxo-LCA, iso-LCA, 3-oxo-DCA (duplicate!), 3-oxo-CDCA, 12-oxo-LCA, 7-oxo-LCA, isoallo-LCA, allo-LCA — Lamichhane references these as general secondary BAs but did not quantify or associate them; attribution should be re-sourced
- **Duplicates / likely misattributions** (3 rows): Gly-DCA (Gly is not in Lamichhane's 22-AA panel), Hyodeoxycholic acid × 2 (the free HDCA core, duplicated with and without the "(HDCA)" suffix)

These existing rows are NOT modified by this walk-through — they stay as-is until the deferred cleanup pass.

---

## Xue 2025 ✓ walked through — annotation-enrichment only (0 new compounds)

**Xue, J. et al. Gut microbiota and derived metabolites mediate obstructive sleep apnea induced atherosclerosis. *Gut Microbes* 17:2474142 (2025).**
DOI: [10.1080/19490976.2025.2474142](https://doi.org/10.1080/19490976.2025.2474142)
PMC: [PMC11881840](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11881840/) (open access)
Paper PDF in repo: `literature/papers/Xue2025_Gut microbiota and derived metabolites mediate obstructive sleep apnea induced atherosclerosis.pdf`
Preprint version (same compound set): bioRxiv `10.1101/2024.11.18.624205` (Nov 2024)

### What Xue 2025 contributes

Mouse study. ApoE⁻/⁻ mice fed HFHC (high-fat/high-cholesterol) ± IHC (intermittent hypoxia/hypercapnia — a sleep-apnea mimic), with regular-chow (RC) + room-air controls. The study shows that gut microbiota and derived **bile acids** mediate OSA-induced atherosclerosis, with **Muribaculaceae vs Akkermansiaceae** as the top discriminating bacterial families. **DCA and LCA** are singled out in the abstract as the key microbial BA effectors.

### Methodology (how they detected the compounds)

| Aspect | Value |
|---|---|
| Sample matrix | Mouse fecal |
| LC | Vanquish UHPLC, Kinetex C18 50×2 mm 1.7 µm |
| MS | Q Exactive (Thermo); **untargeted** LC-MS/MS |
| Feature detection | MZmine 4.2 |
| Annotation | **GNPS2 feature-based molecular networking**, MSI Level 2–3 — "annotated or putatively annotated" against public spectral libraries, **no authentic standards for every ID** |
| Readout | Relative peak-area intensities + PLS-DA VIP > 1 (not absolute concentrations) |
| Statistics | `mixOmics` PLS-DA; ANCOM-BC2 for differential bacterial-family abundance |
| Supplementary tables with the values | S3 (HFHC-associated, 3,120 features) and S4 (IHC-associated, 2,034 features) — released as `KGMI_A_2474142_SM0621.zip` on PMC |

### Cross-check vs existing atlas

Every one of the 23 bile-acid / MCBA compounds Xue reports is **already in the atlas** — mostly via Gentry 2023 (MCBAs) and the standard-BA entries from the original corpus. One name clarification: Xue's "12-ketodeoxycholic acid" is the same compound as the atlas's `12-dehydrodeoxycholate` (CHEBI:233653) — canonical names differ, SMILES byte-identical (verified via RDKit).

NOT attributed: `muricholic acid`, `nor cholic acid`, `isolithocholic acid` — Xue cites these via ref 75 (Ruminococcaceae background) rather than detecting them in their own data.

### Evidence tier for Xue 2025

| Tier | Meaning | UI banner |
|---|---|---|
| **X-2** | Untargeted LC-MS/MS annotation (MSI Level 2–3) + relative-intensity PLS-DA VIP > 1 in Xue's HFHC/IHC mouse model | sky-blue "✓ Detected in HFHC/IHC mouse model" |

Direction of effect appended to `evidence_note` per compound (lifted from Fig 3 + Supp Tables S3/S4):

| Direction | Compounds |
|---|---|
| **Elevated by both HFHC (vs RC-Air) and IHC (vs HFHC-Air)** | `taurocholate`, `taurodeoxycholate`, `12-dehydrodeoxycholate` (= 12-ketodeoxycholic acid) |
| Elevated by HFHC vs RC-Air | `glycocholate`, `Ile-CA`, `Leu-CA`, `Phe-CA` (+ `L-phenylalanocholate` duplicate), `Tyr-CA`, `Val-CA` |
| Elevated by IHC (vs HFHC-Air) | `taurolithocholate` |
| Elevated by HFHC **but** suppressed by IHC | `cholate`, `Lys-CA`, `Phe-DCA`, `Glu-DCA` |
| Suppressed by IHC | `Arg-CA` (+ `L-arginocholate`), `Glu-CA`, `Gln-CA`, `Tyr-DCA` |
| Called out in abstract as key atherosclerosis effectors | `deoxycholate`, `lithocholate` |

Caveats recorded in `evidence_note`:
- Ile-CA and Leu-CA both tagged with "reported as Ile/Leu-CA — MS cannot distinguish the two isomers (same mass)"
- 3 Ile-CA/Leu-CA/Val-CA duplicate pairs exist in the atlas (one row from original Rimal/Guzior ingestion, one from Gentry 2023 library) — both rows updated with Xue annotation; dedup is deferred to end-of-audit cleanup

### Files + scripts

- `scratchpad/xue_parse/xue2025.txt` — pdftotext extraction of the paper
- `scratchpad/enrich_xue_annotations.py` — the enrichment script (no new rows, only `paper_url` + `evidence_note` on 25 rows representing 22 unique compounds)
- `data/molecules.csv.bak_before_xue` — pre-change backup (844 rows, same as post-change)

### What the user sees in the UI

Click any compound in the HFHC-IHC-affected set (e.g., `taurodeoxycholate`) and the evidence-note banner now includes a line like:
> "Xue 2025 (X-2): elevated by both HFHC diet (vs RC-Air) AND IHC (vs HFHC-Air)."

Combined with any pre-existing Gentry or Rimal evidence, the molecule's provenance is now:
`[Gentry 2023 Tier X] | Xue 2025 (X-2): …`

---

## Agongo 2026 ✓ walked through — NEW chemistry class (alkamines)

**Agongo, J. et al. Alkamines reveal a hidden layer of steroid and drug metabolism. *bioRxiv* (2026-05-15).**
DOI: [10.64898/2026.05.13.724743](https://doi.org/10.64898/2026.05.13.724743)
Paper PDF in repo: `literature/papers/Agongo2026_alkemines_steroid.pdf`

### What Agongo 2026 contributes — a NEW chemistry class

**Alkamines = reductive amination products.** Agongo introduces a fundamentally new class of steroid/BA modifications to the atlas:

| | MCBAs / amidates (Gentry, Lamichhane, Mohanty, Lucas) | **Alkamines (Agongo)** |
|---|---|---|
| Linkage | C24 carboxyl → amide | **C3 ketone → secondary amine** |
| Chemistry | `R-COOH + H₂N-R' → R-C(=O)-NH-R' + H₂O` | `R-C(=O)-R' + H₂N-R'' → R-CH(NHR'')-R' + H₂O` (ketone + amine → imine + H₂O, then +H₂ → secondary amine) |
| Starting material | Free BA (CA, CDCA, DCA, UDCA, HDCA) | **Ketone-containing BA** (3-oxo-LCA, 3-oxo-DCA) or steroid (cholestan-3-one) |
| Stereochem outcome | Preserves BA scaffold stereochem | **Creates new stereocenter at C3** (α or β); Agongo tested both R and S of N-acetylcadaverine-LCA alkamine |

Diagnostic MS fragments (for identification): m/z 359.294 (LCA-alkamine), 357.278 (DCA-alkamine), 355.260 (CA-alkamine) — consistent with secondary amine formation at C3.

Library scope: 8,475 synthesized alkamines (96 aldehydes/ketones × amines); 1,626 matched to public LC-MS/MS data; **56 confirmed to Level 1** in biological samples (including **27 bile-acid + steroidal alkamines** and 12 drug-derived); 14 of these further supported by ion mobility drift time. 77% of the synthesized library is absent from PubChem.

### Chemistry provenance (how we built the 18 atlas rows)

- **Ketone cores from atlas:**
  - `3-oxo-LCA` (C₂₄H₃₈O₃, [M] 374.28) and `3-oxo-DCA` (C₂₄H₃₈O₄, [M] 390.28) — pre-existing rows (already CHEBI-provenanced)
  - `5β-cholestan-3-one` (C₂₇H₄₆O, CHEBI:16074, [M] 386.35) — pre-existing row
- **Amine library (12 amines used):** Phe, Ala, Gly, Leu, Val, Tyr, Lys (as free α-amine forms); Tyramine, Putrescine, N-acetylcadaverine, N-acetylputrescine (biogenic polyamines); 5-ASA (5-aminosalicylic acid / mesalamine — drug)
- **Reductive-amination SMARTS** (RDKit): `[CX3:1](=[OX1])(-[CH2;R:2])(-[CH2;R:3]).[NX3H2:4][*:5] >> [C:1](-[CH2;R:2])(-[CH2;R:3])-[NX3H1:4][*:5]`  
  Ring-based ketone carbons (`CH2;R`) ensure the reaction targets C3 and NOT the C24 carboxyl.
- **C3 stereochem left unspecified** — the paper reports both R and S enantiomers of N-acetylcadaverine-LCA alkamine in biological samples. Molecular formula and MW are unaffected.
- **Mass validation:** all 18 formulas and exact masses computed; C3 ketone → sp³ C–NH–R bookkeeping verified via RDKit `ExactMolWt` and `CalcMolFormula`.
- **Ingestion scripts:** `scratchpad/agongo_alkamines.py` (build + validate), `scratchpad/add_agongo_alkamines.py` (append to `molecules.csv`), `scratchpad/render_agongo_alkamines.py` (PNGs)
- **Backup:** `data/molecules.csv.bak_before_agongo`

### 18 compounds added

**Bile acid alkamines (8):**

| Compound | Formula | [M+H]⁺ | Tier | Highlight |
|---|---|---:|---|---|
| Phe-3-oxo-LCA | C₃₃H₄₉NO₄ | 524.37 | **AG-1** | Detected in human feces (Fig 4F); **microbially produced** (B. caccae, B. luhongzhouii, C. innocuum, E. ramosum, cCOM — Fig 4G); 1.9–14.2 µM in animal fecal samples |
| N-acetylcadaverine-3-oxo-LCA | C₃₁H₅₄N₂O₃ | 503.42 | **AG-1** | **Strongest GPCR / nuclear receptor response** in HL-60 Ca²⁺ assay (Fig 4H-I); S-enantiomer is the active form |
| N-acetylcadaverine-3-oxo-DCA | C₃₁H₅₄N₂O₄ | 519.42 | AG-1 | Confirmed Fig 4C |
| Tyramine-3-oxo-DCA | C₃₂H₄₉NO₄ | 512.37 | AG-2 | Confirmed Fig 4C |
| Putrescine-3-oxo-DCA | C₂₈H₅₀N₂O₃ | 463.39 | AG-2 | Confirmed Fig 4C |
| Tyr-3-oxo-DCA | C₃₃H₄₉NO₆ | 556.36 | AG-2 | Confirmed Fig 4C |
| Phe-3-oxo-DCA | C₃₃H₄₉NO₅ | 540.37 | AG-1 | Quantified in microbial culture (Table S6) |
| Lys-3-oxo-DCA | C₃₀H₅₂N₂O₅ | 521.39 | AG-2 | Confirmed Fig 4C |

**Drug-BA alkamines (2):**

| Compound | Formula | [M+H]⁺ | Tier | Highlight |
|---|---|---:|---|---|
| 5-ASA-3-oxo-LCA | C₃₁H₄₅NO₅ | 512.34 | **AG-1** | Detected in **IBD patient stool**; 97% of 5-ASA alkamine MS/MS matches occur in IBD samples (Fig 5A-F); previously unreported |
| 5-ASA-3-oxo-DCA | C₃₁H₄₅NO₆ | 528.33 | AG-1 | IBD patient stool; associated with mesalamine prescribing patterns |

**Cholestan-3-one alkamines (8, non-BA steroidal):**

| Compound | Formula | [M+H]⁺ | Tier |
|---|---|---:|---|
| Phe-cholestan-3-one | C₃₆H₅₇NO₂ | 536.45 | AG-2 |
| Ala-cholestan-3-one | C₃₀H₅₃NO₂ | 460.41 | AG-2 |
| Gly-cholestan-3-one | C₂₉H₅₁NO₂ | 446.40 | AG-2 |
| Leu-cholestan-3-one | C₃₃H₅₉NO₂ | 502.46 | AG-2 |
| Val-cholestan-3-one | C₃₂H₅₇NO₂ | 488.45 | AG-2 |
| Tyramine-cholestan-3-one | C₃₅H₅₇NO | 508.45 | AG-2 |
| N-acetylcadaverine-cholestan-3-one | C₃₄H₆₂N₂O | 515.49 | AG-2 |
| N-acetylputrescine-cholestan-3-one | C₃₃H₆₀N₂O | 501.48 | AG-2 |

### Evidence tiers for Agongo 2026

| Tier | Count | Meaning | UI banner |
|---|---:|---|---|
| **AG-1** | 6 | Level 1 identification = MS/MS + retention time match to authentic synthetic standard in biological samples **+ ion mobility drift time match** (strongest; Agongo's 14-of-27 highest-confidence set) | emerald "✓ Level 1 + ion mobility confirmed" |
| **AG-2** | 12 | Level 1 identification = MS/MS + retention time match to authentic synthetic standard in biological samples (no ion mobility) | sky-blue "✓ Level 1 MS/MS + RT confirmed" |

Not used (no additions at this tier):
- **AG-3**: Library match only (MSI Level 2, no biological sample confirmation) — only biologically confirmed compounds added here.

### UMAP placement

Each alkamine placed by Tanimoto-nearest-neighbor projection + Gaussian jitter (σ=0.15):
- **BA alkamines** → NN is the corresponding MCBA/BA (e.g., `Phe-3-oxo-LCA` → NN `Phe-DCA` sim 0.61; `Tyr-3-oxo-DCA` → NN `Tyr-DCA` sim 0.75). Chemically sensible — they share scaffold despite different attachment chemistry.
- **Cholestane alkamines** → NN is `cholestanol` (sim 0.57–0.68). Consistent with the C27 cholesterol-derived scaffold.
- All inherit **cluster 0** (the main amidate/steroid-conjugate cluster) pending full UMAP recompute.

### Not added (noted in paper but out of atlas scope)

Agongo reports **8,475 synthesized alkamines**, of which 1,626 had at least one MS/MS match in public data. We added only the **18 steroid/bile-acid alkamines with Level 1 biological confirmation** in the paper's main figures. The remaining 27 - 18 = 9 Level 1 confirmed steroid alkamines (from Table S5, not fully enumerated in the main text) could be a follow-up add if the SI is accessed. Drug alkamines beyond 5-ASA-BA (e.g., sulfapyridine conjugates) are out of steroid-atlas scope.

---

## Reilly 2026 ✓ walked through — annotation-enrichment only (0 new compounds)

**Reilly, E.R. et al. Systemic rhythmicity of host and bacterial bile acid amidates in the mouse. *Cell Systems* 17:101541 (April 15, 2026).**
DOI: [10.1016/j.cels.2026.101541](https://doi.org/10.1016/j.cels.2026.101541)
Paper PDF in repo: `literature/papers/Reilly2026.pdf`
Authors include Mohanty and Dorrestein (same ecosystem as Lamichhane, Mohanty, Gentry, Aatsinki, Agongo).

### What Reilly 2026 contributes

**First systematic rhythmicity map of bacterial bile acid amidates (BBAAs = MCBAs).** 690 samples from time-restricted-fed mice, profiled across **14 organs + digestive tract contents + biofluids** at 6 timepoints (ZT0/4/8/12/16/20), to map the temporal and spatial distribution of BAs and their microbial modifications. Introduces the **"encoder/decoder" hypothesis**: microbes "encode" messages via metabolite timing and structure that travel to distant cells/tissues and are "decoded" there.

### Methodology

| Aspect | Value |
|---|---|
| Sample matrix | Mouse: 14 organs (stomach, duodenum, jejunum, ileum, cecum, colon, liver, heart, lungs, spleen, brain, eye, skin, kidney) + GI contents + feces + serum + urine |
| N | 690 samples, 5 mice/timepoint × 6 timepoints × multiple tissues |
| LC-MS | Targeted LC-MS/MS with 460 MS/MS matches total |
| Standards | 13 isotope-labeled internal standards (LCA-d4, TCA-d4, GUDCA-d4, GCA-d4, CA-d4, UDCA-d4, GCDCA-d4, CDCA-d4, HDCA-d5, GLCA-d4, TLCA-d4, TUDCA-d5, GDCA-d4) → **MSI Level 1** for primary BAs and labeled conjugates |
| MCBA library | **Gentry 2023** combinatorial library (168 MCBAs, same Dorrestein-lab source) |
| Coverage | 102 bile acids identified total; 40-102 per organ (skin/brain/eyes fewest; GI contents most) |
| Statistics | **Cosinor-based rhythmometry** for periodic fitting |

### Cross-check vs existing atlas

**All 22 Reilly-named MCBAs are already in the atlas** via Gentry's 168 library + Lamichhane's 5 Cys additions. For "MCA" without α/β/γ specification (Reilly's convention), each generic name maps to **3 atlas rows** (α, β, γ muricholic isomers). Net effect: **45 unique atlas compound names covering 49 rows** (4 pre-existing duplicate pairs for Ile/Leu/Thr/Val-CA picked up).

**0 new compounds added.**

### Evidence tiers for Reilly 2026

| Tier | Direction category | Compounds (atlas names) |
|---|---|---|
| **R-1 / stable** | Hydrophilic or hydrophobic side chain → remained relatively constant across time and tissue | `Gln-CA`, `Ile-CA`, `Leu-CA` (Ile/Leu lumped by MS) |
| **R-1 / aromatic_variable** | Aromatic AA conjugates → exhibited the greatest temporal variability; production or host uptake may be context-dependent | `Phe-CA`, `Tyr-CA`, `Trp-CA`, plus `Phe/Tyr/Trp-{α,β,γ}MCA` (9 more) |
| **R-1 / structural_feature** | Shown structurally in Figure 1 | `Ser-CA`, `Arg-DCA`, `Lys-{α,β,γ}MCA` |
| **R-1 / rhythmic_named** | Named in Figure 4 rhythmicity panel (gut contents peak ZT4-ZT12 = fasted; peripheral organs peak ZT12-ZT20 = fed) | Ala-{CA, α/β/γMCA}, Arg-{CA, α/β/γMCA}, Asn-CA, His-CA, Lys-CA, Met-{CA, α/β/γMCA}, Orn-CA, Thr-CA, Val-CA, Ile-{α/β/γ}MCA, Leu-{α/β/γ}MCA, Lys-DCA |

Additionally **R-2** (1 row — `Hyodeoxycholic acid` free HDCA) applied to pre-existing Reilly-DOI stub rows that are in the 102-BA panel but not individually named in Results. Per-compound data for these lives in Figures S1-S4 (not accessible from the main PDF).

### Key findings recorded in evidence_note

- **Ileal contents** show the greatest temporal synchronization (peak ZT8); synchronization decreases moving distally and peripherally
- **Peripheral organs** (skin, brain, heart, kidney, spleen, lungs, eye) have asynchronous, tissue-specific patterns — first direct evidence of BBAAs **beyond the enterohepatic circulation**
- **Aromatic AA amidates (Phe, Tyr, Trp)** show greatest rhythmic variation across the body
- **Hydrophilic/hydrophobic amidates (Gln-CA, Ile/Leu-CA)** stay relatively constant — cited in paper as representative stable class
- BBAAs appear in gut/serum at ZT8-ZT12 (fasted) and shift to distant organs at ZT12-ZT20 (fed)

### Files + scripts

- `scratchpad/reilly_parse/reilly.txt` — pdftotext extraction of the paper
- `scratchpad/enrich_reilly_annotations.py` — the enrichment script (0 new rows; 49 updates representing 45 unique MCBAs)
- `data/molecules.csv.bak_before_reilly` — pre-change backup

---

## Aatsinki 2025 ✓ walked through — annotation-enrichment only (0 new compounds)

**Aatsinki, S. et al. Early life gut microbial and metabolic shifts reflected in stool consistency — a gut transit time proxy. *bioRxiv* (2025, latest ver 2025-09-16).**
DOI: [10.1101/2024.06.05.597641](https://doi.org/10.1101/2024.06.05.597641)
Paper PDF in repo: `literature/papers/Aatsinki2025_gutmicrobiome.pdf`

### What Aatsinki 2025 contributes

Longitudinal fecal metabolomics of the **FinnBrain Birth Cohort** — 618 stool samples from infants at 2.5, 6, 14, and 30 months. Measures **33 standard BAs + 110 MCBAs** (identical panel to Lamichhane 2025: 5 BA cores × 22 AAs, Gentry co-author). The paper's framing is **stool consistency and water content as a gut-transit-time proxy** — constipation-, breastfeeding-, and polygenic-score-associated MCBAs.

### Methodology (how they detected the compounds)

| Aspect | Value |
|---|---|
| Sample matrix | Infant fecal (freeze-dried + homogenized) |
| LC-MS | Targeted LC-MS/MS |
| Standards | Authentic reference + 10 isotope-labeled internal standards: LCA-d4, TCA-d4, GUDCA-d4, GCA-d4, CA-d4, UDCA-d4, GCDCA-d4, CDCA-d4, DCA-d4, GLCA-d4 |
| MSI confidence | **Level 1** for the 33 standard BAs (highest tier) |
| Per-compound data | Supplementary Tables S9 (237 water-content associations) and S10 (203 breastfeeding-adjusted); S12-S13 for mediation/PGS |

### Cross-check vs existing atlas

**MCBAs (110-panel):** 100% overlap with existing atlas via Gentry 2023 + Lamichhane 2025 (5 Cys additions). No new MCBAs added.

**Standard BAs (33-panel):** 22 of 33 already in the atlas. 11 missing from the atlas that are **listed in Aatsinki's Methods targeted panel** but NOT individually confirmed as detected in the paper text:

- `DHCA` (dihydroxycholestanoic acid), `THCA` (trihydroxycholestanoic acid) — C27 intermediates
- `7-oxo-HCA`, `ω-MCA` — free BA oxidation / mouse-type
- `GHDCA`, `GDHCA`, `GHCA` — glyco conjugates
- `TDHCA`, `TaMCA`, `TbMCA`, `TwMCA` — tauro conjugates

**None of these 11 are added** — per the atlas's detection-evidence criterion (same as Lucas 2026's G-12-oxoLCA/T-12-oxoLCA exclusion), panel-listed compounds without explicit biological-detection evidence in the paper text are not recruited. Aatsinki's Results discusses only "tauroconjugated" / "glycoconjugated" BAs as group categories, never naming these 11 individually. Per-compound data lives in Supplementary Tables S9-S13 which were not accessible from the main PDF. If the SI becomes available, these can be filtered + added in a follow-up.

### MCBA annotation enrichment (what Aatsinki does add)

**34 MCBAs** named individually in Results with specific clinical associations get Aatsinki DOI + direction-of-effect notes appended to existing atlas rows (same pattern as Xue 2025 enrichment):

| Direction | MCBAs |
|---|---|
| **Lower in constipated children** | `Asn-CA`, `Thr-CA`, `Arg-CA` |
| **Higher in constipated children** | `Val-CDCA`, `Lys-UDCA`, `Lys-HDCA`, `Ile-UDCA`, `Leu-UDCA` (UDCA/HDCA and Ile/Leu caveats on isomer indistinguishability) |
| Positively associated with water at 6 & 30 months | `Asn-CDCA`, `Trp-CDCA`, `Met-UDCA` |
| Negatively associated with water across ages | `Val-DCA`, `Val-HDCA`, `Val-UDCA`, `Ala-{CA,CDCA,UDCA,HDCA}`, `Trp-UDCA`, `Ile-DCA`, `Leu-DCA`, `Ile-HDCA`, `Leu-HDCA` |
| Elevated via mediation by constipation | `Orn-UDCA` (reported as "UDCA-Ornithine (1)" — isomer-specific), `Orn-HDCA` |
| Breastfeeding-associated CA conjugates | `His-CA`, `Met-CA`, `Orn-CA`, `Ser-CA` |
| PGS-associated (polygenic score) | `His-UDCA`, `Glu-UDCA`, `Glu-HDCA`, `Lys-CDCA`, `Gln-UDCA` |

### Evidence tiers for Aatsinki 2025

| Tier | Count | Meaning | UI banner |
|---|---:|---|---|
| **A-1** | 34 (40 rows w/ dupes) | Targeted LC-MS/MS detection in FinnBrain 618-sample cohort (MSI Level 1) + **statistically significant clinical association** (constipation / stool water / breastfeeding / PGS) | deep-blue "✓ Associated with stool phenotype / breastfeeding / PGS (FinnBrain cohort)" |
| **A-2** | 7 | Part of the targeted BA panel; measured in the 618-sample cohort but **not individually named in Results** — per-compound data lives in Supplementary Tables S9/S10. Applied to 7 pre-existing Aatsinki-DOI stub rows (free BAs: HCA, α-/β-MCA, 7-oxo-DCA, HDCA; plus a Lamichhane-era misattribution Gly-DCA flagged for cleanup) | gray "Aatsinki panel-measured (no individual Results discussion)" |

### Files + scripts

- `scratchpad/aatsinki_parse/aatsinki.txt` — pdftotext extraction of the paper
- `scratchpad/enrich_aatsinki_annotations.py` — the enrichment script (0 new rows; 40 atlas rows updated representing 34 unique MCBAs; duplicates covered)
- `data/molecules.csv.bak_before_aatsinki` — pre-change backup

---

## Lucas 2026 ✓ walked through

**Lucas, L.N. et al. Investigation of bile salt hydrolase activity in human gut bacteria reveals production of conjugated secondary bile acids. *Nature Communications* 17, 3077 (2026).**
DOI: [10.1038/s41467-026-68556-4](https://doi.org/10.1038/s41467-026-68556-4)
Paper PDF in repo: `literature/papers/Lucas2026_BSH.pdf`

### What Lucas 2026 contributes

The first systematic characterization of **conjugated secondary bile acids (CSBAs)** — bile acids where the steroid core is modified (7-oxo, 3-oxo) **while the glycine or taurine conjugate stays intact**. Lucas profiled 77 bacterial strains across 7 phyla, fed each strain conjugated BA substrates (TCA, TCDCA, TDCA, GCA, GCDCA) at 100 µM and 500 µM, and used uHPLC-MS/MS to track every BA species produced over time (Figs 1, 4B, 5).

The key mechanistic finding: some HSDs (e.g., *C. intestinalis* 3α-HSD, *B. ovatus* 7α-HSD) are active on **intact glyco- or tauro-conjugated primary BAs**, producing conjugated secondary BAs **without first requiring deconjugation by BSH**. This is new biology — previously, the assumption was that secondary-BA production required BSH-mediated deconjugation first.

### Chemistry provenance

- **Free 3-oxo-BA cores:** reused from the atlas's pre-existing `3-oxo-CA` and `3-oxo-CDCA` rows (CHEBI-sourced; mass + stereochem verified)
- **Glycine SMILES:** `NCC(=O)O` — α-amine attacks BA C24 carboxyl to form the amide
- **Amide coupling:** RDKit reaction SMARTS targeting primary α-amine on alpha-carbon
- **Mass validation:**
  - **G-3-oxoCA** = C₂₆H₄₁NO₆, [M+H]⁺ = 464.3007 (expected 464.30, Δ +1.4 ppm) ✓
  - **G-3-oxoCDCA** = C₂₆H₄₁NO₅, [M+H]⁺ = 448.3058 (expected 448.31, Δ −9.5 ppm) ✓
- **Ingestion scripts:** `scratchpad/add_lucas_g3oxo_mcbas.py` and `scratchpad/render_lucas_g3oxo_structures.py`

### Cross-check vs existing atlas (full Lucas CSBA inventory)

| Lucas 2026 CSBA | In atlas? | Action taken |
|---|---|---|
| G-7-oxoLCA | ✓ as `Glycine-7-oxoLCA` (also duplicates as CHEBI:137818 `7-oxoglycolithocholate` — flagged for dedup) | Tier backfilled to **LU-1** |
| G-7-oxoDCA | ✓ as `Glycine-7-oxoDCA` (also CHEBI:137824 `7-oxo-glycodeoxycholate` dup) | Tier backfilled to **LU-2** |
| T-7-oxoLCA | ✓ as `Taurine-7-oxoLCA` (also CHEBI:137724 `7-oxotaurolithocholate` dup) | Tier backfilled to **LU-2** |
| T-7-oxoDCA | ✓ as `Taurine-7-oxoDCA` (also CHEBI:137820 `7-oxo-taurodeoxycholate` dup) | Tier backfilled to **LU-2** |
| **G-3-oxoCA** | ✗ missing | **ADDED** (tier LU-2) |
| **G-3-oxoCDCA** | ✗ missing | **ADDED** (tier LU-2) |
| G-12-oxoLCA | ✗ — **not biologically detected** by Lucas (used only as commercial calibration standard for surrogate quantification) | **NOT added** per user decision: atlas requires biological-detection evidence, not standards-only inclusion |
| T-12-oxoLCA | ✗ — same (commercial calibration standard only) | **NOT added** |

### Evidence tiers for Lucas 2026

| Tier | Count | Meaning | UI banner |
|---|---:|---|---|
| **LU-1** | 1 | Authentic synthesized standard (Lucas synthesized their own reference material) + biological detection in bacterial culture — highest Lucas confidence | emerald "✓ Lucas-synthesized authentic standard + detected" |
| **LU-2** | 5 | LC-MS/MS detection in bacterial culture; exact mass + RT consistent with chemistry; **concentration estimated via G-12-oxoLCA surrogate calibration curve** — no authentic standard for this specific compound | sky-blue "✓ Detected in bacterial culture (surrogate-standard quantification)" |

Per-compound assignments:

| Compound | Tier | Detection context |
|---|---|---|
| Glycine-7-oxoLCA | **LU-1** | Lucas synthesized this standard (Methods: "Conjugated secondary bile acid synthesis"). Detected in *B. ovatus* and *F. varium* cultures (Fig 4B) — CDCA → 7-oxoLCA on intact glyco-conjugate |
| Glycine-7-oxoDCA | LU-2 | *B. ovatus* (from GCA after CA accumulated) and *F. varium* (from GCA directly) |
| Taurine-7-oxoLCA | LU-2 | Product of 7α-HSD activity on T-conjugated substrates |
| Taurine-7-oxoDCA | LU-2 | Product of 7α-HSD activity on T-conjugated substrates |
| **G-3-oxoCA** | **LU-2** | *C. intestinalis* monoculture (Fig 4B) — 3α-HSD action on intact GCA; glyco-conjugate retained because *C. intestinalis* BSH does not cleave glycine conjugates; peaked ~15 µM at 6h |
| **G-3-oxoCDCA** | **LU-2** | *C. intestinalis* monoculture (Fig 4B) — 3α-HSD action on intact GCDCA; also observed in *B. ovatus* and *F. varium* |

### UMAP placement

Both new G-3-oxo CSBAs placed by Tanimoto-nearest-neighbor projection + Gaussian jitter (σ=0.15). NN in both cases is the corresponding **free 3-oxo-BA** (sim ≈ 0.82) — chemistry-consistent. Both inherit cluster 0 (main amidate / conjugated-BA cluster).

### Protein-side (NO backlinks — Lucas did not purify any enzymes; Q8A6H3 removed from atlas)

Lucas 2026 performed **no protein purification** anywhere in the paper. All in vitro experiments are **whole-bacterial-culture** monocultures/cocultures (Figs 1-5). The strongest evidence they provide for the proposed *B. theta* HsdA / 7α-HSD → G-7-oxoLCA link is a **Δhsd knockout in gnotobiotic mice** (Fig 7): WT mice produced 4.79 pmol/mg fecal G-7-oxoLCA on GCDCA water; Δhsd mice undetectable (p=0.0289).

This is **genetic evidence** (the *hsdA* gene is required for in vivo G-7-oxoLCA production), **not biochemical evidence** (the purified protein was never incubated with GCDCA to show it catalyzes the reaction directly). The gene being required doesn't prove the enzyme does the chemistry itself — the effect could be indirect. For the atlas's protein-to-molecule backlink to mean "this enzyme catalyzes this transformation," biochemical purification + substrate assay is needed; Lucas doesn't provide that.

The other 5 Lucas CSBAs have weaker evidence still: paper extrapolation (G-7-oxoDCA, T-7-oxoLCA, T-7-oxoDCA) or whole-cell species-level only with no accession given (G-3-oxoCA, G-3-oxoCDCA).

**No `interacting_protein_accessions` populated for any Lucas CSBA.** Each of the 6 molecule rows has a detailed `evidence_note` entry titled "Lucas 2026 enzyme attribution" recording what the paper actually shows (and does not show). Attribution survives as searchable text metadata only.

**Q8A6H3 removed from the atlas entirely.** Previously the atlas carried Q8A6H3 (B. theta HsdA / BT_2086) as a base Rhea/UniProt corpus row, and an earlier pre-session audit appended the clause *"— directly transforms conjugated bile acids"* to its `protein_names` field based on Lucas's Δhsd result. That annotation claim fails our purified-enzyme rule, so the row was deleted from `data/proteins.csv` during the Lucas 2026 conservative-backlink pass. Atlas protein count: 14088 → 14087. The Glycine-7-oxoLCA `evidence_note` keeps the generic "B. theta HsdA / 7α-HSD" descriptor without naming any UniProt accession.

---

## Global reproducibility

All ingestion scripts live under the scratchpad for the session they were run in:

- Gentry: `scratchpad/build_smiles_table.py`, `scratchpad/add_gentry_molecules.py`, `scratchpad/render_gentry_structures.py`
- Mohanty: `scratchpad/mohanty_amidates.py`, `scratchpad/add_mohanty_amidates.py`
- Lamichhane: `scratchpad/lamichhane_cys_mcbas.py`, `scratchpad/add_lamichhane_cys_mcbas.py`, `scratchpad/render_lamichhane_cys_structures.py`
- Xue (annotation-only, no new rows): `scratchpad/enrich_xue_annotations.py`
- Lucas 2026 completion: `scratchpad/add_lucas_g3oxo_mcbas.py`, `scratchpad/render_lucas_g3oxo_structures.py`
- Aatsinki 2025 (annotation-only, no new rows): `scratchpad/enrich_aatsinki_annotations.py`
- Agongo 2026 (NEW class — alkamines): `scratchpad/agongo_alkamines.py`, `scratchpad/add_agongo_alkamines.py`, `scratchpad/render_agongo_alkamines.py`
- Reilly 2026 (annotation-only, no new rows): `scratchpad/enrich_reilly_annotations.py`

Each script:
1. Loads authoritative BA core SMILES from `scratchpad/final_ba_cores.json`
2. Applies amide coupling via RDKit SMARTS
3. Validates each product's mass against the source paper's reported m/z
4. Appends to `data/molecules.csv` with evidence tier + note
5. Renders 2D structure PNGs into `web/public/atlas/structures/` and updates `structures_index.json`
6. Rebuilds the atlas JSONs via `scripts/build_atlas_data.py`

Backups of molecules.csv before each ingestion are saved at `data/molecules.csv.bak_before_<paper>`.

---

## Next pass (deferred to end-of-audit cleanup)

1. **Dedupe stereochem/ionization duplicates** — rows representing the same molecule with different SMILES (e.g., Gentry L-vs-D amino acid issue, cholate anion vs acid form). Includes the 3 Lamichhane-attributed duplicates/misattributions (Gly-DCA, two Hyodeoxycholic-acid rows).

### Molecule↔protein backlink rule (session-wide policy)

After the Lucas 2026 audit, the atlas uses this strict rule for every paper's `interacting_protein_accessions` population:

- **Backlink requires PURIFIED ENZYME + SUBSTRATE biochemical confirmation.** The purified protein incubated with the candidate substrate produces the candidate product, as measured directly in vitro.
- **Genetic knockout / mutant evidence** (e.g., Δgene strain fails to make product X) shows the gene is required but does NOT prove the enzyme directly catalyzes the reaction — the effect could be indirect, or depend on a cascade of enzymes. Treated as text-level attribution only.
- **Whole-cell / whole-culture evidence** (growing a strain with a BA and detecting a product) is species-level only — does not identify the specific enzyme responsible.
- **Paper-side extrapolation** (e.g., "same enzyme probably makes these other products") does not count as biochemical confirmation.

All currently-unlinked recruited molecules keep an `evidence_note` entry titled "`{Paper} enzyme attribution`" recording which paper attributes production to which species/enzyme/gene (so the information is not lost and is text-searchable), but `interacting_protein_accessions` stays empty until biochemical confirmation surfaces.

Protein rows likewise stay unmodified by literature-recruitment passes unless the paper is the authoritative source for the protein's **existence or biochemical characterization** in the atlas, not merely a user of the already-characterized enzyme.

2. **Re-source the 10 secondary-BA rows currently attributed to Lamichhane** (3-oxo-*, iso-*, allo-*, 12-oxo-LCA, 7-oxo-LCA) — likely came in from a different section of the paper or earlier work; need to pin down the real source
3. **Enrich 105 Gentry overlap rows with Lamichhane clinical context** — for the subset also measured in Lamichhane's DIABIMMUNE cohort, append the Lamichhane DOI to `paper_url` and add an evidence_note line citing cohort detection + any significant association (CDCA-Tyr, CDCA-Ala, DCA-Val, DCA-Ile/Leu/Pro, CDCA-Ser, UDCA-Asn, CA-Pro — per MOESM1 Tables S1–S4)
4. **Recompute UMAP + clusters** across the full 844-row molecule set so newly-added compounds influence the global projection
5. **Backlink molecules to producer proteins** — for compounds where a specific enzyme is known from the source paper (e.g., Rimal's P0DXD2 as producer of MCBAs), populate `interacting_protein_accessions`
6. **Walk through remaining papers:** single-compound cites (Kvitne 2025, Kostenko 2024, Mohanty 2025 preprint)
8. **Agongo SI follow-up:** 9 additional Level 1 confirmed steroid alkamines are in Table S5 but not fully enumerated in main figures; add if SI accessed
7. **Aatsinki SI follow-up:** if Supplementary Tables S9-S13 become available, check whether any of the 11 panel-only standard BAs (DHCA, THCA, 7-oxo-HCA, ω-MCA, GHDCA, GDHCA, GHCA, TDHCA, TaMCA, TbMCA, TwMCA) have detection evidence worth recruiting
7. **Dedupe Lucas-7-oxo rows:** 4 pairs of Lucas-tagged rows (`Glycine-7-oxoLCA` etc) have ChEBI-sourced duplicates (`7-oxoglycolithocholate` CHEBI:137818 etc) with different stereochem/ionization conventions — merge in end-of-audit cleanup
