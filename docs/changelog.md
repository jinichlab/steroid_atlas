# Changelog

## v0.1.0 — 2026-07-23 — Initial public release

### Literature audit pass
- Verified paper attribution for every `is_literature_recruited=1` entry against source-paper text via RAG search
- Removed `A0A8F5DVT9` (wrong Guzior 2024 attribution for *C. minuta* — paper doesn't cover this species)
- Upgraded `P0DXD2` (Rimal 2024 BlBSH) annotation from garbage `"5"` to full citation + role description
- Cleared stale "should be dropped" note from `Q5LF84` (Rimal 2024 B. fragilis BSH) — verified as legitimate KO/complementation host
- Added `WP_243289361` (Guzior 2024 CpBSH/T) — the only Guzior 2024 protein with direct biochemical evidence; ProtT5-embedded and UMAP-projected via cosine-nearest-neighbor (nearest = P54965, cosine sim = 0.9994)

### Strict EC-number curation policy applied
- `P0DXD2`: `3.5.1.24;3.5.1.74` → `2.3.1.-;3.5.1.-` (only Rimal-experimentally-confirmed ECs)
- `Q5LF84`: blank → `2.3.1.-` (KO-evidence-supported)
- `WP_243289361`: `EC 3.5.1.24; EC 2.3.1.-` → `3.5.1.24;2.3.1.-` (Guzior direct kinetics; format normalized)

### ChEBI ID cleanup
- Fixed 6 wrong ChEBI IDs in `molecules.csv`:
  - `CHEBI:34964` → `11909` (isoallopregnanolone; was S-guanyl-cysteine adduct)
  - `CHEBI:1156` → `34461` (THDOC; was 2-hydroxyestrone)
  - `CHEBI:27725` → `30154` (5β-DHP; was butin)
  - `CHEBI:34958` → `16229` (epipregnanolone; was glutathione derivative)
  - `CHEBI:16718` → `1712` (pregnanolone; was reticuline alkaloid)
  - `CHEBI:34979` → `2150` (5β-DHT; was stanolone benzoate)
- Normalized format across atlas: no `CHEBI:` prefix, no `.0` float leak — pure numeric IDs
- Full audit: 594 ChEBI IDs verified against the ChEBI ontology, 0 real mismatches remaining

### Repo structure for public release
- Clean data/ folder with snake_case column names + data dictionary
- app/ with launcher script + README
- literature/ with full audit pipeline preserved
- docs/ with methodology + changelog
- MIT LICENSE, CITATION.cff, requirements.txt
