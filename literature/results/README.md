# literature/results

Per-compound SMILES + mass-validation audit tables regenerated from the build
scripts in `../scripts/`. Everything here is the output of a non-destructive
regeneration — nothing in `data/molecules.csv` is read or written by the
build scripts.

| File | Rows | Produced by |
|---|---:|---|
| `gentry_smiles_validation.tsv` | 168 | `../scripts/build_gentry_168_mcbas.py` |
| `mohanty_smiles_validation.tsv` | 17 | `../scripts/build_mohanty_17_polyamines.py` |
| `lamichhane_smiles_validation.tsv` | 5 | `../scripts/build_lamichhane_5_cys_mcbas.py` |
| `lucas_smiles_validation.tsv` | 6 | `../scripts/build_lucas_g3oxo_mcbas.py` |
| `agongo_smiles_validation.tsv` | 18 | `../scripts/build_agongo_18_alkamines.py` |
| `stereo_audit.tsv` | 826 | `../scripts/stereo_audit_molecules.py` |

Columns (per-paper validation TSVs): `compound_name`, `smiles`,
`molecular_formula`, `computed_mz_*`, `ref_mz_*`, `delta_ppm`, and a
`mass_pass_*` flag. Columns (stereo audit): `compound_name`, `inchikey`,
`flat_inchikey` (stereo-stripped), `formula`, `stereo_defined`,
`stereo_undefined`, `canon_smiles`.

To regenerate any one of these:

```bash
cd literature/scripts
python build_gentry_168_mcbas.py          # or any other build_*.py
```
