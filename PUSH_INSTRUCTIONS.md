# One-time push to `github.com/jinichlab/steroid_atlas`

Before running these:
1. On GitHub, create the empty repo: <https://github.com/organizations/jinichlab/repositories/new> → name `steroid_atlas` → **do NOT initialize with README, license, or .gitignore** (we already have those).
2. Make sure you can authenticate to jinichlab — either an SSH key added to your GitHub account, or a personal-access token (PAT) with `repo` scope.

## Push (SSH — recommended)

```bash
cd /home/adsiordia/marimo_visualizer/steroid-atlas

# Init + commit
git init
git branch -M main
git add .
git commit -m "Initial release: Steroid Atlas v0.1.0

- 35,834 protein sequences (deduplicated) with ProtT5-derived UMAP + clusters
- 681 small molecules with verified ChEBI IDs and SMILES
- 15 literature-recruited proteins from 5 audited papers (Rimal 2024, Guzior 2024,
  McCurry 2024, Jacoby 2025, Arp 2025) with full evidence-based provenance
- Interactive marimo visualizer (app/visualizer.py)
- Reproducible audit pipeline under literature/scripts/
- Curation policies: no fabricated identifiers, EC numbers only when
  experimentally confirmed by the cited paper for that specific protein"

# Add remote + push
git remote add origin git@github.com:jinichlab/steroid_atlas.git
git push -u origin main
```

## Push (HTTPS + token, if SSH not set up)

Same commands as above, except:

```bash
git remote add origin https://github.com/jinichlab/steroid_atlas.git
git push -u origin main
# You'll be prompted for username + password.
# Username: your GitHub username
# Password: a personal-access token (NOT your GitHub password)
#   Generate at https://github.com/settings/tokens with 'repo' scope.
```

## Smoke test from a fresh clone (recommended before announcing)

On another machine (or a fresh directory here), verify a cloner can actually run it:

```bash
git clone https://github.com/jinichlab/steroid_atlas.git test_clone
cd test_clone
pip install -r requirements.txt
./app/run.sh
```

If the app comes up on port 2730 and shows the three views (Protein / Molecule / Nat+Synth), you're good to announce.

## Directory name note

The local staging dir is `steroid-atlas` (hyphen). The GitHub repo will be `steroid_atlas` (underscore). Different naming conventions — GitHub uses whichever you name the remote repo. If you want the local directory renamed to match:

```bash
cd /home/adsiordia/marimo_visualizer
mv steroid-atlas steroid_atlas
```

Do that BEFORE `git init` if you want the whole history to be under the new name. After push, `git clone` will always give you a directory called `steroid_atlas/`, so it doesn't affect anything downstream.

## What's about to be pushed

Total: ~18 MB across ~54 files.

```
steroid_atlas/
├── README.md · LICENSE · CITATION.cff · requirements.txt · .gitignore
├── data/                                    # 4 CSVs + 1 markdown + data dictionary
│   ├── proteins.csv                          (~16 MB — 35,834 proteins × 21 cols)
│   ├── molecules.csv                         (~3 MB — 681 molecules × 10 cols)
│   ├── natural_synthetic_steroids.csv        (~3.5 MB — 2,889 rows × 7 cols)
│   ├── literature_recruited_proteins.csv     (7 KB — 15 rows × 12 cols)
│   ├── literature_recruited_proteins.md
│   └── README.md
├── app/                                     # visualizer + launcher
│   ├── visualizer.py
│   ├── run.sh
│   └── README.md
├── literature/                              # reproducible audit pipeline
│   ├── scripts/                              (11 Python scripts)
│   ├── supplementary/                        (Guzior SI xlsx)
│   ├── sequences/                            (19 FASTAs)
│   ├── embeddings/                           (projection provenance JSON)
│   └── README.md
└── docs/                                    # methodology + changelog
```

Nothing sensitive: no API keys, no credentials, no unreleased data. Everything traces to public databases or open-access papers.

## After the first push

For subsequent changes:

```bash
cd /home/adsiordia/marimo_visualizer/steroid-atlas   # (or _atlas if renamed)
git add <files>
git commit -m "Describe the change"
git push
```
