"""Molecule + protein search prototype — matching styling for both views.

Landing page gates the app: nothing below renders until "Enter the Atlas"
is clicked. "Help" returns to the landing page.

Two toggles once inside:
  · View:   Protein centric | Small molecule centric | Natural + synthetic
  · Method: Multi-column search | Class filter

Same UMAP styling for both views (cluster-colored dots, STAR SVG for
newly-recruited entries, red-ring highlights, gray-out non-matches).

Launch:
    cd /home/adsiordia/marimo_visualizer/MarimoSteroidVisualizer
    ./demo/run_search_prototype.sh
"""

import marimo

__generated_with = "0.23.15"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import altair as alt
    import ast
    import re
    import os, json
    import numpy as np
    import faiss
    from openai import OpenAI
    from dotenv import load_dotenv
    alt.data_transformers.disable_max_rows()
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
        rdkit_ok = True
    except ImportError:
        Chem = None
        Draw = None
        rdkit_ok = False
    return Chem, Draw, alt, ast, faiss, json, mo, np, pd, rdkit_ok, re


@app.cell
def _(pd):
    from pathlib import Path as _Path
    _DATA_DIR = _Path(__file__).resolve().parent.parent / "data"
    molecule_df = pd.read_csv(_DATA_DIR / "molecules.csv", low_memory=False)
    molecule_df = molecule_df.rename(columns={
        "compound_name": "Compound Name",
        "chebi_id": "ChEBI ID",
        "smiles": "SMILES",
        "umap_1": "UMAP_1",
        "umap_2": "UMAP_2",
        "cluster": "clusters",
        "is_literature_recruited": "is_new",
        "paper_url": "Paper",
        "interacting_protein_accessions": "Entry",
    })
    for _col in ("SMILES", "ChEBI ID", "Compound Name", "Paper"):
        if _col in molecule_df.columns:
            molecule_df[_col] = molecule_df[_col].fillna("").astype(str)
    if "is_new" not in molecule_df.columns:
        molecule_df["is_new"] = 0
    molecule_df["clusters"] = molecule_df["clusters"].astype(str)
    return (molecule_df,)


@app.cell
def _(pd):
    from pathlib import Path as _Path
    _DATA_DIR = _Path(__file__).resolve().parent.parent / "data"
    natsyn_df = pd.read_csv(_DATA_DIR / "natural_synthetic_steroids.csv", low_memory=False)
    for _col in ("SMILES", "ChEBI ID", "Compound Name", "Paper"):
        if _col in natsyn_df.columns:
            natsyn_df[_col] = natsyn_df[_col].fillna("").astype(str)
    if "is_new" not in natsyn_df.columns:
        natsyn_df["is_new"] = 0
    if "clusters" not in natsyn_df.columns and "cluster" in natsyn_df.columns:
        natsyn_df["clusters"] = natsyn_df["cluster"]
    if "clusters" not in natsyn_df.columns:
        natsyn_df["clusters"] = "0"
    natsyn_df["clusters"] = natsyn_df["clusters"].astype(str)
    return (natsyn_df,)


@app.cell
def _(pd):
    from pathlib import Path as _Path
    _DATA_DIR = _Path(__file__).resolve().parent.parent / "data"
    protein_df = pd.read_csv(
        _DATA_DIR / "proteins.csv",
        low_memory=False,
        dtype={"interacting_chebi_ids": str, "rhea_reactions": str},
    )
    protein_df = protein_df.rename(columns={
        "accession": "Entry",
        "entry_name": "Entry Name",
        "protein_names": "Protein names",
        "gene_names": "Gene Names",
        "organism": "Organism",
        "length_aa": "Length",
        "sequence": "Sequence",
        "ec_numbers": "reaction_ecs",
        "rhea_reactions": "Rhea ID",
        "interacting_chebi_ids": "ChEBI ID",
        "interacting_compounds": "Compound Name",
        "umap_1": "UMAP_1",
        "umap_2": "UMAP_2",
        "cluster": "clusters",
        "is_literature_recruited": "is_new",
        "paper_url": "Paper",
        "annotation": "Annotation",
        "sequence_source": "Sequence_Source",
        "identifier_type": "Identifier_Type",
    })
    for _col in ("ChEBI ID", "Rhea ID", "SMILES", "Compound Name",
                 "Protein names", "Entry", "Entry Name", "Gene Names",
                 "Organism", "reaction_descriptions", "reaction_ecs",
                 "Annotation", "Paper",
                 "go_ids", "go_labels", "keyword_ids", "keyword_labels",
                 "binder_evidence", "audit_decision", "audit_reason"):
        if _col in protein_df.columns:
            protein_df[_col] = protein_df[_col].fillna("").astype(str)
    if "is_new" not in protein_df.columns:
        protein_df["is_new"] = 0
    protein_df["clusters"] = protein_df["clusters"].astype(str)
    return (protein_df,)


@app.cell
def _(Chem, Draw, molecule_df, rdkit_ok):
    structure_cache = {}
    if rdkit_ok:
        import io, base64
        for _, _row in molecule_df.iterrows():
            _name = str(_row.get("Compound Name", "")).strip()
            _chebi = str(_row.get("ChEBI ID", "")).strip()
            _sm = str(_row.get("SMILES", "")).strip()
            if not _sm:
                continue
            _mol = Chem.MolFromSmiles(_sm)
            if _mol is None:
                continue
            try:
                _img = Draw.MolToImage(_mol, size=(320, 320))
                _buf = io.BytesIO(); _img.save(_buf, format="PNG", optimize=True)
                _b64 = base64.b64encode(_buf.getvalue()).decode()
                # Store under BOTH Compound Name (if any) and ChEBI ID
                if _name:
                    structure_cache[_name.lower()] = _b64
                if _chebi:
                    structure_cache[f"chebi:{_chebi}"] = _b64
            except Exception:
                continue
    return (structure_cache,)


@app.cell
def _(mo):
    # Persistent gate. mo.state survives reruns, unlike run_button.value
    # which is a one-shot trigger.
    get_entered, set_entered = mo.state(False)
    return get_entered, set_entered


@app.cell
def _(mo, set_entered):
    # Depends only on set_entered, so clicking never rebuilds these buttons.
    enter_button = mo.ui.button(
        label="Enter the Atlas", full_width=True,
        on_change=lambda _: set_entered(True),
    )
    help_button = mo.ui.button(
        label="Help", on_change=lambda _: set_entered(False),
    )
    return enter_button, help_button


@app.cell(hide_code=True)
def _(mo, molecule_df, natsyn_df, protein_df):
    intro_text = mo.Html(f"""
    <div style="text-align:center; margin-top:10px;">

    <span style="font-size:34px; font-weight:600;">
    🧭 Welcome to Nature's Steroid Atlas
    </span>

    <p style="max-width:750px; margin:12px auto 0; font-size:15px;">
    Explore Nature's steroids and the proteins they interact with — all
    through an interactive, unified interface.
    </p>

    <p style="font-size:13px; color:#6B7280; margin:6px auto 0;">
    {len(molecule_df):,} molecules · {len(protein_df):,} proteins ·
    {len(natsyn_df):,} natural + synthetic entries
    </p>

    <hr style="width:60%; margin:16px auto;">

    <h3 style="margin:0; text-align:left; text-decoration:underline;">Views Available</h3>

    <div style="text-align:left; max-width:1000px; margin:10px auto; font-size:15px;">

      <b>• Protein-centric view →</b> Explore the protein landscape.<br>
      <span style="margin-left:18px; display:inline-block;">Each point is a
      steroid-binding protein; selecting one reveals associated steroids.</span>
      <br><br>

      <b>• Steroid-centric view →</b> Explore the chemical space of steroids.<br>
      <span style="margin-left:18px; display:inline-block;">Each point is a
      steroid molecule; selecting one highlights interacting proteins.</span>
      <br><br>

      <b>• Natural and Synthetic steroid view →</b> Compare known natural
      steroids against synthetic ones.<br>
      <span style="margin-left:18px; display:inline-block;">Each point is a
      steroid; selecting one reveals associated proteins where available.</span>

      <h3 style="margin-top:20px; text-decoration:underline;">How to Use It</h3>

      <p style="font-size:15px;">
      Pick a view, then narrow the map either by free-text search (name, ChEBI
      ID, SMILES fragment, UniProt accession, gene, sequence) or by structural
      or EC class. Matches gain a <b style="color:#e0144c;">red ring</b>;
      everything else greys out. <b>Stars</b> mark entries newly recruited from
      2024–2026 literature. Click a point, or drag a rectangle for many, and a
      results table appears. Tick rows there to see structures, identifiers
      (ChEBI, UniProt, Rhea) and AlphaFold models. A built-in language model
      companion, trained on scraped information for small molecules, explains 
      any questions asked, once an API key is entered.
      </p>

      <p style="font-size:13px; opacity:0.7;">
      Tip: scroll to zoom, double-click to reset, click pan to move around.
      </p>

    </div>
    </div>
    """)
    return (intro_text,)


@app.cell
def _(mo):
    # Own cell, ungated: the landing cell's mo.stop() would otherwise destroy
    # it on entry and break every cell downstream.
    # .form() holds .value at None until Submit, so nothing reruns per keystroke.
    key_form = mo.ui.text(
        label="OpenAI API key (enables the chat)",
        placeholder="sk-...",
        kind="password",
        full_width=True,
    ).form(submit_button_label="Save & enable chat")
    return (key_form,)


@app.cell
def _(enter_button, get_entered, intro_text, key_form, key_status, mo):
    mo.stop(get_entered())          # already inside → render nothing
    mo.vstack([intro_text, key_form, mo.md(key_status), enter_button], align="center")
    return


@app.cell
def _(index_ok, key_form):
    import os as _os
    from pathlib import Path as _P

    try:
        _ROOT = _P(__file__).resolve().parent.parent
    except NameError:
        _ROOT = _P.cwd()
    _ENV = _ROOT / ".env"

    try:
        from dotenv import load_dotenv as _load
        _load(dotenv_path=_ENV)
    except ImportError:
        pass

    def _looks_like_key(k):
        """Cheap shape check — catches typos before spending a network call."""
        return k.startswith("sk-") and len(k) >= 20 and " " not in k

    def _verify(k):
        """Return (client, status_word). Distinguishes bad key from no network."""
        from openai import OpenAI as _OpenAI
        _c = _OpenAI(api_key=k, max_retries=0, timeout=10.0)
        try:
            _c.models.list()                    # free, instant
            return _c, "ok"
        except Exception as _e:
            _name = type(_e).__name__
            if "Authentication" in _name or "PermissionDenied" in _name:
                return None, "rejected"
            if "Connection" in _name or "Timeout" in _name or "APIStatus" in _name:
                return _c, "offline"            # key may be fine; network isn't
            return None, "rejected"

    def _write_env(k):
        """Set OPENAI_API_KEY in .env, preserving every other line."""
        _lines = []
        if _ENV.exists():
            _lines = [l for l in _ENV.read_text(encoding="utf-8").splitlines()
                      if not l.strip().startswith("OPENAI_API_KEY")]
        _lines.append(f"OPENAI_API_KEY={k}")
        _ENV.write_text("\n".join(_lines).strip() + "\n", encoding="utf-8")
        try:
            _ENV.chmod(0o600)                   # owner-only
        except OSError:
            pass

    _typed = (key_form.value or "").strip()
    _stored = (_os.getenv("OPENAI_API_KEY") or "").strip()

    oai = None
    key_status = ("*Paste an OpenAI key above to enable the chat "
                  "([get one](https://platform.openai.com/api-keys)). "
                  "Everything else in the atlas works without it.*")

    if _typed:
        if not _looks_like_key(_typed):
            key_status = "❌ *That doesn't look like an OpenAI key — they start with `sk-`.*"
        else:
            oai, _res = _verify(_typed)
            if _res == "ok":
                try:
                    _write_env(_typed)
                    _os.environ["OPENAI_API_KEY"] = _typed
                    key_status = "✅ *Chat enabled — key saved, it'll load automatically next time.*"
                except OSError as _e:
                    key_status = (f"✅ *Chat enabled, but couldn't save to .env "
                                  f"({_e.strerror}) — you'll need to re-enter it next launch.*")
            elif _res == "offline":
                key_status = "⚠️ *Couldn't reach OpenAI to check that key — enabling it anyway.*"
            else:
                key_status = "❌ *That key was rejected. Check it and submit again.*"
    elif _stored:
        oai, _res = _verify(_stored)
        if _res == "ok":
            key_status = "✅ *Chat enabled — key loaded from `.env`.*"
        elif _res == "offline":
            key_status = "⚠️ *Couldn't reach OpenAI to check the saved key — enabling it anyway.*"
        else:
            key_status = "❌ *The key saved in `.env` was rejected. Paste a new one above.*"

    rag_ok = bool(oai) and index_ok
    if not index_ok:
        key_status = ("⚠️ *No search index found in `data/rag_store`, so the chat "
                      "is unavailable. Everything else works normally.*")
    return key_status, oai, rag_ok


@app.cell
def _(get_entered, mo, molecule_df, natsyn_df, protein_df):
    mo.stop(not get_entered())
    mo.md(f"""
    # 🧭 Nature's Steroid Atlas

    **{len(molecule_df):,} molecules** · **{len(protein_df):,} proteins** · **{len(natsyn_df):,} natural + synthetic entries**
    · newly recruited from 2024-2026 literature (stars).

    Pick a view and search method below.
    """)
    return


@app.cell
def _(get_entered, help_button, mo):
    mo.stop(not get_entered())
    view = mo.ui.radio(
        options=["Protein centric", "Small molecule centric", "Natural and Synthetic steroids"],
        value="Protein centric",
        label="View",
        inline=True,
    )
    mo.vstack([mo.hstack([help_button], justify="end"), view])
    return (view,)


@app.cell
def _(molecule_df, natsyn_df, protein_df, view):
    # Gated transitively via `view`.
    if view.value == "Protein centric":
        df = protein_df
        view_kind = "protein"
    elif view.value == "Natural and Synthetic steroids":
        df = natsyn_df
        view_kind = "molecule"
    else:
        df = molecule_df
        view_kind = "molecule"
    return df, view_kind


@app.cell
def _(df, mo):
    # Cluster picker: multiselect with all cluster ids + counts. Selecting one
    # or more clusters red-rings their points on the plot, same as name search.
    def _canon(x):
        s = str(x).strip()
        if not s or s.lower() == "nan":
            return ""
        try:
            return str(int(float(s)))
        except (ValueError, TypeError):
            return s
    if "clusters" in df.columns:
        _clean = df["clusters"].apply(_canon)
        _counts = _clean[_clean != ""].value_counts()
        try:
            _sorted_ids = sorted(_counts.index, key=lambda x: int(x))
        except (ValueError, TypeError):
            _sorted_ids = sorted(_counts.index)
        _labels = [f"cluster {c} (n={_counts[c]:,})" for c in _sorted_ids]
        cluster_pick_map = {label: c for label, c in zip(_labels, _sorted_ids)}
    else:
        _labels = []
        cluster_pick_map = {}
    cluster_pick = mo.ui.multiselect(
        options=_labels,
        label=f"Highlight cluster(s) — pick from all {len(_labels)}",
        value=[],
    )
    return cluster_pick, cluster_pick_map


@app.cell
def _(get_entered, mo):
    mo.stop(not get_entered())      # depends only on mo → needs its own gate
    method = mo.ui.radio(
        options=["① Multi-column search", "② Class filter"],
        value="① Multi-column search",
        label="Search method",
        inline=True,
    )
    method
    return (method,)


@app.cell
def _(mo):
    # Definition-only cell — renders nothing, so no gate needed.
    mol_search = mo.ui.text(
        placeholder="Try: estradiol · CHEBI:16469 · C[C@]12CC (SMILES fragment)",
        label="Free-text molecule search",
        full_width=True,
    )
    _mol_classes = [
        "None (show all)",
        "Bile acids (all)",
        "  · Primary bile acids",
        "  · Secondary bile acids",
        "  · Muricholic family",
        "  · Tauro-conjugated bile acids",
        "  · Glyco-conjugated bile acids",
        "  · Sulfate / sulfo bile acids",
        "  · Amino-acid conjugated (MCBAs)",
        "  · 3-oxo / oxidized bile acids",
        "  · 5α (allo) bile acids",
        "  · Ester / methyl bile acids",
        "Corticoids", "Estrogens", "Androgens",
        "Progestins", "Sterols", "Conjugated", "Ecdysteroids",
        "Brassinosteroids", "Backbone / substructure classes",
        "Newly recruited",
    ]
    mol_class = mo.ui.dropdown(_mol_classes, value="None (show all)", label="Structural class")
    return mol_class, mol_search


@app.cell
def _(mo):
    # Definition-only cell — renders nothing, so no gate needed.
    prot_search = mo.ui.text(
        placeholder=(
            "Try: P19410 (UniProt) · Bile salt hydrolase (name) · baiCD (gene) · "
            "MKATVL... (sequence) · GO:0005496 (GO id) · steroid binding (GO label) · "
            "KW-0754 (keyword) · CHEBI:15366 (ligand)"
        ),
        label="Free-text protein search (name · gene · UniProt · GO term · keyword · ChEBI · Rhea · sequence)",
        full_width=True,
    )
    _ec_classes = [
        "All EC classes",
        "EC 1 — Oxidoreductases",
        "EC 2 — Transferases",
        "EC 3 — Hydrolases",
        "EC 4 — Lyases",
        "EC 5 — Isomerases",
        "EC 6 — Ligases",
        "EC 7 — Translocases",
        "Newly recruited",
    ]
    ec_class = mo.ui.dropdown(_ec_classes, value="All EC classes", label="EC top-level class")
    return ec_class, prot_search


@app.cell
def _(cluster_pick, ec_class, method, mo, mol_class, mol_search, prot_search, view_kind):
    # Gated transitively via view_kind.
    if view_kind == "molecule":
        active = mol_search if method.value == "① Multi-column search" else mol_class
        display = mo.vstack([active, cluster_pick])
    else:
        active = prot_search if method.value == "① Multi-column search" else ec_class
        if method.value == "① Multi-column search":
            legend = mo.md(
                "*Searched fields (each result row shows which one matched in the "
                "**matched_in** column):* "
                "protein name · UniProt accession · entry name · gene name · organism · "
                "EC number · sequence · ChEBI ligand · Rhea reaction · compound · "
                "reaction description · GO id · GO label · UniProt keyword id · "
                "UniProt keyword · binder evidence · audit note"
            )
            display = mo.vstack([active, legend, cluster_pick])
        else:
            display = mo.vstack([active, cluster_pick])
    display
    return


@app.cell
def _(
    cluster_pick,
    cluster_pick_map,
    df,
    ec_class,
    method,
    mol_class,
    mol_search,
    pd,
    prot_search,
    re,
    view_kind,
):
    plot_df = df.copy()
    plot_df["_match"] = False
    plot_df["_selected"] = False

    if view_kind == "molecule":
        if method.value == "① Multi-column search":
            _q = (mol_search.value or "").strip().lower()
            if _q.startswith("chebi:"):
                _q = _q.replace("chebi:", "").strip()
            if _q:
                _cols = [c for c in ("Compound Name", "SMILES", "ChEBI ID", "Entry") if c in df.columns]
                _mask = pd.Series(False, index=df.index)
                for _c in _cols:
                    _mask = _mask | df[_c].astype(str).str.lower().str.contains(_q, regex=False, na=False)
                plot_df["_match"] = _mask
        else:
            # Class filter (all 20+ sub-classes)
            _SECONDARY_STEMS = (r"(?i)(?:(?<!cheno)deoxychol(?:ate|ic)|lithochol(?:ate|ic)|"
                r"ursodeoxychol(?:ate|ic)|hyodeoxychol(?:ate|ic)|"
                r"(?:^|[^a-z])(?:uro|urso)chol(?:ate|ic))"
            )
            _PRIMARY_STEMS = (
                r"(?i)(?:(?<![a-z])chol(?:ate|ic|oyl)|chenodeoxychol|"
                r"muricholi[cs]|muricholate|hyocholi[cs]|hyocholate|"
                r"(?:alpha|beta|α|β)-?muricholi[cs]|"
                r"(?:glyco|tauro)chol(?:ate|ic)|(?:glyco|tauro)chenodeoxychol)"
            )
            # _BILE_ANY = r"(?i)chol|muricho|hyocho|urocho"
            _BILE_ANY = (r"(?i)(?:chol(?:ate|ic|oyl|an)|muricho|hyocho|urocho|"
                 r"cholan[-\s]?[0-9]|chol-[0-9])")
            _patterns = {
                "Corticoids":       r"cort|dexamet|prednis|aldos|betameth",
                "Estrogens":        r"estr|estradiol|estrone|estriol|estetrol",
                "Androgens":        r"andr|testost|androst|dihydrotestos|dehydroepi",
                "Progestins":       r"progest|preg|pregnan",
                "Sterols":          r"sterol|cholesterol|ergost|sitosterol|lanost|desmoster|campester|stigmasterol",
                "Conjugated":       r"tauro|glyco|sulfate|sulfat|glucuronide|acyl-|amido",
                "Ecdysteroids":     r"ecdys|ponaster",
                "Brassinosteroids": r"brassin|teasteron|castasteron|typhaster|katasteron|cathasteron",
            }
            _v = mol_class.value
            _lc = plot_df["Compound Name"].fillna("").astype(str).str.strip().str.lower()
            _is_bile = _lc.str.contains(_BILE_ANY, regex=True, na=False)
            _is_secondary = _lc.str.contains(_SECONDARY_STEMS, regex=True, na=False)
            _is_primary = (
                _lc.str.contains(_PRIMARY_STEMS, regex=True, na=False)
                & ~_is_secondary & _is_bile
            )
            if _v == "None (show all)":
                plot_df["_match"] = False
            elif _v == "Newly recruited":
                plot_df["_match"] = (plot_df["is_new"] == 1)
            elif _v == "Backbone / substructure classes":
                plot_df["_match"] = _lc.str.startswith("a ") | _lc.str.startswith("an ")
            elif _v == "Bile acids (all)":
                plot_df["_match"] = _is_bile
            elif _v == "  · Primary bile acids":
                plot_df["_match"] = _is_primary
            elif _v == "  · Secondary bile acids":
                plot_df["_match"] = _is_secondary
            elif _v == "  · Muricholic family":
                plot_df["_match"] = _lc.str.contains(r"muricho", regex=True, na=False)
            elif _v == "  · Tauro-conjugated bile acids":
                plot_df["_match"] = _is_bile & _lc.str.contains(r"tauro", regex=True, na=False)
            elif _v == "  · Glyco-conjugated bile acids":
                plot_df["_match"] = _is_bile & _lc.str.contains(r"glyco", regex=True, na=False)
            elif _v == "  · Sulfate / sulfo bile acids":
                plot_df["_match"] = _is_bile & _lc.str.contains(r"sulfate|sulfo|sulfat", regex=True, na=False)
            elif _v == "  · Amino-acid conjugated (MCBAs)":
                plot_df["_match"] = _is_bile & _lc.str.contains(r"amido|amidate|acyl-", regex=True, na=False)
            elif _v == "  · 3-oxo / oxidized bile acids":
                plot_df["_match"] = _is_bile & _lc.str.contains(r"3-?oxo|7-?oxo|12-?oxo|dehydrocholic|dioxo", regex=True, na=False)
            elif _v == "  · 5α (allo) bile acids":
                plot_df["_match"] = _is_bile & _lc.str.contains(r"5alpha|allocholic|allochol", regex=True, na=False)
            elif _v == "  · Ester / methyl bile acids":
                plot_df["_match"] = _is_bile & _lc.str.contains(r"methyl.*cholate|-cholate.*ester|-cholic.*ester|-cholate$|-CoA", regex=True, na=False)
            else:
                plot_df["_match"] = _lc.str.contains(_patterns[_v], regex=True, na=False)
    else:
        # PROTEIN view
        if method.value == "① Multi-column search":
            _q = (prot_search.value or "").strip()
            if _q:
                # EC-aware search: "EC 1" → anchored regex "1."
                _m_ec = re.match(r"^ec[\s.]*([\d.]+)$", _q, re.IGNORECASE)
                if _m_ec:
                    _eq = _m_ec.group(1).rstrip(".")
                    _esc = re.escape(_eq)
                    if len(_eq.split(".")) >= 4:
                        _re_ec = r"(?<![\d.])" + _esc + r"(?![\d])"
                    else:
                        _re_ec = r"(?<![\d.])" + _esc + r"\."
                    _cols_ec = [c for c in ("Protein names", "reaction_ecs") if c in df.columns]
                    _mask = pd.Series(False, index=df.index)
                    for _c in _cols_ec:
                        _mask = _mask | df[_c].astype(str).str.contains(_re_ec, regex=True, na=False, case=False)
                    plot_df["_match"] = _mask
                else:
                    _cols = [c for c in ("Protein names", "Entry", "Entry Name",
                                          "Gene Names", "Organism", "reaction_ecs",
                                          "Sequence", "ChEBI ID", "Rhea ID",
                                          "Compound Name", "reaction_descriptions",
                                          "go_ids", "go_labels",
                                          "keyword_ids", "keyword_labels",
                                          "binder_evidence", "audit_reason")
                             if c in df.columns]
                    # Human-friendly labels shown in the `matched_in` results column
                    _display = {
                        "Protein names": "protein name",
                        "Entry": "UniProt accession",
                        "Entry Name": "entry name",
                        "Gene Names": "gene name",
                        "Organism": "organism",
                        "reaction_ecs": "EC number",
                        "Sequence": "sequence",
                        "ChEBI ID": "ChEBI ligand",
                        "Rhea ID": "Rhea reaction",
                        "Compound Name": "compound",
                        "reaction_descriptions": "reaction description",
                        "go_ids": "GO id",
                        "go_labels": "GO label",
                        "keyword_ids": "UniProt keyword id",
                        "keyword_labels": "UniProt keyword",
                        "binder_evidence": "binder evidence",
                        "audit_reason": "audit note",
                    }
                    _qlow = _q.lower()
                    _mask = pd.Series(False, index=df.index)
                    _col_masks = {}
                    for _c in _cols:
                        _m = df[_c].astype(str).str.lower().str.contains(_qlow, regex=False, na=False)
                        _col_masks[_c] = _m
                        _mask = _mask | _m
                    plot_df["_match"] = _mask
                    # Build a per-row "matched_in" annotation for rows that hit
                    _matched_in = pd.Series("", index=df.index)
                    if _mask.any():
                        _idx = df.index[_mask]
                        for i in _idx:
                            _hits = [_display.get(c, c) for c, m in _col_masks.items() if m.loc[i]]
                            _matched_in.loc[i] = ", ".join(_hits)
                    plot_df["matched_in"] = _matched_in
        else:
            # EC top-level class filter
            _v = ec_class.value
            if _v == "All EC classes":
                plot_df["_match"] = False
            elif _v == "Newly recruited":
                plot_df["_match"] = (plot_df["is_new"] == 1)
            else:
                _digit = _v.split(" ")[1]  # "EC 1 — ..." → "1"
                _pat = r"(?<![\d.])" + _digit + r"\."
                plot_df["_match"] = df["reaction_ecs"].astype(str).str.contains(_pat, regex=True, na=False)

    # Cluster picker adds red-ring highlight to every point in the picked cluster(s),
    # regardless of which text search / class filter mode is active.
    if cluster_pick.value:
        _picked_ids = {str(cluster_pick_map[label]) for label in cluster_pick.value
                       if label in cluster_pick_map}
        if _picked_ids and "clusters" in plot_df.columns:
            def _canon_cluster(x):
                s = str(x).strip()
                if not s or s.lower() == "nan":
                    return ""
                try:
                    return str(int(float(s)))
                except (ValueError, TypeError):
                    return s
            _plot_clean = plot_df["clusters"].apply(_canon_cluster)
            _cluster_mask = _plot_clean.isin(_picked_ids)
            plot_df["_match"] = plot_df["_match"] | _cluster_mask
            # For rows that only matched via cluster (no prior matched_in text), tag them
            if "matched_in" not in plot_df.columns:
                plot_df["matched_in"] = ""
            plot_df["matched_in"] = plot_df["matched_in"].fillna("").astype(str)
            _only_cluster = _cluster_mask & (plot_df["matched_in"] == "")
            plot_df.loc[_only_cluster, "matched_in"] = "cluster pick"
            _both = _cluster_mask & (plot_df["matched_in"] != "") & (plot_df["matched_in"] != "cluster pick")
            plot_df.loc[_both, "matched_in"] = plot_df.loc[_both, "matched_in"] + ", cluster pick"
    return (plot_df,)


@app.cell
def _(alt, mo, pan_toggle, plot_df, view_kind):
    n_hit = int(plot_df["_match"].sum())
    n_new = int((plot_df["is_new"] == 1).sum())
    n_total = len(plot_df)
    _search_active = n_hit > 0

    STAR = ("M0,-1 L0.22,-0.31 L0.95,-0.31 L0.36,0.12 L0.59,0.81 "
            "L0,0.4 L-0.59,0.81 L-0.36,0.12 L-0.95,-0.31 L-0.22,-0.31 Z")

    # Slim plot for the chart — drop heavy columns and truncate strings so the
    # Vega-Lite spec stays under the WebSocket size limit for large datasets.
    _keep = ["UMAP_1", "UMAP_2", "clusters", "is_new", "_match", "_selected"]
    if view_kind == "molecule":
        _keep += ["Compound Name", "ChEBI ID"]
    else:
        # For proteins, only keep Entry + a short protein name for the tooltip;
        # everything else recovered from plot_df in the detail panel by Entry.
        _keep += ["Entry", "Protein names", "Gene Names"]
    _chart_df = plot_df[[c for c in _keep if c in plot_df.columns]].copy()
    # Aggressively truncate strings — matters most for the 35k-row protein view.
    if view_kind == "protein":
        for _c in ("Protein names", "Gene Names"):
            if _c in _chart_df.columns:
                _chart_df[_c] = _chart_df[_c].astype(str).str.slice(0, 60)

    _chart_df["_kind"] = "Existing"
    _chart_df.loc[_chart_df["is_new"] == 1, "_kind"] = "Newly recruited"
    # Drop _selected to save spec bytes (not used in current logic)
    _chart_df = _chart_df.drop(columns=[c for c in ("_selected",) if c in _chart_df.columns])
    # Precompute size/stroke/strokeWidth — few unique values, Vega compresses well
    _chart_df["_size"] = 30
    _chart_df.loc[_chart_df["is_new"] == 1, "_size"] = 220
    _chart_df["_stroke"] = "white"
    _chart_df["_stroke_w"] = 0.3
    _chart_df.loc[_chart_df["is_new"] == 1, "_stroke"] = "black"
    _chart_df.loc[_chart_df["is_new"] == 1, "_stroke_w"] = 1.0
    _chart_df.loc[_chart_df["_match"], "_stroke"] = "#e0144c"
    _chart_df.loc[_chart_df["_match"], "_stroke_w"] = 2.0

    _shape_scale = alt.Scale(
        domain=["Existing", "Newly recruited"],
        range=["circle", STAR],
    )
    # Pick a color scale that gives high contrast for the number of clusters present.
    # For >20 clusters we generate a unique color per cluster via golden-ratio hue
    # sampling with alternating saturation/value (Altair's tableau20 recycles beyond 20).
    import colorsys as _colorsys
    def _n_distinct_colors(n):
        out = []
        for i in range(n):
            h = (i * 0.6180339887498949) % 1.0
            s = 0.55 + 0.30 * ((i % 3) / 2.0)
            v = 0.60 + 0.28 * ((i + 1) % 2)
            r, g, b = _colorsys.hsv_to_rgb(h, s, v)
            out.append(f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}")
        return out

    # Canonicalize cluster ids to clean integer strings ("0", "1", ..., "94")
    # so the color scale domain and the chart data match exactly regardless of
    # how pandas/CSV serializes floats.
    def _clean_cluster(x):
        s = str(x).strip()
        if not s or s.lower() == "nan":
            return ""
        try:
            return str(int(float(s)))
        except (ValueError, TypeError):
            return s
    _chart_df["clusters"] = _chart_df["clusters"].apply(_clean_cluster).astype(str)

    _n_clusters = int(_chart_df["clusters"].nunique())
    if _n_clusters <= 3:
        _color_scale = alt.Scale(range=["#0E7490", "#F59E0B", "#B91C1C"])
    elif _n_clusters <= 10:
        _color_scale = alt.Scale(scheme="tableau10")
    elif _n_clusters <= 20:
        _color_scale = alt.Scale(scheme="tableau20")
    else:
        # No explicit domain: Vega auto-assigns each cluster the next color from
        # `range`, giving 95 unique colors for 95 clusters. Order-of-first-occurrence
        # is stable within a single render, which is enough for interactive use.
        _color_scale = alt.Scale(range=_n_distinct_colors(_n_clusters))

    # Legend: show every cluster (no ellipsis) in a compact multi-column layout.
    _legend = alt.Legend(
        title="Cluster",
        orient="right",
        columns=4 if _n_clusters > 40 else 2,
        symbolLimit=max(200, _n_clusters + 10),
    )

    if _search_active:
        _color_enc = alt.condition(
            "datum._match",
            alt.Color("clusters:N", scale=_color_scale, legend=_legend),
            alt.value("#D1D5DB"),
        )
    else:
        _color_enc = alt.Color("clusters:N", scale=_color_scale, legend=_legend)
    _size_enc = alt.Size("_size:Q", scale=None, legend=None)
    _stroke_enc = alt.Stroke("_stroke:N", scale=None, legend=None)
    _stroke_w_enc = alt.StrokeWidth("_stroke_w:Q", scale=None, legend=None)

    _tooltip = (["Compound Name:N", "ChEBI ID:N", "clusters:N", "is_new:N"]
                if view_kind == "molecule"
                else ["Entry:N", "Protein names:N", "Gene Names:N", "clusters:N", "is_new:N"])
    _tooltip = [t for t in _tooltip if t.split(":")[0] in _chart_df.columns]

    chart_raw = (
        alt.Chart(_chart_df)
        .mark_point(filled=True, opacity=0.85)
        .encode(
            x=alt.X("UMAP_1:Q", axis=None, scale=alt.Scale(zero=False)),
            y=alt.Y("UMAP_2:Q", axis=None, scale=alt.Scale(zero=False)),
            color=_color_enc,
            shape=alt.Shape("_kind:N", scale=_shape_scale,
                            legend=alt.Legend(title="Type", orient="right")),
            size=_size_enc,
            stroke=_stroke_enc,
            strokeWidth=_stroke_w_enc,
            tooltip=_tooltip,
        )
        .add_params(
            alt.selection_point(name="pt", on="click", clear="dblclick"),
            alt.selection_interval(name="rng"),
            alt.selection_interval(name="zoom", bind="scales",
                                    translate=pan_toggle.value, zoom=True),
        )
        .properties(width=920, height=560)
    )

    chart_widget = mo.ui.altair_chart(chart_raw, chart_selection=True, legend_selection=True)
    return (chart_widget,)


@app.cell
def _(chart_widget, mo, plot_df):
    n_hit_disp = int(plot_df["_match"].sum())
    n_total_disp = len(plot_df)
    n_new_disp = int((plot_df["is_new"] == 1).sum())
    mo.vstack([
        mo.md(
            f"### **{n_hit_disp:,} highlighted** — {n_total_disp:,} total · "
            f"**{n_new_disp} newly-recruited (stars)**"
        ),
        mo.md("*Circles = existing entries · Stars = newly recruited · "
              "**Red ring** = search / class match. "
              "**Click** for one · **drag rectangle** for many · "
              "**scroll to zoom** · **double-click** to reset.*"),
        chart_widget,
    ])
    return


@app.cell
def _(chart_widget, mo, plot_df, view_kind):
    _sel = chart_widget.value
    # Build the pool of rows the user has narrowed down to
    if _sel is not None and hasattr(_sel, "__len__") and 0 < len(_sel) < len(plot_df):
        if view_kind == "molecule" and "Compound Name" in _sel.columns:
            _pool = plot_df[plot_df["Compound Name"].isin(_sel["Compound Name"].tolist())]
        elif view_kind == "protein" and "Entry" in _sel.columns:
            _pool = plot_df[plot_df["Entry"].isin(_sel["Entry"].tolist())]
        else:
            _pool = _sel
    else:
        # No chart selection → fall back to whatever the search/class filter matched
        _pool = plot_df[plot_df["_match"]] if plot_df["_match"].any() else plot_df.head(0)

    # Pick a slim set of columns for the table
    if view_kind == "molecule":
        _cols = [c for c in ("Compound Name", "ChEBI ID", "clusters", "is_new", "Paper") if c in _pool.columns]
    else:
        _cols = [c for c in ("Entry", "Protein names", "Gene Names", "Organism",
                              "reaction_ecs", "matched_in", "go_labels",
                              "keyword_labels", "audit_decision",
                              "clusters", "is_new")
                 if c in _pool.columns]
    _tbl = _pool[_cols].head(500).reset_index(drop=True) if len(_pool) else _pool[_cols].reset_index(drop=True)

    if len(_tbl) == 0:
        selection_table = None
        table_out = mo.md("")
    else:
        selection_table = mo.ui.table(_tbl, page_size=8, selection="multi")
        table_out = mo.vstack([
            mo.md(f"---\n### Results table — {len(_tbl):,} candidates"),
            selection_table,
        ])
    table_out
    return (selection_table,)


@app.cell
def _(ast, mo, plot_df, rdkit_ok, selection_table, structure_cache, view_kind):
    # Priority 1: use the TABLE selection (rows the user ticked)
    # Priority 2: if table has no selection but exists → show nothing (wait for user)
    # Priority 3: if no table at all (empty pool) → show placeholder
    _rows = None
    if selection_table is not None:
        _tbl_sel = selection_table.value
        if _tbl_sel is not None and hasattr(_tbl_sel, "__len__") and len(_tbl_sel) > 0:
            # Recover full plot_df rows by key
            if view_kind == "molecule" and "Compound Name" in _tbl_sel.columns:
                _rows = plot_df[plot_df["Compound Name"].isin(_tbl_sel["Compound Name"].tolist())]
            elif view_kind == "protein" and "Entry" in _tbl_sel.columns:
                _rows = plot_df[plot_df["Entry"].isin(_tbl_sel["Entry"].tolist())]
            else:
                _rows = _tbl_sel

    def _parse_list_string(s):
        s = str(s).strip() if s is not None else ""
        if not s or s.lower() == "nan":
            return []
        if s.startswith("[") and s.endswith("]"):
            try:
                v = ast.literal_eval(s)
                return [str(x).strip() for x in v if str(x).strip()]

            except (ValueError, SyntaxError):
                pass
        if ";" in s:
            return [x.strip() for x in s.split(";") if x.strip()]
        return [s]

    _output = None
    if _rows is None or len(_rows) == 0:
        _output = mo.md("")
    elif view_kind == "molecule":
        # ─── MOLECULE DETAIL — structure + protein list ───
        if not rdkit_ok:
            _output = mo.md("---\n⚠️ RDKit not installed.")
        else:
            _cards = []
            _shown = _rows.head(24)
            for _, _row in _shown.iterrows():
                _name_raw = str(_row.get("Compound Name", "")).strip()
                _chebi = str(_row.get("ChEBI ID", "")).strip()
                _name = _name_raw if _name_raw else (f"[ChEBI:{_chebi}]" if _chebi else "(unnamed)")
                _cluster = str(_row.get("clusters", "?"))
                _is_new = int(_row.get("is_new", 0))
                _paper = str(_row.get("Paper", ""))
                # Try name first, fall back to ChEBI lookup
                _b64 = structure_cache.get(_name_raw.lower())
                if not _b64 and _chebi:
                    _b64 = structure_cache.get(f"chebi:{_chebi}")
                _prot_entries = _parse_list_string(_row.get("Entry", ""))
                _prot_names = _parse_list_string(_row.get("Protein names", ""))
                _n = max(len(_prot_entries), len(_prot_names))
                _prot_entries += [""] * (_n - len(_prot_entries))
                _prot_names += [""] * (_n - len(_prot_names))
                _seen = set(); _paired = []
                for _pe, _pn in zip(_prot_entries, _prot_names):
                    _key = _pe or _pn
                    if _key and _key not in _seen and (_pe or _pn):
                        _seen.add(_key); _paired.append((_pe, _pn))

                _safe_name = _name.replace("<", "&lt;").replace(">", "&gt;")
                _new_badge = ('<span style="background:#FEF3C7; color:#92400E; padding:2px 8px; '
                              'border-radius:12px; font-size:10px; font-weight:700; margin-left:6px;">'
                              '★ NEW</span>') if _is_new else ""
                _paper_html = ""
                if _paper.startswith("http"):
                    _paper_short = _paper[:60] + "…" if len(_paper) > 60 else _paper
                    _paper_html = (
                        f'<div style="margin-top:8px; font-size:11px;">'
                        f'<strong>Source:</strong> <a href="{_paper}" target="_blank" '
                        f'style="color:#0E7490; text-decoration:underline;">{_paper_short}</a></div>'
                    )
                if _b64:
                    _img_html = f'<img src="data:image/png;base64,{_b64}" width="180" style="display:block; margin:auto;"/>'
                else:
                    _img_html = ('<div style="width:180px; height:180px; background:#F3F4F6; '
                                 'display:flex; align-items:center; justify-content:center; '
                                 'color:#9CA3AF; font-size:11px; margin:auto;">(no structure)</div>')
                if _paired:
                    _prot_items = ""
                    import re as _re_id2
                    _UP = _re_id2.compile(
                        r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
                    )
                    for _i, (_pe, _pn) in enumerate(_paired, 1):
                        _pn_safe = _pn[:80].replace("<", "&lt;").replace(">", "&gt;") + ("…" if len(_pn) > 80 else "")
                        if _pe:
                            if _UP.match(_pe):
                                _link = (f'<a href="https://www.uniprot.org/uniprotkb/{_pe}/entry" target="_blank" '
                                         f'style="color:#0E7490; text-decoration:underline; font-family:monospace; font-weight:600;">{_pe}</a>')
                                _af = (f' <a href="https://alphafold.ebi.ac.uk/entry/{_pe}" target="_blank" '
                                       f'style="color:#0E7490; text-decoration:underline; font-size:9px; '
                                       f'padding:1px 4px; border:1px solid #0E7490; border-radius:3px;">AF</a>')
                            else:
                                _link = (f'<span style="font-family:monospace; font-weight:600; color:#111827;">{_pe}</span> '
                                         f'<span style="background:#FEF3C7; color:#92400E; padding:1px 5px; '
                                         f'border-radius:8px; font-size:9px; font-weight:700;">locus tag</span>')
                                _af = ""
                        else:
                            _link = "—"; _af = ""
                        _prot_items += (
                            f'<li style="margin-bottom:8px; font-size:11px; line-height:1.35;">'
                            f'<strong>{_i}.</strong> {_link}{_af}<br>'
                            f'<span style="color:#4B5563; font-size:10.5px;">{_pn_safe}</span>'
                            f'</li>'
                        )
                    _protein_html = (
                        f'<div style="margin-top:12px; padding-top:10px; border-top:1px solid #E5E7EB;">'
                        f'<div style="font-size:12px; font-weight:700; margin-bottom:8px;">'
                        f'Interacts with {len(_paired)} protein{"s" if len(_paired) != 1 else ""}:</div>'
                        f'<div style="max-height:280px; overflow-y:auto; padding:8px 10px; '
                        f'background:#FAFAF7; border:1px solid #E5E7EB; border-radius:4px;">'
                        f'<ol style="margin:0; padding-left:20px;">{_prot_items}</ol>'
                        f'</div></div>'
                    )
                else:
                    _protein_html = ""

                # Vertical layout: BIG steroid on top, proteins list underneath
                _big_img_html = (
                    (f'<img src="data:image/png;base64,{_b64}" width="280" '
                     f'style="display:block; margin:auto;"/>')
                    if _b64 else
                    ('<div style="width:280px; height:280px; background:#F3F4F6; '
                     'display:flex; align-items:center; justify-content:center; '
                     'color:#9CA3AF; font-size:11px; margin:auto;">(no structure)</div>')
                )
                _top_pane = (
                    f'<div style="padding:20px 16px 12px; text-align:center; '
                    f'border-bottom:1px solid #E5E7EB;">'
                    f'{_big_img_html}'
                    f'<div style="font-size:20px; font-weight:600; margin-top:14px; '
                    f'line-height:1.3;">{_safe_name}{_new_badge}</div>'
                    f'<div style="font-size:12px; color:#6B7280; font-family:monospace; margin-top:4px;">'
                    f'CHEBI:{_chebi} · cluster {_cluster}</div>'
                    f'{_paper_html}'
                    f'</div>'
                )
                _bottom_body = _protein_html if _protein_html else (
                    '<div style="color:#9CA3AF; font-size:12px; text-align:center; padding:24px 0;">'
                    '(no interacting proteins linked)</div>'
                )
                _bottom_pane = (
                    f'<div style="padding:14px 18px;">'
                    f'{_bottom_body}'
                    f'</div>'
                )
                _cards.append(
                    f'<div style="margin:12px 0; '
                    f'border:1px solid #E5E7EB; border-radius:10px; background:white; '
                    f'box-shadow:0 2px 6px rgba(0,0,0,0.06); max-width:520px; overflow:hidden;">'
                    f'{_top_pane}{_bottom_pane}'
                    f'</div>'
                )
            _extra = f" · (+{len(_rows) - len(_shown)} more not shown)" if len(_rows) > len(_shown) else ""
            _grid = f'<div style="display:flex; flex-direction:column;">{"".join(_cards)}</div>'
            _output = mo.md(
                f"---\n### Selected molecules — **{len(_rows)}**{_extra}\n\n{_grid}"
            )
    else:
        # ─── PROTEIN DETAIL — protein info + reaction + interacting steroid structures ───
        _cards = []
        _shown = _rows.head(12)  # fewer for protein view because content is denser
        for _, _row in _shown.iterrows():
            _entry = str(_row.get("Entry", ""))
            _entry_name = str(_row.get("Entry Name", ""))
            _prot_names = str(_row.get("Protein names", ""))
            _gene = str(_row.get("Gene Names", ""))
            _org = str(_row.get("Organism", ""))
            _length = str(_row.get("Length", ""))
            _cluster = str(_row.get("clusters", "?"))
            _is_new = int(_row.get("is_new", 0))
            _paper = str(_row.get("Paper", ""))
            _rhea = str(_row.get("Rhea ID", ""))
            _ecs = str(_row.get("reaction_ecs", ""))
            _reactions = str(_row.get("reaction_descriptions", ""))

            _mol_names = _parse_list_string(_row.get("Compound Name", ""))
            _mol_chebis = _parse_list_string(_row.get("ChEBI ID", ""))
            _n = max(len(_mol_names), len(_mol_chebis))
            _mol_names += [""] * (_n - len(_mol_names))
            _mol_chebis += [""] * (_n - len(_mol_chebis))
            _seen = set(); _molecules = []
            for _mn, _mc in zip(_mol_names, _mol_chebis):
                _key = _mn or _mc
                if _key and _key not in _seen:
                    _seen.add(_key); _molecules.append((_mn, _mc))

            _rhea_list = _parse_list_string(_rhea)
            _ec_list = [e.strip() for e in _ecs.replace(",", ";").split(";") if e.strip()]
            _reaction_list = _parse_list_string(_reactions)

            _new_badge = ('<span style="background:#FEF3C7; color:#92400E; padding:2px 8px; '
                          'border-radius:12px; font-size:10px; font-weight:700; margin-left:6px;">'
                          '★ NEW</span>') if _is_new else ""

            _paper_html = ""
            if _paper.startswith("http"):
                _paper_short = _paper[:60] + "…" if len(_paper) > 60 else _paper
                _paper_html = (
                    f'<div style="margin-top:8px; font-size:11px;">'
                    f'<strong>Source:</strong> <a href="{_paper}" target="_blank" '
                    f'style="color:#0E7490; text-decoration:underline;">{_paper_short}</a></div>'
                )

            # Links section — UniProt IDs get external links, locus tags get a
            # paper-source badge so we don't send users to a broken UniProt page.
            import re as _re_id
            _is_uniprot = bool(_re_id.match(
                r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$",
                _entry
            ))
            if _is_uniprot:
                _links_html = (
                    f'<div style="margin-top:8px;">'
                    f'<a href="https://www.uniprot.org/uniprotkb/{_entry}/entry" target="_blank" '
                    f'style="color:#0E7490; text-decoration:underline; font-family:monospace; font-weight:600;">'
                    f'{_entry}</a> · '
                    f'<a href="https://alphafold.ebi.ac.uk/entry/{_entry}" target="_blank" '
                    f'style="color:#0E7490; text-decoration:underline; font-size:11px; '
                    f'padding:2px 6px; border:1px solid #0E7490; border-radius:3px;">AlphaFold structure</a>'
                    f'</div>'
                )
            else:
                _links_html = (
                    f'<div style="margin-top:8px; display:flex; gap:8px; align-items:center; flex-wrap:wrap;">'
                    f'<span style="font-family:monospace; font-weight:600; color:#111827;">{_entry}</span>'
                    f'<span style="background:#FEF3C7; color:#92400E; padding:2px 8px; '
                    f'border-radius:12px; font-size:10px; font-weight:700; letter-spacing:0.04em;">'
                    f'LOCUS TAG · NOT IN UNIPROT</span>'
                    f'</div>'
                    f'<div style="font-size:11px; color:#6B7280; margin-top:4px;">'
                    f'Sequence sourced from the paper\'s supplementary data (see Paper link).'
                    f'</div>'
                )

            # EC + Rhea + reaction descriptions
            _ec_badges = ""
            if _ec_list:
                _ec_badges = ('<div style="margin-top:8px; font-size:11px;"><strong>EC:</strong> ' +
                              " ".join(f'<span style="font-family:monospace; padding:2px 6px; background:#EFF6FF; '
                                       f'color:#1E40AF; border-radius:3px; margin-right:4px;">{e}</span>'
                                       for e in _ec_list[:6]) + '</div>')
            _rhea_html = ""
            if _rhea_list:
                _rhea_links = " · ".join(
                    f'<a href="https://www.rhea-db.org/rhea/{r}" target="_blank" '
                    f'style="color:#0E7490; text-decoration:underline; font-family:monospace;">Rhea:{r}</a>'
                    for r in _rhea_list[:6]
                )
                _rhea_html = f'<div style="margin-top:6px; font-size:11px;"><strong>Rhea:</strong> {_rhea_links}</div>'
            _reaction_html = ""
            if _reaction_list:
                _shown_rxn = _reaction_list[:3]
                _rxn_items = "".join(
                    f'<li style="margin-bottom:4px; line-height:1.35;">{r[:200].replace("<","&lt;")}{"…" if len(r) > 200 else ""}</li>'
                    for r in _shown_rxn
                )
                _more_rxn = f" · (+{len(_reaction_list) - len(_shown_rxn)} more)" if len(_reaction_list) > len(_shown_rxn) else ""
                _reaction_html = (
                    f'<div style="margin-top:8px; font-size:11px;">'
                    f'<strong>Reactions{_more_rxn}:</strong>'
                    f'<ul style="margin:4px 0; padding-left:16px; font-size:10.5px; color:#4B5563;">{_rxn_items}</ul>'
                    f'</div>'
                )

            # Interacting molecule structures — bigger, matching molecule-view style
            _mol_thumbs = ""
            _shown_mols = _molecules[:12]
            _thumb_htmls = []
            for _mn, _mc in _shown_mols:
                _b64_mol = structure_cache.get(_mn.strip().lower()) if _mn else None
                if not _b64_mol and _mc:
                    _b64_mol = structure_cache.get(f"chebi:{_mc.strip()}")
                if _b64_mol:
                    _img = f'<img src="data:image/png;base64,{_b64_mol}" width="180" style="display:block; margin:auto;"/>'
                else:
                    _img = ('<div style="width:180px; height:180px; background:#F3F4F6; '
                            'display:flex; align-items:center; justify-content:center; '
                            'color:#9CA3AF; font-size:11px; margin:auto;">no structure</div>')
                _mn_safe = _mn[:42].replace("<", "&lt;").replace(">", "&gt;")
                _chebi_html = f'<div style="font-size:10px; color:#6B7280; font-family:monospace; margin-top:2px;">CHEBI:{_mc}</div>' if _mc else ""
                _thumb_htmls.append(
                    f'<div style="display:inline-block; margin:8px; padding:12px; text-align:center; '
                    f'border:1px solid #E5E7EB; border-radius:8px; background:white; vertical-align:top; width:200px;">'
                    f'{_img}'
                    f'<div style="font-size:12px; color:#111827; font-weight:600; margin-top:8px; line-height:1.25;">{_mn_safe}</div>'
                    f'{_chebi_html}'
                    f'</div>'
                )
            _more_mols = f' · +{len(_molecules) - len(_shown_mols)} more' if len(_molecules) > len(_shown_mols) else ""
            if _thumb_htmls:
                _mol_thumbs = (
                    f'<div style="padding-top:14px;">'
                    f'<div style="font-size:14px; font-weight:700; margin-bottom:10px; color:#0E7490;">'
                    f'⚡ Acts on {len(_molecules)} steroid{"s" if len(_molecules) != 1 else ""}{_more_mols}:</div>'
                    f'<div style="max-height:520px; overflow-y:auto; padding:4px;">{"".join(_thumb_htmls)}</div>'
                    f'</div>'
                )
            else:
                _mol_thumbs = ('<div style="color:#9CA3AF; font-size:12px; text-align:center; padding:24px 0;">'
                               '(no interacting steroids linked)</div>')

            _prot_names_safe = _prot_names.replace("<", "&lt;").replace(">", "&gt;")[:140]
            _org_safe = _org.replace("<", "&lt;").replace(">", "&gt;")[:100]

            # Top pane — protein header/info (mirror of molecule top pane)
            _top_pane_prot = (
                f'<div style="padding:20px 22px 14px; text-align:center; '
                f'border-bottom:1px solid #E5E7EB;">'
                f'<div style="font-size:20px; font-weight:600; line-height:1.3;">{_prot_names_safe}{_new_badge}</div>'
                f'<div style="font-size:12px; color:#6B7280; margin-top:6px;">'
                f'{_gene} · <em>{_org_safe}</em> · cluster {_cluster} · {_length} aa</div>'
                f'{_links_html}{_ec_badges}{_rhea_html}{_reaction_html}{_paper_html}'
                f'</div>'
            )
            _bottom_pane_prot = f'<div style="padding:14px 22px;">{_mol_thumbs}</div>'

            _cards.append(
                f'<div style="margin:12px 0; '
                f'border:1px solid #E5E7EB; border-radius:10px; background:white; '
                f'box-shadow:0 2px 6px rgba(0,0,0,0.06); max-width:640px; overflow:hidden;">'
                f'{_top_pane_prot}{_bottom_pane_prot}'
                f'</div>'
            )
        _extra = f" · (+{len(_rows) - len(_shown)} more not shown)" if len(_rows) > len(_shown) else ""
        _grid = f'<div style="display:flex; flex-wrap:wrap;">{"".join(_cards)}</div>'
        _output = mo.md(
            f"---\n### Selected proteins — **{len(_rows)}**{_extra}\n\n{_grid}"
        )
    _output
    return


@app.cell
def _(get_entered, mo):
    mo.stop(not get_entered())      # depends only on mo → needs its own gate
    pan_toggle = mo.ui.checkbox(label="Toggle Pan")
    mo.vstack([
        pan_toggle,
        mo.md("""
    ---
    **Legend:** Colored circles = existing entries · Stars = newly recruited from 2024-2026 literature · **Red ring** = search / class match.
    """),
    ])
    return (pan_toggle,)


@app.cell
def _():
    # try:
    #     _ROOT = _P(__file__).resolve().parent.parent
    #     load_dotenv(dotenv_path=_ROOT / ".env")
    #     _STORE = _ROOT / "data" / "RAG_train"
    #     rag_index = faiss.read_index(str(_STORE.parent / "rag_store/index.faiss"))
    #     rag_catalog = [json.loads(l) for l in open(_STORE.parent / "rag_store/catalog.jsonl", encoding="utf-8")]
    #     oai = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))
    #     rag_ok = True
    # except Exception as _e:
    #     rag_index = rag_catalog = oai = None
    #     rag_ok = False
    #     print("RAG disabled:", _e)
    return


@app.cell
def _(faiss, json):
    from pathlib import Path as _P

    rag_index = rag_catalog = None
    index_ok = False

    try:
        _ROOT = _P(__file__).resolve().parent.parent
        _STORE = _ROOT / "data" / "rag_store"
        rag_index = faiss.read_index(str(_STORE / "index.faiss"))
        rag_catalog = [json.loads(l) for l in
                       open(_STORE / "catalog.jsonl", encoding="utf-8")]
        index_ok = True
    except Exception as _e:
        print("Index disabled:", _e)
    return index_ok, rag_catalog, rag_index


@app.cell
def _(np, oai, rag_catalog, rag_index, rag_ok):
    MIN_SCORE = 0.30
    # CITE_SCORE = 0.60

    def retrieve(query, k=6):
        if not rag_ok or not query.strip():
            return []
        _r = oai.embeddings.create(model="text-embedding-3-large", input=[query])
        _v = np.array(_r.data[0].embedding, dtype="float32")
        _v /= np.linalg.norm(_v) + 1e-12
        _scores, _idxs = rag_index.search(_v.reshape(1, -1), k)
        hits = []
        for _s, _i in zip(_scores[0], _idxs[0]):
            if _i == -1 or _s < MIN_SCORE:
                continue
            _row = dict(rag_catalog[_i]); _row["score"] = float(_s)
            hits.append(_row)
        return hits

    return (retrieve,)


@app.cell
def _():
    # Created once, never rebuilt — this is what keeps the chat widget (and
    # your conversation history) alive as you click around the table.
    chat_ctx = {"kind": "molecule", "compounds": [], "labels": []}
    return (chat_ctx,)


@app.cell
def _(chat_ctx, plot_df, selection_table, view_kind):
    # Mutates chat_ctx in place; must NOT redefine it. Gated transitively.
    _sel = selection_table.value if selection_table is not None else None
    chat_ctx["kind"] = view_kind
    chat_ctx["compounds"] = []
    chat_ctx["labels"] = []
    if _sel is not None and len(_sel):
        if view_kind == "molecule" and "Compound Name" in _sel.columns:
            chat_ctx["compounds"] = [str(x) for x in _sel["Compound Name"].head(15)]
            chat_ctx["labels"] = chat_ctx["compounds"]
        elif view_kind == "protein" and "Entry" in _sel.columns:
            _full = plot_df[plot_df["Entry"].isin(_sel["Entry"].tolist())]
            chat_ctx["labels"] = [str(x) for x in _full["Protein names"].head(10)]
            # interacting steroids → what we actually search the corpus with
            _c = []
            for _s in _full["Compound Name"].head(10):
                _c += [p.strip(" '\"[]") for p in str(_s).split(",") if p.strip(" '\"[]")]
            chat_ctx["compounds"] = _c[:15]
    return


@app.cell
def _(chat_ctx, get_entered, mo, oai, rag_ok, retrieve):
    mo.stop(not get_entered())

    # Above this, the answer counts as corpus-backed and gets no footer.
    # (MIN_SCORE, in the retrieve cell, decides what reaches the model at all.)
    CITE_SCORE = 0.60

    _WORDS = ("select", "selected", "selection", "highlighted",
              "chosen", "picked", "these", "this one")

    _SYS_MOLECULE = (
        "You are a precise steroid-chemistry assistant. The context below comes "
        "from per-compound records (ChEBI definitions, identifiers, literature "
        "abstracts). When it genuinely answers the question, ground your answer "
        "in it and name the compound records you used. Never stretch an unrelated "
        "compound record into an answer — if the records don't address the "
        "question, ignore them and answer from your own knowledge instead. "
        "Do not add your own disclaimer about sources; that is handled for you."
    )

    _SYS_PROTEIN = (
        "You are a precise protein-biochemistry assistant for a steroid atlas. "
        "The retrieval corpus covers small molecules only, so it will rarely help "
        "with protein questions — answer those from your own knowledge, directly "
        "and substantively. If the context happens to describe a steroid the "
        "protein acts on, use that detail and name the record. Be concrete about "
        "enzyme mechanism, family, and reaction chemistry where you can, and say "
        "plainly when something is uncertain or organism-dependent. "
        "Do not add your own disclaimer about sources; that is handled for you."
    )

    def steroid_chat(messages, config=None):
        if not rag_ok:
            return ("Chat is offline — add an OpenAI key to `.env` at the repo "
                    "root and restart. Everything else in the atlas works without it.")

        _hist = [{"role": "assistant" if m.role == "assistant" else "user",
                  "content": m.content} for m in messages]
        _q = _hist[-1]["content"] if _hist else "Hello!"
        _kind = chat_ctx["kind"]

        _asks_sel = any(w in _q.lower() for w in _WORDS)
        if _asks_sel and not chat_ctx["labels"]:
            return "Tick some rows in the results table first, then ask me again."

        _sel_txt = ""
        if chat_ctx["labels"]:
            _noun = "proteins" if _kind == "protein" else "molecules"
            _sel_txt = f"\n\nSelected {_noun}: " + "; ".join(chat_ctx["labels"])
            if _kind == "protein" and chat_ctx["compounds"]:
                _sel_txt += "\nSteroids they act on: " + "; ".join(chat_ctx["compounds"])

        # Corpus is compound-centric, so search it with compound words.
        _search_q = _q + (" " + " ".join(chat_ctx["compounds"]) if _asks_sel else "")
        _hits = retrieve(_search_q, k=20)

        # Lookalike IUPAC names make similarity unreliable — if the user ticked
        # rows, float those compounds' chunks to the top.
        if chat_ctx["labels"]:
            _named = [h for h in _hits if h["paper"] in chat_ctx["labels"]]
            _rest = [h for h in _hits if h not in _named]
            _hits = (_named + _rest)[:6]
        else:
            _hits = _hits[:6]

        _top = _hits[0]["score"] if _hits else 0.0

        # Build the context whenever there ARE hits. Weak ones still go in —
        # odd phrasing can depress the score on a perfectly relevant record.
        if _hits:
            _ctx = "\n\n".join(
                f"[{i}] {h['paper']} ({h['score']:.2f}) — "
                f"{h.get('section') or '(no section)'}\n{h['text']}"
                for i, h in enumerate(_hits, 1)
            )
            if _top < CITE_SCORE:
                _ctx = ("NOTE: these records are only loosely related to the "
                        "question. Use anything genuinely relevant, ignore the "
                        "rest, and fill the gaps from your own knowledge.\n\n") + _ctx
        else:
            _ctx = "(nothing in the corpus matched)"

        _sys = _SYS_PROTEIN if _kind == "protein" else _SYS_MOLECULE

        _resp = oai.chat.completions.create(
            model="gpt-4o-mini",
            messages=[{"role": "system", "content": _sys}, *_hist[:-1],
                      {"role": "user",
                       "content": f"Context:\n{_ctx}\n\nQuestion: {_q}{_sel_txt}"}],
            temperature=0.2,
        )
        _answer = _resp.choices[0].message.content

        # Deterministic footer — three states, decided in code.
        if not _hits:
            _answer += ("\n\n---\n*From general knowledge — nothing in the "
                        "compound corpus matched.*")
        elif _top < CITE_SCORE:
            _answer += ("\n\n---\n*Loosely related corpus records were consulted; "
                        "this answer is largely general knowledge.*")
        return _answer

    chat_ui = mo.ui.chat(steroid_chat, prompts=[
        "Tell me about the selected entries",
        "What reaction does this enzyme catalyse?",
        "How does conjugation change this bile acid?",
    ])
    mo.vstack([mo.md("---\n### 💬 Atlas companion"), chat_ui])
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
