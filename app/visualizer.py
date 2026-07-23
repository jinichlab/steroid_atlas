"""Molecule + protein search prototype — matching styling for both views.

Two toggles at the top:
  · View:   Small molecule (681)   |   Protein (35k+)
  · Method: Multi-column search    |   Class filter

Same UMAP styling for both views (cluster-colored dots, STAR SVG for
newly-recruited entries, amber-ring highlights, gray-out non-matches).

Launch:
    cd /home/adsiordia/marimo_visualizer/MarimoSteroidVisualizer
    ./demo/run_search_prototype.sh
"""
import marimo

__generated_with = "0.20.1"
app = marimo.App(width="full")


# ─── Cell 1: imports + RDKit ─────────────────────────────────────────────
@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import altair as alt
    import ast
    import re
    alt.data_transformers.disable_max_rows()
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
        rdkit_ok = True
    except ImportError:
        Chem = None
        Draw = None
        rdkit_ok = False
    return Chem, Draw, alt, ast, mo, pd, rdkit_ok, re


# ─── Cell 2: load molecule dataframe ─────────────────────────────────────
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


# ─── Cell 3: load natural + synthetic steroid catalog ────────────────────
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


# ─── Cell 4: load protein dataframe ──────────────────────────────────────
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
                 "Annotation", "Paper"):
        if _col in protein_df.columns:
            protein_df[_col] = protein_df[_col].fillna("").astype(str)
    if "is_new" not in protein_df.columns:
        protein_df["is_new"] = 0
    protein_df["clusters"] = protein_df["clusters"].astype(str)
    return (protein_df,)


# ─── Cell 4: precompute molecule 2D structure cache ──────────────────────
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


# ─── Cell 5: simple header ──────────────────────────────────────────────
@app.cell
def _(mo, molecule_df, natsyn_df, protein_df):
    mo.md(
        f"""
# 🧭 Nature's Steroid Atlas

**{len(molecule_df):,} molecules** · **{len(protein_df):,} proteins** · **{len(natsyn_df):,} natural + synthetic entries**
· newly recruited from 2024-2026 literature (stars).

Pick a view and search method below.
"""
    )
    return


# ─── Cell 6: view selector — simple radio ───────────────────────────────
@app.cell
def _(mo):
    view = mo.ui.radio(
        options=["Protein centric", "Small molecule centric", "Natural and Synthetic steroids"],
        value="Protein centric",
        label="View",
        inline=True,
    )
    view
    return (view,)


# ─── Cell 7: active dataframe + view kind ────────────────────────────────
@app.cell
def _(molecule_df, natsyn_df, protein_df, view):
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


# ─── Cell 8: method selector ─────────────────────────────────────────────
@app.cell
def _(mo):
    method = mo.ui.radio(
        options=["① Multi-column search", "② Class filter"],
        value="① Multi-column search",
        label="Search method",
        inline=True,
    )
    method
    return (method,)


# ─── Cell 9: search widgets — molecule ───────────────────────────────────
@app.cell
def _(mo):
    mol_search = mo.ui.text(
        placeholder="Try: estradiol · CHEBI:16469 · C[C@]12CC (SMILES fragment)",
        label="Free-text molecule search",
        full_width=True,
    )
    _mol_classes = [
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
    mol_class = mo.ui.dropdown(_mol_classes, value="Bile acids (all)", label="Structural class")
    return mol_class, mol_search


# ─── Cell 10: search widgets — protein ───────────────────────────────────
@app.cell
def _(mo):
    prot_search = mo.ui.text(
        placeholder="Try: P19410 (UniProt) · Bile salt hydrolase (name) · baiCD (gene) · MKATVL... (sequence)",
        label="Free-text protein search",
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


# ─── Cell 11: display active widget ──────────────────────────────────────
@app.cell
def _(ec_class, method, mo, mol_class, mol_search, prot_search, view_kind):
    if view_kind == "molecule":
        active = mol_search if method.value == "① Multi-column search" else mol_class
    else:
        active = prot_search if method.value == "① Multi-column search" else ec_class
    active
    return


# ─── Cell 12: compute match state ────────────────────────────────────────
@app.cell
def _(df, ec_class, method, mol_class, mol_search, pd, prot_search, re, view_kind):
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
            _SECONDARY_STEMS = (
                r"(?i)(?:(?<!cheno)deoxychol(?:ate|ic)|lithochol(?:ate|ic)|"
                r"ursodeoxychol(?:ate|ic)|hyodeoxychol(?:ate|ic)|"
                r"(?:^|[^a-z])(?:uro|urso)chol(?:ate|ic))"
            )
            _PRIMARY_STEMS = (
                r"(?i)(?:(?<![a-z])chol(?:ate|ic|oyl)|chenodeoxychol|"
                r"muricholi[cs]|muricholate|hyocholi[cs]|hyocholate|"
                r"(?:alpha|beta|α|β)-?muricholi[cs]|"
                r"(?:glyco|tauro)chol(?:ate|ic)|(?:glyco|tauro)chenodeoxychol)"
            )
            _BILE_ANY = r"(?i)chol|muricho|hyocho|urocho"
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
            if _v == "Newly recruited":
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
                                          "Sequence") if c in df.columns]
                    _qlow = _q.lower()
                    _mask = pd.Series(False, index=df.index)
                    for _c in _cols:
                        _mask = _mask | df[_c].astype(str).str.lower().str.contains(_qlow, regex=False, na=False)
                    plot_df["_match"] = _mask
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

    return (plot_df,)


# ─── Cell 13: BIG UMAP (view-aware) ──────────────────────────────────────
@app.cell
def _(alt, mo, plot_df, view_kind):
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
    _n_clusters = int(_chart_df["clusters"].nunique())
    if _n_clusters <= 3:
        # Hand-picked contrasting scheme for 2-3 categories (used by natural vs synthetic view)
        _color_scale = alt.Scale(range=["#0E7490", "#F59E0B", "#B91C1C"])
    elif _n_clusters <= 10:
        _color_scale = alt.Scale(scheme="tableau10")
    else:
        _color_scale = alt.Scale(scheme="tableau20")

    if _search_active:
        _color_enc = alt.condition(
            "datum._match",
            alt.Color("clusters:N",
                      scale=_color_scale,
                      legend=alt.Legend(title="Cluster", orient="right", columns=2, symbolLimit=40)),
            alt.value("#D1D5DB"),
        )
    else:
        _color_enc = alt.Color(
            "clusters:N",
            scale=_color_scale,
            legend=alt.Legend(title="Cluster", orient="right", columns=2, symbolLimit=40),
        )
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
                                    translate=False, zoom=True),
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


# ─── Cell 14: SELECTED-ROWS TABLE (drag → table → pick rows to expand) ──
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
                              "reaction_ecs", "clusters", "is_new") if c in _pool.columns]
    _tbl = _pool[_cols].head(500).reset_index(drop=True) if len(_pool) else _pool[_cols].reset_index(drop=True)

    if len(_tbl) == 0:
        selection_table = None
        table_out = mo.md("")
    else:
        selection_table = mo.ui.table(_tbl, page_size=15, selection="multi")
        table_out = mo.vstack([
            mo.md(f"---\n### Results table — {len(_tbl):,} candidates"),
            selection_table,
        ])
    table_out
    return (selection_table,)


# ─── Cell 15: DETAIL PANEL — view-aware ──────────────────────────────────
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


# ─── Cell 16: legend/footer ──────────────────────────────────────────────
@app.cell
def _(mo):
    mo.md(
        """
        ---
        **Legend:** Colored circles = existing entries · Stars = newly recruited from 2024-2026 literature · **Red ring** = search / class match.
        """
    )
    return
