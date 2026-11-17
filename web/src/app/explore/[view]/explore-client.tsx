"use client";

import { useEffect, useMemo, useState } from "react";
import Link from "next/link";
import dynamic from "next/dynamic";
import type {
  ViewKind,
  Protein,
  Molecule,
  NatsynEntry,
  ClusterMeta,
} from "@/lib/types";
import type { UmapPoint, ClusterCentroidLabel } from "@/components/UmapPlot";
import { clusterPalette, distinctPaletteForCentroids } from "@/lib/palette";

const UmapPlot = dynamic(() => import("@/components/UmapPlot"), { ssr: false });

type Row = Protein | Molecule | NatsynEntry;

interface Props {
  kind: ViewKind;
  clusters: ClusterMeta[];
}

const DATA_URL: Record<ViewKind, string> = {
  protein: "/atlas/proteins.json",
  molecule: "/atlas/molecules.json",
  natsyn: "/atlas/natsyn.json",
};
const VIEW_LABELS: Record<ViewKind, string> = {
  protein: "Protein centric",
  molecule: "Steroid centric",
  natsyn: "Natural + synthetic",
};

function rowKey(kind: ViewKind, r: Row): string {
  if (kind === "protein") return (r as Protein).accession;
  if (kind === "molecule") return (r as Molecule).compound_name;
  return `${(r as NatsynEntry).compound_name}::${(r as NatsynEntry).chebi_id}`;
}
function rowLabel(kind: ViewKind, r: Row): string {
  if (kind === "protein") {
    const p = r as Protein;
    const nm = (p.protein_names ?? "").slice(0, 80);
    const org = (p.organism ?? "").slice(0, 40);
    const acc = p.accession ? ` [${p.accession}]` : "";
    return org ? `${nm}${acc} · ${org}` : `${nm}${acc}`;
  }
  return ((r as Molecule).compound_name ?? "").slice(0, 80);
}

function esc(s: string): string {
  return String(s ?? "")
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;");
}

/** Marimo-style multi-line HTML tooltip: name / accession · gene · organism /
 *  EC + ChEBI counts + cluster + newly-recruited badge. */
function tooltipHtmlFor(kind: ViewKind, r: Row): string {
  if (kind === "protein") {
    const p = r as Protein;
    const name = esc((p.protein_names ?? "").slice(0, 90));
    const gene = esc(p.gene_names ?? "");
    const org = esc(p.organism ?? "");
    const acc = esc(p.accession ?? "");
    const ecs = esc((p.ec_numbers ?? "").split(";").slice(0, 2).join(" · "));
    const cluster = ((p.cluster ?? 0) + 1).toString();
    const star = p.is_literature_recruited
      ? `<span style="color:#f59e0b;font-weight:700;">★ NEW</span> `
      : "";
    return `
      <div style="font-size:12px;line-height:1.35;max-width:340px;">
        <div style="font-weight:600;color:currentColor;">${name}</div>
        <div style="margin-top:3px;font-size:11px;opacity:.85;">
          ${star}<span style="font-family:ui-monospace,monospace;">${acc}</span>
          ${gene ? ` · ${gene}` : ""}
        </div>
        <div style="margin-top:2px;font-size:11px;opacity:.75;">
          <em>${org}</em>
        </div>
        <div style="margin-top:2px;font-size:11px;opacity:.75;">
          Cluster ${cluster}${ecs ? ` · EC ${ecs}` : ""}
        </div>
      </div>`;
  }
  const m = r as Molecule;
  const name = esc((m.compound_name ?? "").slice(0, 90));
  const chebi = esc(m.chebi_id ?? "");
  return `
    <div style="font-size:12px;line-height:1.35;max-width:320px;">
      <div style="font-weight:600;">${name}</div>
      ${chebi ? `<div style="margin-top:3px;font-family:ui-monospace,monospace;font-size:11px;opacity:.75;">CHEBI:${chebi}</div>` : ""}
    </div>`;
}
function splitList(s: string): string[] {
  if (!s) return [];
  let txt = String(s).trim();
  // Defense against Python-list leakage ("['A', 'B', ...]"): strip brackets +
  // quotes, convert commas to semicolons, THEN apply the normal split.
  if (txt.startsWith("[") && txt.endsWith("]")) {
    txt = txt.slice(1, -1);
    if (txt.includes("'") || txt.includes('"') || txt.includes(",")) {
      txt = txt.replace(/['"]/g, "").replace(/,/g, ";");
    }
  }
  return Array.from(
    new Set(
      txt
        .replace(/;/g, "\n")
        .split("\n")
        .map((x) => x.trim())
        .filter(Boolean),
    ),
  );
}

type StructIndex = Record<string, string>;

export default function ExploreClient({ kind, clusters }: Props) {
  const [rows, setRows] = useState<Row[]>([]);
  const [loading, setLoading] = useState(true);
  const [structIdx, setStructIdx] = useState<StructIndex>({});
  // Multi-select cluster picker — click any cluster row to toggle it in
  // or out of the highlight set. Empty set = no highlight (show all).
  const [highlightSet, setHighlightSet] = useState<Set<number>>(new Set());
  const toggleClusterHighlight = (cid: number) =>
    setHighlightSet((cur) => {
      const next = new Set(cur);
      if (next.has(cid)) next.delete(cid);
      else next.add(cid);
      return next;
    });
  const clearClusterHighlight = () => setHighlightSet(new Set());
  // Convenience: single scalar for the old code paths (first selection,
  // used for cluster-label centroid display). Empty set → null.
  const highlight = highlightSet.size ? Array.from(highlightSet)[0] : null;
  const setHighlight = (v: number | null) => {
    if (v == null) setHighlightSet(new Set());
    else setHighlightSet(new Set([v]));
  };
  const [selectionSet, setSelectionSet] = useState<Set<string> | null>(null);
  const [focusedPoint, setFocusedPoint] = useState<Row | null>(null);
  const [focusedCompound, setFocusedCompound] = useState<string | null>(null);
  // Multi-select for the bottom detail cards — clicking a row in the
  // per-compound list toggles that protein in the set; the "See full
  // detail" button on the focused-point card adds that one protein too.
  const [focusedProteins, setFocusedProteins] = useState<Protein[]>([]);
  const toggleProtein = (p: Protein) =>
    setFocusedProteins((cur) => {
      const i = cur.findIndex((x) => x.accession === p.accession);
      if (i >= 0) return [...cur.slice(0, i), ...cur.slice(i + 1)];
      return [...cur, p];
    });
  const isProteinFocused = (acc: string) =>
    focusedProteins.some((x) => x.accession === acc);
  const clearFocusedProteins = () => setFocusedProteins([]);

  // Same multi-select pattern for molecules in the steroid-centric view.
  const [focusedMolecules, setFocusedMolecules] = useState<Molecule[]>([]);
  // About-the-atlas collapsible legend (molecule / natsyn views only).
  const [aboutOpen, setAboutOpen] = useState(false);
  const molKey = (m: Molecule) => `${m.compound_name}|${m.chebi_id ?? ""}`;
  const toggleMolecule = (m: Molecule) =>
    setFocusedMolecules((cur) => {
      const k = molKey(m);
      const i = cur.findIndex((x) => molKey(x) === k);
      if (i >= 0) return [...cur.slice(0, i), ...cur.slice(i + 1)];
      return [...cur, m];
    });
  const isMoleculeFocused = (m: Molecule) =>
    focusedMolecules.some((x) => molKey(x) === molKey(m));
  const clearFocusedMolecules = () => setFocusedMolecules([]);
  const [legendOpen, setLegendOpen] = useState(true);
  // Free-text search across names / accessions / gene / GO / EC / keywords /
  // ChEBI / compound / sequence.
  const [search, setSearch] = useState("");
  // A bumping counter that tells UmapPlot to erase its sticky rectangle
  // when the user hits the top-of-plot "clear" button.
  const [lassoClearSignal, setLassoClearSignal] = useState(0);
  // Steroid filter — protein view only. When a steroid is picked from
  // the autocomplete, we highlight proteins whose interacting_chebi_ids
  // contain that ChEBI id. Populated from molecules.json.
  const [steroidList, setSteroidList] = useState<
    { name: string; chebi: string; smiles: string }[]
  >([]);
  // Compound-name → all its synonyms (lowercased). Used by the protein-view
  // search to expand queries: typing "Trp-CA" finds proteins whose
  // interacting_compounds contains "L-tryptophanocholate".
  const [compoundSyn, setCompoundSyn] = useState<Map<string, string[]>>(
    new Map(),
  );
  const [steroidQuery, setSteroidQuery] = useState("");
  const [pickedSteroid, setPickedSteroid] = useState<{
    name: string;
    chebi: string;
    smiles: string;
  } | null>(null);
  const [steroidDropdownOpen, setSteroidDropdownOpen] = useState(false);

  useEffect(() => {
    let alive = true;
    setLoading(true);
    fetch(DATA_URL[kind])
      .then((r) => r.json())
      .then((d: Row[]) => {
        if (!alive) return;
        setRows(d);
        setLoading(false);
      })
      .catch(() => setLoading(false));
    return () => {
      alive = false;
    };
  }, [kind]);

  useEffect(() => {
    let alive = true;
    fetch("/atlas/structures_index.json")
      .then((r) => r.json())
      .then((d: StructIndex) => alive && setStructIdx(d))
      .catch(() => {});
    return () => {
      alive = false;
    };
  }, []);

  // Load molecule list once — used to power the steroid filter dropdown
  // in the protein view. Cheap: ~3 MB JSON, cached by the browser.
  useEffect(() => {
    if (kind !== "protein") return;
    let alive = true;
    fetch("/atlas/molecules.json")
      .then((r) => r.json())
      .then((d: Molecule[]) => {
        if (!alive) return;
        const list = d
          .map((m) => ({
            name: (m.compound_name || "").trim(),
            chebi: (m.chebi_id || "").toString().trim(),
            smiles: (m.smiles || "").trim(),
          }))
          .filter((x) => x.name && x.chebi)
          .sort((a, b) => a.name.localeCompare(b.name));
        setSteroidList(list);
        // Build synonym map: canonical compound_name (lowercase) → [all synonyms lowercased].
        // Splits on '|' AND ',' (atlas synonyms often list multiple shorthand forms
        // comma-separated within a single '|'-delimited entry, e.g. "Trp-CA, tryptophan-CA, tryptophan-cholate").
        const synMap = new Map<string, string[]>();
        for (const m of d) {
          const canon = (m.compound_name || "").toLowerCase().trim();
          if (!canon) continue;
          const raw = (m as Molecule & { synonyms?: string }).synonyms ?? "";
          const expanded: string[] = [];
          for (const chunk of String(raw).split("|")) {
            for (const sub of chunk.split(",")) {
              const t = sub.replace(/^\s*aka\s+/i, "").trim().toLowerCase();
              if (t && t !== canon) expanded.push(t);
            }
          }
          if (expanded.length) synMap.set(canon, expanded);
        }
        setCompoundSyn(synMap);
      })
      .catch(() => {});
    return () => {
      alive = false;
    };
  }, [kind]);

  const structUrl = (name: string, chebi?: string): string | null => {
    const nm = (name || "").trim().toLowerCase();
    if (nm && structIdx[nm]) return `/atlas/structures/${structIdx[nm]}`;
    const cb = (chebi || "").toString().trim();
    if (cb && structIdx[`chebi:${cb}`])
      return `/atlas/structures/${structIdx[`chebi:${cb}`]}`;
    return null;
  };

  const points: UmapPoint[] = useMemo(
    () =>
      rows.map((r) => {
        const raw = (r as any).cluster;
        const clusterNum =
          typeof raw === "number"
            ? raw
            : typeof raw === "string"
              ? raw === "natural"
                ? 0
                : 1
              : 0;
        return {
          x: (r as any).umap_1,
          y: (r as any).umap_2,
          cluster: clusterNum,
          is_star: (r as any).is_literature_recruited ?? 0,
          label: rowLabel(kind, r),
          key: rowKey(kind, r),
          tooltipHtml: tooltipHtmlFor(kind, r),
        };
      }),
    [rows, kind],
  );

  // ── Search — filters into a Set<key> that feeds the same "selection" slot
  //  the lasso uses, so search matches highlight on the map + populate the
  //  catalogue below the same way.
  const searchMatchKeys = useMemo<Set<string> | null>(() => {
    const q = search.trim().toLowerCase();
    if (!q) return null;
    // Multi-column contains-any-token match
    const tokens = q.split(/\s+/).filter(Boolean);
    const match = (r: Row): boolean => {
      const parts: string[] = [];
      if (kind === "protein") {
        const p = r as Protein;
        parts.push(
          p.accession, p.entry_name, p.protein_names, p.gene_names,
          p.organism, p.ec_numbers, p.rhea_reactions,
          p.interacting_chebi_ids, p.interacting_compounds,
          p.go_ids, p.go_labels, p.keyword_ids, p.keyword_labels,
        );
        // Expand each compound in interacting_compounds through the synonym map
        // so searches for the paper shorthand (e.g. "Trp-CA") find proteins
        // whose interacting_compounds field uses the canonical name
        // (e.g. "L-tryptophanocholate").
        if (compoundSyn.size) {
          for (const c of splitList(p.interacting_compounds ?? "")) {
            const syns = compoundSyn.get(c.toLowerCase().trim());
            if (syns) parts.push(syns.join(" "));
          }
        }
      } else if (kind === "molecule") {
        const m = r as Molecule;
        parts.push(m.compound_name, m.chebi_id, m.smiles,
                   m.interacting_protein_accessions,
                   (m as Molecule & { synonyms?: string }).synonyms);
      } else {
        const n = r as NatsynEntry;
        parts.push(n.compound_name, n.chebi_id, n.smiles, n.protein_entries,
                   (n as NatsynEntry & { synonyms?: string }).synonyms);
      }
      const blob = parts.join("\n").toLowerCase();
      return tokens.every((t) => blob.includes(t));
    };
    const s = new Set<string>();
    for (const r of rows) if (match(r)) s.add(rowKey(kind, r));
    return s;
  }, [rows, search, kind, compoundSyn]);

  const focusedPointKey = focusedPoint ? rowKey(kind, focusedPoint) : null;

  // Cluster centroids drive on-plot labels AND the distinct-color palette
  // (so neighbors on the map don't end up as similar shades of green).
  const centroids = useMemo(() => {
    const byCid = new Map<number, { sx: number; sy: number; n: number }>();
    for (const p of points) {
      const e = byCid.get(p.cluster);
      if (e) {
        e.sx += p.x;
        e.sy += p.y;
        e.n += 1;
      } else {
        byCid.set(p.cluster, { sx: p.x, sy: p.y, n: 1 });
      }
    }
    return Array.from(byCid.entries()).map(([id, v]) => ({
      id,
      cx: v.sx / v.n,
      cy: v.sy / v.n,
    }));
  }, [points]);

  const paletteMap = useMemo(() => {
    if (!centroids.length) return new Map<number, [number, number, number]>();
    return distinctPaletteForCentroids(centroids);
  }, [centroids]);

  const paletteFallback = useMemo(
    () => clusterPalette(Math.max(clusters.length, 2)),
    [clusters],
  );
  const swatchFor = (cid: number): string => {
    const m = paletteMap.get(cid);
    if (m) return `rgb(${m[0]},${m[1]},${m[2]})`;
    return paletteFallback[cid] ?? "#999";
  };

  const clusterLabels: ClusterCentroidLabel[] = useMemo(
    () =>
      clusters.map((c) => ({
        cluster: c.id,
        display_id: c.display_id,
        name: c.name || "",
        n: c.n,
      })),
    [clusters],
  );

  // ── Per-cluster protein-type classification (runtime) ──────────────────
  // For each cluster we count how many of its proteins are Enzyme /
  // Transporter / Receptor / Other, then assign a dominant type when one
  // class holds ≥60% share. Powers the sidebar sort + filter controls.
  type ProteinType =
    | "Enzyme"
    | "Transporter"
    | "Receptor"
    | "Other"
    | "Mixed";
  const clusterTypeMap = useMemo<
    Map<number, { dominant: ProteinType; counts: Record<ProteinType, number> }>
  >(() => {
    const m = new Map<
      number,
      { dominant: ProteinType; counts: Record<ProteinType, number> }
    >();
    if (kind !== "protein" || !rows.length) return m;

    const classify = (r: Protein): ProteinType => {
      const ec = (r.ec_numbers ?? "").trim();
      if (ec && /\d/.test(ec)) return "Enzyme";
      const kw = (r.keyword_labels ?? "").toLowerCase();
      if (kw.includes("transport") || kw.includes("ion channel") ||
          kw.includes("symporter") || kw.includes("antiporter")) {
        return "Transporter";
      }
      if (kw.includes("receptor")) return "Receptor";
      return "Other";
    };

    const perCluster = new Map<number, Record<ProteinType, number>>();
    for (const r of rows as Protein[]) {
      const c = Number((r as any).cluster);
      const t = classify(r);
      let e = perCluster.get(c);
      if (!e) {
        e = { Enzyme: 0, Transporter: 0, Receptor: 0, Other: 0, Mixed: 0 };
        perCluster.set(c, e);
      }
      e[t] += 1;
    }

    for (const [cid, counts] of perCluster) {
      const total = counts.Enzyme + counts.Transporter + counts.Receptor +
        counts.Other;
      let dominant: ProteinType = "Mixed";
      for (const t of ["Enzyme", "Transporter", "Receptor", "Other"] as
          ProteinType[]) {
        if (counts[t] / total >= 0.6) {
          dominant = t;
          break;
        }
      }
      m.set(cid, { dominant, counts });
    }
    return m;
  }, [rows, kind]);

  // Sidebar controls
  const [clusterSortMode, setClusterSortMode] = useState<
    "id" | "size" | "type" | "name"
  >("id");
  const [typeFilter, setTypeFilter] = useState<Set<ProteinType>>(
    new Set(["Enzyme", "Transporter", "Receptor", "Other", "Mixed"]),
  );

  const TYPE_COLOR: Record<ProteinType, string> = {
    Enzyme:      "#2563eb",
    Transporter: "#16a34a",
    Receptor:    "#ea580c",
    Other:       "#a855f7",
    Mixed:       "#64748b",
  };

  // Compute the sorted + filtered cluster list for the sidebar.
  const clustersForSidebar = useMemo(() => {
    let arr = [...clusters];
    if (kind === "protein" && clusterTypeMap.size) {
      arr = arr.filter((c) => {
        const t = clusterTypeMap.get(c.id)?.dominant ?? "Mixed";
        return typeFilter.has(t);
      });
    }
    if (clusterSortMode === "size") {
      arr.sort((a, b) => b.n - a.n);
    } else if (clusterSortMode === "name") {
      arr.sort((a, b) => (a.name || "").localeCompare(b.name || ""));
    } else if (clusterSortMode === "type") {
      const order: Record<ProteinType, number> = {
        Enzyme: 0, Transporter: 1, Receptor: 2, Other: 3, Mixed: 4,
      };
      arr.sort((a, b) => {
        const ta = clusterTypeMap.get(a.id)?.dominant ?? "Mixed";
        const tb = clusterTypeMap.get(b.id)?.dominant ?? "Mixed";
        if (order[ta] !== order[tb]) return order[ta] - order[tb];
        return b.n - a.n;
      });
    } else {
      arr.sort((a, b) => a.id - b.id);
    }
    return arr;
  }, [clusters, clusterSortMode, typeFilter, clusterTypeMap, kind]);

  // Steroid filter — for each protein, does its interacting_chebi_ids or
  // interacting_compounds mention the picked steroid? Match is by ChEBI
  // id first (precise) with a compound-name fallback (for entries that
  // list the compound by name but not chebi id).
  const steroidMatchKeys = useMemo<Set<string> | null>(() => {
    if (kind !== "protein" || !pickedSteroid) return null;
    const chebi = pickedSteroid.chebi;
    const nameLower = pickedSteroid.name.toLowerCase();
    const s = new Set<string>();
    for (const r of rows as Protein[]) {
      const chebiField = (r.interacting_chebi_ids ?? "").toString();
      const compField = (r.interacting_compounds ?? "").toLowerCase();
      // Match ChEBI id as a stand-alone token (avoid substring collisions
      // like 46 matching 46000). Look for it between non-digit boundaries.
      const chebiRe = new RegExp(`(^|[^0-9])${chebi}([^0-9]|$)`);
      if (chebiRe.test(chebiField) || compField.includes(nameLower)) {
        s.add(rowKey(kind, r));
      }
    }
    return s;
  }, [kind, pickedSteroid, rows]);

  const pool = useMemo(() => {
    if (steroidMatchKeys && steroidMatchKeys.size) {
      return rows.filter((r) => steroidMatchKeys.has(rowKey(kind, r)));
    }
    if (searchMatchKeys && searchMatchKeys.size) {
      return rows.filter((r) => searchMatchKeys.has(rowKey(kind, r)));
    }
    if (selectionSet && selectionSet.size) {
      return rows.filter((r) => selectionSet.has(rowKey(kind, r)));
    }
    if (highlightSet.size > 0) {
      return rows.filter((r) => highlightSet.has(Number((r as any).cluster)));
    }
    return [];
  }, [rows, highlightSet, selectionSet, searchMatchKeys, steroidMatchKeys, kind]);

  const compoundTiles = useMemo(() => {
    if (kind !== "protein" || !pool.length) return [];
    const byCompound = new Map<string, { chebis: Set<string>; prots: Protein[] }>();
    for (const r of pool as Protein[]) {
      const names = splitList(r.interacting_compounds ?? "");
      const chebis = splitList(r.interacting_chebi_ids ?? "");
      names.forEach((nm, i) => {
        const cb = chebis[i] ?? "";
        const entry = byCompound.get(nm);
        if (entry) {
          entry.prots.push(r);
          if (cb) entry.chebis.add(cb);
        } else {
          byCompound.set(nm, {
            prots: [r],
            chebis: new Set(cb ? [cb] : []),
          });
        }
      });
    }
    return Array.from(byCompound.entries())
      .map(([name, v]) => ({
        name,
        chebi: v.chebis.values().next().value ?? "",
        prots: v.prots,
      }))
      .sort((a, b) => b.prots.length - a.prots.length);
  }, [pool, kind]);

  const focusedCompoundProts = useMemo(() => {
    if (!focusedCompound) return [];
    const t = compoundTiles.find((t) => t.name === focusedCompound);
    return t?.prots ?? [];
  }, [compoundTiles, focusedCompound]);

  useEffect(() => {
    clearFocusedProteins();
  }, [focusedCompound, highlight]);

  return (
    <div className="min-h-screen bg-slate-50 text-slate-900">
      <div className="mx-auto max-w-[1500px] px-4 py-6">
        {/* Header */}
        <div className="mb-3 flex flex-wrap items-center gap-3">
          <Link href="/" className="text-sm text-slate-500 hover:text-slate-900">
            ← Atlas home
          </Link>
          <div className="ml-2 text-xl font-semibold tracking-tight text-slate-900">
            {VIEW_LABELS[kind]}
          </div>
          <div className="text-sm text-slate-500">
            {loading ? "loading…" : `${rows.length.toLocaleString()} entries`}
          </div>
          <div className="ml-auto flex gap-2 text-xs">
            {(["protein", "molecule", "natsyn"] as ViewKind[]).map((v) => (
              <Link
                key={v}
                href={`/explore/${v}`}
                className={`rounded-md border px-3 py-1 transition ${
                  v === kind
                    ? "border-sky-600 bg-sky-600 text-white"
                    : "border-slate-300 bg-white text-slate-700 hover:bg-slate-100"
                }`}
              >
                {VIEW_LABELS[v]}
              </Link>
            ))}
          </div>
        </div>

        {/* Search bar */}
        <div className="mb-3 flex flex-wrap items-center gap-2">
          <input
            type="text"
            placeholder={
              kind === "protein"
                ? "Search — name, gene, accession, GO, EC, ChEBI, Rhea, keyword, sequence…"
                : kind === "molecule"
                  ? "Search — compound name, ChEBI, SMILES…"
                  : "Search — compound name, ChEBI, SMILES…"
            }
            value={search}
            onChange={(e) => setSearch(e.target.value)}
            className="w-full max-w-2xl rounded-md border border-slate-300 bg-white px-3 py-1.5 text-sm shadow-sm focus:border-sky-500 focus:outline-none focus:ring-2 focus:ring-sky-100"
          />
          {search && (
            <div className="text-xs text-slate-500">
              {searchMatchKeys
                ? `${searchMatchKeys.size.toLocaleString()} matches`
                : ""}
              <button
                className="ml-2 rounded border border-slate-300 bg-white px-2 py-0.5 text-slate-700 hover:bg-slate-100"
                onClick={() => setSearch("")}
              >
                ✕ clear
              </button>
            </div>
          )}
        </div>

        {/* About the atlas — molecule / natsyn views only. Collapsible. */}
        {(kind === "molecule" || kind === "natsyn") && (
          <div className="mb-3 rounded-md border border-slate-200 bg-white shadow-sm">
            <button
              onClick={() => setAboutOpen((v) => !v)}
              className="flex w-full items-center justify-between gap-2 px-3 py-2 text-left text-sm"
            >
              <span className="flex flex-wrap items-center gap-2 text-slate-700">
                <span className="rounded bg-slate-100 px-1.5 py-0.5 text-xs font-semibold uppercase tracking-wide text-slate-600">
                  About
                </span>
                <span className="text-xs text-slate-500">
                  Tier 1 = standard · Tier 2+ = biologically detected · <span className="font-mono">Trp-CA</span> = tryptophan amide at C24 of cholate
                </span>
              </span>
              <span className="text-xs text-slate-400">{aboutOpen ? "▲ hide" : "▼ show details"}</span>
            </button>
            {aboutOpen && (
              <div className="grid grid-cols-1 gap-4 border-t border-slate-200 px-3 py-3 text-xs text-slate-700 lg:grid-cols-2">
                <div>
                  <div className="mb-1 text-sm font-semibold text-slate-800">Evidence tiers</div>
                  <p className="mb-2 italic leading-relaxed text-slate-600">
                    Evidence tiers apply only to literature-recruited rows (molecules added from
                    reverse-metabolomics papers). Pre-existing ChEBI / UniProt entries, which form the
                    base atlas, carry no tier.
                  </p>
                  <ul className="space-y-2 leading-relaxed">
                    <li>
                      <div className="font-semibold text-slate-800">Tier 1 — chemistry-only anchor</div>
                      <div>
                        Synthetic standard made for MS/MS libraries but not yet detected in any biological
                        sample. Kept as a negative-space marker for the MCBA chemistry landscape.
                      </div>
                    </li>
                    <li>
                      <div className="font-semibold text-slate-800">Tiers 2 and higher — biologically detected</div>
                      <div>
                        LC-MS/MS retention-time match to an authentic standard in at least one biological
                        sample. Higher tiers indicate stronger detection (public-data match → cohort
                        identification → cross-cohort replication → bacterial-culture production). Each
                        compound card shows the full tier label and source.
                      </div>
                    </li>
                    <li>
                      <div className="font-semibold text-slate-800">Canonical reference — no tier</div>
                      <div>
                        Pre-existing ChEBI / UniProt entries; no evidence tier applies.
                      </div>
                    </li>
                  </ul>
                </div>
                <div>
                  <div className="mb-1 text-sm font-semibold text-slate-800">Conjugate naming</div>
                  <p className="mb-2 leading-relaxed">
                    Microbially-conjugated bile acids (MCBAs) are named{" "}
                    <span className="font-mono">Residue-BileAcid</span> — the residue is the amine
                    nucleophile, the bile-acid is the steroid core.
                  </p>
                  <div className="mb-2 rounded border border-slate-200 bg-slate-50 px-3 py-2 font-mono text-[11px] leading-5 text-slate-800">
                    <div className="mb-1 text-[10px] uppercase tracking-wide text-slate-500">
                      C24 amide conjugate (most MCBAs)
                    </div>
                    <div>R-NH₂ + HO-C(=O)-[core] → R-NH-C(=O)-[core] + H₂O</div>
                    <div className="mt-2 mb-1 text-[10px] uppercase tracking-wide text-slate-500">
                      C3 alkamine conjugate (reductive amination)
                    </div>
                    <div>R-NH₂ + O=C(C3)-[core] → R-NH-CH(C3)-[core] + H₂O</div>
                  </div>
                  <ul className="space-y-1.5 leading-relaxed">
                    <li>
                      <span className="font-mono text-slate-900">Trp-CA</span> — R = tryptophan α-amine;
                      C24 amide on cholate.{" "}
                      <span className="italic">cholate + L-Trp → N-tryptophanoyl-cholate + H₂O</span>.
                    </li>
                    <li>
                      <span className="font-mono text-slate-900">Put-CA</span> — R = putrescine (polyamine
                      primary amine); C24 amide on cholate.
                    </li>
                    <li>
                      <span className="font-mono text-slate-900">Cys-CA</span> — R = cysteine α-amine; C24
                      amide on cholate (thiol stays free).
                    </li>
                    <li>
                      <span className="font-mono text-slate-900">G-3-oxoCA</span> — R = glycine; C24 amide
                      on a 3-dehydro (3-oxo) cholate core.
                    </li>
                    <li>
                      <span className="font-mono text-slate-900">Phe-C3-DCA</span> — R = phenylalanine
                      α-amine; C–N bond at the C3 ring carbon of deoxycholate (reductive amination). A
                      different bond and position from the C24 amide family above.
                    </li>
                  </ul>
                  <div className="mt-2 leading-relaxed">
                    <span className="font-semibold text-slate-800">BA cores — </span>
                    <span className="font-mono">CA</span> cholate (3α,7α,12α),{" "}
                    <span className="font-mono">CDCA</span> chenodeoxycholate (3α,7α),{" "}
                    <span className="font-mono">DCA</span> deoxycholate (3α,12α),{" "}
                    <span className="font-mono">UDCA</span> ursodeoxycholate (3α,7β),{" "}
                    <span className="font-mono">HDCA</span> hyodeoxycholate (3α,6α),{" "}
                    <span className="font-mono">α/β/γ/ω-MCA</span> muricholate isomers,{" "}
                    <span className="font-mono">LCA</span> lithocholate (3α).
                  </div>
                </div>
              </div>
            )}
          </div>
        )}

        {/* Steroid filter — protein view only. Autocomplete over the 677
            steroid molecules; picking one highlights every protein that
            interacts with it via ChEBI id or compound-name match. */}
        {kind === "protein" && (
          <div className="mb-3 flex flex-wrap items-center gap-2">
            <div className="relative">
              <input
                type="text"
                placeholder={
                  pickedSteroid
                    ? pickedSteroid.name
                    : "🧬 Filter by steroid — name, ChEBI id, or SMILES…"
                }
                value={steroidQuery}
                onFocus={() => setSteroidDropdownOpen(true)}
                onChange={(e) => {
                  setSteroidQuery(e.target.value);
                  setSteroidDropdownOpen(true);
                }}
                className="w-[520px] max-w-full rounded-md border border-amber-300 bg-amber-50 px-3 py-1.5 text-sm shadow-sm focus:border-amber-500 focus:outline-none focus:ring-2 focus:ring-amber-100"
              />
              {steroidDropdownOpen && steroidQuery.trim().length >= 2 && (
                <div className="absolute left-0 top-full z-30 mt-1 max-h-[320px] w-[520px] max-w-full overflow-y-auto rounded-md border border-slate-300 bg-white shadow-lg">
                  {(() => {
                    const raw = steroidQuery.trim();
                    const q = raw.toLowerCase();
                    const isChebi = /^chebi[:_ -]?\d+$/i.test(raw)
                      || /^\d{2,}$/.test(raw);
                    const chebiNum = isChebi
                      ? raw.replace(/[^0-9]/g, "")
                      : "";
                    const looksLikeSmiles = /[()=\[\]#\\/]/.test(raw)
                      || /[A-Z]/.test(raw);

                    // Score each match — lower is better. Sorts so:
                    // 1) exact name match  (score 0)
                    // 2) name starts with query at a word boundary (1)
                    // 3) name starts with query mid-word (2)
                    // 4) query appears as a whole word (3)
                    // 5) substring match anywhere (4)
                    // 6) chebi / smiles matches (5)
                    // Ties broken by shorter name (parent compound wins).
                    const score = (s: { name: string; smiles: string; chebi: string }) => {
                      const lower = s.name.toLowerCase();
                      if (lower === q) return 0;
                      // followed by common boundary chars
                      if (lower.startsWith(q + " ") ||
                          lower.startsWith(q + "-") ||
                          lower.startsWith(q + ",") ||
                          lower.startsWith(q + "(") ||
                          lower.startsWith(q + "/") ||
                          lower.startsWith(q + ":")) return 1;
                      if (lower.startsWith(q)) return 2;
                      // whole-word match somewhere in the name
                      const wordRe = new RegExp(
                        `(?:^|[\\s\\-,/():])${q.replace(
                          /[-\/\\^$*+?.()|[\]{}]/g, "\\$&")}(?:$|[\\s\\-,/():])`);
                      if (wordRe.test(lower)) return 3;
                      if (lower.includes(q)) return 4;
                      // ChEBI id or SMILES fallback
                      if (chebiNum && s.chebi === chebiNum) return 5;
                      if (looksLikeSmiles && s.smiles &&
                          s.smiles.includes(raw)) return 5;
                      return 999;
                    };

                    const matches = steroidList
                      .map((s) => ({ s, sc: score(s) }))
                      .filter((x) => x.sc < 999)
                      .sort((a, b) => {
                        if (a.sc !== b.sc) return a.sc - b.sc;
                        // tie-break: shorter names first (parent compound
                        // usually wins over its acetate/glucuronate)
                        return a.s.name.length - b.s.name.length;
                      })
                      .slice(0, 30)
                      .map((x) => x.s);
                    if (matches.length === 0) {
                      return (
                        <div className="px-3 py-2 text-xs text-slate-500">
                          No steroid matches "{steroidQuery}"{" "}
                          <span className="text-slate-400">
                            (tried name, ChEBI id, and SMILES substring)
                          </span>
                        </div>
                      );
                    }
                    return matches.map((s) => (
                      <button
                        key={s.chebi}
                        onClick={() => {
                          setPickedSteroid(s);
                          setSteroidQuery("");
                          setSteroidDropdownOpen(false);
                          setHighlight(null);
                          setFocusedPoint(null);
                          setFocusedCompound(null);
                          clearFocusedProteins();
                        }}
                        className="flex w-full flex-col gap-0.5 border-b border-slate-100 px-3 py-1.5 text-left text-xs hover:bg-amber-50"
                      >
                        <div className="flex items-baseline gap-2">
                          <span className="text-slate-900">{s.name}</span>
                          <span className="ml-auto font-mono text-[10px] text-slate-500">
                            CHEBI:{s.chebi}
                          </span>
                        </div>
                        {s.smiles && (
                          <div className="truncate font-mono text-[10px] text-slate-500">
                            {s.smiles}
                          </div>
                        )}
                      </button>
                    ));
                  })()}
                </div>
              )}
            </div>
            {pickedSteroid && (
              <>
                <div className="rounded-md border border-amber-500 bg-amber-100 px-2 py-1 text-xs text-amber-800">
                  🧬 <b>{pickedSteroid.name}</b>
                  <span className="ml-2 font-mono text-[10px]">
                    CHEBI:{pickedSteroid.chebi}
                  </span>
                  <span className="ml-2 font-semibold">
                    → {steroidMatchKeys?.size.toLocaleString() ?? 0}{" "}
                    interacting proteins
                  </span>
                </div>
                <button
                  className="rounded border border-slate-300 bg-white px-2 py-1 text-xs text-slate-700 hover:bg-slate-100"
                  onClick={() => setPickedSteroid(null)}
                >
                  ✕ clear steroid
                </button>
              </>
            )}
            {steroidList.length === 0 && !pickedSteroid && (
              <div className="text-xs text-slate-400">
                loading steroid list…
              </div>
            )}
          </div>
        )}

        {/* Close the steroid dropdown when clicking elsewhere on the page */}
        {steroidDropdownOpen && (
          <div
            className="fixed inset-0 z-20"
            onClick={() => setSteroidDropdownOpen(false)}
          />
        )}

        {/* Above-plot legend */}
        <div className="mb-2 flex flex-wrap items-center gap-4 text-xs text-slate-600">
          <span className="inline-flex items-center gap-1">
            <span className="inline-block h-3 w-3 rounded-full bg-cyan-500" />{" "}
            existing
          </span>
          <span className="inline-flex items-center gap-1">
            <svg width="14" height="14" viewBox="-1 -1 2 2" className="inline-block">
              <path
                d="M0,-0.95 L0.22,-0.31 L0.9,-0.31 L0.36,0.12 L0.56,0.78 L0,0.38 L-0.56,0.78 L-0.36,0.12 L-0.9,-0.31 L-0.22,-0.31 Z"
                fill="currentColor"
                className="text-cyan-500"
              />
            </svg>
            newly recruited (colored by cluster)
          </span>
          <span className="inline-flex items-center gap-1">
            <span className="inline-block h-3 w-3 rounded-full bg-cyan-500 ring-2 ring-rose-500" />{" "}
            selected dot
          </span>
          {(highlight != null || focusedPoint || selectionSet || search ||
            pickedSteroid) && (
            <button
              className="ml-auto rounded-md border border-slate-300 bg-white px-2 py-1 text-slate-700 hover:bg-slate-100"
              onClick={() => {
                setHighlight(null);
                setSelectionSet(null);
                setFocusedPoint(null);
                setFocusedCompound(null);
                clearFocusedProteins();
                setSearch("");
                setPickedSteroid(null);
                setLassoClearSignal((n) => n + 1);
              }}
            >
              ✕ clear
            </button>
          )}
        </div>

        <div className="grid gap-4 md:grid-cols-[minmax(0,1fr)_260px]">
          <div>
            <UmapPlot
              points={points}
              clusterLabels={clusterLabels}
              highlightCluster={highlightSet.size === 1 ? highlight : null}
              matchKeys={
                steroidMatchKeys && steroidMatchKeys.size
                  ? steroidMatchKeys
                  : searchMatchKeys && searchMatchKeys.size
                    ? searchMatchKeys
                    : selectionSet && selectionSet.size
                      ? selectionSet
                      : highlightSet.size > 1
                        ? new Set(
                            points
                              .filter((p) => highlightSet.has(p.cluster))
                              .map((p) => p.key!)
                              .filter(Boolean),
                          )
                        : null
              }
              focusedKey={focusedPointKey}
              darkMode={false}
              paletteOverride={paletteMap}
              clearLassoSignal={lassoClearSignal}
              onClickPoint={(p) => {
                const r = rows.find((r) => rowKey(kind, r) === p.key);
                if (r) setFocusedPoint(r);
              }}
              onLassoSelect={(pts) => {
                const s = new Set(pts.map((p) => p.key!).filter(Boolean));
                setSelectionSet(s.size ? s : null);
                setHighlight(null);
                setFocusedPoint(null);
                setFocusedCompound(null);
                clearFocusedProteins();
              }}
              height={640}
            />
            {focusedPoint && (
              <div className="mt-3 rounded-xl border-2 border-sky-500 bg-sky-50 p-4">
                <div className="text-xs uppercase tracking-wide text-sky-700">
                  {kind === "protein" ? "Selected protein" : "Selected"}
                </div>
                <div className="mt-1 text-lg font-semibold text-slate-900">
                  {kind === "protein"
                    ? (focusedPoint as Protein).protein_names
                    : (focusedPoint as Molecule).compound_name}
                </div>
                {kind === "protein" && (
                  <div className="mt-1 text-sm text-slate-600">
                    {(focusedPoint as Protein).gene_names} ·{" "}
                    <em>{(focusedPoint as Protein).organism}</em> ·{" "}
                    {(focusedPoint as Protein).length_aa} aa · Cluster{" "}
                    {((focusedPoint as Protein).cluster ?? 0) + 1}
                    {(focusedPoint as Protein).is_literature_recruited
                      ? " · ★ newly recruited"
                      : ""}
                  </div>
                )}
                {kind === "protein" && (focusedPoint as Protein).putative_role && (
                  <div className="mt-2 inline-flex items-start gap-2 rounded border border-amber-400 bg-amber-50 px-2 py-1 text-xs text-amber-900">
                    <span className="font-semibold">⚠ Putative:</span>
                    <span>{(focusedPoint as Protein).putative_role}</span>
                  </div>
                )}
                {kind === "molecule" && (() => {
                  const m = focusedPoint as Molecule;
                  const src = structUrl(m.compound_name, m.chebi_id);
                  return (
                    <>
                      <div className="mt-1 text-sm text-slate-600">
                        {m.chebi_id
                          ? <span>CHEBI:{m.chebi_id}</span>
                          : <span className="italic text-slate-500">no canonical ChEBI</span>}
                        {m.is_literature_recruited
                          ? <span className="ml-2 rounded-full bg-amber-100 px-2 py-0.5 text-[10px] font-semibold text-amber-800">★ newly recruited</span>
                          : null}
                      </div>
                      {src ? (
                        <div className="mt-2 flex items-center justify-center rounded-md border border-slate-200 bg-white p-2" style={{maxHeight: 220}}>
                          {/* eslint-disable-next-line @next/next/no-img-element */}
                          <img src={src} alt={m.compound_name} className="max-h-[200px] max-w-full object-contain" />
                        </div>
                      ) : (
                        <div className="mt-2 rounded-md border border-dashed border-slate-300 bg-slate-50 p-3 text-center text-xs italic text-slate-500">
                          (no 2D structure image available — SMILES: <span className="font-mono">{(m.smiles || "").slice(0,60)}…</span>)
                        </div>
                      )}
                    </>
                  );
                })()}
                {kind === "molecule" && (focusedPoint as Molecule).evidence_tier && (() => {
                  const tier = (focusedPoint as Molecule).evidence_tier;
                  const label =
                    tier === "1" ? "⚠ Chemistry reference only — no biological detection"
                    : tier === "2" ? "⚠ MS/MS match in public data (isomer uncertain)"
                    : tier === "3" ? "✓ Detected in human IBD stool (iHMP2, isomer-resolved)"
                    : tier === "4-5" ? "✓ Cross-cohort validated in human IBD (iHMP2 + Boland)"
                    : tier === "5" ? "✓ Produced by gut bacteria in culture"
                    : "";
                  const style =
                    tier === "1" ? "border-red-400 bg-red-50 text-red-900"
                    : tier === "2" ? "border-amber-400 bg-amber-50 text-amber-900"
                    : tier === "3" ? "border-sky-400 bg-sky-50 text-sky-900"
                    : "border-emerald-500 bg-emerald-50 text-emerald-900";
                  return (
                    <div className={`mt-2 rounded border px-2 py-1 text-xs ${style}`}>
                      <div className="font-semibold">{label}</div>
                      <div className="mt-0.5 italic text-[10px] opacity-75">Suggested classification — see source paper for primary data and methodology.</div>
                    </div>
                  );
                })()}
                {kind === "molecule" && (focusedPoint as Molecule).gentry_producers && (
                  <div className="mt-2 rounded border border-slate-200 bg-slate-50 p-2 text-xs text-slate-700">
                    <div className="font-semibold text-slate-900">Suggested bacterial producers (not definitive):</div>
                    <div className="mt-0.5">{(focusedPoint as Molecule).gentry_producers?.replace(/;/g, "; ")}</div>
                    <div className="mt-1 italic text-slate-500">From HMP isolate culture screen — the specific enzyme in each strain is not identified. See source paper for details.</div>
                  </div>
                )}
                {kind === "molecule" && (focusedPoint as Molecule).paper_url && (
                  <div className="mt-2 text-xs text-slate-600">
                    <span className="font-semibold">Source paper(s): </span>
                    {(focusedPoint as Molecule).paper_url.split("|").map((url, i, arr) => {
                      const u = url.trim();
                      return u ? (
                        <a key={i} href={u} target="_blank" className="mr-2 underline text-sky-600 hover:text-sky-800">
                          paper {arr.length > 1 ? i + 1 : ""} ↗
                        </a>
                      ) : null;
                    })}
                    <span className="italic text-slate-400"> (see paper for full methodology)</span>
                  </div>
                )}
                {kind === "molecule" && (() => {
                  const accStr = (focusedPoint as Molecule).interacting_protein_accessions ?? "";
                  const accs = splitList(accStr);
                  if (accs.length === 0) {
                    return (
                      <div className="mt-2 text-xs">
                        <span className="italic text-amber-700">⚠ No proteins associated in the atlas — this compound is listed as a literature-recruited structure without a direct enzyme attribution. For bacterial producer suggestions, see above and consult the source paper.</span>
                      </div>
                    );
                  }
                  return (
                    <div className="mt-3">
                      <div className="mb-1 text-xs font-semibold uppercase tracking-wide text-slate-500">
                        Linked proteins ({accs.length})
                      </div>
                      <div className="flex flex-wrap gap-1.5">
                        {accs.slice(0, 40).map((acc) => (
                          <a
                            key={acc}
                            href={`https://www.uniprot.org/uniprotkb/${acc}/entry`}
                            target="_blank"
                            rel="noopener noreferrer"
                            className="rounded border border-sky-300 bg-sky-50 px-2 py-0.5 font-mono text-[11px] text-sky-800 hover:border-sky-500 hover:bg-sky-100"
                            title={`Open ${acc} on UniProt`}
                          >
                            {acc}
                          </a>
                        ))}
                        {accs.length > 40 && (
                          <span className="rounded border border-slate-200 px-2 py-0.5 text-[11px] text-slate-500">
                            …and {accs.length - 40} more
                          </span>
                        )}
                      </div>
                    </div>
                  );
                })()}
                <div className="mt-2 flex flex-wrap gap-2 text-sm">
                  {kind === "protein" && (
                    <>
                      {(focusedPoint as Protein).identifier_type === "GenBank_Protein" ? (
                        <a
                          className="rounded border border-sky-500 px-2 py-1 text-sky-700 hover:bg-sky-100"
                          href={`https://www.ncbi.nlm.nih.gov/protein/${(focusedPoint as Protein).accession}.1`}
                          target="_blank"
                        >
                          NCBI Protein
                        </a>
                      ) : (focusedPoint as Protein).identifier_type === "Locus_Tag" ? (
                        <span className="rounded border border-slate-300 px-2 py-1 text-xs text-slate-500">
                          no external DB (locus tag)
                        </span>
                      ) : (
                        <>
                          <a
                            className="rounded border border-sky-500 px-2 py-1 text-sky-700 hover:bg-sky-100"
                            href={`https://www.uniprot.org/uniprotkb/${(focusedPoint as Protein).accession}/entry`}
                            target="_blank"
                          >
                            UniProt
                          </a>
                          <a
                            className="rounded border border-emerald-500 px-2 py-1 text-emerald-700 hover:bg-emerald-100"
                            href={`https://alphafold.ebi.ac.uk/entry/${(focusedPoint as Protein).accession}`}
                            target="_blank"
                          >
                            AlphaFold
                          </a>
                        </>
                      )}
                      <button
                        className="rounded border border-slate-400 px-2 py-1 text-slate-700 hover:bg-slate-100"
                        onClick={() =>
                          toggleProtein(focusedPoint as Protein)
                        }
                      >
                        {isProteinFocused((focusedPoint as Protein).accession)
                          ? "✓ In details"
                          : "＋ Add to details ↓"}
                      </button>
                    </>
                  )}
                  <button
                    className="ml-auto text-xs text-slate-500 hover:text-slate-900"
                    onClick={() => setFocusedPoint(null)}
                  >
                    clear
                  </button>
                </div>
              </div>
            )}
          </div>

          {/* Cluster sidebar with sort + filter controls */}
          <aside className="rounded-xl border border-slate-200 bg-white shadow-sm">
            <div className="flex items-center justify-between border-b border-slate-100 p-3">
              <button
                onClick={() => setLegendOpen((v) => !v)}
                className="flex-1 text-left"
              >
                <div className="text-sm font-semibold text-slate-900">
                  Clusters ({clustersForSidebar.length}
                  {clustersForSidebar.length !== clusters.length
                    ? ` / ${clusters.length}`
                    : ""})
                  {highlightSet.size > 0 && (
                    <span className="ml-2 rounded-full bg-sky-100 px-2 py-0.5 text-[10px] font-semibold text-sky-800">
                      {highlightSet.size} selected
                    </span>
                  )}
                </div>
                <div className="text-xs text-slate-500">
                  {legendOpen
                    ? highlightSet.size > 0
                      ? "click a row to add/remove from selection"
                      : "click any row to highlight (click more to add)"
                    : "click to expand"}
                </div>
              </button>
              {highlightSet.size > 0 && (
                <button
                  onClick={clearClusterHighlight}
                  title="Clear all selected clusters"
                  className="ml-2 rounded border border-slate-300 bg-white px-2 py-1 text-[11px] font-semibold text-slate-600 hover:bg-slate-100"
                >
                  ✕ clear
                </button>
              )}
              <button
                onClick={() => setLegendOpen((v) => !v)}
                className="ml-1 text-slate-500"
              >
                {legendOpen ? "▾" : "▸"}
              </button>
            </div>

            {/* Sort + filter controls (protein view only — the smaller views
                 don't have a meaningful type split) */}
            {legendOpen && kind === "protein" && (
              <div className="border-b border-slate-100 bg-slate-50 p-2 text-[11px]">
                <div className="mb-1.5 flex items-center gap-1">
                  <span className="text-slate-500">Sort:</span>
                  {(
                    [
                      ["id", "cluster id"],
                      ["size", "size"],
                      ["type", "function"],
                      ["name", "name"],
                    ] as [typeof clusterSortMode, string][]
                  ).map(([v, label]) => (
                    <button
                      key={v}
                      onClick={() => setClusterSortMode(v)}
                      className={`rounded border px-1.5 py-0.5 transition ${
                        clusterSortMode === v
                          ? "border-sky-500 bg-sky-100 text-sky-800"
                          : "border-slate-300 bg-white text-slate-600 hover:bg-slate-100"
                      }`}
                    >
                      {label}
                    </button>
                  ))}
                </div>
                <div className="flex flex-wrap items-center gap-1">
                  <span className="text-slate-500">Show:</span>
                  {(
                    ["Enzyme", "Transporter", "Receptor", "Other", "Mixed"] as
                      ProteinType[]
                  ).map((t) => {
                    const on = typeFilter.has(t);
                    return (
                      <button
                        key={t}
                        onClick={() => {
                          const next = new Set(typeFilter);
                          if (on) next.delete(t);
                          else next.add(t);
                          setTypeFilter(next);
                        }}
                        className={`flex items-center gap-1 rounded border px-1.5 py-0.5 transition ${
                          on
                            ? "border-slate-800 text-slate-900"
                            : "border-slate-200 bg-white text-slate-400 line-through"
                        }`}
                        style={{
                          background: on ? TYPE_COLOR[t] + "22" : undefined,
                        }}
                      >
                        <span
                          className="inline-block h-2 w-2 rounded-full"
                          style={{ background: TYPE_COLOR[t] }}
                        />
                        {t}
                      </button>
                    );
                  })}
                </div>
              </div>
            )}

            {legendOpen && (
              <div className="max-h-[540px] overflow-y-auto p-2">
                {clustersForSidebar.map((c) => {
                  const isPick = highlightSet.has(c.id);
                  const clr = swatchFor(c.id);
                  const typeInfo = clusterTypeMap.get(c.id);
                  const dom = typeInfo?.dominant;
                  return (
                    <button
                      key={c.id}
                      onClick={() => {
                        toggleClusterHighlight(c.id);
                        setFocusedCompound(null);
                        clearFocusedProteins();
                      }}
                      className={`mb-1 flex w-full items-center gap-2 rounded-md border px-2 py-1.5 text-left text-[12px] transition ${
                        isPick
                          ? "border-sky-500 bg-sky-50 text-slate-900"
                          : "border-slate-200 bg-white text-slate-700 hover:bg-slate-50"
                      }`}
                    >
                      <span
                        className="inline-block h-4 w-4 flex-shrink-0 rounded-sm ring-1 ring-black/10"
                        style={{ background: clr }}
                      />
                      <span className="w-6 flex-shrink-0 font-mono text-[11px] text-slate-500">
                        {c.display_id}
                      </span>
                      <span className="flex-1 truncate">
                        {c.name || "(unnamed)"}
                      </span>
                      {dom && kind === "protein" && (
                        <span
                          title={`${dom} · ${
                            Object.entries(typeInfo!.counts)
                              .filter(([, n]) => n > 0)
                              .map(([t, n]) => `${t}: ${n}`)
                              .join(", ")
                          }`}
                          className="flex-shrink-0 rounded-full px-1.5 py-0.5 text-[9px] font-semibold"
                          style={{
                            background: TYPE_COLOR[dom] + "22",
                            color: TYPE_COLOR[dom],
                          }}
                        >
                          {dom.slice(0, 3).toUpperCase()}
                        </span>
                      )}
                      <span className="flex-shrink-0 text-[10px] text-slate-500">
                        n={c.n.toLocaleString()}
                      </span>
                    </button>
                  );
                })}
                {clustersForSidebar.length === 0 && (
                  <div className="px-2 py-4 text-center text-xs text-slate-400">
                    No clusters match the current filter
                  </div>
                )}
              </div>
            )}
          </aside>
        </div>

        {/* Flat entry table — one row per pooled entry (from lasso,
            cluster, or search). Same-shape data as the marimo table.
            Includes a Select-all-N-proteins button. */}
        {pool.length > 0 && (
          <div className="mt-6">
            <div className="mb-2 flex flex-wrap items-baseline gap-3">
              <div className="text-lg font-semibold text-slate-900">
                📋 Selected entries
              </div>
              <div className="text-sm text-slate-500">
                {pool.length.toLocaleString()} row
                {pool.length !== 1 ? "s" : ""}
                {(kind === "protein" || kind === "molecule")
                  ? " · tick any row to add its detail card below"
                  : ""}
              </div>
              {kind === "protein" && (
                <div className="ml-auto flex items-center gap-2">
                  {(() => {
                    // Cap at 100 to keep the browser responsive when
                    // rendering many detail cards (each has structures).
                    const MAX = 100;
                    const targets = (pool as Protein[]).slice(0, MAX);
                    const allPicked = targets.every((p) =>
                      isProteinFocused(p.accession),
                    );
                    return (
                      <>
                        <button
                          onClick={() => {
                            if (allPicked) {
                              setFocusedProteins((cur) =>
                                cur.filter(
                                  (x) =>
                                    !targets.some(
                                      (t) => t.accession === x.accession,
                                    ),
                                ),
                              );
                            } else {
                              setFocusedProteins((cur) => {
                                const have = new Set(cur.map((x) => x.accession));
                                const additions = targets.filter(
                                  (t) => !have.has(t.accession),
                                );
                                return [...cur, ...additions];
                              });
                            }
                          }}
                          className={`rounded-md border px-3 py-1 text-xs font-semibold transition ${
                            allPicked
                              ? "border-rose-500 bg-rose-50 text-rose-700 hover:bg-rose-100"
                              : "border-sky-500 bg-sky-50 text-sky-700 hover:bg-sky-100"
                          }`}
                        >
                          {allPicked
                            ? `☒ Deselect all ${targets.length}`
                            : `☑ Select all ${targets.length}${
                                pool.length > MAX ? ` (of ${pool.length})` : ""
                              }`}
                        </button>
                        {pool.length > MAX && !allPicked && (
                          <span
                            title="Rendering more than 100 detail cards at once slows the browser; the button caps at the first 100 rows."
                            className="cursor-help self-center text-xs text-slate-500"
                          >
                            (capped at 100)
                          </span>
                        )}
                      </>
                    );
                  })()}
                </div>
              )}
              {kind === "molecule" && (
                <div className="ml-auto flex items-center gap-2">
                  {(() => {
                    const MAX = 100;
                    const targets = (pool as Molecule[]).slice(0, MAX);
                    const allPicked = targets.every((m) => isMoleculeFocused(m));
                    return (
                      <>
                        <button
                          onClick={() => {
                            if (allPicked) {
                              setFocusedMolecules((cur) =>
                                cur.filter(
                                  (x) => !targets.some((t) => molKey(t) === molKey(x)),
                                ),
                              );
                            } else {
                              setFocusedMolecules((cur) => {
                                const have = new Set(cur.map(molKey));
                                return [...cur, ...targets.filter((t) => !have.has(molKey(t)))];
                              });
                            }
                          }}
                          className={`rounded-md border px-3 py-1 text-xs font-semibold transition ${
                            allPicked
                              ? "border-rose-500 bg-rose-50 text-rose-700 hover:bg-rose-100"
                              : "border-sky-500 bg-sky-50 text-sky-700 hover:bg-sky-100"
                          }`}
                        >
                          {allPicked
                            ? `☒ Deselect all ${targets.length}`
                            : `☑ Select all ${targets.length}${
                                pool.length > MAX ? ` (of ${pool.length})` : ""
                              }`}
                        </button>
                        {pool.length > MAX && !allPicked && (
                          <span
                            title="Capped at 100 to keep the browser responsive."
                            className="cursor-help self-center text-xs text-slate-500"
                          >
                            (capped at 100)
                          </span>
                        )}
                      </>
                    );
                  })()}
                </div>
              )}
            </div>
            <div className="max-h-[360px] overflow-auto rounded-lg border border-slate-200 bg-white">
              <table className="w-full text-[12px]">
                <thead className="sticky top-0 z-10 bg-slate-100 text-left text-[11px] uppercase tracking-wide text-slate-500">
                  {kind === "protein" ? (
                    <tr>
                      <th className="px-2 py-1.5"></th>
                      <th className="px-2 py-1.5">Accession</th>
                      <th className="px-2 py-1.5">Name</th>
                      <th className="px-2 py-1.5">Gene</th>
                      <th className="px-2 py-1.5">Organism</th>
                      <th className="px-2 py-1.5 text-right">Cluster</th>
                      <th className="px-2 py-1.5">EC</th>
                      <th className="px-2 py-1.5 text-right">Refs</th>
                      <th className="px-2 py-1.5">★</th>
                    </tr>
                  ) : kind === "molecule" ? (
                    <tr>
                      <th className="px-2 py-1.5"></th>
                      <th className="px-2 py-1.5">Compound</th>
                      <th className="px-2 py-1.5">ChEBI</th>
                      <th className="px-2 py-1.5 text-right">Cluster</th>
                      <th className="px-2 py-1.5 text-right">Proteins</th>
                      <th className="px-2 py-1.5">★</th>
                    </tr>
                  ) : (
                    <tr>
                      <th className="px-2 py-1.5">Compound</th>
                      <th className="px-2 py-1.5">ChEBI</th>
                      <th className="px-2 py-1.5 text-right">Cluster</th>
                    </tr>
                  )}
                </thead>
                <tbody>
                  {pool.slice(0, 500).map((r) => {
                    if (kind === "protein") {
                      const p = r as Protein;
                      const isPick = isProteinFocused(p.accession);
                      return (
                        <tr
                          key={p.accession}
                          onClick={() => toggleProtein(p)}
                          className={`cursor-pointer border-b border-slate-100 transition ${
                            isPick ? "bg-sky-100" : "hover:bg-slate-50"
                          }`}
                        >
                          <td className="px-2 py-1.5">
                            <span
                              className={`inline-block h-3.5 w-3.5 rounded-sm border transition ${
                                isPick
                                  ? "border-sky-600 bg-sky-600"
                                  : "border-slate-400 bg-white"
                              }`}
                            >
                              {isPick && (
                                <span className="block text-center text-[9px] leading-[13px] text-white">
                                  ✓
                                </span>
                              )}
                            </span>
                          </td>
                          <td className="px-2 py-1.5 font-mono font-semibold text-sky-700">
                            {p.accession}
                          </td>
                          <td className="px-2 py-1.5 text-slate-900">
                            {(p.protein_names || "").slice(0, 60)}
                          </td>
                          <td className="px-2 py-1.5 text-slate-600">
                            {(p.gene_names || "").slice(0, 24)}
                          </td>
                          <td className="px-2 py-1.5 italic text-slate-600">
                            {(p.organism || "").slice(0, 32)}
                          </td>
                          <td className="px-2 py-1.5 text-right text-slate-600">
                            {((p.cluster ?? 0) + 1).toString()}
                          </td>
                          <td className="px-2 py-1.5 font-mono text-[11px] text-slate-500">
                            {(p.ec_numbers || "")
                              .split(/[,;]/)
                              .slice(0, 2)
                              .join(" ")
                              .slice(0, 24)}
                          </td>
                          <td className="px-2 py-1.5 text-right text-slate-500">
                            {p.pubmed_count || ""}
                          </td>
                          <td className="px-2 py-1.5">
                            {p.is_literature_recruited ? (
                              <span
                                className="text-amber-500"
                                title="Newly recruited"
                              >
                                ★
                              </span>
                            ) : (
                              ""
                            )}
                          </td>
                        </tr>
                      );
                    }
                    if (kind === "molecule") {
                      const m = r as Molecule;
                      const isPick = isMoleculeFocused(m);
                      const nProts = splitList(m.interacting_protein_accessions ?? "").length;
                      return (
                        <tr
                          key={rowKey(kind, r)}
                          onClick={() => toggleMolecule(m)}
                          className={`cursor-pointer border-b border-slate-100 transition ${
                            isPick ? "bg-sky-100" : "hover:bg-slate-50"
                          }`}
                        >
                          <td className="px-2 py-1.5">
                            <span
                              className={`inline-block h-3.5 w-3.5 rounded-sm border transition ${
                                isPick ? "border-sky-600 bg-sky-600" : "border-slate-400 bg-white"
                              }`}
                            >
                              {isPick && (
                                <span className="block text-center text-[9px] leading-[13px] text-white">✓</span>
                              )}
                            </span>
                          </td>
                          <td className="px-2 py-1.5 text-slate-900">
                            {(m.compound_name || "").slice(0, 60)}
                          </td>
                          <td className="px-2 py-1.5 font-mono text-[11px] text-slate-500">
                            {m.chebi_id ? `CHEBI:${String(m.chebi_id).replace(/\.0$/, "")}` : ""}
                          </td>
                          <td className="px-2 py-1.5 text-right text-slate-600">
                            {String((m as any).cluster ?? "")}
                          </td>
                          <td className="px-2 py-1.5 text-right text-slate-600">
                            {nProts > 0 ? nProts : ""}
                          </td>
                          <td className="px-2 py-1.5">
                            {m.is_literature_recruited ? (
                              <span className="text-amber-500" title="Newly recruited">★</span>
                            ) : ""}
                          </td>
                        </tr>
                      );
                    }
                    const m = r as Molecule | NatsynEntry;
                    return (
                      <tr
                        key={rowKey(kind, r)}
                        onClick={() => setFocusedPoint(r)}
                        className="cursor-pointer border-b border-slate-100 hover:bg-slate-50"
                      >
                        <td className="px-2 py-1.5 text-slate-900">
                          {(m.compound_name || "").slice(0, 60)}
                        </td>
                        <td className="px-2 py-1.5 font-mono text-slate-500">
                          {m.chebi_id ? `CHEBI:${m.chebi_id}` : ""}
                        </td>
                        <td className="px-2 py-1.5 text-right text-slate-600">
                          {String((m as any).cluster ?? "")}
                        </td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
              {pool.length > 500 && (
                <div className="border-t border-slate-200 bg-slate-50 px-3 py-1.5 text-xs text-slate-500">
                  …showing first 500 of {pool.length.toLocaleString()} rows
                </div>
              )}
            </div>
          </div>
        )}

        {/* Steroid catalogue */}
        {kind === "protein" && compoundTiles.length > 0 && (
          <div className="mt-6">
            <div className="mb-3 flex flex-wrap items-baseline gap-3">
              <div className="text-lg font-semibold text-slate-900">
                🧬 Steroids in this pool
              </div>
              <div className="text-sm text-slate-500">
                {compoundTiles.length.toLocaleString()} compounds across{" "}
                {pool.length.toLocaleString()} proteins
              </div>
              {focusedCompound && (
                <button
                  className="ml-auto rounded-md border border-slate-300 bg-white px-2 py-1 text-xs text-slate-700 hover:bg-slate-100"
                  onClick={() => setFocusedCompound(null)}
                >
                  ✕ clear steroid focus
                </button>
              )}
            </div>
            <div className="grid gap-3 sm:grid-cols-2 md:grid-cols-3 xl:grid-cols-4">
              {compoundTiles.slice(0, 60).map((t) => {
                const isPick = focusedCompound === t.name;
                const src = structUrl(t.name, t.chebi);
                return (
                  <button
                    key={t.name}
                    onClick={() =>
                      setFocusedCompound(isPick ? null : t.name)
                    }
                    className={`overflow-hidden rounded-xl border bg-white p-3 text-left text-xs transition ${
                      isPick
                        ? "border-sky-500 ring-2 ring-sky-400"
                        : "border-slate-200 hover:border-slate-300 hover:shadow-md"
                    }`}
                  >
                    <div className="flex items-center justify-center rounded-lg bg-white"
                         style={{ height: 210, padding: 6 }}>
                      {src ? (
                        // eslint-disable-next-line @next/next/no-img-element
                        <img
                          src={src}
                          alt={t.name}
                          className="h-full w-full object-contain"
                          style={{ imageRendering: "auto" }}
                        />
                      ) : (
                        <span className="text-[11px] text-slate-400">
                          no structure
                        </span>
                      )}
                    </div>
                    <div className="mt-2 line-clamp-2 min-h-[2.4em] text-sm font-medium text-slate-900 text-center">
                      {t.name}
                    </div>
                    {t.chebi && (
                      <div className="mt-0.5 text-center font-mono text-[10px] text-slate-500">
                        CHEBI:{t.chebi}
                      </div>
                    )}
                    <div className="mt-1 flex justify-center">
                      <span className="rounded-full bg-indigo-50 px-2 py-0.5 text-[11px] font-medium text-indigo-700">
                        {t.prots.length.toLocaleString()} protein
                        {t.prots.length !== 1 ? "s" : ""}
                      </span>
                    </div>
                    {isPick && (
                      <div className="mt-1 text-center text-[11px] font-semibold text-sky-700">
                        ◉ FOCUSED
                      </div>
                    )}
                  </button>
                );
              })}
            </div>
            {compoundTiles.length > 60 && (
              <div className="mt-2 text-xs text-slate-500">
                …and {compoundTiles.length - 60} more compounds
              </div>
            )}

            {focusedCompound && focusedCompoundProts.length > 0 && (
              <div className="mt-4 rounded-xl border-2 border-sky-400 bg-sky-50 p-4">
                <div className="mb-2 flex flex-wrap items-baseline gap-2 text-sm text-sky-800">
                  <b className="text-slate-900">{focusedCompound}</b>
                  <span className="text-slate-500">
                    ({focusedCompoundProts.length.toLocaleString()})
                  </span>
                  <span className="text-slate-600">
                    — tick any protein below to add its full detail card
                  </span>
                  {focusedProteins.length > 0 && (
                    <button
                      className="ml-auto rounded border border-slate-300 bg-white px-2 py-0.5 text-xs text-slate-700 hover:bg-slate-100"
                      onClick={clearFocusedProteins}
                    >
                      ✕ clear all {focusedProteins.length} details
                    </button>
                  )}
                </div>
                <div className="max-h-[380px] overflow-y-auto rounded-md border border-slate-200 bg-white">
                  {focusedCompoundProts.slice(0, 200).map((p) => {
                    const isPick = isProteinFocused(p.accession);
                    return (
                      <button
                        key={p.accession}
                        onClick={() => toggleProtein(p)}
                        className={`flex w-full flex-wrap items-baseline gap-2 border-b border-slate-100 px-3 py-2 text-left text-[12px] transition ${
                          isPick ? "bg-sky-100" : "hover:bg-slate-50"
                        }`}
                      >
                        <span
                          className={`inline-block h-3.5 w-3.5 flex-shrink-0 rounded-sm border transition ${
                            isPick
                              ? "border-sky-600 bg-sky-600"
                              : "border-slate-400 bg-white"
                          }`}
                        >
                          {isPick && (
                            <span className="block text-center text-[9px] leading-[13px] text-white">
                              ✓
                            </span>
                          )}
                        </span>
                        <span className="font-mono font-semibold text-sky-700">
                          {p.accession}
                        </span>
                        <span className="text-slate-900">
                          {(p.protein_names || "").slice(0, 60)}
                        </span>
                        <span className="text-slate-500">
                          {p.gene_names ? `· ${p.gene_names.slice(0, 20)}` : ""}{" "}
                          · {(p.organism || "").slice(0, 40)}
                          {p.pubmed_count ? ` · ${p.pubmed_count} refs` : ""}
                        </span>
                      </button>
                    );
                  })}
                </div>
                {focusedCompoundProts.length > 200 && (
                  <div className="mt-2 text-xs text-slate-500">
                    …showing first 200 of {focusedCompoundProts.length}
                  </div>
                )}
              </div>
            )}
          </div>
        )}

        {/* Molecule / natsyn view — steroid tile grid (each is its OWN
            structure). Shows up when a cluster is picked or the user
            drags to select a region. */}
        {kind === "molecule" && focusedMolecules.length > 0 && (
          <div className="mt-6 space-y-4">
            <div className="flex flex-wrap items-baseline gap-3">
              <div className="text-lg font-semibold text-slate-900">
                📄 Detail cards ({focusedMolecules.length})
              </div>
              <button
                className="ml-auto rounded-md border border-slate-300 bg-white px-2 py-1 text-xs text-slate-700 hover:bg-slate-100"
                onClick={clearFocusedMolecules}
              >
                ✕ close all
              </button>
            </div>
            {focusedMolecules.map((m) => (
              <MoleculeDetailCard
                key={molKey(m)}
                molecule={m}
                structUrl={structUrl}
                onClose={() => toggleMolecule(m)}
              />
            ))}
          </div>
        )}

        {(kind === "molecule" || kind === "natsyn") && pool.length > 0 && (
          <div className="mt-6">
            <div className="mb-3 flex flex-wrap items-baseline gap-3">
              <div className="text-lg font-semibold text-slate-900">
                🧬 Steroids in this selection
              </div>
              <div className="text-sm text-slate-500">
                {pool.length.toLocaleString()} compounds
              </div>
            </div>
            <div className="grid gap-3 sm:grid-cols-2 md:grid-cols-3 xl:grid-cols-4">
              {(pool as (Molecule | NatsynEntry)[])
                .slice(0, 60)
                .map((m) => {
                  const nm = m.compound_name;
                  const cb = m.chebi_id;
                  const src = structUrl(nm, cb);
                  const isPick = focusedPoint
                    ? rowKey(kind, focusedPoint) ===
                      rowKey(kind, m as any)
                    : false;
                  return (
                    <button
                      key={rowKey(kind, m as any)}
                      onClick={() => setFocusedPoint(m as any)}
                      className={`overflow-hidden rounded-xl border bg-white p-3 text-left text-xs transition ${
                        isPick
                          ? "border-sky-500 ring-2 ring-sky-400"
                          : "border-slate-200 hover:border-slate-300 hover:shadow-md"
                      }`}
                    >
                      <div
                        className="flex items-center justify-center rounded-lg bg-white"
                        style={{ height: 200, padding: 6 }}
                      >
                        {src ? (
                          // eslint-disable-next-line @next/next/no-img-element
                          <img
                            src={src}
                            alt={nm}
                            className="h-full w-full object-contain"
                          />
                        ) : (
                          <span className="text-[11px] text-slate-400">
                            no structure
                          </span>
                        )}
                      </div>
                      <div className="mt-2 line-clamp-2 min-h-[2.4em] text-center text-sm font-medium text-slate-900">
                        {nm}
                      </div>
                      {cb && (
                        <div className="mt-0.5 text-center font-mono text-[10px] text-slate-500">
                          CHEBI:{cb}
                        </div>
                      )}
                    </button>
                  );
                })}
            </div>
            {pool.length > 60 && (
              <div className="mt-2 text-xs text-slate-500">
                …and {pool.length - 60} more compounds
              </div>
            )}
          </div>
        )}

        {focusedProteins.length > 0 && (
          <div className="mt-6 space-y-4">
            <div className="flex flex-wrap items-baseline gap-3">
              <div className="text-lg font-semibold text-slate-900">
                📄 Detail cards ({focusedProteins.length})
              </div>
              <button
                className="ml-auto rounded-md border border-slate-300 bg-white px-2 py-1 text-xs text-slate-700 hover:bg-slate-100"
                onClick={clearFocusedProteins}
              >
                ✕ close all
              </button>
            </div>
            {focusedProteins.map((p) => (
              <ProteinDetailCard
                key={p.accession}
                protein={p}
                structUrl={structUrl}
                onClose={() => toggleProtein(p)}
              />
            ))}
          </div>
        )}

      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
function ProteinDetailCard({
  protein,
  structUrl,
  onClose,
}: {
  protein: Protein;
  structUrl: (n: string, c?: string) => string | null;
  onClose: () => void;
}) {
  const compounds = splitList(protein.interacting_compounds ?? "");
  const chebis = splitList(protein.interacting_chebi_ids ?? "");
  const rheas = splitList(protein.rhea_reactions ?? "");
  const gos = splitList(protein.go_labels ?? "").slice(0, 8);
  const kws = splitList(protein.keyword_labels ?? "").slice(0, 8);
  const pubmeds = splitList(protein.pubmed_ids ?? "").slice(0, 5);

  return (
    <div className="mt-6 rounded-xl border border-slate-200 bg-white p-5 shadow-md">
      <div className="mb-3 flex flex-wrap items-baseline gap-3">
        <div className="text-lg font-bold text-slate-900">
          {protein.protein_names}
        </div>
        {protein.is_literature_recruited ? (
          <span className="rounded-full bg-amber-100 px-2 py-0.5 text-[10px] font-semibold text-amber-800">
            ★ NEW
          </span>
        ) : null}
        <button
          className="ml-auto rounded-md border border-slate-300 bg-white px-2 py-1 text-xs text-slate-700 hover:bg-slate-100"
          onClick={onClose}
        >
          ✕ close
        </button>
      </div>
      <div className="mb-3 text-sm text-slate-600">
        <span className="font-mono text-sky-700">{protein.accession}</span> ·{" "}
        {protein.gene_names || "—"} · <em>{protein.organism}</em> ·{" "}
        {protein.length_aa} aa · Cluster {(protein.cluster ?? 0) + 1}
      </div>

      {protein.putative_role && (
        <div className="mb-3 flex items-start gap-2 rounded-md border border-amber-400 bg-amber-50 px-3 py-2 text-sm text-amber-900">
          <span className="text-base leading-none">⚠</span>
          <div>
            <div className="font-semibold">Putative annotation</div>
            <div className="text-xs">
              {protein.putative_role}
              {protein.cluster_level_evidence === "yes" && (
                <span className="ml-1 italic">— individual catalytic activity not biochemically confirmed; see paper for details.</span>
              )}
            </div>
          </div>
        </div>
      )}

      <div className="mb-4 flex flex-wrap gap-2 text-sm">
        {protein.identifier_type === "GenBank_Protein" ? (
          <a
            href={`https://www.ncbi.nlm.nih.gov/protein/${protein.accession}.1`}
            target="_blank"
            className="rounded border border-sky-500 px-2 py-1 text-sky-700 hover:bg-sky-50"
          >
            NCBI Protein
          </a>
        ) : protein.identifier_type === "Locus_Tag" ? (
          <span className="rounded border border-slate-300 px-2 py-1 text-xs text-slate-500">
            no external DB (locus tag — see paper)
          </span>
        ) : (
          <>
            <a
              href={`https://www.uniprot.org/uniprotkb/${protein.accession}/entry`}
              target="_blank"
              className="rounded border border-sky-500 px-2 py-1 text-sky-700 hover:bg-sky-50"
            >
              UniProt
            </a>
            <a
              href={`https://alphafold.ebi.ac.uk/entry/${protein.accession}`}
              target="_blank"
              className="rounded border border-emerald-500 px-2 py-1 text-emerald-700 hover:bg-emerald-50"
            >
              AlphaFold
            </a>
          </>
        )}
        {protein.paper_url && protein.paper_url.startsWith("http") && (
          <a
            href={protein.paper_url}
            target="_blank"
            className="rounded border border-amber-500 px-2 py-1 text-amber-700 hover:bg-amber-50"
          >
            Source paper
          </a>
        )}
      </div>

      {protein.annotation && (
        <div className="mb-3 inline-flex items-center gap-2 rounded border border-slate-200 bg-slate-50 px-2 py-1 text-xs text-slate-700">
          <span className="font-semibold uppercase tracking-wide text-slate-500">UniProt annotation score</span>
          <span className="font-mono text-slate-900">{protein.annotation}/5</span>
        </div>
      )}
      {protein.sequence_source && (
        <div className="mb-3 text-xs text-slate-500">
          <span className="font-semibold">Sequence source:</span> {protein.sequence_source}
        </div>
      )}

      <div className="grid gap-4 md:grid-cols-2">
        <div>
          {protein.ec_numbers && (
            <div className="mb-2 text-sm text-slate-700">
              <span className="font-semibold text-slate-900">EC:</span>{" "}
              {splitList(protein.ec_numbers.replace(/,/g, ";"))
                .slice(0, 8)
                .map((e) => (
                  <span
                    key={e}
                    className="mr-1 inline-block rounded bg-slate-100 px-1.5 py-0.5 font-mono text-[11px] text-sky-700"
                  >
                    {e}
                  </span>
                ))}
            </div>
          )}
          {rheas.length > 0 && (
            <div className="mb-2 text-sm text-slate-700">
              <span className="font-semibold text-slate-900">Rhea:</span>{" "}
              {rheas.slice(0, 6).map((r) => (
                <a
                  key={r}
                  href={`https://www.rhea-db.org/rhea/${r}`}
                  target="_blank"
                  className="mr-1 inline-block rounded bg-slate-100 px-1.5 py-0.5 font-mono text-[11px] text-sky-700 hover:bg-slate-200"
                >
                  Rhea:{r}
                </a>
              ))}
            </div>
          )}
          {gos.length > 0 && (
            <div className="mb-2 text-sm text-slate-700">
              <span className="font-semibold text-slate-900">GO:</span>{" "}
              {gos.map((g) => (
                <span
                  key={g}
                  className="mr-1 inline-block rounded bg-indigo-50 px-1.5 py-0.5 text-[11px] text-indigo-800"
                >
                  {g.slice(0, 40)}
                </span>
              ))}
            </div>
          )}
          {kws.length > 0 && (
            <div className="mb-2 text-sm text-slate-700">
              <span className="font-semibold text-slate-900">Keywords:</span>{" "}
              {kws.map((k) => (
                <span
                  key={k}
                  className="mr-1 inline-block rounded bg-amber-50 px-1.5 py-0.5 text-[11px] text-amber-800"
                >
                  {k.slice(0, 30)}
                </span>
              ))}
            </div>
          )}
          {pubmeds.length > 0 && (
            <div className="mb-2 text-sm text-slate-700">
              <span className="font-semibold text-slate-900">
                PubMed ({protein.pubmed_count}):
              </span>{" "}
              {pubmeds.map((pm) => (
                <a
                  key={pm}
                  href={`https://pubmed.ncbi.nlm.nih.gov/${pm}`}
                  target="_blank"
                  className="mr-1 inline-block rounded bg-emerald-50 px-1.5 py-0.5 text-[11px] text-emerald-800 hover:bg-emerald-100"
                >
                  PMID {pm}
                </a>
              ))}
            </div>
          )}
        </div>
        <div>
          <div className="mb-2 text-sm font-semibold text-slate-900">
            Interacting steroids ({compounds.length})
          </div>
          {compounds.length === 0 && (
            <div className="mb-3 rounded-md border border-dashed border-amber-400 bg-amber-50 p-3 text-xs text-amber-900">
              <div className="font-semibold">⚠ No steroids associated with this protein in the atlas.</div>
              <div className="mt-1">
                This is unusual for a protein in an atlas whose premise is protein ↔ steroid linkage. Possible reasons: (1) the protein has no documented steroid substrate in Rhea / UniProt, (2) a cited paper has not yet been audited to link it to specific steroid rows, or (3) this is a data gap worth reporting. See <span className="font-semibold">Source paper</span> (above) for the primary activity claim.
              </div>
            </div>
          )}
          <div className="grid grid-cols-3 gap-2">
            {compounds.slice(0, 12).map((nm, i) => {
              const cb = chebis[i] ?? "";
              const src = structUrl(nm, cb);
              return (
                <div
                  key={nm + i}
                  className="overflow-hidden rounded-md border border-slate-200 bg-white p-2 text-center text-[10px]"
                >
                  <div
                    className="flex items-center justify-center"
                    style={{ height: 110 }}
                  >
                    {src ? (
                      // eslint-disable-next-line @next/next/no-img-element
                      <img
                        src={src}
                        alt={nm}
                        className="max-h-full max-w-full object-contain"
                      />
                    ) : (
                      <span className="text-slate-400">no img</span>
                    )}
                  </div>
                  <div className="mt-1 line-clamp-2 min-h-[2.2em] text-slate-700">
                    {nm.slice(0, 40)}
                  </div>
                </div>
              );
            })}
            {compounds.length > 12 && (
              <div className="col-span-3 text-center text-[11px] text-slate-500">
                …and {compounds.length - 12} more
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  );
}

function MoleculeDetailCard({
  molecule,
  structUrl,
  onClose,
}: {
  molecule: Molecule;
  structUrl: (n: string, c?: string) => string | null;
  onClose: () => void;
}) {
  const accs = splitList(molecule.interacting_protein_accessions ?? "");
  const src = structUrl(molecule.compound_name, molecule.chebi_id);
  const synonyms = splitList(molecule.synonyms ?? "")
    .map((s) => s.replace(/^aka\s+/i, "").trim())
    .filter(Boolean)
    .slice(0, 8);
  const tier = molecule.evidence_tier;
  const tierLabel =
    tier === "1" ? "⚠ Chemistry reference only — no biological detection"
    : tier === "2" ? "⚠ MS/MS match in public data (isomer uncertain)"
    : tier === "3" ? "✓ Detected in human IBD stool (iHMP2, isomer-resolved)"
    : tier === "4-5" ? "✓ Cross-cohort validated in human IBD (iHMP2 + Boland)"
    : tier === "5" ? "✓ Produced by gut bacteria in culture"
    : "";

  return (
    <div className="mt-6 rounded-xl border border-slate-200 bg-white p-5 shadow-md">
      <div className="mb-3 flex flex-wrap items-baseline gap-3">
        <div className="text-lg font-bold text-slate-900">
          {molecule.compound_name}
        </div>
        {molecule.is_literature_recruited ? (
          <span className="rounded-full bg-amber-100 px-2 py-0.5 text-[10px] font-semibold text-amber-800">
            ★ NEW
          </span>
        ) : null}
        <button
          className="ml-auto rounded-md border border-slate-300 bg-white px-2 py-1 text-xs text-slate-700 hover:bg-slate-100"
          onClick={onClose}
        >
          ✕ close
        </button>
      </div>
      <div className="mb-3 text-sm text-slate-600">
        {molecule.chebi_id ? (
          <a
            href={`https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:${String(molecule.chebi_id).replace(/\.0$/, "")}`}
            target="_blank"
            rel="noopener noreferrer"
            className="font-mono text-sky-700 hover:underline"
          >
            CHEBI:{String(molecule.chebi_id).replace(/\.0$/, "")}
          </a>
        ) : (
          <span className="italic text-slate-500">no canonical ChEBI</span>
        )}
        {molecule.cluster !== undefined && molecule.cluster !== null && (
          <span className="ml-3">Cluster {((molecule.cluster as number) ?? 0) + 1}</span>
        )}
      </div>

      {tierLabel && (
        <div className="mb-3 inline-flex items-center gap-2 rounded-md border border-emerald-300 bg-emerald-50 px-2 py-1 text-xs text-emerald-900">
          <span className="font-semibold">Evidence tier {tier}:</span>
          <span>{tierLabel}</span>
        </div>
      )}

      {molecule.evidence_note && (
        <div className="mb-3 rounded border border-slate-200 bg-slate-50 p-3 text-sm text-slate-800">
          <div className="mb-1 text-xs font-semibold uppercase tracking-wide text-slate-500">Note</div>
          <div className="whitespace-pre-wrap">{molecule.evidence_note}</div>
        </div>
      )}

      {molecule.paper_url && (
        <div className="mb-3">
          <a
            href={molecule.paper_url}
            target="_blank"
            rel="noopener noreferrer"
            className="inline-block rounded border border-sky-500 px-2 py-1 text-xs text-sky-700 hover:bg-sky-50"
          >
            📄 Source paper
          </a>
        </div>
      )}

      <div className="grid gap-4 md:grid-cols-2">
        <div>
          {src ? (
            <div className="mb-3 flex items-center justify-center rounded-lg border border-slate-200 bg-white p-3" style={{ minHeight: 260 }}>
              {/* eslint-disable-next-line @next/next/no-img-element */}
              <img
                src={src}
                alt={molecule.compound_name}
                className="max-h-[260px] max-w-full object-contain"
              />
            </div>
          ) : (
            <div className="mb-3 rounded-lg border border-dashed border-slate-300 bg-slate-50 p-5 text-center text-xs italic text-slate-500">
              no 2D structure image available
            </div>
          )}
          {molecule.smiles && (
            <div className="mb-2 text-xs text-slate-700">
              <div className="mb-1 font-semibold uppercase tracking-wide text-slate-500">SMILES</div>
              <div className="break-all rounded bg-slate-50 px-2 py-1 font-mono text-[11px] text-slate-800">
                {molecule.smiles}
              </div>
            </div>
          )}
          {synonyms.length > 0 && (
            <div className="mb-2 text-sm text-slate-700">
              <span className="font-semibold text-slate-900">Synonyms:</span>{" "}
              {synonyms.map((s) => (
                <span
                  key={s}
                  className="mr-1 inline-block rounded bg-slate-100 px-1.5 py-0.5 text-[11px] text-slate-700"
                >
                  {s.slice(0, 40)}
                </span>
              ))}
            </div>
          )}
        </div>
        <div>
          <div className="mb-2 text-sm font-semibold text-slate-900">
            Interacting proteins ({accs.length})
          </div>
          {accs.length === 0 ? (
            <div className="rounded-md border border-dashed border-amber-400 bg-amber-50 p-3 text-xs text-amber-900">
              <div className="font-semibold">⚠ No proteins associated in the atlas.</div>
              <div className="mt-1">
                This compound is in the atlas as a literature-recruited structure but no enzyme in the atlas has a direct substrate attribution to it.
              </div>
            </div>
          ) : (
            <div className="max-h-[360px] overflow-auto rounded-lg border border-slate-200 bg-white">
              <table className="w-full text-[12px]">
                <thead className="sticky top-0 z-10 bg-slate-100 text-left text-[11px] uppercase tracking-wide text-slate-500">
                  <tr>
                    <th className="px-2 py-1.5">Accession</th>
                    <th className="px-2 py-1.5 text-right">Open</th>
                  </tr>
                </thead>
                <tbody>
                  {accs.slice(0, 200).map((acc) => (
                    <tr key={acc} className="border-b border-slate-100 hover:bg-slate-50">
                      <td className="px-2 py-1.5 font-mono text-[12px] font-semibold text-sky-800">
                        {acc}
                      </td>
                      <td className="px-2 py-1.5 text-right">
                        <a
                          href={`https://www.uniprot.org/uniprotkb/${acc}/entry`}
                          target="_blank"
                          rel="noopener noreferrer"
                          className="inline-block rounded border border-sky-400 bg-white px-1.5 py-0.5 text-[10px] text-sky-700 hover:bg-sky-50"
                          title={`Open ${acc} on UniProt`}
                        >
                          UniProt ↗
                        </a>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
              {accs.length > 200 && (
                <div className="border-t border-slate-200 bg-slate-50 px-3 py-1.5 text-xs text-slate-500">
                  …showing first 200 of {accs.length} proteins
                </div>
              )}
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
