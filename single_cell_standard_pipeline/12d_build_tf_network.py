#!/usr/bin/env python
# =============================================================================
# 12d_build_tf_network.py  --  TF -> FOCUS -> targets GRN around a focus gene
# =============================================================================
# Builds a 3-layer, directed regulatory subnetwork centered on one gene (default
# Nr4a1) from the pySCENIC outputs of 12b, per label:
#
#     upstream TFs  -->  FOCUS (Nr4a1)  -->  downstream target genes
#
# Layer 1 (upstream): TFs that predict FOCUS   = adjacency rows with target==FOCUS
# Layer 2          : FOCUS itself
# Layer 3 (downstream): genes FOCUS predicts   = adjacency rows with TF==FOCUS
# Edges drawn = every adjacency edge whose SOURCE is in {upstream TFs, FOCUS} and
# whose TARGET is in the node set, so you also see the "cross" wiring (an upstream
# TF that also hits a downstream target) and TF->TF links -- the real induced GRN,
# not just a star. Each edge is flagged `confirmed` if it is a motif-pruned regulon
# edge (from regulons.p), so the overlay separates discovery vs confident edges.
#
# IMPORTANT: GRNBoost2 `importance` is gradient-boosted feature importance (a
# predictive weight), NOT mutual information and NOT causal. Read edges as
# "coordinated/predictive regulation", same caveat as the AUC.
#
# Reads regulons.p (a pickle of pyscenic Regulon objects) -> RUN IN THE pyscenic
# conda env (the class defs must be importable to unpickle). If regulons.p is
# missing/unreadable, it degrades to adjacency-only (all edges `confirmed=False`).
#
# OUTPUT  (per label, under SCENIC/<label>/network/):
#   <FOCUS>_network_nodes.csv   gene, role, layer, is_tf, imp_to_focus, imp_from_focus
#   <FOCUS>_network_edges.csv   source, target, importance, confirmed, edge_type
#   <FOCUS>_network.png         networkx multipartite layout (solid=regulon, dashed=adj)
#
# Portable: same NR4A1_* env vars as config.R / 12b.
#   NR4A1_SCENIC_DIR   SCENIC I/O dir     (default: <NR4A1_OUTPUT>/SCENIC)
#   NR4A1_FOCUS_GENE   focus gene         (default: Nr4a1)
#   NR4A1_NET_TOP_UP   # upstream TFs     (default: 15)
#   NR4A1_NET_TOP_DOWN # downstream genes (default: 40)
#
# USAGE:  python 12d_build_tf_network.py                 # all labels
#         python 12d_build_tf_network.py --labels A B    # a subset
#         NR4A1_FOCUS_GENE=Foxp3 python 12d_build_tf_network.py
# =============================================================================
import os
import sys
import glob
import pickle
import argparse

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx

# --- CONFIG (portable) -------------------------------------------------------
ROOT_PATH  = os.environ.get(
    "NR4A1_ROOT",
    "/run/media/external/scdc_projects/2026_nr4a1_ack/r_process")
OUTPUT_DIR = os.environ.get("NR4A1_OUTPUT", os.path.join(ROOT_PATH, "seurat_output"))
SCENIC_DIR = os.environ.get("NR4A1_SCENIC_DIR", os.path.join(OUTPUT_DIR, "SCENIC"))
FOCUS_GENE = os.environ.get("NR4A1_FOCUS_GENE", "Nr4a1-cust")
TOP_UP     = int(os.environ.get("NR4A1_NET_TOP_UP", "15"))
TOP_DOWN   = int(os.environ.get("NR4A1_NET_TOP_DOWN", "40"))
# Genes removed from EVERY network (nodes + edges) before building/plotting. The
# focus is the custom Nr4a1-cust feature; the plain "Nr4a1" is the same gene with
# a different exon count, so it is dropped to avoid a duplicate hub.
DROP_GENES = set(g.strip() for g in
                 os.environ.get("NR4A1_NET_DROP_GENES", "Nr4a1").split(",") if g.strip())
# cap on target<->target ring edges (strongest by importance) to avoid a hairball
RING_MAX   = int(os.environ.get("NR4A1_NET_RING_MAX", "60"))
# Edge trimming by IMPORTANCE PERCENTILE to suppress GRNBoost2 false positives:
# start from ALL candidate edges, then keep only the top X% strongest within each
# edge class. One plot + CSV set is written per percentile (default: 5% and 10%).
EDGE_PCTLS = [float(x) for x in
              os.environ.get("NR4A1_NET_EDGE_PCTLS", "5,10").split(",") if x.strip()]


def find_adj(d):
    """Path to the GRN adjacency file: Python-API 'adjacencies.tsv' OR CLI 'adj.tsv'."""
    for fn in ("adjacencies.tsv", "adj.tsv"):
        p = os.path.join(d, fn)
        if os.path.exists(p) and os.path.getsize(p) > 0:
            return p
    return None


def _regulon_edges_from_csv(p):
    """Parse the pyscenic CLI ctx output (regulons.csv) -> set of (TF, target).
    Format: multi-header CSV indexed by (TF, MotifID); the 'TargetGenes' column
    holds a string list of (gene, weight) tuples."""
    import ast
    try:
        df = pd.read_csv(p, header=[0, 1], index_col=[0, 1], skipinitialspace=True)
    except Exception as e:
        print(f"    [note] could not read regulons.csv ({e}); adjacency-only.", flush=True)
        return None
    tg_cols = [c for c in df.columns if str(c[1]) == "TargetGenes"]
    if not tg_cols:
        return None
    tg = tg_cols[0]
    edges = set()
    for idx, row in df.iterrows():
        tf = idx[0] if isinstance(idx, tuple) else idx
        try:
            for gene, _w in ast.literal_eval(row[tg]):
                edges.add((str(tf), str(gene)))
        except Exception:
            continue
    return edges if edges else None


def load_regulon_edges(label_dir):
    """Set of (TF, target) motif-pruned edges. Tries regulons.p (Python 12b pickle)
    then regulons.csv (CLI ctx output). Returns None if neither is usable."""
    p = os.path.join(label_dir, "regulons.p")
    if os.path.exists(p):
        try:
            with open(p, "rb") as fh:
                regs = pickle.load(fh)
            edges = set()
            for r in regs:
                tf = getattr(r, "name", "").split("(")[0]
                genes = getattr(r, "genes", None)
                if genes is None:
                    g2w = getattr(r, "gene2weight", {})
                    genes = list(g2w.keys()) if g2w else []
                for g in genes:
                    edges.add((tf, g))
            if edges:
                return edges
        except Exception as e:
            print(f"    [note] could not unpickle regulons.p ({e}); trying regulons.csv.", flush=True)
    csvp = os.path.join(label_dir, "regulons.csv")
    if os.path.exists(csvp):
        return _regulon_edges_from_csv(csvp)
    return None


def build_for_label(label):
    d = os.path.join(SCENIC_DIR, label)
    adj_path = find_adj(d)
    if adj_path is None:
        print(f"  [SKIP] {label}: no adjacencies.tsv / adj.tsv", flush=True)
        return
    adj = pd.read_csv(adj_path, sep="\t")
    # GRNBoost2 columns are TF, target, importance.
    cols = {c.lower(): c for c in adj.columns}
    tf_c  = cols.get("tf", adj.columns[0])
    tg_c  = cols.get("target", adj.columns[1])
    imp_c = cols.get("importance", adj.columns[2])
    adj = adj.rename(columns={tf_c: "TF", tg_c: "target", imp_c: "importance"})

    if FOCUS_GENE not in set(adj["TF"]).union(adj["target"]):
        print(f"  [SKIP] {label}: focus gene '{FOCUS_GENE}' not in adjacencies.", flush=True)
        return

    # drop unwanted genes (e.g. the plain "Nr4a1") from the adjacency table up front
    drop = DROP_GENES - {FOCUS_GENE}
    if drop:
        adj = adj[~adj["TF"].isin(drop) & ~adj["target"].isin(drop)].copy()

    tf_universe = set(adj["TF"].unique())
    reg_edges   = load_regulon_edges(d)
    out = os.path.join(d, "network")
    os.makedirs(out, exist_ok=True)

    # ALL candidate edges (untrimmed). Percentile trimming happens per class below.
    up_all = adj[adj["target"] == FOCUS_GENE].copy()    # upstream TF -> focus
    dn_all = adj[adj["TF"] == FOCUS_GENE].copy()         # focus -> target
    if up_all.empty and dn_all.empty:
        print(f"  [SKIP] {label}: focus '{FOCUS_GENE}' has no edges.", flush=True)
        return

    def _top_pct(df, pct):
        if df.empty: return df
        k = max(1, int(np.ceil(len(df) * pct / 100.0)))
        return df.sort_values("importance", ascending=False).head(k)

    def _assemble(pct):
        """Keep the top `pct`% strongest edges within each class, build node/edge
        tables. Returns (nodes_df, sub) or None."""
        up   = _top_pct(up_all, pct)
        down = _top_pct(dn_all, pct)
        U = [g for g in dict.fromkeys(up["TF"].tolist()) if g != FOCUS_GENE]
        # Roles are MUTUALLY EXCLUSIVE: a gene that both regulates and is regulated
        # by the focus is shown as an upstream regulator (blue) only, and dropped
        # from the target set. This guarantees blue nodes touch ONLY the pink hub,
        # the hub touches only orange targets, and orange<->orange stays among
        # pure targets -- no blue<->orange edges.
        Uset = set(U)
        D = [g for g in dict.fromkeys(down["target"].tolist())
             if g != FOCUS_GENE and g not in Uset]
        Dset = set(D)
        # edges, restricted so each class only connects its allowed roles:
        up_e = up[up["TF"].isin(Uset)].copy();            up_e["edge_type"] = "upstream"    # blue -> pink
        dn_e = down[down["target"].isin(Dset)].copy();    dn_e["edge_type"] = "downstream"  # pink -> orange
        ti = adj[adj["TF"].isin(Dset) & adj["target"].isin(Dset) &
                 (adj["TF"] != adj["target"])].copy()                                        # orange <-> orange
        ti = _top_pct(ti, pct)
        ti_e = ti.copy();   ti_e["edge_type"] = "target"
        sub = pd.concat([up_e, dn_e, ti_e], ignore_index=True)
        if sub.empty: return None
        sub["confirmed"] = False if reg_edges is None else \
            [ (r.TF, r.target) in reg_edges for r in sub.itertuples() ]
        nodes = set(U) | {FOCUS_GENE} | set(D)
        imp_to   = dict(zip(up["TF"], up["importance"]))
        imp_from = dict(zip(down["target"], down["importance"]))
        rows = []
        for g in nodes:
            if g == FOCUS_GENE:  role, layer = "focus", 1
            elif g in set(U):    role, layer = "upstream", 0
            else:                role, layer = "target", 2
            rows.append(dict(gene=g, role=role, layer=layer, is_tf=g in tf_universe,
                             imp_to_focus=imp_to.get(g, np.nan),
                             imp_from_focus=imp_from.get(g, np.nan)))
        nodes_df = pd.DataFrame(rows).sort_values(["layer", "role", "gene"])
        return nodes_df, sub

    def _write(nodes_df, sub, suffix, title_sfx):
        nodes_df.to_csv(os.path.join(out, f"{FOCUS_GENE}_network_nodes{suffix}.csv"), index=False)
        sub[["TF", "target", "importance", "confirmed", "edge_type"]].rename(
            columns={"TF": "source"}).to_csv(
            os.path.join(out, f"{FOCUS_GENE}_network_edges{suffix}.csv"), index=False)
        nU = int((nodes_df.role == "upstream").sum()); nD = int((nodes_df.role == "target").sum())
        print(f"  {label}{suffix}: {len(nodes_df)} nodes ({nU} up, {nD} down) | "
              f"{len(sub)} edges ({int(sub['confirmed'].sum())} regulon-confirmed)", flush=True)
        _plot(nodes_df, sub, label + title_sfx, os.path.join(out, f"{FOCUS_GENE}_network{suffix}.png"))

    # one plot + CSV set per percentile (top 5%, top 10%, ...)
    results = {}
    for pct in EDGE_PCTLS:
        r = _assemble(pct)
        if r is None:
            print(f"  [skip] {label}: no edges at top {pct}%", flush=True); continue
        pstr = (f"{pct:g}").replace(".", "p")
        _write(r[0], r[1], f"_top{pstr}pct", f"  (top {pct:g}% edges)")
        results[pct] = r
    # unsuffixed default (= most inclusive percentile) for the R per-condition plot
    if results:
        nodes_df, sub = results[max(results)]
        nodes_df.to_csv(os.path.join(out, f"{FOCUS_GENE}_network_nodes.csv"), index=False)
        sub[["TF", "target", "importance", "confirmed", "edge_type"]].rename(
            columns={"TF": "source"}).to_csv(
            os.path.join(out, f"{FOCUS_GENE}_network_edges.csv"), index=False)


def _layered_pos(G, y_gap=0.9):
    """Multipartite-style layout with CONTROLLED vertical spacing (y_gap between
    adjacent nodes in a layer) so labels don't collide. Returns (pos, max_nodes)."""
    layers = {}
    for n, d in G.nodes(data=True):
        layers.setdefault(int(d.get("layer", 0)), []).append(n)
    pos = {}
    for li, lyr in enumerate(sorted(layers)):
        nodes = sorted(layers[lyr])
        k = len(nodes)
        for i, n in enumerate(nodes):
            pos[n] = (float(li), (i - (k - 1) / 2.0) * y_gap)
    return pos, max((len(v) for v in layers.values()), default=1)


def _draw_labels(G, pos, ax, font_size=7):
    """Label EVERY node, offset slightly right so text sits beside the node."""
    lab_pos = {n: (x + 0.07, y) for n, (x, y) in pos.items()}
    nx.draw_networkx_labels(G, lab_pos, labels={n: n for n in G.nodes},
                            font_size=font_size, horizontalalignment="left", ax=ax)


def _plot(nodes_df, edges_df, label, png):
    """Hub-and-circle layout: upstream TFs (left column) -> FOCUS hub -> targets
    arranged in a RING on the right that interacts with itself (target<->target
    GRN from the adjacency importances). Only edges touching FOCUS, plus the ring
    interactions, are drawn (no upstream->target cross clutter)."""
    G = nx.DiGraph()
    for r in nodes_df.itertuples():
        G.add_node(r.gene, role=r.role)
    for r in edges_df.itertuples():
        G.add_edge(r.TF, r.target, importance=float(r.importance),
                   confirmed=bool(r.confirmed), etype=getattr(r, "edge_type", "star"))

    focus = FOCUS_GENE
    up  = sorted(n for n in G if G.nodes[n].get("role") == "upstream")
    tgt = [g for g in nodes_df.gene[nodes_df.role == "target"]]

    # --- layout -----------------------------------------------------------------
    # upstream column pushed well LEFT with wide vertical spacing so left-anchored
    # labels never collide with the hub spokes.
    pos = {}
    up_x, up_gap = -5.2, 1.4
    for i, n in enumerate(up):
        pos[n] = (up_x, (i - (len(up) - 1) / 2.0) * up_gap)
    pos[focus] = (0.0, 0.0)
    nt = max(1, len(tgt))
    r  = max(1.8, nt * 0.17)
    cx = r + 1.8                       # ring kept ENTIRELY right of the hub (x=0):
    for i, n in enumerate(tgt):        # leftmost ring point = cx - r = 1.8 > 0, so
        ang = 2 * np.pi * i / nt + np.pi / 2.0   # it never wraps into the blue column
        pos[n] = (cx + r * np.cos(ang), r * np.sin(ang))

    role_color = {"upstream": "#2171b5", "focus": "#f768a1", "target": "#fd8d3c"}
    ncolors = [role_color.get(G.nodes[n].get("role"), "#fd8d3c") for n in G.nodes]
    nsize   = [1200 if G.nodes[n].get("role") == "focus" else
               440 if G.nodes[n].get("role") == "upstream" else 380 for n in G.nodes]

    imps = np.array([d["importance"] for *_ , d in G.edges(data=True)]) if G.number_of_edges() else np.array([1.0])
    wmax = imps.max() if imps.size else 1.0
    # wider edges (visible), scaled by GRNBoost2 importance; floor keeps weak links drawn
    def W(elist): return [1.1 + 5.5 * (G[u][v]["importance"] / wmax) for u, v in elist]

    # Edge groups. ORANGE = the GRN: it only ever touches orange targets and/or the
    # focus hub -> (a) focus -> target (hub spokes) and (b) target <-> target (ring).
    # GREY = the regulatory INPUT: upstream TFs -> focus only (never orange).
    up_star = [(u, v) for u, v, d in G.edges(data=True) if d.get("etype") == "upstream"]
    up_s    = [(u, v) for u, v in up_star if G[u][v]["confirmed"]]
    up_d    = [(u, v) for u, v in up_star if not G[u][v]["confirmed"]]
    dn_star = [(u, v) for u, v, d in G.edges(data=True) if d.get("etype") == "downstream"]
    ring    = [(u, v) for u, v, d in G.edges(data=True) if d.get("etype") == "target"]

    # size to the actual extent so the ring is circular (equal aspect) and labels fit
    xspan = (cx + r) - up_x + 6.0
    yspan = max(2 * r + 2.0, len(up) * up_gap + 2.0)
    fig_w = min(42, max(16, xspan * 0.62))
    fig_h = min(42, max(11, yspan * 0.62))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_aspect("equal")
    # (a) focus -> target hub spokes (orange, directed)
    nx.draw_networkx_edges(G, pos, edgelist=dn_star, edge_color="#d94801", alpha=0.75,
                           width=W(dn_star), arrows=True, arrowsize=9, ax=ax)
    # (b) target <-> target GRN inside the ring (orange, curved arcs)
    nx.draw_networkx_edges(G, pos, edgelist=ring, edge_color="#fd8d3c", alpha=0.55,
                           width=W(ring), arrows=True, arrowsize=6,
                           connectionstyle="arc3,rad=0.18", ax=ax)
    # upstream input: dashed = adjacency only, solid = motif-confirmed regulon (grey)
    nx.draw_networkx_edges(G, pos, edgelist=up_d, style="dashed", edge_color="#9e9e9e",
                           alpha=0.5, width=W(up_d), arrows=True, arrowsize=9, ax=ax)
    nx.draw_networkx_edges(G, pos, edgelist=up_s, style="solid", edge_color="#404040",
                           alpha=0.8, width=W(up_s), arrows=True, arrowsize=11, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_color=ncolors, node_size=nsize,
                           edgecolors="white", linewidths=0.6, alpha=0.85, ax=ax)

    # labels: upstream to the LEFT, focus to the right (bold), targets RADIALLY out
    for n, (x, y) in pos.items():
        role = G.nodes[n].get("role")
        if role == "upstream":
            ax.text(x - 0.22, y, n, fontsize=13, ha="right", va="center")
        elif role == "focus":
            ax.text(x + 0.22, y, n, fontsize=16, ha="left", va="center", fontweight="bold")
        else:
            vx, vy = x - cx, y
            nrm = (vx * vx + vy * vy) ** 0.5 or 1.0
            ax.text(x + 0.30 * vx / nrm, y + 0.30 * vy / nrm, n, fontsize=12,
                    ha=("left" if vx >= 0 else "right"), va="center")

    ax.set_title(f"{FOCUS_GENE} regulatory neighborhood — {label}\n"
                 f"grey = upstream TF input → {FOCUS_GENE} (solid = motif-confirmed regulon, dashed = adjacency)   ·   "
                 f"orange = {FOCUS_GENE} target GRN (hub spokes + target↔target ring)",
                 fontsize=13)
    ax.margins(x=0.24, y=0.08)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(png, dpi=200)
    plt.close(fig)


def plot_from_csv(nodes_csv, edges_csv, out_png):
    """Render an augmented network (backbone + co-expression) from CSVs, e.g. the
    *_augmented.csv files written by 12d_plot_tf_network.R, with networkx."""
    nodes = pd.read_csv(nodes_csv)
    edges = pd.read_csv(edges_csv)
    G = nx.DiGraph()
    mods = sorted({int(m) for m in nodes.get("module", pd.Series(dtype=float)).dropna().tolist()})
    cmap = plt.get_cmap("tab10")
    mod_color = {m: cmap(i % 10) for i, m in enumerate(mods)}
    for r in nodes.itertuples():
        G.add_node(r.gene, layer=int(r.layer), role=r.role,
                   module=(int(r.module) if ("module" in nodes.columns and pd.notna(r.module)) else None))
    for r in edges.itertuples():
        G.add_edge(r.source, r.target, importance=float(r.importance),
                   confirmed=(bool(r.confirmed) if pd.notna(r.confirmed) else False),
                   directed=(bool(r.directed) if "directed" in edges.columns else True),
                   edge_type=getattr(r, "edge_type", "backbone"))
    pos, maxk = _layered_pos(G, y_gap=0.9)

    def ncolor(n):
        d = G.nodes[n]
        if d["role"] == "focus":    return "#d73027"
        if d["role"] == "upstream": return "#4575b4"
        return mod_color.get(d["module"], "#bdbdbd")
    ncolors = [ncolor(n) for n in G.nodes]
    nsize   = [900 if G.nodes[n]["role"] == "focus" else 320 for n in G.nodes]

    be = [(u, v, d) for u, v, d in G.edges(data=True) if d["directed"]]
    ce = [(u, v, d) for u, v, d in G.edges(data=True) if not d["directed"]]
    wmax = max([d["importance"] for *_ , d in G.edges(data=True)] + [1.0])
    def W(elist): return [0.4 + 2.5 * (d["importance"] / wmax) for *_ , d in elist]

    fig, ax = plt.subplots(figsize=(15, max(9, maxk * 0.42)))
    # co-expression mesh (undirected, no arrows)
    nx.draw_networkx_edges(G, pos, edgelist=[(u, v) for u, v, _ in ce],
                           width=W(ce), edge_color="#1b9e77", alpha=0.30,
                           arrows=False, ax=ax)
    # backbone: dashed = adjacency-only, solid = regulon-confirmed
    bd = [(u, v) for u, v, d in be if not d["confirmed"]]
    bs = [(u, v) for u, v, d in be if d["confirmed"]]
    nx.draw_networkx_edges(G, pos, edgelist=bd, style="dashed", edge_color="#bdbdbd",
                           alpha=0.35, arrows=True, arrowsize=8, ax=ax)
    nx.draw_networkx_edges(G, pos, edgelist=bs, style="solid", edge_color="#444444",
                           alpha=0.75, arrows=True, arrowsize=10, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_color=ncolors, node_size=nsize,
                           edgecolors="white", linewidths=0.5, alpha=0.85, ax=ax)
    _draw_labels(G, pos, ax, font_size=7)
    ax.margins(x=0.18, y=0.03)
    ax.set_title("Augmented regulatory network (SCENIC backbone + Spearman co-expression)\n"
                 "directed = TF/focus edges (solid=regulon, dashed=adjacency); "
                 "green = target co-expression", fontsize=11)
    ax.axis("off"); fig.tight_layout(); fig.savefig(out_png, dpi=200); plt.close(fig)
    print(f"  plotted augmented network -> {out_png}", flush=True)


def discover_labels():
    if not os.path.isdir(SCENIC_DIR):
        return []
    out = []
    for d in sorted(os.listdir(SCENIC_DIR)):
        p = os.path.join(SCENIC_DIR, d)
        if os.path.isdir(p) and find_adj(p) is not None:
            out.append(d)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", nargs="+", default=None,
                    help="labels (subfolders under SCENIC/) to build; default = auto-discover")
    ap.add_argument("--nodes", help="plot-only mode: an *_nodes_augmented.csv to render")
    ap.add_argument("--edges", help="plot-only mode: an *_edges_augmented.csv to render")
    ap.add_argument("--out", default="augmented_network.png", help="plot-only output PNG")
    args = ap.parse_args()

    # plot-only mode: render an augmented graph (e.g. from 12d_plot_tf_network.R)
    if args.nodes and args.edges:
        plot_from_csv(args.nodes, args.edges, args.out)
        return

    labels = args.labels if args.labels else discover_labels()
    if not labels:
        sys.exit(f"[ERROR] no labels with adjacencies.tsv / adj.tsv under {SCENIC_DIR}. Run 12b first.")
    print(f"Focus gene: {FOCUS_GENE} | top up={TOP_UP} down={TOP_DOWN} | labels: {labels}",
          flush=True)
    for lab in labels:
        try:
            build_for_label(lab)
        except Exception as e:
            print(f"  [FAIL] {lab}: {e}", flush=True)
    print("\n=== 12d done === node/edge CSVs + networkx PNG under SCENIC/<label>/network/",
          flush=True)


if __name__ == "__main__":
    main()
