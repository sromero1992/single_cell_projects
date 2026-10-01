#!/usr/bin/env python
# =============================================================================
# 12b_run_pyscenic.py  --  pySCENIC (grn -> ctx -> aucell) on SELECTED looms
# =============================================================================
# Runs the full pySCENIC regulon workflow on one or more cell-type looms produced
# by 12a_export_for_pyscenic.R. Single-cell scRNA only -- no ATAC/peaks; regulons
# come from the pre-built cisTarget motif databases.
#
# Per selected label it writes, into SCENIC/<label>/:
#   adjacencies.tsv   (GRNBoost2 TF->target importances)
#   regulons.p        (motif-pruned regulons, pickled)
#   aucell.csv        (per-cell regulon activity matrix; cells x regulons)
#   <label>_pyscenic.loom (expression + regulons + AUCell, for SCope/AnnData)
#
# RUN (from the activated env; use the libstdc++ preload if you hit GLIBCXX):
#   conda activate scanpy_env_311
#   LD_PRELOAD=$CONDA_PREFIX/lib/libstdc++.so.6 python 12b_run_pyscenic.py
#   # or a subset:   python 12b_run_pyscenic.py --labels Tumor_epithelium Stem_cells
# =============================================================================
import os
import sys
import glob
import pickle
import argparse

import numpy as np
import pandas as pd
import loompy as lp

# =============================================================================
# --- CONFIG (portable; nr4a1 defaults, override via env vars) ----------------
# =============================================================================
# These env vars are the SAME ones config.R reads, so one `export` drives both the
# R and Python halves of the pipeline. Defaults target the Nr4a1 study:
#   NR4A1_ROOT       project r_process dir       (default: the optimus nr4a1 path)
#   NR4A1_OUTPUT     output dir                   (default: <ROOT>/seurat_output)
#   NR4A1_SCENIC_DIR SCENIC I/O dir               (default: <OUTPUT>/SCENIC)
#   NR4A1_CISTARGET  cisTarget database dir       (default: ~ssromerogon/cisTarget_databases)
#   NR4A1_SCENIC_WORKERS  dask/aucell workers     (default: 4)
ROOT_PATH   = os.environ.get(
    "NR4A1_ROOT",
    "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_nr4a1_ack/r_process")
OUTPUT_DIR  = os.environ.get("NR4A1_OUTPUT", os.path.join(ROOT_PATH, "seurat_output"))
SCENIC_DIR  = os.environ.get("NR4A1_SCENIC_DIR", os.path.join(OUTPUT_DIR, "SCENIC"))
DB_DIR      = os.environ.get("NR4A1_CISTARGET", "/home/ssromerogon/cisTarget_databases")

TF_FILE     = os.path.join(DB_DIR, "allTFs_mm.txt")
DB_FEATHERS = [
    os.path.join(DB_DIR, "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"),
    os.path.join(DB_DIR, "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"),
]
MOTIF_TBL   = os.path.join(DB_DIR, "motifs-v9-nr.mgi-m0.001-o0.0.tbl")

N_WORKERS   = int(os.environ.get("NR4A1_SCENIC_WORKERS", "4"))   # fewer = less RAM/teardown mess

# EXPORT_LOOM: also write <label>_pyscenic.loom via export2loom. This builds a
# t-SNE embedding (sklearn) and is the most segfault-prone / warning-noisy step;
# it is NOT needed downstream (aucell.csv + regulons.p are the essentials). Off.
EXPORT_LOOM = False

# Labels to run. Default None = AUTO-DISCOVER every SCENIC/<label>/<label>.loom
# (or *_expr_cells_x_genes.csv.gz) that 12a wrote -- so it always matches the
# actual folders (e.g. sex-split "Stem_cells__Female"). Override with --labels A B.
SELECTION = None

# =============================================================================
# --- Helpers -----------------------------------------------------------------
# =============================================================================
def load_expression(label):
    """Return a cells x genes DataFrame from SCENIC/<label>/<label>.loom (or csv.gz)."""
    d = os.path.join(SCENIC_DIR, label)
    loom = os.path.join(d, label + ".loom")
    if os.path.exists(loom):
        lf = lp.connect(loom, mode="r", validate=False)
        gk = next((k for k in ("Gene", "var_names", "gene_names", "GeneName", "Accession")
                   if k in lf.ra.keys()), list(lf.ra.keys())[0])
        ck = next((k for k in ("CellID", "obs_names", "cell_names", "CellName")
                   if k in lf.ca.keys()), list(lf.ca.keys())[0])
        genes = lf.ra[gk]; cells = lf.ca[ck]
        expr = pd.DataFrame(lf[:, :], index=np.asarray(genes), columns=np.asarray(cells)).T
        lf.close()
        # Force plain python str gene names (avoids numpy.str_ vs str mixed-type
        # index errors in ctx/df2regulons with some pandas/pyarrow combos).
        expr.columns = [str(g) for g in expr.columns]
        expr.index   = [str(c) for c in expr.index]
        return expr
    csvs = glob.glob(os.path.join(d, "*_expr_cells_x_genes.csv.gz"))
    if csvs:
        return pd.read_csv(csvs[0], index_col=0)  # already cells x genes
    raise FileNotFoundError(f"No loom or csv.gz for '{label}' in {d}")


def run_label(label, tf_names, dbs, client):
    from arboreto.algo import grnboost2
    from pyscenic.utils import modules_from_adjacencies
    from pyscenic.prune import prune2df, df2regulons
    from pyscenic.aucell import aucell

    d = os.path.join(SCENIC_DIR, label)
    print(f"\n########## {label} ##########", flush=True)
    expr = load_expression(label)                       # cells x genes
    print(f"  expression: {expr.shape[0]} cells x {expr.shape[1]} genes", flush=True)

    # Only regulators present in the matrix -> avoids arboreto's
    # "Must supply at least one delayed object" when TF names don't overlap.
    genes_in = set(map(str, expr.columns))
    tfs_in   = [t for t in tf_names if t in genes_in]
    print(f"  TFs present in matrix: {len(tfs_in)} / {len(tf_names)}", flush=True)
    if len(tfs_in) == 0:
        print("  [SKIP] No TF names overlap the gene columns -- symbol case/format mismatch.")
        print("    sample genes:", list(map(str, expr.columns[:8])), flush=True)
        print("    sample TFs:  ", tf_names[:8], flush=True)
        return

    # 1) GRN: TF -> target importances (GRNBoost2)
    adj_path = os.path.join(d, "adjacencies.tsv")
    adjacencies = grnboost2(expression_data=expr, tf_names=tfs_in,
                            client_or_address=client, verbose=True)
    adjacencies.to_csv(adj_path, sep="\t", index=False)
    print(f"  [grn]    {len(adjacencies)} TF-target links -> {adj_path}", flush=True)

    # 2) CTX: build co-expression modules, prune to motif-supported regulons
    modules = list(modules_from_adjacencies(adjacencies, expr))
    df = prune2df(dbs, modules, MOTIF_TBL, num_workers=N_WORKERS)
    regulons = df2regulons(df)
    with open(os.path.join(d, "regulons.p"), "wb") as fh:
        pickle.dump(regulons, fh)
    print(f"  [ctx]    {len(regulons)} regulons -> regulons.p", flush=True)

    # 3) AUCell: per-cell regulon activity
    auc_mtx = aucell(expr, regulons, num_workers=N_WORKERS)
    auc_mtx.to_csv(os.path.join(d, "aucell.csv"))
    print(f"  [aucell] {auc_mtx.shape[1]} regulons x {auc_mtx.shape[0]} cells -> aucell.csv",
          flush=True)

    # Optional combined loom (expression + regulons + AUCell). Off by default:
    # this is the t-SNE/loompy step that is slow, warning-noisy and segfault-prone,
    # and it is not needed downstream (aucell.csv + regulons.p are the essentials).
    if EXPORT_LOOM:
        try:
            from pyscenic.export import export2loom
            export2loom(ex_mtx=expr, regulons=regulons,
                        out_fname=os.path.join(d, f"{label}_pyscenic.loom"),
                        title=label, num_workers=N_WORKERS)
            print(f"  [loom]   {label}_pyscenic.loom", flush=True)
        except Exception as e:
            print(f"  [loom]   skipped ({e})", flush=True)

    # Quick peek: is Nr4a1 a recovered regulon?
    reg_names = [r.name for r in regulons]
    hit = [r for r in reg_names if r.lower().startswith("nr4a1")]
    print("  Nr4a1 regulon(s):", hit if hit else "none recovered here", flush=True)


def discover_labels():
    """Every SCENIC/<label>/ that contains <label>.loom or a *_expr_cells_x_genes.csv.gz."""
    out = []
    if not os.path.isdir(SCENIC_DIR):
        return out
    for d in sorted(os.listdir(SCENIC_DIR)):
        p = os.path.join(SCENIC_DIR, d)
        if os.path.isdir(p) and (
                os.path.exists(os.path.join(p, d + ".loom")) or
                glob.glob(os.path.join(p, "*_expr_cells_x_genes.csv.gz"))):
            out.append(d)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", nargs="+", default=SELECTION,
                    help="labels (subfolders under SCENIC/) to run; default = auto-discover")
    args = ap.parse_args()

    labels = args.labels if args.labels else discover_labels()
    if not labels:
        sys.exit(f"[ERROR] No looms found under {SCENIC_DIR}. Run 12a first (check SPLIT_BY).")

    # Sanity on inputs
    for f in [TF_FILE, MOTIF_TBL, *DB_FEATHERS]:
        if not os.path.exists(f):
            sys.exit(f"[ERROR] Missing database file: {f}")

    tf_names = [ln.strip() for ln in open(TF_FILE) if ln.strip()]
    print(f"TFs: {len(tf_names)} | databases: {len(DB_FEATHERS)} | labels: {labels}",
          flush=True)

    from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
    dbs = [RankingDatabase(fname=f, name=os.path.basename(f).replace(".feather", ""))
           for f in DB_FEATHERS]

    from distributed import Client, LocalCluster
    cluster = LocalCluster(n_workers=N_WORKERS, threads_per_worker=1)
    client = Client(cluster)
    try:
        for label in labels:
            try:
                run_label(label, tf_names, dbs, client)
            except Exception as e:
                print(f"  [FAIL] {label}: {e}", flush=True)
    finally:
        client.close()
        cluster.close()
    print("\n=== pySCENIC done ===", flush=True)


if __name__ == "__main__":
    main()
