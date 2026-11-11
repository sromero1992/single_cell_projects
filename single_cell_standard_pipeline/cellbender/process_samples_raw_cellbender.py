#!/usr/bin/env python3
# =============================================================================
# process_samples_raw_cellbender.py
#   Per-sample ambient-RNA removal with CellBender `remove-background`.
# =============================================================================
# PORTABLE: nothing project-specific is hard-coded. Every path/column/parameter
# comes from an ENV VAR (with a sensible default) or a --flag (flag > env >
# default). Reads a metadata .xlsx that has at least a sample-id column, runs
# CellBender on each sample's RAW 10x matrix, and writes the corrected output
# (optionally also an .h5ad, and/or a 10x-layout filtered .h5 that Script 01 can
# ingest). Already-finished samples are skipped so re-runs resume.
#
# LAYOUT EXPECTED (same as copy_sc_data.sh / Script 01):
#   <H5_DIR>/<SampleID>/<RAW_H5>        e.g. h5_files/<id>/sample_raw_feature_bc_matrix.h5
#
# CONFIG (env var  |  flag  |  default):
#   CB_SOURCE          --source            working dir                 (cwd)
#   CB_METADATA        --metadata          metadata .xlsx              (<source>/metadata.xlsx)
#   CB_SAMPLE_COL      --sample-col        sample-id column            (SampleID)
#   CB_H5_DIR          --h5-dir            per-sample 10x dir          (h5_files)
#   CB_RAW_H5          --raw-h5            raw matrix filename         (sample_raw_feature_bc_matrix.h5)
#   CB_OUT             --out              output root                 (<source>/cellbender_samples)
#   CB_BIN             --cellbender        cellbender executable       (cellbender)
#   CB_CUDA            --cuda / --cpu      GPU on/off                  (auto: GPU if torch sees one)
#   CB_EXPECTED_CELLS  --expected-cells    expected real cells         (20000; "" -> CellBender estimates)
#   CB_TOTAL_DROPLETS  --total-droplets    total droplets included     (50000; "" -> CellBender estimates)
#   CB_EPOCHS          --epochs                                         (150)
#   CB_FPR             --fpr                                            (0.01)
#   CB_LR              --learning-rate                                  (1e-4)
#   CB_EXTRA_ARGS      --extra-args        extra raw args for cellbender ("")
#   CB_WRITE_H5AD      --h5ad / --no-h5ad  also write .h5ad             (on if scanpy present)
#   CB_PIPELINE_H5     --pipeline-h5       also drop a 10x filtered h5 into
#                                          <H5_DIR>/<id>/sample_cellbender_filtered.h5  (off)
#
# PER-SAMPLE OVERRIDES: if the metadata has columns `expected_cells` and/or
# `total_droplets`, those values are used per sample (blank -> global default).
# Any extra metadata columns are carried onto the .h5ad `.obs` automatically.
#
# ENV TO RUN IT:  see environment_cellbender.yml in this folder (CellBender +
# PyTorch + scanpy). The scanpy_env_311 env does NOT contain CellBender.
# =============================================================================
import argparse
import os
import shutil
import subprocess
import sys
from datetime import datetime

import pandas as pd


def _b(env, default):
    v = os.environ.get(env)
    if v is None:
        return default
    return str(v).strip().lower() in ("1", "true", "yes", "on")


def parse_args():
    p = argparse.ArgumentParser(description="Per-sample CellBender remove-background.")
    p.add_argument("--source",         default=os.environ.get("CB_SOURCE", os.getcwd()))
    p.add_argument("--metadata",       default=os.environ.get("CB_METADATA"))
    p.add_argument("--sample-col",     default=os.environ.get("CB_SAMPLE_COL", "SampleID"))
    p.add_argument("--h5-dir",         default=os.environ.get("CB_H5_DIR", "h5_files"))
    p.add_argument("--raw-h5",         default=os.environ.get("CB_RAW_H5", "sample_raw_feature_bc_matrix.h5"))
    p.add_argument("--out",            default=os.environ.get("CB_OUT"))
    p.add_argument("--cellbender",     default=os.environ.get("CB_BIN", "cellbender"))
    p.add_argument("--expected-cells", default=os.environ.get("CB_EXPECTED_CELLS", "20000"))
    p.add_argument("--total-droplets", default=os.environ.get("CB_TOTAL_DROPLETS", "50000"))
    p.add_argument("--epochs",         default=os.environ.get("CB_EPOCHS", "150"))
    p.add_argument("--fpr",            default=os.environ.get("CB_FPR", "0.01"))
    p.add_argument("--learning-rate",  default=os.environ.get("CB_LR", "1e-4"))
    p.add_argument("--extra-args",     default=os.environ.get("CB_EXTRA_ARGS", ""))
    cuda = p.add_mutually_exclusive_group()
    cuda.add_argument("--cuda", dest="cuda", action="store_true", default=None)
    cuda.add_argument("--cpu",  dest="cuda", action="store_false")
    h5ad = p.add_mutually_exclusive_group()
    h5ad.add_argument("--h5ad",    dest="h5ad", action="store_true", default=None)
    h5ad.add_argument("--no-h5ad", dest="h5ad", action="store_false")
    p.add_argument("--pipeline-h5", dest="pipeline_h5", action="store_true",
                   default=_b("CB_PIPELINE_H5", False))
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()
    if a.cuda is None:
        a.cuda = _b("CB_CUDA", _detect_gpu())
    if a.h5ad is None:
        a.h5ad = _b("CB_WRITE_H5AD", _have_scanpy())
    return a


def _detect_gpu():
    try:
        import torch
        return bool(torch.cuda.is_available())
    except Exception:
        return False


def _have_scanpy():
    try:
        import scanpy  # noqa: F401
        import anndata  # noqa: F401
        return True
    except Exception:
        return False


def main():
    a = parse_args()
    source = os.path.abspath(a.source)
    os.chdir(source)
    meta_file = a.metadata or os.path.join(source, "metadata.xlsx")
    out_root = a.out or os.path.join(source, "cellbender_samples")
    os.makedirs(out_root, exist_ok=True)
    log_path = os.path.join(source, "cellbender_run_log.txt")

    def log(msg):
        line = f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
        print(line, flush=True)
        with open(log_path, "a") as fh:
            fh.write(line + "\n")

    if not os.path.exists(meta_file):
        sys.exit(f"[ERROR] metadata not found: {meta_file} (set --metadata / CB_METADATA)")
    D = pd.read_excel(meta_file)
    if a.sample_col not in D.columns:
        sys.exit(f"[ERROR] sample column '{a.sample_col}' not in {meta_file}. "
                 f"Columns: {list(D.columns)}")
    D[a.sample_col] = D[a.sample_col].astype(str)

    log(f"CellBender run | source={source} | samples={len(D)} | "
        f"GPU={'on' if a.cuda else 'off (CPU, slow)'} | h5ad={a.h5ad} | out={out_root}")

    ok, failed, skipped = [], [], []
    for _, row in D.iterrows():
        sid = str(row[a.sample_col])
        fbase = "".join(c if c.isalnum() else "_" for c in sid)
        sdir = os.path.join(out_root, fbase)
        os.makedirs(sdir, exist_ok=True)

        raw_in = os.path.abspath(os.path.join(a.h5_dir, sid, a.raw_h5))
        cb_out = os.path.join(sdir, f"{fbase}_cellbender.h5")                 # cellbender writes this
        cb_filt = os.path.join(sdir, f"{fbase}_cellbender_filtered.h5")        # + a _filtered.h5
        h5ad_out = os.path.join(sdir, f"{fbase}_cellbender.h5ad")

        done = os.path.exists(h5ad_out) if a.h5ad else os.path.exists(cb_out)
        if done:
            log(f"[skip] {sid}: output exists"); skipped.append(sid); continue
        if not os.path.exists(raw_in):
            log(f"[skip] {sid}: raw matrix not found ({raw_in})"); failed.append(sid); continue

        # per-sample params (metadata columns override globals; blank -> default)
        exp_cells = _pick(row, "expected_cells", a.expected_cells)
        tot_drop  = _pick(row, "total_droplets", a.total_droplets)

        cmd = [a.cellbender, "remove-background", "--input", raw_in, "--output", cb_out,
               "--epochs", str(a.epochs), "--fpr", str(a.fpr),
               "--learning-rate", str(a.learning_rate), "--checkpoint-mins", "60"]
        if a.cuda:
            cmd.append("--cuda")
        if exp_cells:
            cmd += ["--expected-cells", str(int(float(exp_cells)))]
        if tot_drop:
            cmd += ["--total-droplets-included", str(int(float(tot_drop)))]
        if a.extra_args:
            cmd += a.extra_args.split()

        log(f"[run] {sid}: {' '.join(cmd)}")
        if a.dry_run:
            continue

        try:
            subprocess.run(cmd, check=True, cwd=sdir)
        except subprocess.CalledProcessError as e:
            log(f"[FAIL] {sid}: cellbender exited {e.returncode}"); failed.append(sid); continue
        except FileNotFoundError:
            sys.exit(f"[ERROR] '{a.cellbender}' not found. Activate the CellBender env "
                     f"or set --cellbender / CB_BIN to its path.")

        # rename the default checkpoint so it doesn't clash between samples
        ckpt = os.path.join(sdir, "ckpt.tar.gz")
        if os.path.exists(ckpt):
            os.replace(ckpt, os.path.join(sdir, f"{fbase}_ckpt.tar.gz"))

        # optional: drop a 10x filtered h5 into the pipeline layout for Script 01
        if a.pipeline_h5 and os.path.exists(cb_filt):
            dest = os.path.join(a.h5_dir, sid, "sample_cellbender_filtered.h5")
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            shutil.copyfile(cb_filt, dest)
            log(f"       pipeline h5 -> {dest}")

        # optional: convert to .h5ad, carrying ALL metadata columns onto .obs
        if a.h5ad:
            try:
                import numpy as np
                import scanpy as sc
                try:
                    import hdf5plugin
                    comp = dict(compression=hdf5plugin.FILTERS["zstd"])
                except Exception:
                    comp = dict(compression="gzip")
                adata = sc.read_10x_h5(cb_out)
                adata.X = adata.X.astype(np.float32)
                adata.var_names_make_unique()
                for col in D.columns:
                    adata.obs[col] = row[col]
                adata.write_h5ad(h5ad_out, **comp)
                log(f"       h5ad -> {h5ad_out}")
            except Exception as e:
                log(f"[warn] {sid}: h5ad conversion failed ({e}); CellBender .h5 kept.")

        ok.append(sid)

    log(f"DONE. ok={len(ok)} skipped={len(skipped)} failed={len(failed)}")
    if failed:
        log(f"  failed/missing: {', '.join(failed)}")


def _pick(row, col, default):
    if col in row and pd.notna(row[col]) and str(row[col]).strip() != "":
        return row[col]
    return default


if __name__ == "__main__":
    main()
