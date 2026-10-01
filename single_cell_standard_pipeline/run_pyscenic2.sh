#!/usr/bin/env bash
# =============================================================================
# run_pyscenic2.sh  --  stable runner for 12b_run_pyscenic.py
# =============================================================================
# Activates the pySCENIC conda env and applies the stability settings we learned
# the hard way, then runs 12b. Any extra args pass THROUGH to 12b (e.g. --labels
# Stem_cells__Female). Portable: honours the SAME env vars as config.R / 12b;
# defaults target the Nr4a1 study.
#
#   NR4A1_PY_SCENIC        conda env name          (default: pyscenic)
#   NR4A1_ROOT / NR4A1_OUTPUT / NR4A1_SCENIC_DIR / NR4A1_CISTARGET   (read by 12b)
#   NR4A1_SCENIC_WORKERS   dask / aucell workers   (default: 4)
#   NR4A1_PYSCENIC_SCRIPT  path to 12b             (default: alongside this file)
#
# USAGE:
#   ./run_pyscenic2.sh                 # auto-discover + run every label
#   ./run_pyscenic2.sh --labels A B    # run a subset
# =============================================================================
set -uo pipefail

ENV_NAME="${NR4A1_PY_SCENIC:-pyscenic}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="${NR4A1_PYSCENIC_SCRIPT:-$HERE/12b_run_pyscenic.py}"

if [ ! -f "$PY_SCRIPT" ]; then
  echo "ERROR: 12b script not found at '$PY_SCRIPT' (set NR4A1_PYSCENIC_SCRIPT)."; exit 1
fi

# Single-threaded BLAS/OMP: the Python API's num_workers does NOT cap these, and
# oversubscription alongside the dask workers is what triggered the OOM/segfaults.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export NR4A1_SCENIC_WORKERS="${NR4A1_SCENIC_WORKERS:-4}"

# Activate conda in a non-interactive shell.
if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not on PATH. Run from a shell where 'conda' works."; exit 1
fi
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME" || { echo "ERROR: cannot activate conda env '$ENV_NAME'."; exit 1; }

# libstdc++ GLIBCXX_3.4.30 fix: the env's newer libstdc++ must be preloaded
# because the system one is too old for scipy/arboreto.
if [ -f "$CONDA_PREFIX/lib/libstdc++.so.6" ]; then
  export LD_PRELOAD="$CONDA_PREFIX/lib/libstdc++.so.6${LD_PRELOAD:+:$LD_PRELOAD}"
fi

echo "[run_pyscenic2] env=$ENV_NAME  workers=$NR4A1_SCENIC_WORKERS"
echo "[run_pyscenic2] script=$PY_SCRIPT"
echo "[run_pyscenic2] $(python --version 2>&1) @ $(command -v python)"
exec python "$PY_SCRIPT" "$@"
