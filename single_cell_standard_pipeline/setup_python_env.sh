#!/bin/bash
# =============================================================================
# scRNA-seq PIPELINE - PYTHON ENVIRONMENT SETUP (Script 10 / CellRank)
# Version: 1.0
#
# PURPOSE:
#   Add CellRank and its companions to an EXISTING conda environment, then
#   verify the result and print the exact reticulate lines for R.
#
#   Only Script 10 (trajectory inference) needs Python. Scripts 00-09 are
#   pure R - including CytoTRACE 2 in Script 09, which is an R package.
#
#   CPU only. Nothing here touches CUDA, and nothing in this pipeline uses a
#   GPU: CellRank runs on pyGPCCA/scipy and CytoTRACE 2 runs in R.
#
# USAGE:
#   chmod +x setup_python_env.sh
#   ./setup_python_env.sh                 # install into the env named below
#   ./setup_python_env.sh --check         # verify only, install nothing
#   ./setup_python_env.sh --env other_env # target a different environment
#
#   Safe to re-run. pip skips anything already satisfied.
# =============================================================================

set -uo pipefail

# --- CONFIGURATION -----------------------------------------------------------
ENV_NAME="scanpy_env_311"

# Core requirement for Script 10.
PKGS_CORE=(
  "cellrank>=2.0"      # trajectory inference; provides CytoTRACEKernel
  "scvelo"             # supplies scv.pp.moments, the default CT_IMPUTE_METHOD
)

# Usually already present in a scanpy env; listed so a bare env also works.
PKGS_BASE=(
  "scanpy"
  "anndata"
  "igraph"
  "leidenalg"
)

# Optional: only needed if you set CT_IMPUTE_METHOD <- "cellmapper" in Script 10.
PKGS_OPTIONAL=(
  "cellmapper"
)

INSTALL_OPTIONAL="no"    # "yes" to also install PKGS_OPTIONAL

# =============================================================================
# --- ARGUMENT PARSING --------------------------------------------------------
# =============================================================================
CHECK_ONLY="no"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --check)     CHECK_ONLY="yes"; shift ;;
    --optional)  INSTALL_OPTIONAL="yes"; shift ;;
    --env)       ENV_NAME="$2"; shift 2 ;;
    -h|--help)
      grep '^#' "$0" | sed 's/^# \{0,1\}//' | head -25
      exit 0 ;;
    *) echo "Unknown option: $1  (try --help)"; exit 1 ;;
  esac
done

echo "============================================================"
echo " Python environment setup for Script 10 (CellRank)"
echo "============================================================"
echo " Target env : ${ENV_NAME}"
echo " Mode       : $([ "$CHECK_ONLY" = "yes" ] && echo 'CHECK ONLY' || echo 'INSTALL')"
echo " GPU/CUDA   : not used by this pipeline"
echo "============================================================"
echo

# =============================================================================
# --- PRE-FLIGHT --------------------------------------------------------------
# =============================================================================
if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found on PATH."
  echo "       Run this from a shell where 'conda env list' works."
  exit 1
fi

if ! conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "ERROR: conda environment '${ENV_NAME}' not found. Available:"
  conda env list | sed 's/^/    /'
  exit 1
fi

ENV_PY="$(conda run -n "${ENV_NAME}" python -c 'import sys; print(sys.executable)' 2>/dev/null)"
if [ -z "${ENV_PY}" ]; then
  echo "ERROR: could not resolve the python interpreter inside '${ENV_NAME}'."
  exit 1
fi
echo "Interpreter : ${ENV_PY}"
conda run -n "${ENV_NAME}" python -c 'import sys; print("Python      :", sys.version.split()[0])'
echo

# =============================================================================
# --- INSTALL -----------------------------------------------------------------
# =============================================================================
if [ "$CHECK_ONLY" = "no" ]; then

  # Cheap insurance: snapshot the env so it can be restored if a dependency
  # resolution goes badly. Costs a second, saves an afternoon.
  BACKUP="env_snapshot_${ENV_NAME}_$(date +%Y%m%d_%H%M%S).txt"
  echo "--- Snapshotting current package list -> ${BACKUP}"
  conda run -n "${ENV_NAME}" python -m pip freeze > "${BACKUP}" 2>/dev/null \
    && echo "    saved ($(wc -l < "${BACKUP}") packages)" \
    || echo "    [WARNING] snapshot failed; continuing"
  echo "    To roll back:  conda run -n ${ENV_NAME} python -m pip install -r ${BACKUP}"
  echo

  ALL_PKGS=("${PKGS_BASE[@]}" "${PKGS_CORE[@]}")
  if [ "$INSTALL_OPTIONAL" = "yes" ]; then
    ALL_PKGS+=("${PKGS_OPTIONAL[@]}")
  fi

  echo "--- Installing: ${ALL_PKGS[*]}"
  echo "    (pip skips anything already satisfied)"
  echo

  if conda run --no-capture-output -n "${ENV_NAME}" \
       python -m pip install --upgrade "${ALL_PKGS[@]}"; then
    echo
    echo "--- pip completed"
  else
    echo
    echo "[ERROR] pip install failed. The environment is unchanged for any"
    echo "        package that did not install. Common causes:"
    echo "          - numpy/scipy ABI conflict with existing packages"
    echo "          - no network access from this node"
    echo "        Restore with: conda run -n ${ENV_NAME} python -m pip install -r ${BACKUP}"
    exit 1
  fi
  echo
fi

# =============================================================================
# --- VERIFY ------------------------------------------------------------------
# =============================================================================
echo "============================================================"
echo " VERIFICATION"
echo "============================================================"

conda run --no-capture-output -n "${ENV_NAME}" python - <<'PYEOF'
import importlib, sys, warnings

# scanpy/anndata deprecate module-level __version__, so reading it emits a
# FutureWarning. Ask importlib.metadata first and only fall back.
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

try:
    from importlib.metadata import version as _dist_version, PackageNotFoundError
except ImportError:
    _dist_version = None

REQUIRED = ["scanpy", "anndata", "cellrank", "numpy", "scipy", "pandas"]
IMPUTE   = ["scvelo"]          # CT_IMPUTE_METHOD = "moments"
OPTIONAL = ["cellmapper"]      # CT_IMPUTE_METHOD = "cellmapper"

def get_version(mod, name):
    if _dist_version is not None:
        try:
            return _dist_version(name)
        except Exception:
            pass
    return getattr(mod, "__version__", "?")

def check(mods, label):
    ok = True
    print(f"\n{label}")
    for m in mods:
        try:
            mod = importlib.import_module(m)
            print(f"  [OK]      {m:<12} {get_version(mod, m)}")
        except Exception as e:
            print(f"  [MISSING] {m:<12} {type(e).__name__}: {e}")
            ok = False
    return ok

req_ok = check(REQUIRED, "Required by Script 10:")
imp_ok = check(IMPUTE,   "Imputation backend (CT_IMPUTE_METHOD = 'moments'):")
check(OPTIONAL, "Optional (CT_IMPUTE_METHOD = 'cellmapper'):")

# The specific thing Script 10 calls. CellRank v1 does not have this kernel,
# so an old install passes the import check above but fails at run time.
print("\nCytoTRACEKernel availability:")
try:
    import cellrank as cr
    if hasattr(cr.kernels, "CytoTRACEKernel"):
        print(f"  [OK]      cellrank {cr.__version__} exposes CytoTRACEKernel")
    else:
        print(f"  [FAIL]    cellrank {cr.__version__} has NO CytoTRACEKernel")
        print("            Upgrade: pip install --upgrade 'cellrank>=2.0'")
        req_ok = False
except Exception as e:
    print(f"  [FAIL]    {type(e).__name__}: {e}")
    req_ok = False

print("\n" + "="*60)
if req_ok and imp_ok:
    print("READY - Script 10 can run against this environment.")
elif req_ok:
    print("MOSTLY READY - core is fine, but the 'moments' imputation backend")
    print("is missing. Install scvelo, or set CT_IMPUTE_METHOD <- \"none\".")
else:
    print("NOT READY - see the [MISSING]/[FAIL] lines above.")
print("="*60)
sys.exit(0 if req_ok else 1)
PYEOF

VERIFY_RC=$?
echo

# =============================================================================
# --- R SIDE ------------------------------------------------------------------
# =============================================================================
echo "============================================================"
echo " NEXT: WIRE THIS INTO R"
echo "============================================================"
cat <<EOF
1. Install reticulate once (R side only, tiny):

     install.packages("reticulate")

2. Confirm R can see this environment. IMPORTANT: reticulate binds an
   interpreter ONCE per R session, so run this in a FRESH session before
   anything else imports Python:

     library(reticulate)
     use_condaenv("${ENV_NAME}", required = TRUE)
     py_config()
     py_module_available("cellrank")   # must be TRUE

3. Script 10 is already configured for it:

     PYTHON_ENV_TYPE <- "conda"
     PYTHON_ENV      <- "${ENV_NAME}"

   Script 10 verifies every module before doing any work, so a
   misconfiguration fails in seconds rather than after loading the object.

4. On an HPC batch job, this override is the most reliable binding:

     export RETICULATE_PYTHON=${ENV_PY}

REMINDER: Scripts 00-09 need none of this. CytoTRACE 2 in Script 09 is an
R package - install it with:

     devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r")
EOF
echo "============================================================"

exit ${VERIFY_RC}
