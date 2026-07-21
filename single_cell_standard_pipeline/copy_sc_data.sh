#!/bin/bash
# =============================================================================
# scRNA-seq PIPELINE - DATA STAGING SCRIPT (UNIFIED)
# Version: 2.0
#
# PURPOSE:
#   Pull the per-sample Cell Ranger outputs (H5 matrices, and optionally probe
#   matrices and web summaries) out of the mapped-data tree and stage them into
#   a flat h5_files/ directory laid out exactly the way Script 01 expects:
#
#     h5_files/<SampleID>/sample_filtered_feature_bc_matrix.h5
#     h5_files/<SampleID>/sample_raw_probe_bc_matrix.h5   (optional)
#     h5_files/<SampleID>/web_summary.html                (optional)
#
#   <SampleID> must match the SampleID column of your metadata .xlsx file.
#
# UNIFIED BUILD PROVENANCE:
#   Merges github_pipeline/copy_sc_data.sh (probe-matrix copying, Safe pools)
#   and eithan_coffee/copy_sc_data.sh (web_summary.html copying, Eitan pools)
#   into a single parameterised script.
#
#   Fixes carried over from those versions:
#     - eithan_coffee/copy_sc_data.sh had `path2` commented out but still ran
#       `cp -rv ${path2}` on the empty variable, producing a cp error on every
#       sample. Optional files are now guarded by explicit toggles.
#     - Added `set -uo pipefail`, existence checks, and a per-sample success/
#       failure tally, so a silently-missing sample cannot slip through
#       unnoticed into the R pipeline.
#
# USAGE:
#   1. Edit the CONFIGURATION block below.
#   2. chmod +x copy_sc_data.sh
#   3. ./copy_sc_data.sh
#   Optionally run with DRY_RUN=1 ./copy_sc_data.sh to preview without copying.
# =============================================================================

set -uo pipefail

# =============================================================================
# --- CONFIGURATION (EDIT THIS BLOCK) -----------------------------------------
# =============================================================================

# SOURCE_ROOT: top-level mapped-data directory containing the pool folders.
SOURCE_ROOT="/mnt/matrix/tigger/2026_Safe_Study17_Sc/mapped-data"

# POOLS: the pool/run folder names to pull from, as an array.
#   Nr4a1 / Study 17 : ("Safe_Pool53" "Safe_Pool54" "Safe_Pool55" "Safe_Pool56")
#   Eitan coffee     : ("Eitan-Pool57" "Eitan-Pool58" "Eitan-Pool59" "Eitan-Pool60")
POOLS=("Safe_Pool53" "Safe_Pool54" "Safe_Pool55" "Safe_Pool56")

# SUBDIR: path inside each pool folder that holds one directory per sample.
SUBDIR="outs/per_sample_outs"

# DEST_DIR: final staging directory name. Point H5_DIR in Script 01 here.
DEST_DIR="h5_files"

# --- What to copy ------------------------------------------------------------
# COPY_FILTERED : the main gene-expression matrix. Required by Script 01.
# COPY_PROBE    : probe matrix. Set "yes" only if this run used a probe assay
#                 (i.e. ADD_PROBE_DATA <- TRUE in Script 01).
# COPY_SUMMARY  : Cell Ranger web_summary.html. Cheap and very useful for QC.
COPY_FILTERED="yes"
COPY_PROBE="yes"
COPY_SUMMARY="yes"

# OVERWRITE: "no" skips files already staged (safe to re-run); "yes" re-copies.
OVERWRITE="no"

# =============================================================================
# --- EXECUTION (DO NOT EDIT BELOW THIS LINE) ---------------------------------
# =============================================================================

FILTERED_NAME="sample_filtered_feature_bc_matrix.h5"
PROBE_NAME="sample_raw_probe_bc_matrix.h5"
SUMMARY_NAME="web_summary.html"

DRY_RUN="${DRY_RUN:-0}"

n_samples=0
n_copied=0
n_skipped=0
declare -a MISSING_REQUIRED=()
declare -a MISSING_OPTIONAL=()

echo "============================================================"
echo " scRNA-seq data staging"
echo "============================================================"
echo " Source root : ${SOURCE_ROOT}"
echo " Pools       : ${POOLS[*]}"
echo " Destination : ${DEST_DIR}"
echo " Copying     : filtered=${COPY_FILTERED} probe=${COPY_PROBE} summary=${COPY_SUMMARY}"
[ "${DRY_RUN}" = "1" ] && echo " MODE        : DRY RUN (nothing will be copied)"
echo "============================================================"
echo

# Fail early if the source root is wrong - the most common misconfiguration.
if [ ! -d "${SOURCE_ROOT}" ]; then
  echo "ERROR: SOURCE_ROOT does not exist or is not readable:"
  echo "       ${SOURCE_ROOT}"
  echo "Check the path and that the storage volume is mounted."
  exit 1
fi

mkdir -p "${DEST_DIR}"

# copy_one <source_file> <target_dir> <required:yes|no> <sample_label>
copy_one() {
  local src="$1" tgt="$2" required="$3" label="$4"
  local base
  base="$(basename "${src}")"

  if [ ! -f "${src}" ]; then
    if [ "${required}" = "yes" ]; then
      echo "        [MISSING-REQUIRED] ${base}"
      MISSING_REQUIRED+=("${label}/${base}")
    else
      echo "        [missing-optional]  ${base}"
      MISSING_OPTIONAL+=("${label}/${base}")
    fi
    return 1
  fi

  if [ -f "${tgt}/${base}" ] && [ "${OVERWRITE}" = "no" ]; then
    echo "        [skip, exists]      ${base}"
    n_skipped=$((n_skipped + 1))
    return 0
  fi

  if [ "${DRY_RUN}" = "1" ]; then
    echo "        [would copy]        ${base}"
  else
    if cp -f "${src}" "${tgt}/"; then
      echo "        [copied]            ${base}"
      n_copied=$((n_copied + 1))
    else
      echo "        [COPY FAILED]       ${base}"
      MISSING_REQUIRED+=("${label}/${base} (copy failed)")
      return 1
    fi
  fi
  return 0
}

for pool in "${POOLS[@]}"; do
  pool_path="${SOURCE_ROOT}/${pool}/${SUBDIR}"
  echo "--- Pool: ${pool}"

  if [ ! -d "${pool_path}" ]; then
    echo "    [WARNING] Not found, skipping pool: ${pool_path}"
    echo
    continue
  fi

  # Iterate sample directories directly rather than parsing `ls` output, so
  # that sample names containing spaces or unusual characters are handled.
  shopt -s nullglob
  sample_dirs=("${pool_path}"/*/)
  shopt -u nullglob

  if [ ${#sample_dirs[@]} -eq 0 ]; then
    echo "    [WARNING] No sample directories inside ${pool_path}"
    echo
    continue
  fi

  for sample_path in "${sample_dirs[@]}"; do
    sample="$(basename "${sample_path}")"
    echo "    Sample: ${sample}"
    n_samples=$((n_samples + 1))

    target="${DEST_DIR}/${sample}"
    [ "${DRY_RUN}" = "1" ] || mkdir -p "${target}"

    [ "${COPY_FILTERED}" = "yes" ] && \
      copy_one "${sample_path}${FILTERED_NAME}" "${target}" "yes" "${sample}"
    [ "${COPY_PROBE}" = "yes" ] && \
      copy_one "${sample_path}${PROBE_NAME}"    "${target}" "no"  "${sample}"
    [ "${COPY_SUMMARY}" = "yes" ] && \
      copy_one "${sample_path}${SUMMARY_NAME}"  "${target}" "no"  "${sample}"
  done
  echo
done

# =============================================================================
# --- SUMMARY -----------------------------------------------------------------
# =============================================================================
echo "============================================================"
echo " STAGING SUMMARY"
echo "============================================================"
echo " Samples found   : ${n_samples}"
echo " Files copied    : ${n_copied}"
echo " Files skipped   : ${n_skipped} (already present; OVERWRITE=no)"
echo " Staged into     : ${DEST_DIR}/"

if [ ${#MISSING_OPTIONAL[@]} -gt 0 ]; then
  echo
  echo " Optional files not found (${#MISSING_OPTIONAL[@]}):"
  printf '   - %s\n' "${MISSING_OPTIONAL[@]}"
  echo "   If probe matrices are listed here, set ADD_PROBE_DATA <- FALSE in Script 01."
fi

if [ ${#MISSING_REQUIRED[@]} -gt 0 ]; then
  echo
  echo " *** REQUIRED FILES MISSING (${#MISSING_REQUIRED[@]}) ***"
  printf '   - %s\n' "${MISSING_REQUIRED[@]}"
  echo
  echo " Script 01 will FAIL on these samples. Resolve before running R."
  echo "============================================================"
  exit 1
fi

echo
echo " All required files staged successfully."
echo
echo " NEXT STEPS:"
echo "   1. Confirm every SampleID in your metadata .xlsx has a matching"
echo "      folder inside ${DEST_DIR}/  (names must match EXACTLY)."
echo "   2. Set H5_DIR in 01_process_data.R to point at ${DEST_DIR}/"
echo "   3. Run 00_rlibs_installation.R once, then 01_process_data.R"
echo "============================================================"
