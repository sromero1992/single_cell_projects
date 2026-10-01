# =============================================================================
# config.R  --  shared, portable configuration for the unified scRNA-seq pipeline
# =============================================================================
# SINGLE source of truth for machine/project-specific settings. Every numbered
# script sources this first (it looks for ./config.R in the working directory,
# or the path in the NR4A1_CONFIG environment variable).
#
# Defaults target the Nr4a1 colon KO study ("Nr4a1_s17_ack"). To run a DIFFERENT
# dataset or on a different machine you have two options:
#   (A) edit the literal defaults below once, OR
#   (B) change nothing here and export environment variables (best for HPC / for
#       sharing the same code across setups):
#
#     export NR4A1_PROJECT=MyProjectName
#     export NR4A1_ROOT=/path/to/project/r_process
#     export NR4A1_OUTPUT=/path/to/output            # default: <ROOT>/seurat_output
#     export NR4A1_CISTARGET=/path/to/cisTarget_databases
#     export NR4A1_PY_SCENIC=pyscenic                # conda env for 12b (pySCENIC)
#     export NR4A1_PY_CELLRANK=scanpy_env_311        # conda env for 10 (CellRank)
#
# Example — the Wu diet project that scripts 02-06/08 previously pointed at:
#     export NR4A1_PROJECT=Wu_Diet_project2
#     export NR4A1_ROOT=/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_wu_project2/r_process
#
# NOTE: this file sets ONLY the machine/project knobs (project name, root/output
# dirs, external DB dir, python env names). Which .rds each script READS (e.g.
# _unified_annotated vs _with_cell_scores) stays in that script, since it is a
# per-step choice, not a machine setting.
# =============================================================================

# --- Project + root (nr4a1 defaults) -----------------------------------------
PROJECT_NAME  <- Sys.getenv("NR4A1_PROJECT", "Nr4a1_s17_ack")
ROOT_PATH     <- Sys.getenv(
  "NR4A1_ROOT",
  "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_nr4a1_ack/r_process"
)
OUTPUT_DIR    <- Sys.getenv("NR4A1_OUTPUT", file.path(ROOT_PATH, "seurat_output"))

# --- External databases / python environments --------------------------------
CISTARGET_DIR   <- Sys.getenv("NR4A1_CISTARGET", "/home/ssromerogon/cisTarget_databases")
PY_ENV_SCENIC   <- Sys.getenv("NR4A1_PY_SCENIC",   "pyscenic")        # used by 12b
PY_ENV_CELLRANK <- Sys.getenv("NR4A1_PY_CELLRANK", "scanpy_env_311")  # used by 10

# --- Guards + banner ----------------------------------------------------------
if (!nzchar(PROJECT_NAME)) stop("config.R: PROJECT_NAME is empty.", call. = FALSE)
if (!nzchar(ROOT_PATH))    stop("config.R: ROOT_PATH is empty.",    call. = FALSE)
if (!dir.exists(OUTPUT_DIR))
  message("[config] NOTE: OUTPUT_DIR does not exist yet (will be created on write): ",
          OUTPUT_DIR)

.NR4A1_CONFIG_LOADED <- TRUE
message("[config] PROJECT_NAME = ", PROJECT_NAME)
message("[config] ROOT_PATH    = ", ROOT_PATH)
message("[config] OUTPUT_DIR   = ", OUTPUT_DIR)
