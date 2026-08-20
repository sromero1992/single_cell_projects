# =============================================================================
# smoke_test.R — verify TamuScDSC works on a fresh machine, with NO real data.
#
# Run AFTER 00_rlibs_installation.R (which installs deps + TamuScDSC). From R:
#   source("TamuScDSC/smoke_test.R")        # cwd = repo root
#
# It builds a tiny synthetic 2-sample dataset in memory and walks the pipeline.
# Each stage is guarded: heavy optional packages (scDblFinder / DoubletFinder /
# harmony) are exercised only if installed, and skipped with a note otherwise,
# so the core plumbing is validated even on a partial install.
# =============================================================================

ok   <- function(msg) cat(sprintf("  [PASS] %s\n", msg))
skip <- function(msg) cat(sprintf("  [SKIP] %s\n", msg))
step <- function(msg) cat(sprintf("\n== %s ==\n", msg))
stopifnot_msg <- function(cond, msg) { if (!isTRUE(cond)) stop("[FAIL] ", msg, call. = FALSE); ok(msg) }

step("0. Package loads")
suppressPackageStartupMessages(library(TamuScDSC))
cat("  TamuScDSC version:", as.character(utils::packageVersion("TamuScDSC")), "\n")

# --- 1. No-data contracts (never need Seurat/Bioconductor) -------------------
step("1. Schema + markers (no data required)")
s <- TamuScDSC_schema()
stopifnot_msg(s$df_score == "DF_score" && s$sc_score == "scDblFinder_score",
              "schema column names match 01b expectations")

m <- get_markers()
stopifnot_msg(nrow(m) > 0 && all(c("cell_type", "gene") %in% colnames(m)),
              "get_markers() reads the bundled CSV")
stopifnot_msg(!any(grepl("\\|", m$gene)),
              "pipe-delimited gene lists are exploded to one gene per row")
m2 <- add_markers(m, "Tuft", c("Dclk1", "Trpm5"))
stopifnot_msg(all(c("Dclk1", "Trpm5") %in% markers_as_list(m2)[["Tuft"]]),
              "add_markers + markers_as_list round-trip")

# --- 2. Synthetic 2-sample dataset -------------------------------------------
step("2. Build a tiny synthetic dataset")
if (!requireNamespace("Seurat", quietly = TRUE)) {
  skip("Seurat not installed — cannot run the data-dependent steps. Fix Script 00 first.")
} else {
  suppressPackageStartupMessages(library(Seurat))
  set.seed(1)
  make_sample <- function(sid, ncell = 120, ngene = 300) {
    counts <- matrix(rpois(ncell * ngene, lambda = 0.5), nrow = ngene, ncol = ncell)
    rownames(counts) <- c(paste0("Gene", seq_len(ngene - 5)),
                          paste0("mt-", c("Nd1", "Co1", "Cytb", "Atp6", "Nd4")))  # 5 mito genes
    colnames(counts) <- paste0(sid, "_cell", seq_len(ncell))
    so <- CreateSeuratObject(counts = counts, project = sid, min.cells = 0)
    so$SampleID <- sid
    so
  }
  samples <- list(A = make_sample("A"), B = make_sample("B"))
  ok(sprintf("built %d samples, %d cells each", length(samples), ncol(samples[[1]])))

  step("3. Ingest generic (all three scenarios reach the same list)")
  l_from_list   <- as_sample_list(samples)
  merged        <- merge_samples(samples)
  l_from_merged <- as_sample_list(merged, sample_col = "SampleID")
  stopifnot_msg(length(l_from_list) == 2 && length(l_from_merged) == 2,
                "as_sample_list() handles a list AND a merged object")

  step("4. Per-sample QC (labels, does not delete)")
  samples <- apply_qc(samples, p = list(min_features = 1, max_mt = 100), apply = FALSE)
  stopifnot_msg("qc_pass" %in% colnames(samples[[1]]@meta.data),
                "apply_qc writes qc_pass + percent_mt")

  step("5. Doublet detection (uses whichever caller is installed)")
  if (requireNamespace("scDblFinder", quietly = TRUE)) {
    samples <- detect_doublets(samples, method = "scDblFinder",
                               p = doublet_params(method = "scDblFinder"))
    stopifnot_msg(all(c("scDblFinder_score", "scDblFinder_class", "doublet_consensus")
                      %in% colnames(samples[[1]]@meta.data)),
                  "detect_doublets writes score/class/consensus columns")
  } else {
    skip("scDblFinder not installed — doublet step not exercised (BiocManager::install('scDblFinder')).")
  }

  step("6. Merge + dataset QC + integrate")
  data <- merge_samples(samples)
  data <- apply_qc(data, p = list(min_features = 1, max_mt = 100, min_cells_per_gene = 0), apply = TRUE)
  if (requireNamespace("harmony", quietly = TRUE)) {
    data <- integrate_data(data, method = "RunHarmony", group_by = "SampleID",
                           n_pcs = 10, n_features = 100, resolution = 0.8)
    stopifnot_msg("umap_harmony" %in% names(data@reductions),
                  "integrate_data (RunHarmony) produced a UMAP")
  } else {
    skip("harmony not installed — integration not exercised (install.packages('harmony')).")
  }

  step("7. Subset + graft-back (barcode-keyed)")
  subA <- subset_samples(data, ids = "A")
  subA$test_col <- 99
  data <- graft_meta(data, subA, cols = "test_col")
  stopifnot_msg(all(data$test_col[data$SampleID == "A"] == 99, na.rm = TRUE) &&
                  all(is.na(data$test_col[data$SampleID == "B"])),
                "graft_meta wrote only the child's cells, by barcode")

  step("8. Provenance recorded")
  provenance(data)
}

cat("\n=== smoke_test.R finished. Any [FAIL] above would have stopped the run. ===\n")
