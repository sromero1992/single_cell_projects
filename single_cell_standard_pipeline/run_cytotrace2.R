## =============================================================================
## run_cytotrace2.R
## Standalone analysis: CytoTRACE2 + AUCell cell-cycle + CellRank trajectories
## =============================================================================
##
## USAGE:
##   source("run_cytotrace2.R")
##   run_cytotrace2("my_data.rds", out_dir = "ct2_output")
##
## REQUIRED R PACKAGES:
##   CytoTRACE2   remotes::install_github("digitalcytometry/cytotrace2", subdir="cytotrace2_r")
##   Seurat       install.packages("Seurat")
##   AUCell       BiocManager::install("AUCell")
##   ggplot2      install.packages(c("ggplot2","patchwork","viridis","scales"))
##   org.Mm.eg.db BiocManager::install(c("org.Mm.eg.db","AnnotationDbi"))  # if mouse
##   reticulate   install.packages("reticulate")   # only for CellRank section
##
## PYTHON ENV (for CellRank section only):
##   conda create -n scanpy_env_311 python=3.11
##   pip install scanpy cellrank
## =============================================================================

## ── CONFIG ────────────────────────────────────────────────────────────────────
PYTHON_ENV  <- "scanpy_env_311"   # conda env with scanpy + cellrank installed
SPECIES     <- "mouse"            # "mouse" or "human"
N_NEIGHBORS <- 30L                # kNN for CellRank graph
SEED        <- 42L

## =============================================================================
## MAIN FUNCTION
## =============================================================================

run_cytotrace2 <- function(
    rds_path,
    out_dir      = "cytotrace2_output",
    species      = SPECIES,
    python_env   = PYTHON_ENV,
    run_cellrank = TRUE,
    seed         = SEED
) {

  ## ── 0. Setup ----------------------------------------------------------------
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  set.seed(seed)

  .check_pkg <- function(pkg, bioc = FALSE, github = NULL) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (!is.null(github)) {
        stop(sprintf("Package '%s' not found. Install with:\n  remotes::install_github('%s')", pkg, github))
      } else if (bioc) {
        stop(sprintf("Package '%s' not found. Install with:\n  BiocManager::install('%s')", pkg, pkg))
      } else {
        stop(sprintf("Package '%s' not found. Install with:\n  install.packages('%s')", pkg, pkg))
      }
    }
  }

  .check_pkg("Seurat")
  .check_pkg("CytoTRACE2", github = "digitalcytometry/cytotrace2/cytotrace2_r")
  .check_pkg("AUCell",     bioc = TRUE)
  .check_pkg("ggplot2")
  .check_pkg("patchwork")
  .check_pkg("viridis")

  ## ── 1. Load Seurat object ---------------------------------------------------
  cat("[1/6] Loading Seurat object from:", rds_path, "\n")
  so <- readRDS(rds_path)
  cat(sprintf("  %d cells | %d genes\n", ncol(so), nrow(so)))

  ## Check for UMAP (needed for plots)
  has_umap <- "umap" %in% names(so@reductions)
  if (!has_umap) {
    cat("  No UMAP found — running PCA + UMAP (default settings)...\n")
    so <- Seurat::NormalizeData(so, verbose = FALSE)
    so <- Seurat::FindVariableFeatures(so, verbose = FALSE)
    so <- Seurat::ScaleData(so, verbose = FALSE)
    so <- Seurat::RunPCA(so, verbose = FALSE)
    so <- Seurat::RunUMAP(so, dims = 1:30, verbose = FALSE)
  }

  ## ── 2. CytoTRACE2 -----------------------------------------------------------
  cat("[2/6] Running CytoTRACE2...\n")

  ## Extract normalised log expression matrix (genes × cells) with gene symbols
  expr_mat <- as.matrix(Seurat::GetAssayData(so, layer = "data"))

  ## If row names look like Ensembl IDs, try to map to gene symbols
  if (all(grepl("^ENS", rownames(expr_mat)))) {
    cat("  Row names look like Ensembl IDs — mapping to gene symbols...\n")
    if (species == "mouse" &&
        requireNamespace("org.Mm.eg.db", quietly = TRUE) &&
        requireNamespace("AnnotationDbi", quietly = TRUE)) {
      sym <- suppressMessages(AnnotationDbi::mapIds(
        org.Mm.eg.db::org.Mm.eg.db,
        keys    = rownames(expr_mat),
        column  = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"))
      keep <- !is.na(sym) & sym != "" & !duplicated(sym)
      expr_mat <- expr_mat[keep, ]
      rownames(expr_mat) <- sym[keep]
      cat(sprintf("  Mapped %d Ensembl IDs → gene symbols.\n", sum(keep)))
    } else if (species == "human" &&
               requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
               requireNamespace("AnnotationDbi", quietly = TRUE)) {
      sym <- suppressMessages(AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys    = rownames(expr_mat),
        column  = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"))
      keep <- !is.na(sym) & sym != "" & !duplicated(sym)
      expr_mat <- expr_mat[keep, ]
      rownames(expr_mat) <- sym[keep]
    }
  }

  ## Deduplicate gene symbols (keep highest-mean row)
  if (any(duplicated(rownames(expr_mat)))) {
    expr_mat <- expr_mat[order(rowMeans(expr_mat), decreasing = TRUE), ]
    expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]
    cat(sprintf("  After dedup: %d unique gene symbols.\n", nrow(expr_mat)))
  }

  ## Convert to dense data.frame (CytoTRACE2 requirement)
  ct2_input <- as.data.frame(expr_mat)
  cat(sprintf("  Input: %d genes × %d cells | max value: %.2f\n",
              nrow(ct2_input), ncol(ct2_input), max(ct2_input, na.rm = TRUE)))

  ## ── THE EXACT WORKING CALL ─────────────────────────────────────────────────
  ##   Do NOT add verbose=TRUE — it is not a valid argument.
  ##   parallelize_smoothing=FALSE avoids OS-level issues (safe on Linux too).
  ct2_result <- CytoTRACE2::cytotrace2(
    ct2_input,
    species               = species,
    is_seurat             = FALSE,
    parallelize_models    = FALSE,
    parallelize_smoothing = FALSE,
    seed                  = seed
  )

  ## Extract and store CytoTRACE2_Score (continuous, 0=differentiated, 1=stem)
  ct2_scores                <- ct2_result[["CytoTRACE2_Score"]]
  names(ct2_scores)         <- rownames(ct2_result)
  ct2_scores                <- ct2_scores[colnames(so)]          # align to Seurat
  ct2_scores                <- pmin(pmax(ct2_scores, 0), 1)
  so$CytoTRACE2_Score       <- ct2_scores
  so$CytoTRACE2_Potency     <- ct2_result[["CytoTRACE2_Potency"]][colnames(so)]
  cat(sprintf("  Done. Score range: [%.3f, %.3f]\n",
              min(ct2_scores, na.rm = TRUE), max(ct2_scores, na.rm = TRUE)))

  ## ── 3. Cell-cycle scoring (AUCell — G2M and S kept SEPARATE) ---------------
  ##  Decision rationale: G2M (mitotic spindle) and S (DNA replication) capture
  ##  distinct biological phases and should not be collapsed into one score.
  ##  Two separate AUCell scores let us ask: does high potency co-occur more
  ##  with S-phase (progenitors replicating) or G2M (active division)?
  cat("[3/6] Computing AUCell cell-cycle scores (G2M and S phase separately)...\n")

  ## Built-in Seurat cycling gene sets (Tirosh et al. 2016, human symbols).
  ## Convert to mouse Title Case if needed (e.g. CDK1 → Cdk1).
  data("cc.genes.updated.2019", package = "Seurat", envir = environment())
  s_genes   <- cc.genes.updated.2019$s.genes
  g2m_genes <- cc.genes.updated.2019$g2m.genes
  if (species == "mouse") {
    s_genes   <- paste0(toupper(substr(s_genes,   1, 1)),
                        tolower(substring(s_genes,   2)))
    g2m_genes <- paste0(toupper(substr(g2m_genes, 1, 1)),
                        tolower(substring(g2m_genes, 2)))
  }

  ## Filter to genes present in data
  s_genes   <- intersect(s_genes,   rownames(expr_mat))
  g2m_genes <- intersect(g2m_genes, rownames(expr_mat))
  cat(sprintf("  S-phase genes found: %d / %d\n",   length(s_genes),   length(cc.genes.updated.2019$s.genes)))
  cat(sprintf("  G2M-phase genes found: %d / %d\n", length(g2m_genes), length(cc.genes.updated.2019$g2m.genes)))

  gene_sets <- list(S_phase = s_genes, G2M_phase = g2m_genes)

  ## Build gene rankings once, compute AUC for both sets
  cat("  Building AUCell rankings (this takes ~30s for large datasets)...\n")
  expr_mat_int <- round(expr_mat)   # AUCell expects integer-like counts
  rankings     <- AUCell::AUCell_buildRankings(expr_mat_int,
                                                plotStats = FALSE,
                                                verbose   = FALSE)
  auc_res <- AUCell::AUCell_calcAUC(gene_sets, rankings, verbose = FALSE)
  auc_mat <- AUCell::getAUC(auc_res)   # 2 × n_cells matrix

  so$AUCell_S_Score   <- auc_mat["S_phase",   colnames(so)]
  so$AUCell_G2M_Score <- auc_mat["G2M_phase", colnames(so)]
  cat(sprintf("  S-phase AUC:   [%.3f, %.3f]\n",
              min(so$AUCell_S_Score),   max(so$AUCell_S_Score)))
  cat(sprintf("  G2M-phase AUC: [%.3f, %.3f]\n",
              min(so$AUCell_G2M_Score), max(so$AUCell_G2M_Score)))

  ## ── 4. Figures --------------------------------------------------------------
  cat("[4/6] Generating figures...\n")

  ## Build a unified data.frame for plotting
  umap_coords <- as.data.frame(Seurat::Embeddings(so, "umap"))
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  df_plot <- data.frame(
    umap_coords,
    CytoTRACE2  = so$CytoTRACE2_Score,
    S_AUC       = so$AUCell_S_Score,
    G2M_AUC     = so$AUCell_G2M_Score,
    CellType    = if (!is.null(so$celltype)) so$celltype else
                  if (!is.null(so$seurat_clusters)) paste0("C", so$seurat_clusters) else
                  "unknown",
    row.names   = colnames(so)
  )

  ## ── 4a. UMAP panels ──
  .umap_col <- function(var, label, palette = "viridis") {
    ggplot2::ggplot(df_plot, ggplot2::aes(x = UMAP_1, y = UMAP_2, colour = .data[[var]])) +
      ggplot2::geom_point(size = 0.4, alpha = 0.7) +
      viridis::scale_colour_viridis_c(option = palette, name = label) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "right",
                     axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank()) +
      ggplot2::labs(x = "UMAP 1", y = "UMAP 2", title = label)
  }

  p_ct2     <- .umap_col("CytoTRACE2", "CytoTRACE2 Score", "plasma")
  p_s_auc   <- .umap_col("S_AUC",      "S-phase AUC",      "viridis")
  p_g2m_auc <- .umap_col("G2M_AUC",    "G2M-phase AUC",    "viridis")

  fig_umap <- p_ct2 + p_s_auc + p_g2m_auc +
    patchwork::plot_layout(ncol = 3) +
    patchwork::plot_annotation(title = "CytoTRACE2 and cell-cycle AUCell scores")

  ggplot2::ggsave(file.path(out_dir, "fig_umap_ct2_cycle.pdf"),
                  fig_umap, width = 16, height = 5)
  ggplot2::ggsave(file.path(out_dir, "fig_umap_ct2_cycle.png"),
                  fig_umap, width = 16, height = 5, dpi = 200)

  ## ── 4b. 2-D KDE contour: CytoTRACE2 vs AUCell (G2M and S separate) ──
  ##  We keep G2M and S as two panels because they mark distinct phases and
  ##  can have different relationships to potency (e.g. S may correlate with
  ##  progenitor expansion while G2M marks active mitosis).
  .contour_panel <- function(y_var, y_lab, colour_var = "CellType") {
    ggplot2::ggplot(df_plot,
           ggplot2::aes(x = CytoTRACE2, y = .data[[y_var]])) +
      ## filled density contours (background)
      ggplot2::stat_density_2d(
        ggplot2::aes(fill = ggplot2::after_stat(level)),
        geom   = "polygon",
        alpha  = 0.35,
        colour = NA,
        bins   = 12
      ) +
      viridis::scale_fill_viridis_c(option = "mako", name = "Density") +
      ## individual points coloured by cell type
      ggplot2::geom_point(
        ggplot2::aes(colour = .data[[colour_var]]),
        size  = 0.6,
        alpha = 0.6
      ) +
      ## regression line
      ggplot2::geom_smooth(method = "loess", se = TRUE,
                           colour = "firebrick", linewidth = 0.8, alpha = 0.2) +
      ggplot2::labs(
        x     = "CytoTRACE2 Score (0 = differentiated → 1 = stem)",
        y     = y_lab,
        title = paste("CytoTRACE2 vs", y_lab),
        colour = "Cell type"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "right")
  }

  p_contour_s   <- .contour_panel("S_AUC",   "S-phase AUC (AUCell)")
  p_contour_g2m <- .contour_panel("G2M_AUC", "G2M-phase AUC (AUCell)")

  fig_contour <- p_contour_s + p_contour_g2m +
    patchwork::plot_layout(ncol = 2, guides = "collect") +
    patchwork::plot_annotation(
      title    = "Lineage potency × cell-cycle activity (2-D KDE contours)",
      subtitle = "Contour fill = cell density  |  Points = individual cells (coloured by lineage)  |  Line = LOESS"
    )

  ggplot2::ggsave(file.path(out_dir, "fig_contour_ct2_vs_cycle.pdf"),
                  fig_contour, width = 14, height = 6)
  ggplot2::ggsave(file.path(out_dir, "fig_contour_ct2_vs_cycle.png"),
                  fig_contour, width = 14, height = 6, dpi = 200)
  cat("  Saved: fig_umap_ct2_cycle + fig_contour_ct2_vs_cycle\n")

  ## ── 5. CellRank trajectories (reticulate) -----------------------------------
  if (!run_cellrank) {
    cat("[5/6] CellRank skipped (run_cellrank = FALSE).\n")
  } else {
    cat("[5/6] CellRank trajectory (CytoTRACEKernel + PseudotimeKernel via reticulate)...\n")

    ## Check / initialise Python
    if (!requireNamespace("reticulate", quietly = TRUE))
      stop("reticulate not installed. Run: install.packages('reticulate')")

    if (!reticulate::py_available(initialize = FALSE)) {
      tryCatch(
        reticulate::use_condaenv(python_env, required = FALSE),
        error = function(e)
          cat(sprintf("  Warning: could not init env '%s': %s\n", python_env, e$message))
      )
    }
    if (!reticulate::py_available(initialize = TRUE))
      stop("Python not available. Set PYTHON_ENV at top of script to your conda env name.")

    for (.mod in c("scanpy", "cellrank", "scipy")) {
      if (!reticulate::py_module_available(.mod))
        stop(sprintf("Python module '%s' not found. In your env run: pip install %s", .mod, .mod))
    }

    ## Export metadata + expression for Python
    cr_dir  <- file.path(out_dir, "cellrank_output")
    dir.create(cr_dir, showWarnings = FALSE)

    meta_tsv <- file.path(cr_dir, "metadata.tsv")
    mm_mat   <- file.path(cr_dir, "expr_matrix.mtx.gz")  # sparse expression
    mm_genes <- file.path(cr_dir, "genes.txt")
    mm_cells <- file.path(cr_dir, "cells.txt")

    ## Metadata
    umap_df <- as.data.frame(Seurat::Embeddings(so, "umap"))
    pca_df  <- as.data.frame(Seurat::Embeddings(so, "pca")[, 1:min(50, ncol(Seurat::Embeddings(so, "pca")))])
    meta_df <- cbind(
      barcode       = colnames(so),
      umap_df,
      pca_df,
      CytoTRACE2    = so$CytoTRACE2_Score,
      S_AUC         = so$AUCell_S_Score,
      G2M_AUC       = so$AUCell_G2M_Score,
      celltype      = if (!is.null(so$celltype)) so$celltype else NA_character_
    )
    write.table(meta_df, meta_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

    ## Sparse expression matrix (genes × cells) as MatrixMarket
    expr_sparse <- Seurat::GetAssayData(so, layer = "counts")
    Matrix::writeMM(expr_sparse, gzfile(mm_mat))
    writeLines(rownames(expr_sparse), mm_genes)
    writeLines(colnames(expr_sparse), mm_cells)

    cat("  Metadata + expression exported to:", cr_dir, "\n")

    ## Pass paths to Python
    .py_main <- reticulate::import_main()
    .py_main$cr_meta_tsv  <- meta_tsv
    .py_main$cr_expr_mtx  <- mm_mat
    .py_main$cr_genes_txt <- mm_genes
    .py_main$cr_cells_txt <- mm_cells
    .py_main$cr_out_dir   <- cr_dir
    .py_main$cr_n_neighbors <- N_NEIGHBORS

    ## Run CytoTRACEKernel AND PseudotimeKernel (with CytoTRACE2 scores)
    reticulate::py_run_string("
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
from cellrank.kernels import CytoTRACEKernel, PseudotimeKernel, ConnectivityKernel

## ── Build AnnData from exported files ────────────────────────────────────────
meta  = pd.read_csv(cr_meta_tsv, sep='\\t', index_col='barcode')

## Load sparse expression matrix
X = scipy.io.mmread(cr_expr_mtx).T.tocsr()   # cells × genes
genes = open(cr_genes_txt).read().splitlines()
cells = open(cr_cells_txt).read().splitlines()

adata = sc.AnnData(X=X)
adata.obs_names  = cells
adata.var_names  = genes
adata.obs        = meta.reindex(cells)

## Store UMAP and PCA
umap_cols = [c for c in meta.columns if c.startswith('UMAP_')]
pca_cols  = [c for c in meta.columns if c.startswith('PC')]
if umap_cols:
    adata.obsm['X_umap'] = meta.reindex(cells)[umap_cols].values.astype(float)
if pca_cols:
    adata.obsm['X_pca']  = meta.reindex(cells)[pca_cols].values.astype(float)

## Normalise for CytoTRACEKernel (expects library-size-normalised log data)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['logcounts'] = adata.X.copy()

## Build kNN graph (required by all kernels)
sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=cr_n_neighbors,
                n_pcs=min(30, len(pca_cols)))
print(f'AnnData built: {adata.n_obs} cells x {adata.n_vars} genes')

## Connectivity kernel (shared background for all methods)
ck = ConnectivityKernel(adata).compute_transition_matrix()
U  = adata.obsm['X_umap']

## ── A. CytoTRACEKernel (uses its OWN CytoTRACE score from expression) ───────
print('Running CytoTRACEKernel...')
try:
    ctk = CytoTRACEKernel(adata)
    ctk.compute_cytotrace()          ## REQUIRED: sets internal pseudotime state
    ctk.compute_transition_matrix(threshold_scheme='soft', nu=0.5)
    T_ctk = (0.8 * ctk + 0.2 * ck).transition_matrix.toarray()
    V_ctk = T_ctk @ U - U
    pd.DataFrame(V_ctk, columns=['dUMAP_1','dUMAP_2'],
                 index=adata.obs_names).to_csv(cr_out_dir + '/velocity_CytoTRACEKernel.csv')
    ## Also save the internal CytoTRACE score so R can colour the UMAP
    pd.DataFrame({'CytoTRACE_internal': adata.obs['ct_pseudotime']},
                 index=adata.obs_names).to_csv(cr_out_dir + '/cytotrace_internal_scores.csv')
    print('  CytoTRACEKernel: Done')
except Exception as e:
    print(f'  CytoTRACEKernel FAILED: {e}')

## ── B. PseudotimeKernel with CytoTRACE2 scores (R-computed) ─────────────────
print('Running PseudotimeKernel (CytoTRACE2 scores)...')
try:
    ct2_col = 'CytoTRACE2'
    score   = adata.obs[ct2_col].astype(float)
    score   = score.fillna(score.median())
    smin, smax = score.min(), score.max()
    if smax > smin:
        score = (score - smin) / (smax - smin)
    adata.obs['_ct2_pt'] = score
    ptk = PseudotimeKernel(adata, time_key='_ct2_pt').compute_transition_matrix(
        threshold_scheme='soft', nu=0.5)
    T_ptk = (0.8 * ptk + 0.2 * ck).transition_matrix.toarray()
    V_ptk = T_ptk @ U - U
    pd.DataFrame(V_ptk, columns=['dUMAP_1','dUMAP_2'],
                 index=adata.obs_names).to_csv(cr_out_dir + '/velocity_CytoTRACE2_PT.csv')
    print('  PseudotimeKernel (CytoTRACE2): Done')
except Exception as e:
    print(f'  PseudotimeKernel FAILED: {e}')

print('CellRank sweep complete.')
")

    ## ── 5b. Read velocity CSVs and plot in R ──────────────────────────────────
    .plot_velocity <- function(vel_csv, method_label, df_base) {
      if (!file.exists(vel_csv)) return(NULL)
      vel <- read.csv(vel_csv, row.names = 1)
      df  <- df_base
      df$dx <- vel[rownames(df), "dUMAP_1"]
      df$dy <- vel[rownames(df), "dUMAP_2"]
      ## bin into grid for arrow display
      df$gx <- cut(df$UMAP_1, 20, labels = FALSE)
      df$gy <- cut(df$UMAP_2, 20, labels = FALSE)
      grid <- aggregate(cbind(UMAP_1, UMAP_2, dx, dy) ~ gx + gy,
                        data = df, FUN = mean, na.rm = TRUE)
      grid$len <- sqrt(grid$dx^2 + grid$dy^2)
      grid     <- grid[grid$len > quantile(grid$len, 0.25, na.rm = TRUE), ]

      ggplot2::ggplot(df, ggplot2::aes(x = UMAP_1, y = UMAP_2,
                                        colour = CytoTRACE2)) +
        ggplot2::geom_point(size = 0.4, alpha = 0.6) +
        viridis::scale_colour_viridis_c(option = "plasma", name = "CT2 Score") +
        ggplot2::geom_segment(
          data = grid,
          ggplot2::aes(x = UMAP_1, y = UMAP_2,
                       xend = UMAP_1 + dx * 0.5,
                       yend = UMAP_2 + dy * 0.5),
          inherit.aes = FALSE,
          arrow = ggplot2::arrow(length = ggplot2::unit(0.08, "inches"),
                                 type = "closed"),
          colour = "black", linewidth = 0.4, alpha = 0.7
        ) +
        ggplot2::labs(title = method_label,
                      x = "UMAP 1", y = "UMAP 2") +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(axis.text = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank())
    }

    df_base <- data.frame(
      UMAP_1     = umap_df[, 1],
      UMAP_2     = umap_df[, 2],
      CytoTRACE2 = so$CytoTRACE2_Score,
      row.names  = colnames(so)
    )

    p_ctk  <- .plot_velocity(
      file.path(cr_dir, "velocity_CytoTRACEKernel.csv"),
      "CytoTRACEKernel\n(CellRank built-in CytoTRACE)", df_base)
    p_ct2k <- .plot_velocity(
      file.path(cr_dir, "velocity_CytoTRACE2_PT.csv"),
      "PseudotimeKernel\n(CytoTRACE2 scores)", df_base)

    panels <- Filter(Negate(is.null), list(p_ctk, p_ct2k))
    if (length(panels) > 0) {
      fig_traj <- patchwork::wrap_plots(panels, ncol = length(panels)) +
        patchwork::plot_annotation(
          title    = "CellRank trajectory arrows",
          subtitle = "Arrow direction = inferred differentiation flow  |  Colour = CytoTRACE2 Score"
        )
      ggplot2::ggsave(file.path(out_dir, "fig_cellrank_trajectories.pdf"),
                      fig_traj, width = 7 * length(panels), height = 6)
      ggplot2::ggsave(file.path(out_dir, "fig_cellrank_trajectories.png"),
                      fig_traj, width = 7 * length(panels), height = 6, dpi = 200)
      cat("  Saved: fig_cellrank_trajectories\n")
    } else {
      cat("  No velocity CSVs produced — check Python output above.\n")
    }
  }

  ## ── 6. Save updated Seurat object and summary table -----------------------
  cat("[6/6] Saving outputs...\n")
  saveRDS(so, file.path(out_dir, "seurat_with_ct2_cycle.rds"))

  summary_df <- data.frame(
    cell    = colnames(so),
    CT2_Score   = so$CytoTRACE2_Score,
    CT2_Potency = so$CytoTRACE2_Potency,
    S_AUC       = so$AUCell_S_Score,
    G2M_AUC     = so$AUCell_G2M_Score
  )
  write.csv(summary_df, file.path(out_dir, "scores_summary.csv"), row.names = FALSE)

  cat("\n── Analysis complete ──────────────────────────────────────────────────\n")
  cat("Outputs saved to:", out_dir, "\n")
  cat("  seurat_with_ct2_cycle.rds  — Seurat object with all scores\n")
  cat("  scores_summary.csv         — per-cell score table\n")
  cat("  fig_umap_ct2_cycle.*       — UMAP panels (CT2, S-AUC, G2M-AUC)\n")
  cat("  fig_contour_ct2_vs_cycle.* — 2-D KDE contour CT2 vs S-AUC / G2M-AUC\n")
  if (run_cellrank)
    cat("  fig_cellrank_trajectories.* — CytoTRACEKernel + PT-kernel arrows\n")
  cat("────────────────────────────────────────────────────────────────────────\n")

  invisible(list(seurat = so, scores = summary_df))
}

## =============================================================================
## EXAMPLE CALL
## =============================================================================
## source("run_cytotrace2.R")
##
## result <- run_cytotrace2(
##   rds_path     = "my_data.rds",
##   out_dir      = "ct2_output",
##   species      = "mouse",         # or "human"
##   python_env   = "scanpy_env_311",
##   run_cellrank = TRUE
## )
##
## Access scores directly:
##   result$seurat$CytoTRACE2_Score
##   result$seurat$AUCell_S_Score
##   result$seurat$AUCell_G2M_Score
## =============================================================================
