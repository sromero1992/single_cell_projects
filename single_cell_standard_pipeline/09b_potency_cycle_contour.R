# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 09b: POTENCY x CYCLING STATE LANDSCAPE
# =============================================================================
# Runs AFTER Script 09 (cell scores). Places every cell on a 2-D state landscape:
#
#       x-axis  =  POTENCY          (CytoTRACE2 score; 0 = differentiated,
#                                     1 = stem/totipotent)
#       y-axis  =  CYCLING ACTIVITY (1 - H_cycle, pooled Shannon entropy over
#                                     S + G2M genes; HIGH activity = cycling)
#
# and draws KDE-contoured plots classifying four quadrants, named by the two AXES
# (not a specific cell identity) so they stay valid for any cycling axis:
#       high potency + LOW  activity  ->  Stem-like / non-cycling    (reserve stem)
#       high potency + HIGH activity  ->  Stem-like / cycling        (crypt-base stem)
#       low  potency + HIGH activity  ->  Differentiated / cycling   (transit-amplifying)
#       low  potency + LOW  activity  ->  Differentiated / non-cycling (mature)
# The cycling axis is drawn three ways -- combined S+G2M, S-phase only, G2M only --
# each with its OWN median split (phase-specific quadrants; files *_Sphase / *_G2M).
#
# The potency_cycle_contour() function below is SELF-CONTAINED and portable:
# source() it into any project and call it on an rds path or a Seurat object.
# When this file is RUN as a script (Rscript 09b_...R) the driver at the bottom
# wires it to the pipeline via config.R (nr4a1 defaults; override via env vars).
# =============================================================================


## =============================================================================
## 1.  CYCLING GENE SETS  (Tirosh et al. 2016, human symbols)
## -----------------------------------------------------------------------------
## CITATION (cite this for the cell-cycle / S & G2M gene lists):
##   Tirosh I, Izar B, Prakadan SM, et al. "Dissecting the multicellular
##   ecosystem of metastatic melanoma by single-cell RNA-seq."
##   Science. 2016 Apr 8;352(6282):189-196. doi:10.1126/science.aad0501
## These are the canonical 43 S-phase + 54 G2/M marker genes from that paper
## (the same lists distributed as Seurat::cc.genes). Used here for AUCell S/G2M
## scoring and the pooled cycling-entropy metric.
## =============================================================================
.PCC_S_GENES <- c(
  "MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6",
  "CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP",
  "RAD51AP1","GMNN","ZWINT","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM",
  "CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8"
)
.PCC_G2M_GENES <- c(
  "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2",
  "NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2",
  "CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B",
  "HJURP","CDCA3","HN1","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2",
  "DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN",
  "LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA"
)

## Convert a vector of HUMAN gene symbols to mouse Title-case (CDK1 -> Cdk1)
.pcc_to_mouse <- function(g) paste0(toupper(substr(g, 1, 1)), tolower(substring(g, 2)))


## =============================================================================
## 2.  HELPERS
## =============================================================================
.pcc_check_pkg <- function(pkg, how = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Package '%s' not found.%s", pkg,
                 if (is.null(how)) "" else paste0(" Install with:\n  ", how)))
}

## --- Match a gene set against a matrix's rownames, with fallbacks -------------
## Exact match first; then case-insensitive; then (mouse) Ensembl-ID lookup.
.pcc_match_genes <- function(genes, all_genes, species = "mouse") {
  hit <- intersect(genes, all_genes)
  if (length(hit) >= 3) return(hit)

  up_map <- setNames(all_genes, toupper(all_genes))
  hit2   <- up_map[intersect(toupper(genes), toupper(all_genes))]
  hit2   <- unname(hit2[!is.na(hit2)])
  if (length(hit2) >= 3) return(hit2)

  if (species == "mouse" && any(grepl("^ENSMUSG", all_genes)) &&
      requireNamespace("org.Mm.eg.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    ens <- tryCatch(suppressMessages(AnnotationDbi::mapIds(
      org.Mm.eg.db::org.Mm.eg.db, keys = .pcc_to_mouse(genes),
      column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")),
      error = function(e) NULL)
    if (!is.null(ens)) {
      hit3 <- intersect(ens[!is.na(ens)], all_genes)
      if (length(hit3) >= 3) return(hit3)
    }
  }
  hit2
}


## =============================================================================
## 3.  POOLED SHANNON CELL-CYCLE ENTROPY  (ported from scDVEP.R)
## =============================================================================
## H_cycle(c) = normalised Shannon entropy of expression over the pooled S + G2M
## panel, per cell. HIGH H_cycle = spread uniformly -> quiescent; LOW = focused on
## S or G2M -> cycling. cycle_activity = 1 - H_cycle (HIGH = cycling). Cells with
## near-zero total cycling expression are gated out ("Unscored").
compute_pooled_cycle_entropy <- function(expr_mat, s_genes, g2m_genes,
                                         pseudocount = 1e-6, verbose = TRUE) {
  cycle_use <- union(s_genes, g2m_genes)
  cycle_use <- intersect(cycle_use, rownames(expr_mat))
  if (length(cycle_use) < 5)
    stop(sprintf("Only %d cycling genes found in the matrix; need >= 5.",
                 length(cycle_use)))
  if (verbose)
    message(sprintf("[H_cycle] %d cycling genes pooled (S+G2M).", length(cycle_use)))

  mat_cyc <- as.matrix(expr_mat[cycle_use, , drop = FALSE])

  total_cyc <- colSums(mat_cyc)
  gate      <- quantile(total_cyc, 0.05, na.rm = TRUE)
  low_expr  <- total_cyc < gate

  Z     <- colSums(mat_cyc + pseudocount)
  P     <- sweep(mat_cyc + pseudocount, 2, Z, "/")
  H_raw <- -colSums(P * log(P + pseudocount))
  H     <- H_raw / log(length(cycle_use))       # normalise to [0,1]
  H     <- pmin(pmax(H, 0), 1)
  H[low_expr] <- NA_real_
  names(H) <- colnames(expr_mat)

  med   <- median(H, na.rm = TRUE)
  state <- ifelse(is.na(H), "Unscored", ifelse(H >= med, "Quiescent", "Cycling"))
  names(state) <- colnames(expr_mat)

  list(H_cycle        = H,          # entropy: HIGH = quiescent
       cycle_activity = 1 - H,      # HIGH = cycling (plotting axis)
       cycle_state    = state,
       genes_used     = cycle_use)
}


## =============================================================================
## 4.  MAIN FUNCTION
## =============================================================================
potency_cycle_contour <- function(
    object,                              # rds path OR a Seurat object
    out_dir       = "state_landscape",
    species       = "mouse",             # "mouse" or "human"
    assay         = "RNA",
    ct2_score_col = "CytoTRACE2_Score",  # reused if this column already exists
    cycle_axis    = c("auc", "entropy"), # y-axis metric: AUCell S+G2M (default,
                                         # graded/clean) or 1 - H_cycle (banded)
    subset_col    = NULL,                # e.g. "CellType_broad": restrict to a lineage
    subset_val    = NULL,                # e.g. "Colonocytes": the value(s) to keep
    celltype_col  = NULL,                # column for islands/facets (e.g. subtype
                                         # "CellType"); NULL = auto-detect
    condition_col = NULL,                # grouping for the stat figures; NULL =
                                         # auto-detect "Genotype_sex"
    condition_levels = NULL,             # x-axis order of the conditions (NULL = sorted)
    comparisons   = NULL,                # list of c(g1,g2) pairs for Wilcoxon; NULL =
                                         # the 4 Nr4a1 contrasts if those levels exist
    seed          = 42L
) {
  cycle_axis <- match.arg(cycle_axis)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  set.seed(seed)

  .pcc_check_pkg("Seurat")
  .pcc_check_pkg("ggplot2")
  .pcc_check_pkg("AUCell",  "BiocManager::install('AUCell')")
  # (viridis colour scale comes from ggplot2 itself; no separate 'viridis' pkg needed)

  ## -- 4.0  Load object --------------------------------------------------------
  so <- if (is.character(object)) {
    cat("[1/6] Loading Seurat object:", object, "\n"); readRDS(object)
  } else {
    cat("[1/6] Using in-memory Seurat object.\n"); object
  }
  cat(sprintf("      %d cells x %d genes\n", ncol(so), nrow(so)))

  # Optional: restrict to a broad lineage (e.g. CellType_broad == "Colonocytes")
  # so the landscape + median thresholds are LOCAL to that lineage, and the
  # islands/facets split by its subtypes (set celltype_col = "CellType").
  if (!is.null(subset_col) && !is.null(subset_val)) {
    if (!subset_col %in% colnames(so@meta.data))
      stop("subset_col '", subset_col, "' not in metadata.", call. = FALSE)
    keep_cells <- colnames(so)[as.character(so@meta.data[[subset_col]]) %in% subset_val]
    if (length(keep_cells) < 20)
      stop(sprintf("Only %d cells match %s in {%s}.", length(keep_cells),
                   subset_col, paste(subset_val, collapse = ", ")), call. = FALSE)
    cat(sprintf("      Subsetting to %s in {%s}: %d cells\n",
                subset_col, paste(subset_val, collapse = ", "), length(keep_cells)))
    # MEMORY: strip the (large) parent to just the RNA assay FIRST -- drop
    # reductions / graphs / scale.data / other assays -- so the peak during the
    # subset is small, then subset and free the parent. 09b only needs RNA
    # counts/data + meta.data (CytoTRACE2_Score etc.).
    so <- tryCatch(
      Seurat::DietSeurat(so, assays = assay, dimreducs = NULL, graphs = NULL),
      error = function(e) so)
    so <- subset(so, cells = keep_cells)
    gc()
  }

  # Keep SPARSE -- densifying the whole genes x cells matrix was the ~25 GiB
  # allocation. compute_pooled_cycle_entropy() densifies only the ~100 cycling-gene
  # rows, and AUCell handles a sparse matrix, so nothing needs the full dense copy.
  expr_mat <- tryCatch(
    SeuratObject::LayerData(so, assay = assay, layer = "data"),
    error = function(e) Seurat::GetAssayData(so, assay = assay, layer = "data"))

  if (all(grepl("^ENS", rownames(expr_mat)))) {
    cat("      Row names look like Ensembl IDs -- mapping to symbols...\n")
    db  <- if (species == "mouse") "org.Mm.eg.db" else "org.Hs.eg.db"
    if (requireNamespace(db, quietly = TRUE) &&
        requireNamespace("AnnotationDbi", quietly = TRUE)) {
      sym <- suppressMessages(AnnotationDbi::mapIds(
        getExportedValue(db, db), keys = rownames(expr_mat),
        column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
      keep <- !is.na(sym) & sym != "" & !duplicated(sym)
      expr_mat <- expr_mat[keep, ]; rownames(expr_mat) <- sym[keep]
    }
  }
  if (any(duplicated(rownames(expr_mat)))) {
    expr_mat <- expr_mat[order(rowMeans(expr_mat), decreasing = TRUE), ]
    expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]
  }

  s_src   <- if (species == "mouse") .pcc_to_mouse(.PCC_S_GENES)   else .PCC_S_GENES
  g2m_src <- if (species == "mouse") .pcc_to_mouse(.PCC_G2M_GENES) else .PCC_G2M_GENES
  s_use   <- .pcc_match_genes(s_src,   rownames(expr_mat), species)
  g2m_use <- .pcc_match_genes(g2m_src, rownames(expr_mat), species)
  cat(sprintf("      Cycling genes matched: %d S + %d G2M\n",
              length(s_use), length(g2m_use)))

  ## -- 4.1  POTENCY: CytoTRACE2 (reuse if present, else run) -------------------
  ct2_candidates <- unique(c(ct2_score_col, "CytoTRACE2_Score", "CytoTRACE2"))
  ct2_found <- ct2_candidates[
    ct2_candidates %in% colnames(so@meta.data) &
    vapply(ct2_candidates, function(cc)
      cc %in% colnames(so@meta.data) && any(!is.na(so@meta.data[[cc]])),
      logical(1))]
  if (length(ct2_found) > 0) {
    use_col <- ct2_found[1]
    cat("[2/6] Reusing existing CytoTRACE2 scores from column '", use_col, "'.\n", sep = "")
    ct2 <- so@meta.data[[use_col]]; names(ct2) <- colnames(so)
  } else {
    cat("[2/6] No existing potency column -- running CytoTRACE2...\n")
    .pcc_check_pkg("CytoTRACE2",
      "remotes::install_github('digitalcytometry/cytotrace2', subdir='cytotrace2_r')")
    ct2_input  <- as.data.frame(expr_mat)
    ct2_result <- CytoTRACE2::cytotrace2(
      ct2_input, species = species, is_seurat = FALSE,
      parallelize_models = FALSE, parallelize_smoothing = FALSE, seed = seed)
    ct2 <- ct2_result[["CytoTRACE2_Score"]]
    names(ct2) <- rownames(ct2_result)
    ct2 <- pmin(pmax(ct2[colnames(so)], 0), 1)
    so$CytoTRACE2_Score   <- ct2
    so$CytoTRACE2_Potency <- ct2_result[["CytoTRACE2_Potency"]][colnames(so)]
  }
  cat(sprintf("      CytoTRACE2 range: [%.3f, %.3f]\n",
              min(ct2, na.rm = TRUE), max(ct2, na.rm = TRUE)))

  ## -- 4.2  POOLED SHANNON CYCLING ENTROPY -------------------------------------
  cat("[3/6] Computing pooled Shannon cell-cycle entropy (S+G2M together)...\n")
  cyc <- compute_pooled_cycle_entropy(expr_mat, s_use, g2m_use)
  so$H_cycle        <- cyc$H_cycle[colnames(so)]
  so$cycle_activity <- cyc$cycle_activity[colnames(so)]
  so$cycle_state    <- cyc$cycle_state[colnames(so)]

  ## -- 4.3  PER-PHASE SCORES: AUCell (S and G2M kept SEPARATE) -----------------
  cat("[4/6] Scoring S-phase and G2M with AUCell (separate)...\n")
  gene_sets <- list(S_phase = s_use, G2M_phase = g2m_use)
  rankings  <- AUCell::AUCell_buildRankings(round(expr_mat),
                                            plotStats = FALSE, verbose = FALSE)
  auc_mat   <- AUCell::getAUC(AUCell::AUCell_calcAUC(gene_sets, rankings, verbose = FALSE))
  so$AUCell_S_Score   <- auc_mat["S_phase",   colnames(so)]
  so$AUCell_G2M_Score <- auc_mat["G2M_phase", colnames(so)]
  cat(sprintf("      S-phase AUC   [%.3f, %.3f] | G2M AUC [%.3f, %.3f]\n",
              min(so$AUCell_S_Score),   max(so$AUCell_S_Score),
              min(so$AUCell_G2M_Score), max(so$AUCell_G2M_Score)))

  ## -- 4.4  Assign 2-D state quadrants (median splits) -------------------------
  ct_col <- if (!is.null(celltype_col)) celltype_col else
            if (!is.null(so$CellType)) "CellType" else
            if (!is.null(so$celltype)) "celltype" else
            if (!is.null(so$seurat_clusters)) "seurat_clusters" else NULL
  cell_type <- if (!is.null(ct_col) && ct_col %in% colnames(so@meta.data)) {
    v <- as.character(so@meta.data[[ct_col]])
    if (ct_col == "seurat_clusters") paste0("C", v) else v
  } else "unknown"

  cond_col  <- if (!is.null(condition_col)) condition_col else
               if ("Genotype_sex" %in% colnames(so@meta.data)) "Genotype_sex" else NULL
  condition <- if (!is.null(cond_col) && cond_col %in% colnames(so@meta.data))
                 as.character(so@meta.data[[cond_col]]) else NA_character_
  # Seurat cluster per cell -> pseudobulk resolves to cluster-level cell STATES
  # (several per cell type), not one blob per cell type. Falls back to cell type.
  clust_col <- if ("seurat_clusters" %in% colnames(so@meta.data)) "seurat_clusters" else
               if ("SCT_snn_res.0.8" %in% colnames(so@meta.data)) "SCT_snn_res.0.8" else NULL
  cluster   <- if (!is.null(clust_col))
                 paste0("C", as.character(so@meta.data[[clust_col]])) else cell_type
  # Two cycling metrics stored side by side. The plotting/quadrant 'activity'
  # axis is AUCell S+G2M by default (graded, rank-based, no banding); the entropy
  # (1 - H_cycle) version is kept for reference and is banded by gene-detection
  # count (few cycling genes -> artificially high). Pick with cycle_axis=.
  activity_auc     <- (so$AUCell_S_Score + so$AUCell_G2M_Score) / 2
  df <- data.frame(
    barcode          = colnames(so),
    potency          = ct2[colnames(so)],
    H_cycle          = so$H_cycle,
    activity_entropy = so$cycle_activity,              # 1 - H_cycle (banded)
    activity_auc     = activity_auc,                   # mean AUCell S+G2M (clean)
    S_AUC            = so$AUCell_S_Score,
    G2M_AUC          = so$AUCell_G2M_Score,
    celltype         = cell_type,
    cluster          = cluster,
    condition        = condition,
    row.names = colnames(so),
    stringsAsFactors = FALSE
  )
  df$activity  <- if (cycle_axis == "entropy") df$activity_entropy else df$activity_auc
  activity_lab <- if (cycle_axis == "entropy")
    "Cycling activity  (1 - H_cycle,  high = cycling)" else
    "Cycling activity  (mean AUCell S+G2M,  high = cycling)"
  p_thr <- median(df$potency, na.rm = TRUE)   # potency split is IDENTICAL for all axes

  ## -- 4.4b  QUADRANT SCHEME (biology-safe, axis-aware) ------------------------
  ## The cycling axis can be combined S+G2M, S-phase only, or G2M only; the quadrant
  ## names describe the two AXES directly -- potency (Stem-like vs Differentiated)
  ## x activity (high/low ON THAT AXIS) -- instead of asserting a specific cell
  ## identity (the old "Quiescent-Stem (LT-HSC)" tag was hematopoietic, wrong for a
  ## colon epithelium). Colours are fixed by quadrant position; only the activity
  ## wording changes per axis. Biological reading (colon): Stem-like+non-cycling =
  ## reserve/quiescent stem; Stem-like+cycling = active crypt-base stem;
  ## Differentiated+cycling = transit-amplifying; Differentiated+non-cycling = mature.
  .quad_words <- function(axis_key) switch(axis_key,
    combined = c(hi = "cycling",      lo = "non-cycling"),
    S        = c(hi = "S-phase high", lo = "S-phase low"),
    G2M      = c(hi = "G2/M high",    lo = "G2/M low"))
  .quad_scheme <- function(axis_key) {
    w <- .quad_words(axis_key)
    labs4 <- c(stem_lo = paste0("Stem-like\n", w[["lo"]]),       # hi potency, lo activity
               stem_hi = paste0("Stem-like\n", w[["hi"]]),       # hi potency, hi activity
               diff_hi = paste0("Differentiated\n", w[["hi"]]),  # lo potency, hi activity
               diff_lo = paste0("Differentiated\n", w[["lo"]]))  # lo potency, lo activity
    cols4 <- c(stem_lo = "#B2182B", stem_hi = "#F4A582",
               diff_hi = "#92C5DE", diff_lo = "#2166AC")
    list(labels = labs4, cols = cols4, levels = unname(labs4),
         state_cols = setNames(unname(cols4), unname(labs4)))
  }
  .classify_state <- function(potency, activity, p_thr, a_thr, sch) {
    key <- ifelse(is.na(potency) | is.na(activity), NA_character_,
           ifelse(potency >= p_thr & activity <  a_thr, "stem_lo",
           ifelse(potency >= p_thr & activity >= a_thr, "stem_hi",
           ifelse(potency <  p_thr & activity >= a_thr, "diff_hi", "diff_lo"))))
    factor(unname(sch$labels[key]), levels = sch$levels)
  }

  # combined-axis quadrants: used by the inline combined figures + the saved CSV
  .sch       <- .quad_scheme("combined")
  a_thr      <- median(df$activity, na.rm = TRUE)
  df$state   <- .classify_state(df$potency, df$activity, p_thr, a_thr, .sch)
  state_cols <- .sch$state_cols
  .ql        <- .sch$labels   # stem_lo, stem_hi, diff_hi, diff_lo (combined labels)

  ## -- 4.5  FIGURE 1 — potency x cycling-activity contour, state quadrants -----
  cat("[5/6] Drawing contour landscapes...\n")
  library(ggplot2)
  d1 <- df[!is.na(df$potency) & !is.na(df$activity), ]
  xmn <- min(d1$potency); xmx <- max(d1$potency)
  ymn <- min(d1$activity); ymx <- max(d1$activity)

  fig1 <- ggplot(d1, aes(potency, activity)) +
    geom_point(aes(colour = state), size = 0.5, alpha = 0.45) +
    geom_density_2d(colour = "grey25", linewidth = 0.3, alpha = 0.7) +
    geom_vline(xintercept = p_thr, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = a_thr, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    scale_colour_manual(values = state_cols, name = "State", na.value = "grey80") +
    annotate("text", x = xmx, y = ymn, label = .ql[["stem_lo"]],
             size = 3, hjust = 1, vjust = 0, colour = .sch$cols[["stem_lo"]], fontface = "bold") +
    annotate("text", x = xmx, y = ymx, label = .ql[["stem_hi"]],
             size = 3, hjust = 1, vjust = 1, colour = .sch$cols[["stem_hi"]], fontface = "bold") +
    annotate("text", x = xmn, y = ymx, label = .ql[["diff_hi"]],
             size = 3, hjust = 0, vjust = 1, colour = .sch$cols[["diff_hi"]], fontface = "bold") +
    annotate("text", x = xmn, y = ymn, label = .ql[["diff_lo"]],
             size = 3, hjust = 0, vjust = 0, colour = .sch$cols[["diff_lo"]], fontface = "bold") +
    theme_classic(base_size = 11) +
    labs(title    = "Cell-state landscape: potency vs cycling activity",
         subtitle = paste0("x = CytoTRACE2 potency (high = stem) | y = ",
                           if (cycle_axis == "entropy") "1 - H_cycle" else "mean AUCell S+G2M",
                           " (high = cycling)\nPoints = cells | contours = density | dashed = median thresholds"),
         x = "CytoTRACE2 potency  (0 = differentiated -> 1 = stem)",
         y = activity_lab)
  ggplot2::ggsave(file.path(out_dir, "fig_potency_cycle_contour.pdf"), fig1, width = 8, height = 6)
  ggplot2::ggsave(file.path(out_dir, "fig_potency_cycle_contour.png"), fig1, width = 8, height = 6, dpi = 200)

  ## -- 4.6  FIGURE 2 — S-score x G2M-score plane, coloured by potency ----------
  d2 <- df[!is.na(df$S_AUC) & !is.na(df$G2M_AUC) & !is.na(df$potency), ]
  fig2 <- ggplot(d2, aes(S_AUC, G2M_AUC)) +
    geom_point(aes(colour = potency), size = 0.6, alpha = 0.75) +
    geom_density_2d(colour = "grey25", linewidth = 0.3, alpha = 0.7) +
    ggplot2::scale_colour_viridis_c(option = "plasma", name = "CytoTRACE2\npotency") +
    theme_classic(base_size = 11) +
    labs(title    = "S vs G2M phase landscape, coloured by potency",
         subtitle = "x = S-phase AUCell | y = G2M AUCell | colour = CytoTRACE2 | contours = cell density",
         x = "S-phase score (AUCell)", y = "G2M score (AUCell)")
  ggplot2::ggsave(file.path(out_dir, "fig_phase_landscape.pdf"), fig2, width = 8, height = 6)
  ggplot2::ggsave(file.path(out_dir, "fig_phase_landscape.png"), fig2, width = 8, height = 6, dpi = 200)

  ## -- 4.6b FIGURE 3 — cell-type ISLANDS on the landscape ----------------------
  ## Each cell type's own density contours overlaid on the same potency x activity
  ## plane -> which cell types sit in which quadrant.
  n_ct <- length(unique(d1$celltype))
  # Legend driven by the POINTS layer (every cell type has points) not the contour
  # layer (density_2d is dropped for sparse types -> blank legend keys). Points carry
  # the legend with solid dots; contours are show.legend = FALSE.
  fig3 <- ggplot(d1, aes(potency, activity, colour = celltype)) +
    geom_point(size = 0.3, alpha = 0.12) +
    geom_density_2d(linewidth = 0.4, bins = 6, show.legend = FALSE) +
    geom_vline(xintercept = p_thr, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = a_thr, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    labs(title    = "Cell-type islands on the potency x cycling landscape",
         subtitle = "each contour set = one cell type's cell density | dashed = global median thresholds",
         x = "CytoTRACE2 potency  (0 = differentiated -> 1 = stem)",
         y = activity_lab, colour = "cell type") +
    theme_classic(base_size = 11) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
  ggplot2::ggsave(file.path(out_dir, "fig_celltype_islands.pdf"), fig3, width = 10, height = 7)
  ggplot2::ggsave(file.path(out_dir, "fig_celltype_islands.png"), fig3, width = 10, height = 7, dpi = 200)

  ## -- 4.6c FIGURE 4 — per-cell-type facets with LOCAL median potency ----------
  ## One mini-landscape per cell type; the red line is that cell type's own median
  ## potency (its "local" potency), on a shared axis so cell types are comparable.
  med_ct <- aggregate(potency ~ celltype, data = d1, FUN = median, na.rm = TRUE)
  names(med_ct)[2] <- "med_potency"
  fig4 <- ggplot(d1, aes(potency, activity)) +
    geom_point(size = 0.25, alpha = 0.18, colour = "grey45") +
    geom_density_2d(colour = "#08519C", linewidth = 0.3, bins = 6) +
    geom_vline(data = med_ct, aes(xintercept = med_potency),
               colour = "#B2182B", linewidth = 0.6) +
    facet_wrap(~ celltype) +
    labs(title    = "Per-cell-type landscape (local potency)",
         subtitle = "red line = that cell type's median CytoTRACE2 potency",
         x = "CytoTRACE2 potency  (0 = differentiated -> 1 = stem)", y = activity_lab) +
    theme_bw(base_size = 10) +
    theme(strip.text = element_text(face = "bold", size = 8),
          plot.title = element_text(face = "bold"))
  n_rows <- ceiling(n_ct / 4)
  ggplot2::ggsave(file.path(out_dir, "fig_celltype_facets.png"), fig4,
                  width = 14, height = max(6, n_rows * 3), dpi = 200, limitsize = FALSE)

  ## -- 4.6d FIGURE 5 — PSEUDOBULK landscape (per CLUSTER; de-noised) -----------
  ## Each point = ONE Seurat cluster's MEAN (potency, activity) -> a distinct cell
  ## STATE. Several points per cell type reveal within-type state heterogeneity.
  ## Sized by n cells, coloured by cell type (NO text labels). If conditions exist
  ## -> one panel per sex_phenotype group, means recomputed per cluster x group.
  has_cond <- !all(is.na(d1$condition))
  # Canonical sex_phenotype order used across the rest of the pipeline (02/03/05):
  # all FEMALE first, then all MALE, genotype WT -> Polyp -> KO. With facet_wrap
  # ncol = 3 this puts females in the top row, males in the bottom, genotype-aligned
  # columns. Falls back to sorted order for non-nr4a1 datasets.
  .canon_cond <- c("WT_Female", "Polyp_Female", "Polyp_NR4a1_KO_Female",
                   "WT_Male",   "Polyp_Male",   "Polyp_NR4a1_KO_Male")
  cond_lev <- if (has_cond) {
    l <- if (!is.null(condition_levels)) condition_levels
         else if (all(unique(stats::na.omit(d1$condition)) %in% .canon_cond)) .canon_cond
         else sort(unique(d1$condition))
    intersect(l, unique(d1$condition))
  } else NULL
  # Split each Seurat cluster BY CELL TYPE (not a single majority vote): one point
  # per cell type present in a cluster, so minority populations -- STEM CELLS above
  # all -- get their own points wherever they exist, instead of being swallowed by
  # the cluster's dominant type. Tiny sub-groups are dropped to limit clutter.
  min_sub  <- 15   # min cells per (cluster x cell type [x condition]) to plot a point
  grp_vars <- if (has_cond) c("cluster", "celltype", "condition")
              else c("cluster", "celltype")
  dagg <- d1[, c("potency", "activity", grp_vars), drop = FALSE]
  agg  <- aggregate(cbind(potency, activity) ~ ., data = dagg, FUN = mean, na.rm = TRUE)
  cnt  <- aggregate(list(n = d1$potency), by = d1[grp_vars], FUN = length)
  agg  <- merge(agg, cnt, by = grp_vars)
  agg  <- agg[agg$n >= min_sub, , drop = FALSE]
  if (has_cond) {
    d1$condition  <- factor(d1$condition,  levels = cond_lev)
    agg$condition <- factor(agg$condition, levels = cond_lev)
  }
  n_clust <- length(unique(agg$cluster))
  n_pts   <- nrow(agg)
  # draw STEM CELL points last so they sit on top (not hidden behind abundant types)
  agg <- agg[order(agg$celltype == "Stem cells"), , drop = FALSE]
  # zoom window (cluster-mean extent + pad); corners reused for quadrant annotations
  xr <- range(c(agg$potency,  p_thr), na.rm = TRUE)
  yr <- range(c(agg$activity, a_thr), na.rm = TRUE)
  xpad <- diff(xr) * 0.10 + 1e-6; ypad <- diff(yr) * 0.12 + 1e-6
  xlo <- xr[1] - xpad; xhi <- xr[2] + xpad; ylo <- yr[1] - ypad; yhi <- yr[2] + ypad
  fig5 <- ggplot() +
    geom_point(data = d1, aes(potency, activity), colour = "grey85", size = 0.15, alpha = 0.08) +
    geom_vline(xintercept = p_thr, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
    geom_hline(yintercept = a_thr, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
    geom_point(data = agg, aes(potency, activity, colour = celltype, size = n)) +
    # quadrant labels: four SCALAR annotate() calls (scalars recycle safely under
    # facets; a single vector annotate() mis-recycles to 4 corners x N panels).
    annotate("text", x = xhi, y = ylo, label = .ql[["stem_lo"]],
             colour = .sch$cols[["stem_lo"]], hjust = 1, vjust = 0, size = 2.6, fontface = "bold") +
    annotate("text", x = xhi, y = yhi, label = .ql[["stem_hi"]],
             colour = .sch$cols[["stem_hi"]], hjust = 1, vjust = 1, size = 2.6, fontface = "bold") +
    annotate("text", x = xlo, y = yhi, label = .ql[["diff_hi"]],
             colour = .sch$cols[["diff_hi"]], hjust = 0, vjust = 1, size = 2.6, fontface = "bold") +
    annotate("text", x = xlo, y = ylo, label = .ql[["diff_lo"]],
             colour = .sch$cols[["diff_lo"]], hjust = 0, vjust = 0, size = 2.6, fontface = "bold") +
    scale_size(range = c(2, 9), name = "n cells") +
    labs(title    = "Pseudobulk landscape (Seurat-cluster cell states)",
         subtitle = paste0(
           if (has_cond) "one panel per sex_phenotype | " else "",
           "each point = one cell type within a Seurat cluster (mean; ", n_pts,
           " pts over ", n_clust, " clusters, >=", min_sub,
           " cells) | colour = cell type | grey = single cells"),
         x = "CytoTRACE2 potency  (0 = differentiated -> 1 = stem)",
         y = activity_lab, colour = "cell type") +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(face = "bold"))
  if (has_cond) fig5 <- fig5 + facet_wrap(~ condition, ncol = 3)
  # zoom so the crowded centre is readable (corners above frame the quadrant labels)
  fig5 <- fig5 + coord_cartesian(xlim = c(xlo, xhi), ylim = c(ylo, yhi))
  pw <- if (has_cond) 15 else 10
  ph <- if (has_cond) 9  else 7
  ggplot2::ggsave(file.path(out_dir, "fig_pseudobulk_landscape.pdf"), fig5, width = pw, height = ph)
  ggplot2::ggsave(file.path(out_dir, "fig_pseudobulk_landscape.png"), fig5, width = pw, height = ph, dpi = 200)
  utils::write.csv(agg, file.path(out_dir, "pseudobulk_cluster_means.csv"), row.names = FALSE)

  ## -- 4.6e FIGURES 6 — potency, cycling, S, G2M BY CONDITION, Wilcoxon --------
  if (has_cond && length(unique(stats::na.omit(d1$condition))) >= 2) {
    lev <- cond_lev
    cmp <- comparisons
    if (is.null(cmp)) {
      default_cmp <- list(c("Polyp_NR4a1_KO_Female", "Polyp_Female"),
                          c("Polyp_Female", "WT_Female"),
                          c("Polyp_NR4a1_KO_Male", "Polyp_Male"),
                          c("Polyp_Male", "WT_Male"))
      cmp <- Filter(function(x) all(x %in% lev), default_cmp)
    }
    have_ggpubr <- requireNamespace("ggpubr", quietly = TRUE)
    .violin_by_cond <- function(metric, ylab, fname) {
      dd <- d1[!is.na(d1[[metric]]) & !is.na(d1$condition), ]
      dd$condition <- factor(dd$condition, levels = lev)
      p <- ggplot(dd, aes(.data[["condition"]], .data[[metric]], fill = .data[["condition"]])) +
        geom_violin(scale = "width", trim = TRUE, alpha = 0.7, linewidth = 0.3) +
        geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.7, linewidth = 0.3) +
        facet_wrap(~ celltype, scales = "free_y") +
        labs(title = paste0(ylab, " by condition"), x = NULL, y = ylab, fill = "condition") +
        # legend laid out as a 2-row grid filled BY ROW: females top, males bottom
        guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
        theme_classic(base_size = 11) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              legend.position = "bottom", strip.text = element_text(face = "bold"))
      if (have_ggpubr && length(cmp) > 0)
        p <- p + ggpubr::stat_compare_means(comparisons = cmp, method = "wilcox.test",
                   method.args = list(exact = FALSE), label = "p.signif",
                   step.increase = 0.08, size = 3)
      ggplot2::ggsave(file.path(out_dir, fname), p,
                      width = 14, height = max(6, n_rows * 3.2), dpi = 200,
                      bg = "white", limitsize = FALSE)
    }
    if (!have_ggpubr)
      cat("      [note] ggpubr not installed -> condition violins without Wilcoxon brackets.\n")
    .violin_by_cond("potency",  "CytoTRACE2 potency", "fig_violin_potency_by_condition.png")
    .violin_by_cond("activity", activity_lab,         "fig_violin_activity_by_condition.png")
    # S-phase and G2M SEPARATELY -> read cycling changes across cancer stage
    # (WT -> Polyp -> KO) and sex, per cell type (comparisons stay within-sex).
    .violin_by_cond("S_AUC",   "S-phase activity (AUCell)", "fig_violin_Sphase_by_condition.png")
    .violin_by_cond("G2M_AUC", "G2M activity (AUCell)",     "fig_violin_G2M_by_condition.png")

    ## -- 4.6f FIGURE 7 — quadrant-state COMPOSITION per cell type x condition --
    cdf <- d1[!is.na(d1$state) & !is.na(d1$celltype) & !is.na(d1$condition), ]
    if (nrow(cdf) > 0) {
      cdf$condition <- factor(cdf$condition, levels = lev)
      comp <- aggregate(list(n = rep(1L, nrow(cdf))),
                        by = list(celltype = cdf$celltype, condition = cdf$condition,
                                  state = cdf$state), FUN = sum)
      tot  <- aggregate(list(tot = comp$n),
                        by = list(celltype = comp$celltype, condition = comp$condition), FUN = sum)
      comp <- merge(comp, tot, by = c("celltype", "condition"))
      comp$pct <- 100 * comp$n / comp$tot
      pcomp <- ggplot(comp, aes(condition, pct, fill = state)) +
        geom_col(colour = "white", linewidth = 0.2) +
        facet_wrap(~ celltype, scales = "free_x") +
        scale_fill_manual(values = state_cols, name = "State", drop = FALSE) +
        labs(title = "Quadrant-state composition by condition",
             subtitle = "% of each cell type's cells in the 4 states, per condition",
             x = NULL, y = "% of cells") +
        theme_bw(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              strip.text = element_text(face = "bold"),
              plot.title = element_text(face = "bold"))
      ggplot2::ggsave(file.path(out_dir, "fig_quadrant_composition.png"), pcomp,
                      width = 14, height = max(6, n_rows * 3), dpi = 200,
                      bg = "white", limitsize = FALSE)
      utils::write.csv(comp, file.path(out_dir, "quadrant_composition.csv"), row.names = FALSE)
    }
  } else {
    cat("      [note] no usable condition column -> skipping by-condition violins/composition.\n")
  }

  ## -- 4.6g  PHASE-SPECIFIC LANDSCAPES: S-only and G2M-only --------------------
  ## Redraw the landscape figures (contour, islands, facets, pseudobulk, quadrant
  ## composition) using S-phase ONLY and G2M ONLY as the cycling axis -- each with
  ## its OWN median split, so the quadrants are phase-specific (a cell high in S but
  ## low in G2M, or vice versa, lands in a different quadrant than on the combined
  ## axis). Files get a _Sphase / _G2M suffix; the combined-axis figures are above.
  .phase_landscape <- function(act_col, axis_key, act_lab, sfx) {
    sch  <- .quad_scheme(axis_key); ql <- sch$labels; scols <- sch$state_cols
    dloc <- df; dloc$activity <- dloc[[act_col]]
    d1 <- dloc[!is.na(dloc$potency) & !is.na(dloc$activity), ]
    if (nrow(d1) < 10) { message("  [skip ", axis_key, "] too few scored cells"); return(invisible()) }
    a_thr <- median(d1$activity, na.rm = TRUE)
    d1$state <- .classify_state(d1$potency, d1$activity, p_thr, a_thr, sch)
    xmn <- min(d1$potency); xmx <- max(d1$potency)
    ymn <- min(d1$activity); ymx <- max(d1$activity)
    n_ct <- length(unique(d1$celltype)); n_rows <- ceiling(n_ct / 4)

    ## FIGURE 1 — contour + quadrant states
    fig1 <- ggplot(d1, aes(potency, activity)) +
      geom_point(aes(colour = state), size = 0.5, alpha = 0.45) +
      geom_density_2d(colour = "grey25", linewidth = 0.3, alpha = 0.7) +
      geom_vline(xintercept = p_thr, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
      geom_hline(yintercept = a_thr, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
      scale_colour_manual(values = scols, name = "State", na.value = "grey80") +
      annotate("text", x = xmx, y = ymn, label = ql[["stem_lo"]], size = 3, hjust = 1, vjust = 0, colour = sch$cols[["stem_lo"]], fontface = "bold") +
      annotate("text", x = xmx, y = ymx, label = ql[["stem_hi"]], size = 3, hjust = 1, vjust = 1, colour = sch$cols[["stem_hi"]], fontface = "bold") +
      annotate("text", x = xmn, y = ymx, label = ql[["diff_hi"]], size = 3, hjust = 0, vjust = 1, colour = sch$cols[["diff_hi"]], fontface = "bold") +
      annotate("text", x = xmn, y = ymn, label = ql[["diff_lo"]], size = 3, hjust = 0, vjust = 0, colour = sch$cols[["diff_lo"]], fontface = "bold") +
      theme_classic(base_size = 11) +
      labs(title = paste0("Cell-state landscape: potency vs ", act_lab),
           subtitle = paste0("x = CytoTRACE2 potency (high = stem) | y = ", act_lab,
                             "\nPoints = cells | contours = density | dashed = median thresholds"),
           x = "CytoTRACE2 potency  (0 = differentiated -> 1 = stem)", y = act_lab)
    ggplot2::ggsave(file.path(out_dir, paste0("fig_potency_cycle_contour", sfx, ".png")), fig1, width = 8, height = 6, dpi = 200)

    ## FIGURE 3 — cell-type islands
    fig3 <- ggplot(d1, aes(potency, activity, colour = celltype)) +
      geom_point(size = 0.3, alpha = 0.12) +
      geom_density_2d(linewidth = 0.4, bins = 6, show.legend = FALSE) +
      geom_vline(xintercept = p_thr, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
      geom_hline(yintercept = a_thr, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
      labs(title = paste0("Cell-type islands (", act_lab, ")"),
           subtitle = "each contour set = one cell type's cell density | dashed = global median thresholds",
           x = "CytoTRACE2 potency  (0 = differentiated -> 1 = stem)", y = act_lab, colour = "cell type") +
      theme_classic(base_size = 11) +
      guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
    ggplot2::ggsave(file.path(out_dir, paste0("fig_celltype_islands", sfx, ".png")), fig3, width = 10, height = 7, dpi = 200)

    ## FIGURE 4 — per-cell-type facets
    med_ct <- aggregate(potency ~ celltype, data = d1, FUN = median, na.rm = TRUE)
    names(med_ct)[2] <- "med_potency"
    fig4 <- ggplot(d1, aes(potency, activity)) +
      geom_point(size = 0.25, alpha = 0.18, colour = "grey45") +
      geom_density_2d(colour = "#08519C", linewidth = 0.3, bins = 6) +
      geom_vline(data = med_ct, aes(xintercept = med_potency), colour = "#B2182B", linewidth = 0.6) +
      facet_wrap(~ celltype) +
      labs(title = paste0("Per-cell-type landscape (", act_lab, ")"),
           subtitle = "red line = that cell type's median CytoTRACE2 potency",
           x = "CytoTRACE2 potency  (0 = differentiated -> 1 = stem)", y = act_lab) +
      theme_bw(base_size = 10) +
      theme(strip.text = element_text(face = "bold", size = 8), plot.title = element_text(face = "bold"))
    ggplot2::ggsave(file.path(out_dir, paste0("fig_celltype_facets", sfx, ".png")), fig4, width = 14, height = max(6, n_rows * 3), dpi = 200, limitsize = FALSE)

    ## FIGURE 5 — pseudobulk landscape (cluster x cell type)
    has_cond <- !all(is.na(d1$condition)); min_sub <- 15
    cl <- if (has_cond) {
      l <- if (!is.null(condition_levels)) condition_levels
           else if (all(unique(stats::na.omit(d1$condition)) %in% .canon_cond)) .canon_cond
           else sort(unique(stats::na.omit(d1$condition)))
      intersect(l, unique(stats::na.omit(d1$condition)))
    } else NULL
    grp_vars <- if (has_cond) c("cluster", "celltype", "condition") else c("cluster", "celltype")
    dagg <- d1[, c("potency", "activity", grp_vars), drop = FALSE]
    agg  <- aggregate(cbind(potency, activity) ~ ., data = dagg, FUN = mean, na.rm = TRUE)
    cnt  <- aggregate(list(n = d1$potency), by = d1[grp_vars], FUN = length)
    agg  <- merge(agg, cnt, by = grp_vars); agg <- agg[agg$n >= min_sub, , drop = FALSE]
    if (has_cond) { d1$condition <- factor(d1$condition, levels = cl); agg$condition <- factor(agg$condition, levels = cl) }
    n_clust <- length(unique(agg$cluster)); n_pts <- nrow(agg)
    agg <- agg[order(agg$celltype == "Stem cells"), , drop = FALSE]
    xr <- range(c(agg$potency, p_thr), na.rm = TRUE); yr <- range(c(agg$activity, a_thr), na.rm = TRUE)
    xpad <- diff(xr) * 0.10 + 1e-6; ypad <- diff(yr) * 0.12 + 1e-6
    xlo <- xr[1]-xpad; xhi <- xr[2]+xpad; ylo <- yr[1]-ypad; yhi <- yr[2]+ypad
    fig5 <- ggplot() +
      geom_point(data = d1, aes(potency, activity), colour = "grey85", size = 0.15, alpha = 0.08) +
      geom_vline(xintercept = p_thr, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
      geom_hline(yintercept = a_thr, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
      geom_point(data = agg, aes(potency, activity, colour = celltype, size = n)) +
      annotate("text", x = xhi, y = ylo, label = ql[["stem_lo"]], colour = sch$cols[["stem_lo"]], hjust = 1, vjust = 0, size = 2.6, fontface = "bold") +
      annotate("text", x = xhi, y = yhi, label = ql[["stem_hi"]], colour = sch$cols[["stem_hi"]], hjust = 1, vjust = 1, size = 2.6, fontface = "bold") +
      annotate("text", x = xlo, y = yhi, label = ql[["diff_hi"]], colour = sch$cols[["diff_hi"]], hjust = 0, vjust = 1, size = 2.6, fontface = "bold") +
      annotate("text", x = xlo, y = ylo, label = ql[["diff_lo"]], colour = sch$cols[["diff_lo"]], hjust = 0, vjust = 0, size = 2.6, fontface = "bold") +
      scale_size(range = c(2, 9), name = "n cells") +
      labs(title = paste0("Pseudobulk landscape (", act_lab, ")"),
           subtitle = paste0(if (has_cond) "one panel per sex_phenotype | " else "",
             "each point = one cell type within a Seurat cluster (mean; ", n_pts, " pts, >=", min_sub, " cells)"),
           x = "CytoTRACE2 potency  (0 = differentiated -> 1 = stem)", y = act_lab, colour = "cell type") +
      theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold"))
    if (has_cond) fig5 <- fig5 + facet_wrap(~ condition, ncol = 3)
    fig5 <- fig5 + coord_cartesian(xlim = c(xlo, xhi), ylim = c(ylo, yhi))
    pw <- if (has_cond) 15 else 10; ph <- if (has_cond) 9 else 7
    ggplot2::ggsave(file.path(out_dir, paste0("fig_pseudobulk_landscape", sfx, ".png")), fig5, width = pw, height = ph, dpi = 200)
    utils::write.csv(agg, file.path(out_dir, paste0("pseudobulk_cluster_means", sfx, ".csv")), row.names = FALSE)

    ## FIGURE 7 — quadrant-state composition per cell type x condition
    if (has_cond) {
      cdf <- d1[!is.na(d1$state) & !is.na(d1$celltype) & !is.na(d1$condition), ]
      if (nrow(cdf) > 0) {
        cdf$condition <- factor(cdf$condition, levels = cl)
        comp <- aggregate(list(n = rep(1L, nrow(cdf))), by = list(celltype = cdf$celltype, condition = cdf$condition, state = cdf$state), FUN = sum)
        tot  <- aggregate(list(tot = comp$n), by = list(celltype = comp$celltype, condition = comp$condition), FUN = sum)
        comp <- merge(comp, tot, by = c("celltype", "condition")); comp$pct <- 100 * comp$n / comp$tot
        pcomp <- ggplot(comp, aes(condition, pct, fill = state)) +
          geom_col(colour = "white", linewidth = 0.2) +
          facet_wrap(~ celltype, scales = "free_x") +
          scale_fill_manual(values = scols, name = "State", drop = FALSE) +
          labs(title = paste0("Quadrant-state composition (", act_lab, ")"),
               subtitle = "% of each cell type's cells in the 4 states, per condition", x = NULL, y = "% of cells") +
          theme_bw(base_size = 10) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))
        ggplot2::ggsave(file.path(out_dir, paste0("fig_quadrant_composition", sfx, ".png")), pcomp, width = 14, height = max(6, n_rows * 3), dpi = 200, bg = "white", limitsize = FALSE)
        utils::write.csv(comp, file.path(out_dir, paste0("quadrant_composition", sfx, ".csv")), row.names = FALSE)
      }
    }
    message("  phase landscape done: ", axis_key, " (", sfx, ")")
  }
  cat("[5b/6] Drawing S-only and G2M-only landscapes...\n")
  .phase_landscape("S_AUC",   "S",   "S-phase activity (AUCell)", "_Sphase")
  .phase_landscape("G2M_AUC", "G2M", "G2/M activity (AUCell)",    "_G2M")

  ## -- 4.7  Save scores + updated object ---------------------------------------
  cat("[6/6] Saving outputs...\n")
  write.csv(df[, c("barcode","potency","H_cycle","activity_entropy","activity_auc",
                   "activity","S_AUC","G2M_AUC","celltype","cluster","condition","state")],
            file.path(out_dir, "state_scores.csv"), row.names = FALSE)
  saveRDS(so, file.path(out_dir, "seurat_with_state.rds"))

  cat("\n-- State composition -------------------------------------------------\n")
  print(table(df$state, useNA = "ifany"))
  cat("\nOutputs written to: ", normalizePath(out_dir), "\n", sep = "")
  cat("  (cycle axis = ", cycle_axis, ")\n", sep = "")
  cat("  fig_potency_cycle_contour.*  -- potency x cycling-activity, quadrant states\n")
  cat("  fig_phase_landscape.*        -- S vs G2M plane, coloured by potency\n")
  cat("  fig_celltype_islands.*       -- per-cell-type density contours (islands)\n")
  cat("  fig_celltype_facets.png      -- per-cell-type mini-landscapes + local median potency\n")
  cat("  fig_pseudobulk_landscape.*   -- de-noised: one mean point per Seurat cluster (cell state)\n")
  cat("  fig_violin_potency_by_condition.png  -- potency by condition (Wilcoxon)\n")
  cat("  fig_violin_activity_by_condition.png -- cycling by condition (Wilcoxon)\n")
  cat("  fig_violin_Sphase_by_condition.png   -- S-phase by condition (Wilcoxon)\n")
  cat("  fig_violin_G2M_by_condition.png      -- G2M by condition (Wilcoxon)\n")
  cat("  fig_quadrant_composition.png -- %% cells in each state per cell type x condition\n")
  cat("  --- phase-specific landscapes (own median quadrants) ---\n")
  cat("  *_Sphase.png  -- contour/islands/facets/pseudobulk/composition, S-phase axis only\n")
  cat("  *_G2M.png     -- same set, G2/M axis only\n")
  cat("  state_scores.csv / pseudobulk_cluster_means{,_Sphase,_G2M}.csv / quadrant_composition{,_Sphase,_G2M}.csv\n")
  cat("  seurat_with_state.rds        -- Seurat object with all columns added\n")
  cat("----------------------------------------------------------------------\n")

  invisible(list(seurat = so, scores = df,
                 thresholds = c(potency = p_thr, activity = a_thr)))
}


## =============================================================================
## 5.  PIPELINE DRIVER (standalone; hardcoded paths)
## =============================================================================
## Runs only when this file is executed directly (Rscript 09b_...R). When it is
## source()d from another script the function above stays reusable and nothing
## runs. Edit the three paths below for a different project/machine.
if (sys.nframe() == 0L) {
  PROJECT_NAME <- "Nr4a1_s17_ack"
  ROOT_PATH    <- "/run/media/external/scdc_projects/2026_nr4a1_ack/r_process"
  OUTPUT_DIR   <- file.path(ROOT_PATH, "seurat_output")

  # Input = the cell-scores object from Script 09 (carries CytoTRACE2_Score).
  RDS_PATH <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_with_cell_scores.rds"))
  SPECIES  <- "mouse"

  # --- WHICH CELLS to build the landscape on ---------------------------------
  # Restrict to one broad lineage and facet by its subtypes. Set SUBSET_COL and
  # SUBSET_VAL to NULL to use the whole object (then CELLTYPE_COL can stay NULL).
  SUBSET_COL   <- "CellType_broad"   # column holding the broad labels; NULL = all cells
  SUBSET_VAL   <- "Colonocytes"      # the broad type(s) to pull
  CELLTYPE_COL <- "CellType"         # column used for islands/facets (the subtypes)

  # Output dir gets a subfolder per subset so runs don't overwrite each other.
  tag     <- if (!is.null(SUBSET_VAL))
               gsub("[^A-Za-z0-9]+", "_", paste(SUBSET_VAL, collapse = "_")) else "all_cells"
  OUT_DIR <- file.path(OUTPUT_DIR, "state_landscape", tag)

  if (!file.exists(RDS_PATH))
    stop("RDS not found: ", RDS_PATH, "\n  Run Script 09 first.", call. = FALSE)
  potency_cycle_contour(RDS_PATH, out_dir = OUT_DIR, species = SPECIES,
                        subset_col = SUBSET_COL, subset_val = SUBSET_VAL,
                        celltype_col = CELLTYPE_COL)
}
