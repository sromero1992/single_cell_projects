# ============================================================================
# TimeSCape v0.2 — Pathway Circadian Analysis
#
# Pipeline:
#   1. pull_genesets()      — download GOBP / KEGG gene sets via msigdbr
#   2. .deduplicate_genesets() — remove "POSITIVE_REGULATION_OF_X" when "X" exists
#   3. auc_score_cells()    — AUCell scores (one score per cell per pathway)
#   4. pathway_cosinor()    — cosinor fit on AUCell scores → pathway-level circadian
#   5. write_pathway_results() — Excel workbook (all + confident sheets)
#   6. .build_zt_network()  — Pearson correlation GRN at one ZT time point
#   7. plot_grn_timeseries() — 8-panel GRN grid across ZT time points
# ============================================================================


# ── 0. Interactive collection chooser ────────────────────────────────────────

#' Interactively choose an MSigDB collection and subcollection
#'
#' Prints a numbered table of all available collections from msigdbr and
#' prompts the user to pick one.  Returns a list with \code{$collection} and
#' \code{$subcategory} ready to pass directly to \code{pull_genesets()}.
#'
#' In non-interactive sessions (e.g. Rscript batch runs) the function returns
#' the supplied defaults without prompting.
#'
#' @param default_collection  Fallback collection if session is non-interactive.
#' @param default_subcategory Fallback subcategory if session is non-interactive.
#'
#' @return A named list: \code{list(collection = "C2", subcategory = "CP:KEGG_LEGACY")}.
#' @export
choose_collection <- function(default_collection  = "C2",
                               default_subcategory = "CP:KEGG_LEGACY") {
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("Install msigdbr first:  install.packages('msigdbr')")

  cols <- tryCatch(
    as.data.frame(msigdbr::msigdbr_collections()),
    error = function(e) stop("Could not fetch msigdbr collections: ", e$message)
  )

  # Build display label: "C2 | CP:KEGG_LEGACY | KEGG Legacy Pathways (186)"
  labels <- mapply(function(cat, sub, name, n) {
    code <- if (nchar(trimws(sub)) > 0)
               paste0(cat, " | ", sub)
             else
               cat
    sprintf("%-22s  %-45s [%d gene sets]", code, name, n)
  },
  cols$gs_collection,
  cols$gs_subcollection,
  cols$gs_collection_name,
  cols$num_genesets,
  SIMPLIFY = TRUE)

  if (!interactive()) {
    message("Non-interactive session — using default: ",
            default_collection, "/", default_subcategory)
    return(list(collection  = default_collection,
                subcategory = default_subcategory))
  }

  cat("\n Available MSigDB collections (msigdbr ", as.character(packageVersion("msigdbr")), "):\n", sep = "")
  cat(rep("-", 80), "\n", sep = "")
  for (i in seq_along(labels)) cat(sprintf(" %2d. %s\n", i, labels[i]))
  cat(rep("-", 80), "\n", sep = "")

  choice <- utils::menu(labels, title = "\nSelect a collection (0 to cancel):")
  if (choice == 0L) stop("No collection selected.")

  selected_cat <- cols$gs_collection[choice]
  selected_sub <- cols$gs_subcollection[choice]

  cat(sprintf("\n  Selected: %s / %s  (%s)\n",
              selected_cat,
              if (nchar(trimws(selected_sub)) > 0) selected_sub else "(no subcategory)",
              cols$gs_collection_name[choice]))

  list(collection  = selected_cat,
       subcategory = selected_sub)
}


# ── 1. Pull and deduplicate gene sets ─────────────────────────────────────────

#' Download gene sets from MSigDB via msigdbr, then deduplicate
#'
#' Removes "child" terms whose name is a prefix-extended form of a shorter
#' "parent" term — e.g. "GOBP_POSITIVE_REGULATION_OF_T_CELL_ACTIVATION" is
#' dropped when "GOBP_T_CELL_ACTIVATION" already exists in the set.
#'
#' @param collection  MSigDB collection code. Common choices:
#'   \describe{
#'     \item{"C5"}{Gene Ontology (use with subcategory "GO:BP", "GO:MF", "GO:CC")}
#'     \item{"C2"}{Curated (use with subcategory "CP:KEGG", "CP:REACTOME")}
#'   }
#' @param subcategory Subcategory within the collection.
#'   Examples: \code{"GO:BP"}, \code{"CP:KEGG"}, \code{"CP:REACTOME"}.
#' @param species     Species name as in msigdbr. Default \code{"Mus musculus"}.
#' @param min_size    Minimum gene set size (default 10).
#' @param max_size    Maximum gene set size (default 500).
#' @param deduplicate Logical. Remove positive/negative regulation variants
#'   when the base process exists (default TRUE).
#'
#' @return A named list: each element is a character vector of gene symbols.
#' @export
pull_genesets <- function(
    collection  = "C2",
    subcategory = "CP:KEGG_LEGACY",
    species     = "Mus musculus",
    min_size    = 10L,
    max_size    = 500L,
    deduplicate = FALSE
) {
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("Install msigdbr first:  install.packages('msigdbr')")

  # MSigDB renamed CP:KEGG → CP:KEGG_LEGACY in recent releases (v2023+).
  # Try the requested subcategory first; if it fails with "Unknown subcollection",
  # automatically try the known aliases before giving up.
  kegg_aliases <- list(
    "CP:KEGG"        = c("CP:KEGG_LEGACY", "CP:KEGG_MEDICUS"),
    "CP:KEGG_LEGACY" = c("CP:KEGG"),
    "KEGG"           = c("CP:KEGG_LEGACY", "CP:KEGG")
  )

  try_subcats <- c(subcategory, kegg_aliases[[subcategory]])
  try_subcats <- unique(try_subcats[!is.na(try_subcats)])

  db <- NULL
  for (sc in try_subcats) {
    message("Downloading MSigDB gene sets: ", collection, "/", sc,
            " for ", species, " ...")
    db <- tryCatch(
      msigdbr::msigdbr(species     = species,
                       category    = collection,
                       subcategory = sc),
      error = function(e) {
        if (grepl("Unknown subcollection|Unknown sub", e$message, ignore.case = TRUE)) {
          message("  '", sc, "' not found in this msigdbr version — trying next alias...")
          NULL
        } else {
          stop(e)
        }
      }
    )
    if (!is.null(db)) { subcategory <- sc; break }
  }

  if (is.null(db)) {
    avail <- tryCatch(
      paste(unique(msigdbr::msigdbr_collections()$gs_subcat), collapse = ", "),
      error = function(e) "(run msigdbr::msigdbr_collections() to check)"
    )
    stop(sprintf(
      "Could not download '%s/%s'. Available subcategories:\n  %s",
      collection, subcategory, avail))
  }

  # Build named list: pathway_name → unique gene symbols
  gs_list <- split(db$gene_symbol, db$gs_name)
  gs_list <- lapply(gs_list, unique)

  # Size filter
  sz      <- lengths(gs_list)
  gs_list <- gs_list[sz >= min_size & sz <= max_size]
  message(sprintf("  %d gene sets after size filter [%d–%d genes]",
                  length(gs_list), min_size, max_size))

  # Deduplication
  if (deduplicate) {
    gs_list <- .deduplicate_genesets(gs_list)
    message(sprintf("  %d gene sets after deduplication", length(gs_list)))
  }

  gs_list
}


#' Remove redundant positive/negative/regulation variants from a gene-set list
#'
#' Case-insensitive removal of five modifier prefixes when the base process
#' already exists in the list:
#'   POSITIVE_REGULATION_OF_X  → base X
#'   NEGATIVE_REGULATION_OF_X  → base X
#'   REGULATION_OF_X           → base X
#'   POSITIVE_X                → base X
#'   NEGATIVE_X                → base X
#'
#' Only removes a variant if the exact base name (same prefix, same suffix,
#' same case as stored) is present after size filtering.
#'
#' @param gs_list  Named list of gene-set character vectors.
#' @param verbose  Print a table of removed sets (default FALSE).
#' @return Deduplicated named list.
#' @keywords internal
.deduplicate_genesets <- function(gs_list, verbose = FALSE) {
  nms    <- names(gs_list)
  nms_up <- toupper(nms)           # upper-case copy for case-insensitive matching

  remove   <- logical(length(nms))
  base_of  <- character(length(nms))   # records which base each variant maps to
  reason   <- character(length(nms))

  # Five modifier patterns, checked in order of specificity (longest first)
  modifiers <- c(
    "POSITIVE_REGULATION_OF_",
    "NEGATIVE_REGULATION_OF_",
    "REGULATION_OF_",
    "POSITIVE_",
    "NEGATIVE_"
  )

  # Build a regex that matches any of the five modifiers (case-insensitive via toupper)
  mod_alt <- paste(modifiers, collapse = "|")
  rgx     <- paste0("^(.+?)_(", mod_alt, ")(.+)$")

  m <- regmatches(nms_up, regexec(rgx, nms_up, perl = TRUE))

  for (i in seq_along(nms)) {
    parts <- m[[i]]
    if (length(parts) < 4) next          # no modifier match — keep as-is

    prefix_up  <- parts[2]               # e.g. "GOBP"
    mod_up     <- parts[3]               # e.g. "POSITIVE_REGULATION_OF_"
    suffix_up  <- parts[4]               # e.g. "T_CELL_ACTIVATION"
    base_up    <- paste0(prefix_up, "_", suffix_up)

    # Find the actual (original-case) base in the list
    hit <- which(nms_up == base_up)
    if (length(hit) > 0) {
      remove[i]  <- TRUE
      base_of[i] <- nms[hit[1]]
      reason[i]  <- gsub("_$", "", mod_up)   # tidy label for the table
    }
  }

  n_removed <- sum(remove)
  message(sprintf("  Removed %d modifier-variant sets (base process present)",
                  n_removed))

  if (verbose && n_removed > 0) {
    tbl <- data.frame(
      removed  = nms[remove],
      base     = base_of[remove],
      modifier = reason[remove],
      stringsAsFactors = FALSE
    )
    print(tbl, row.names = FALSE)
  }

  gs_list[!remove]
}


#' Inspect gene-set list: show names, sizes, and deduplication status
#'
#' Useful for debugging.  Returns (and optionally writes) a data.frame with
#' every gene set name, its size, and whether it would be removed by
#' \code{.deduplicate_genesets()}.
#'
#' @param gs_list   Named list of gene-set vectors (before or after dedup).
#' @param outfile   Optional CSV path to write the table. NULL = no file.
#' @return A data.frame with columns: name, size, modifier, base, would_remove.
#' @export
inspect_genesets <- function(gs_list, outfile = NULL) {
  nms    <- names(gs_list)
  nms_up <- toupper(nms)

  modifiers <- c(
    "POSITIVE_REGULATION_OF_",
    "NEGATIVE_REGULATION_OF_",
    "REGULATION_OF_",
    "POSITIVE_",
    "NEGATIVE_"
  )
  mod_alt <- paste(modifiers, collapse = "|")
  rgx     <- paste0("^(.+?)_(", mod_alt, ")(.+)$")
  m       <- regmatches(nms_up, regexec(rgx, nms_up, perl = TRUE))

  modifier_col    <- character(length(nms))
  base_col        <- character(length(nms))
  would_remove    <- logical(length(nms))

  for (i in seq_along(nms)) {
    parts <- m[[i]]
    if (length(parts) < 4) next
    prefix_up <- parts[2]; mod_up <- parts[3]; suffix_up <- parts[4]
    base_up   <- paste0(prefix_up, "_", suffix_up)
    hit       <- which(nms_up == base_up)
    modifier_col[i] <- gsub("_$", "", mod_up)
    if (length(hit) > 0) {
      base_col[i]     <- nms[hit[1]]
      would_remove[i] <- TRUE
    }
  }

  tbl <- data.frame(
    name         = nms,
    size         = lengths(gs_list),
    modifier     = modifier_col,
    base         = base_col,
    would_remove = would_remove,
    stringsAsFactors = FALSE
  )

  cat(sprintf("Total sets     : %d\n", nrow(tbl)))
  cat(sprintf("Would remove   : %d\n", sum(tbl$would_remove)))
  cat(sprintf("Would keep     : %d\n", sum(!tbl$would_remove)))
  cat(sprintf("Modifier types removed:\n"))
  print(table(tbl$modifier[tbl$would_remove]))

  if (!is.null(outfile)) {
    write.csv(tbl, outfile, row.names = FALSE)
    message("Written → ", outfile)
  }

  invisible(tbl)
}


# ── 2. Phase-bin enrichment ───────────────────────────────────────────────────

#' Normalize a pathway name for cross-source matching
#'
#' Strips collection prefixes (KEGG_, GOBP_, REACTOME_, etc.), converts to
#' lowercase, and collapses punctuation/whitespace into single spaces.
#' Used internally to align phyper and EnrichR pathway names.
#'
#' @param x Character vector of pathway names.
#' @return Character vector of normalized names (same length as x).
#' @keywords internal
.norm_pw_name <- function(x) {
  # Remove MSigDB collection prefixes
  x <- gsub("^(KEGG|GOBP|GOCC|GOMF|REACTOME|HALLMARK|WP|HP|BIOCARTA|PID)_",
             "", x, ignore.case = TRUE, perl = TRUE)
  # Remove EnrichR database suffixes like " (Homo sapiens)" or trailing digits
  x <- gsub("\\s*\\([^)]+\\)\\s*$", "", x)
  # Normalize separators: _ / - → space
  x <- gsub("[_/\\-]+", " ", x)
  # Collapse whitespace, lowercase
  tolower(trimws(gsub("\\s+", " ", x)))
}


#' Phase-bin enrichment: group circadian genes by acrophase, enrich each bin,
#' and build phase-restricted gene sets for AUCell
#'
#' Genes in the same pathway but different acrophases CANCEL each other when
#' scored together with AUCell (a ZT6-peaking gene and a ZT18-peaking gene
#' in the same set produce a flat average AUC with no circadian signal).
#' This function solves that by:
#'   1. Binning confident circadian genes into acrophase windows.
#'   2. Running over-representation analysis (ORA) per bin:
#'      a. Hypergeometric test (phyper) against the local genesets — uses the
#'         correct tested-gene background (~13,790 genes) instead of all genome.
#'      b. (Optional) EnrichR API query — multi-database, no local download,
#'         but uses genome-wide background (~20,000 genes).
#'      Pathways significant in BOTH tests are labelled "consensus" and are the
#'      highest-confidence hits; those found only in one test are still reported.
#'   3. Building PHASE-RESTRICTED gene sets = pathway genes ∩ bin genes.
#'      All members peak within the same narrow phase window, so their
#'      AUCell scores will oscillate coherently.
#'
#' @param T1          Stats data.frame returned by \code{run_timescape()}$T1.
#' @param conf_mask   Logical vector (length = nrow(T1)) marking confident genes.
#' @param genesets    Named list of gene-symbol vectors (from \code{pull_genesets()}).
#' @param bin_width   Phase window width in hours (default 1). Smaller = tighter
#'                    co-regulation; 3 = one ZT interval for typical 8-ZT designs.
#' @param n_top       Max enriched pathways to keep per bin (default 5).
#' @param min_overlap Minimum gene overlap required for ORA (default 3).
#' @param min_bin_genes Minimum genes in a bin to run ORA (default 5).
#' @param p_thresh    P-value threshold for ORA (default 0.05).
#' @param use_padj    Logical. If TRUE, filter on BH-adjusted p-value; if FALSE
#'                    (default), filter on raw p-value.  Raw is recommended for
#'                    small bins (< ~200 genes) where BH is too conservative.
#' @param use_enrichr Logical. If TRUE (default), also query the EnrichR API for
#'                    each bin and produce a consensus result.  Requires the
#'                    \code{enrichR} package and internet access.
#' @param enrichr_dbs Character vector of EnrichR database names to query.
#'                    Default covers mouse and human KEGG + Reactome.
#'                    Run \code{enrichR::listEnrichrDbs()} to see all options.
#'
#' @return A list with:
#'   \describe{
#'     \item{$bin_table}{data.frame — confident genes + phase_bin column.}
#'     \item{$ora_results}{Named list (one entry per bin) of combined ORA tables
#'       with a \code{source} column: "phyper", "enrichr", or "consensus".}
#'     \item{$phase_gs}{Named list of phase-restricted gene sets ready for
#'       \code{auc_score_cells()}.  Built from consensus hits first, then
#'       phyper-only hits.  Names follow "ZT<a>-<b>__<PATHWAY>".}
#'   }
#' @export
phase_bin_analysis <- function(
    T1,
    conf_mask,
    genesets,
    bin_width       = 1,
    n_top           = 5L,
    min_overlap     = 3L,
    min_bin_genes   = 5L,
    p_thresh        = 0.05,
    use_padj        = FALSE,   # FALSE = filter on raw p-value (recommended for small bins)
                               # TRUE  = filter on BH-adjusted p-value (stricter, needs many genes)
    use_enrichr     = TRUE,    # also query EnrichR API for consensus validation
    enrichr_dbs     = c(       # EnrichR databases to query (edit as needed)
      "KEGG_2026",
      "Reactome_Pathways_2024",
      "GO_Biological_Process_2025"
    ),
    exclude_patterns = NULL    # optional regex vector to remove irrelevant pathways
                               # applied to both phyper and EnrichR hits by pathway name
                               # e.g. c("SARS|DENGUE|CARDIAC|ADIPOGEN|PARKINSON|AXON_GUID")
) {
  T1_conf   <- T1[conf_mask, , drop = FALSE]
  all_genes <- T1$Genes          # background = all genes tested by run_timescape
  n_bg      <- length(all_genes)

  # ── Check EnrichR availability + validate database names ────────────────────
  # enrichR stores its API base URL in an R option (enrichR.base.address) that
  # is set inside .onAttach() — the hook that runs when you call library() or
  # require().  requireNamespace() only loads the namespace and intentionally
  # skips .onAttach(), so the option stays NULL and every enrichR:: call either
  # crashes ("invalid 'pattern' argument" in setEnrichrSite) or builds a
  # malformed URL ("Could not resolve hostname [datasetStatistics]").
  # Fix: attach the package with require() so .onAttach() fires normally.
  enrichr_ok        <- FALSE
  valid_enrichr_dbs <- character(0)

  if (use_enrichr) {
    if (!requireNamespace("enrichR", quietly = TRUE)) {
      message("  [EnrichR] Package 'enrichR' not installed — skipping.",
              "\n  Install with:  install.packages('enrichR')")
    } else {
      # requireNamespace() only loads the namespace — it skips .onAttach(),
      # which is where enrichR sets enrichR.base.address and related options.
      # Attaching the package (require/library) fires .onAttach() and makes
      # all enrichR:: calls work exactly as they do interactively.
      suppressPackageStartupMessages(
        require("enrichR", quietly = TRUE, character.only = TRUE)
      )

      avail_dbs <- tryCatch(
        enrichR::listEnrichrDbs()$libraryName,
        error = function(e) {
          msg <- tryCatch(conditionMessage(e), error = function(.) "(unknown error)")
          message("  [EnrichR] Could not reach Enrichr API: ", msg,
                  "\n  Falling back to phyper-only.",
                  "\n  Check internet access / firewall settings.")
          NULL
        }
      )

      if (!is.null(avail_dbs)) {
        valid_enrichr_dbs <- enrichr_dbs[enrichr_dbs %in% avail_dbs]
        missing_dbs       <- setdiff(enrichr_dbs, avail_dbs)

        if (length(missing_dbs) > 0) {
          # Suggest close matches using partial string matching
          close <- lapply(missing_dbs, function(db) {
            hits <- avail_dbs[grepl(sub("_[0-9]+$", "", db), avail_dbs,
                                    ignore.case = TRUE)]
            if (length(hits) > 0) paste(head(hits, 3), collapse = " / ") else "—"
          })
          for (i in seq_along(missing_dbs))
            message(sprintf("  [EnrichR] DB not found: '%s'  (similar: %s)",
                            missing_dbs[i], close[[i]]))
          message("  Run enrichR::listEnrichrDbs() for the full list.")
        }

        if (length(valid_enrichr_dbs) > 0) {
          enrichr_ok <- TRUE
          message(sprintf("  [EnrichR] Connected — querying %d database(s): %s",
                          length(valid_enrichr_dbs),
                          paste(valid_enrichr_dbs, collapse = ", ")))
        } else {
          message("  [EnrichR] None of the requested databases matched — ",
                  "falling back to phyper-only.")
        }
      }
    }
  }

  # ── 1. Bin confident genes by acrophase ─────────────────────────────────────
  breaks     <- seq(0, 24, by = bin_width)
  bin_labels <- sprintf("ZT%02.0f-%02.0f", head(breaks, -1), tail(breaks, -1))

  T1_conf$phase_bin <- cut(
    T1_conf$Acrophase_24,
    breaks         = breaks,
    labels         = bin_labels,
    include.lowest = TRUE,
    right          = FALSE
  )

  bin_counts <- table(T1_conf$phase_bin)
  cat(sprintf("\nPhase-bin analysis  (bin_width = %g hr, %d bins)\n",
              bin_width, length(breaks) - 1L))
  cat(sprintf("  Confident genes: %d  |  Background (tested): %d\n",
              nrow(T1_conf), n_bg))
  if (enrichr_ok)
    cat(sprintf("  ORA methods: phyper (custom bg) + EnrichR (%d dbs)\n",
                length(valid_enrichr_dbs)))
  else
    cat("  ORA methods: phyper (custom background only)\n")
  cat("  Genes per bin:\n")
  bc_show <- bin_counts[bin_counts > 0]
  for (i in seq_along(bc_show))
    cat(sprintf("    %s : %d genes\n", names(bc_show)[i], bc_show[i]))

  # ── 2. ORA per bin ──────────────────────────────────────────────────────────
  ora_results <- list()
  phase_gs    <- list()
  filter_lbl  <- if (use_padj) "BH-adjusted p" else "raw p"

  active_bins <- names(bin_counts)[bin_counts >= min_bin_genes]
  if (length(active_bins) == 0)
    stop("No phase bin has >= ", min_bin_genes,
         " genes. Lower min_bin_genes or increase bin_width.")

  for (bin in active_bins) {
    bin_genes <- T1_conf$Genes[as.character(T1_conf$phase_bin) == bin]
    n_bin     <- length(bin_genes)
    cat(sprintf("\n  [%s]  %d genes", bin, n_bin))
    flush.console()

    # ── 2a. Hypergeometric ORA (phyper, custom background) ────────────────────
    cat(sprintf(" — phyper vs %d gene sets...", length(genesets)))
    flush.console()

    phyper_rows <- lapply(names(genesets), function(pw) {
      pw_genes <- genesets[[pw]]
      ov_genes <- intersect(bin_genes, pw_genes)
      overlap  <- length(ov_genes)
      if (overlap < min_overlap) return(NULL)

      # P(X >= overlap) under H0: random draw of n_bin genes from background
      p <- phyper(overlap - 1L,
                  length(pw_genes),
                  n_bg - length(pw_genes),
                  n_bin,
                  lower.tail = FALSE)

      data.frame(
        Pathway        = pw,
        Overlap        = overlap,
        Pathway_size   = length(pw_genes),
        Bin_genes      = n_bin,
        Overlap_genes  = paste(ov_genes, collapse = ";"),
        pvalue         = p,
        stringsAsFactors = FALSE
      )
    })

    phyper_df <- do.call(rbind, Filter(Negate(is.null), phyper_rows))
    if (!is.null(phyper_df) && nrow(phyper_df) > 0) {
      phyper_df$pvalue_adj <- p.adjust(phyper_df$pvalue, method = "BH")
      phyper_df <- phyper_df[order(phyper_df$pvalue), ]
      filter_col_ph <- if (use_padj) phyper_df$pvalue_adj else phyper_df$pvalue
      phyper_sig     <- phyper_df[filter_col_ph < p_thresh, ]
    } else {
      phyper_sig <- data.frame()
    }
    cat(sprintf("  %d sig.", nrow(phyper_sig)))
    flush.console()

    # ── 2b. EnrichR ORA (API, genome-wide background) ─────────────────────────
    enrichr_sig <- data.frame()
    if (enrichr_ok && length(bin_genes) > 0) {
      cat("  EnrichR...")
      flush.console()
      tryCatch({
        enr_raw <- enrichR::enrichr(bin_genes, valid_enrichr_dbs)
        enr_rows <- lapply(names(enr_raw), function(db) {
          df <- enr_raw[[db]]
          if (is.null(df) || nrow(df) == 0) return(NULL)
          # Keep significant hits (use Adjusted.P.value from EnrichR directly)
          sig_df <- df[df$Adjusted.P.value < p_thresh, , drop = FALSE]
          if (nrow(sig_df) == 0) return(NULL)
          data.frame(
            Pathway_enrichr = sig_df$Term,
            DB              = db,
            Overlap_enrichr = sig_df$Overlap,
            pvalue_enrichr  = sig_df$P.value,
            pvalue_adj_enrichr = sig_df$Adjusted.P.value,
            stringsAsFactors = FALSE
          )
        })
        enrichr_raw_df <- do.call(rbind, Filter(Negate(is.null), enr_rows))
        if (!is.null(enrichr_raw_df) && nrow(enrichr_raw_df) > 0) {
          enrichr_raw_df <- enrichr_raw_df[order(enrichr_raw_df$pvalue_enrichr), ]
          enrichr_sig    <- enrichr_raw_df
          cat(sprintf(" %d sig.", nrow(enrichr_sig)))
        } else {
          cat(" 0 sig.")
        }
      }, error = function(e) {
        cat(sprintf(" [error: %s]", e$message))
      })
      flush.console()
    }

    # ── 2c. Consensus: match phyper and EnrichR hits by normalized name ────────
    if (nrow(phyper_sig) == 0) {
      cat(sprintf("  → 0 phyper hits, skipping bin."))
      next
    }

    # ── 2c. Optional exclusion filter ─────────────────────────────────────────
    # Remove pathways matching any exclude_patterns regex from BOTH result sets.
    # Useful for cell-type irrelevant hits (viral, cardiac, adipocyte, neuronal)
    # that appear due to shared gene families rather than true biology.
    if (!is.null(exclude_patterns) && length(exclude_patterns) > 0) {
      excl_rgx <- paste(exclude_patterns, collapse = "|")
      if (nrow(phyper_sig) > 0) {
        keep_ph      <- !grepl(excl_rgx, phyper_sig$Pathway, ignore.case = TRUE)
        n_excl_ph    <- sum(!keep_ph)
        phyper_sig   <- phyper_sig[keep_ph, , drop = FALSE]
        if (n_excl_ph > 0)
          cat(sprintf("  [excl] removed %d phyper hits matching exclude_patterns.", n_excl_ph))
      }
      if (nrow(enrichr_sig) > 0) {
        keep_enr     <- !grepl(excl_rgx, enrichr_sig$Pathway_enrichr, ignore.case = TRUE)
        n_excl_enr   <- sum(!keep_enr)
        enrichr_sig  <- enrichr_sig[keep_enr, , drop = FALSE]
        if (n_excl_enr > 0)
          cat(sprintf("  [excl] removed %d EnrichR hits matching exclude_patterns.", n_excl_enr))
      }
    }

    if (nrow(phyper_sig) == 0) {
      cat(sprintf("  → 0 phyper hits after exclusion filter, skipping bin."))
      next
    }

    # Add normalized name columns for matching
    phyper_sig$norm_name <- .norm_pw_name(phyper_sig$Pathway)

    if (!enrichr_ok) {
      # EnrichR was not run (use_enrichr = FALSE or connection failed) —
      # label everything phyper_only, no consensus possible.
      phyper_sig$source <- "phyper_only"
      enr_only_df       <- data.frame()

    } else if (nrow(enrichr_sig) > 0) {
      enrichr_sig$norm_name <- .norm_pw_name(enrichr_sig$Pathway_enrichr)

      # Tag each phyper hit as "consensus" or "phyper_only".
      # NOTE: consensus requires both tests to hit the SAME pathway name after
      # normalization.  Zero consensus is normal when phyper uses one database
      # (e.g. KEGG from msigdbr, mouse symbols) and EnrichR hits a different
      # database (e.g. Reactome or GO:BP) — those names will never match.
      # To maximise consensus, either:
      #   a) Pull all three databases into genesets (pull_genesets Reactome + GO:BP)
      #      so phyper tests the same universe as EnrichR, OR
      #   b) Query EnrichR with only the same collection as genesets (e.g. KEGG only).
      phyper_sig$source <- ifelse(
        phyper_sig$norm_name %in% enrichr_sig$norm_name,
        "consensus", "phyper_only"
      )

      # EnrichR-only hits (significant in EnrichR, not tested or not sig in phyper)
      enr_only_names <- enrichr_sig$norm_name[
        !enrichr_sig$norm_name %in% phyper_sig$norm_name
      ]
      enr_only_df <- enrichr_sig[enrichr_sig$norm_name %in% enr_only_names, , drop = FALSE]

      n_consensus   <- sum(phyper_sig$source == "consensus")
      n_phyper_only <- sum(phyper_sig$source == "phyper_only")
      n_enr_only    <- nrow(enr_only_df)

      cat(sprintf(
        "\n    Consensus: %d | phyper-only: %d | EnrichR-only: %d",
        n_consensus, n_phyper_only, n_enr_only))

    } else {
      # EnrichR ran but returned 0 significant hits for this bin
      phyper_sig$source <- "phyper_only"
      enr_only_df       <- data.frame()
      cat("  EnrichR: 0 sig.")
    }

    # ── 2d. Build the combined ORA table for this bin ──────────────────────────
    # Consensus and phyper-only both carry gene-level information (overlap genes)
    # because phyper is computed against the local genesets.
    # EnrichR-only rows are appended without gene lists (EnrichR doesn't return them).
    combined_rows <- phyper_sig[, c("Pathway","Overlap","Pathway_size","Bin_genes",
                                     "Overlap_genes","pvalue","pvalue_adj","source"),
                                 drop = FALSE]

    if (nrow(enr_only_df) > 0) {
      enr_append <- data.frame(
        Pathway       = enr_only_df$Pathway_enrichr,
        Overlap       = NA_integer_,
        Pathway_size  = NA_integer_,
        Bin_genes     = n_bin,
        Overlap_genes = NA_character_,
        pvalue        = enr_only_df$pvalue_enrichr,
        pvalue_adj    = enr_only_df$pvalue_adj_enrichr,
        source        = "enrichr_only",
        stringsAsFactors = FALSE
      )
      combined_rows <- rbind(combined_rows, enr_append)
    }

    # Sort: consensus first, then by raw p-value
    source_order <- match(combined_rows$source,
                          c("consensus","phyper_only","enrichr_only"))
    combined_rows <- combined_rows[order(source_order, combined_rows$pvalue), ]

    # Keep n_top rows, prioritising consensus
    top_df        <- head(combined_rows, n_top)
    ora_results[[bin]] <- top_df

    n_conf_hit <- sum(top_df$source == "consensus", na.rm = TRUE)
    cat(sprintf("  → Keeping top %d (%d consensus).\n",
                nrow(top_df), n_conf_hit))

    # ── 2e. Build phase-restricted gene sets (Option B: intersect only) ────────
    # Only genes that co-peak in THIS bin AND belong to the pathway are kept.
    # This guarantees all AUCell set members are co-phased → no cancellation.
    for (i in seq_len(nrow(top_df))) {
      pw <- top_df$Pathway[i]
      # phyper rows carry gene lists; enrichr-only rows do not
      if (!is.na(top_df$Overlap_genes[i]) && nchar(top_df$Overlap_genes[i]) > 0) {
        ov_g <- strsplit(top_df$Overlap_genes[i], ";")[[1]]
      } else if (!is.null(genesets[[pw]])) {
        # Fallback: intersect local geneset with bin genes (e.g. enrichr-only hit
        # whose pathway happens to be present under a different name in genesets)
        ov_g <- intersect(genesets[[pw]], bin_genes)
      } else {
        next  # EnrichR-only with no local gene list — can't build a phase set
      }
      if (length(ov_g) < 2L) next
      gs_name             <- paste0(bin, "__", pw)
      phase_gs[[gs_name]] <- ov_g
    }
  }

  cat(sprintf("\n\nSummary:\n"))
  cat(sprintf("  Bins with ≥1 enriched pathway : %d\n", length(ora_results)))
  cat(sprintf("  Phase-restricted gene sets     : %d\n", length(phase_gs)))
  if (length(phase_gs) > 0) {
    sz <- lengths(phase_gs)
    cat(sprintf("  Gene set sizes  min=%d  median=%g  max=%d\n",
                min(sz), median(sz), max(sz)))
    if (enrichr_ok) {
      n_cons_sets <- sum(sapply(names(ora_results), function(b) {
        sum(ora_results[[b]]$source == "consensus", na.rm = TRUE)
      }))
      cat(sprintf("  Consensus hits (phyper ∩ EnrichR): %d\n", n_cons_sets))
    }
  }

  invisible(list(
    bin_table   = T1_conf,
    ora_results = ora_results,
    phase_gs    = phase_gs
  ))
}


# ── 3. AUCell scoring ──────────────────────────────────────────────────────────

#' Score cells for pathway activity using AUCell
#'
#' Builds a ranking of genes per cell from the expression matrix, then
#' computes the Area Under the Curve (AUC) for each gene set.  Returns a
#' matrix of AUC scores (pathways × cells).
#'
#' @param obj         Seurat or SingleCellExperiment object.
#' @param genesets    Named list of gene-symbol vectors (from \code{pull_genesets()}).
#' @param use_norm    Logical. Use normalised slot (TRUE, default) or raw counts.
#' @param auc_max_rank Fraction of genes to use for AUC calculation (default 0.05 = top 5\%).
#' @param n_cores     Number of cores for AUCell (default 1).
#'
#' @return A pathways × cells matrix of AUC scores (values 0–1).
#' @export
auc_score_cells <- function(
    obj,
    genesets,
    use_norm     = TRUE,
    auc_max_rank = 0.05,
    n_cores      = 1L,
    min_gs_size  = 3L    # lower than default 5 to support phase-restricted sets
) {
  for (pkg in c("AUCell", "BiocParallel")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Install ", pkg, " first:\n",
           "BiocManager::install(c('AUCell', 'BiocParallel'))")
  }

  message("Building per-cell gene rankings for AUCell ...")
  expr_mat <- .get_matrix(obj, use_normalized = use_norm)
  ngenes   <- nrow(expr_mat)

  # Convert AUC max rank from fraction to absolute gene count
  auc_rank <- max(1L, as.integer(round(auc_max_rank * ngenes)))
  message(sprintf("  Matrix: %d genes × %d cells | AUC top-%d genes (%.1f%%)",
                  ngenes, ncol(expr_mat), auc_rank, auc_max_rank * 100))

  # Only keep gene sets whose genes overlap with the expression matrix
  gs_filt  <- lapply(genesets, intersect, rownames(expr_mat))
  sz       <- lengths(gs_filt)
  gs_filt  <- gs_filt[sz >= min_gs_size]
  message(sprintf("  %d / %d gene sets have ≥%d genes in the expression matrix",
                  length(gs_filt), length(genesets), min_gs_size))
  if (length(gs_filt) == 0)
    stop("No gene sets overlap the expression matrix. Check species / gene symbols.")

  # Build rankings (cells ranked by expression, one column per cell).
  # nCores / BPPARAM were removed from both AUCell functions in recent versions;
  # passing either causes an error.  Run single-threaded (AUCell is fast enough
  # for most datasets at n_cores = 1 via its own internal chunking).
  rankings <- AUCell::AUCell_buildRankings(
    expr_mat,
    plotStats = FALSE
  )

  # Score gene sets in chunks so progress is visible.
  # AUCell_calcAUC no longer exposes a progress bar in recent versions,
  # so we split gs_filt into chunks and report after each one.
  n_sets     <- length(gs_filt)
  chunk_size <- 1L   # 1 = progress per pathway; increase (e.g. 10) for fewer updates
  chunks     <- split(seq_len(n_sets),
                      ceiling(seq_len(n_sets) / chunk_size))
  n_chunks   <- length(chunks)

  cat(sprintf("Computing AUC scores (%d pathways)...\n", n_sets))

  auc_rows <- vector("list", n_chunks)
  for (i in seq_along(chunks)) {
    idx      <- chunks[[i]]
    done     <- idx[length(idx)]
    pct      <- round(100 * done / n_sets)
    cat(sprintf("\r  [%-20s] %3d%%  (%d / %d pathways)",
                strrep("=", round(20 * done / n_sets)), pct, done, n_sets))
    flush.console()

    auc_chunk   <- AUCell::AUCell_calcAUC(
      gs_filt[idx],
      rankings,
      aucMaxRank = auc_rank,
      verbose    = FALSE
    )
    auc_rows[[i]] <- AUCell::getAUC(auc_chunk)
  }
  cat(sprintf("\r  [%s] 100%%  (%d / %d pathways) — done.\n",
              strrep("=", 20), n_sets, n_sets))

  # Return as plain numeric matrix: pathways × cells
  do.call(rbind, auc_rows)
}


# ── 3. Pathway-level cosinor ──────────────────────────────────────────────────

#' Run cosinor on AUCell pathway scores for one cell type
#'
#' Treats each pathway's AUC score as if it were a "gene" and fits the same
#' cosinor model used by \code{estimate_phaseR}.  Returns a results table
#' comparable to \code{run_timescape}'s T1, plus per-ZT mean scores for heatmap.
#'
#' @param auc_mat      Pathways × cells matrix from \code{auc_score_cells()}.
#' @param meta         data.frame of cell-level metadata (.get_meta(obj)).
#' @param celltype_col Metadata column name for cell type.
#' @param zt_col       Metadata column name for ZT strings.
#' @param tmeta        data.frame from \code{build_tmeta()}.
#' @param target_ct    Character. Cell type to analyse (must match celltype_col values).
#' @param period12     Logical. 12-hr (TRUE) or 24-hr (FALSE, default).
#' @param custom_zt    Optional character vector to restrict ZT time points.
#'
#' @return A list with:
#'   \describe{
#'     \item{stats}{data.frame — one row per pathway; same columns as run_timescape T1.}
#'     \item{zt_means}{data.frame — per-ZT mean AUC; pathways × ZT.}
#'   }
#' @export
pathway_cosinor <- function(
    auc_mat,
    meta,
    celltype_col,
    zt_col,
    tmeta,
    target_ct,
    period12  = FALSE,
    custom_zt = NULL
) {
  # Cell mask
  ct_mask  <- as.character(meta[[celltype_col]]) == target_ct
  zt_strs  <- as.character(meta[[zt_col]])[ct_mask]
  if (!is.null(custom_zt)) {
    keep    <- zt_strs %in% custom_zt
    ct_mask[ct_mask][!keep] <- FALSE
    zt_strs <- zt_strs[keep]
  }

  n_cells <- sum(ct_mask)
  if (n_cells < 10L)
    stop(sprintf("Only %d cells for '%s' — too few for pathway cosinor.", n_cells, target_ct))

  # Subset AUC matrix to this cell type
  auc_ct   <- auc_mat[, ct_mask, drop = FALSE]
  zt_num   <- tmeta$ZT_times[match(zt_strs, tmeta$zt_str)]
  zt_present <- sort(unique(zt_num[!is.na(zt_num)]))
  nzts       <- length(zt_present)

  if (nzts < 4L)
    stop(sprintf("Only %d ZT time points for '%s' — need ≥4.", nzts, target_ct))

  zt_labels  <- tmeta$zt_str[match(zt_present, tmeta$ZT_times)]
  npaths     <- nrow(auc_ct)
  path_names <- rownames(auc_ct)
  message(sprintf("Pathway cosinor: %d pathways × %d cells (%d ZT points) for '%s'",
                  npaths, n_cells, nzts, target_ct))

  period <- if (period12) 12 else 24

  # Fit cosinor for each pathway
  res <- lapply(seq_len(npaths), function(i) {
    Xg_zts <- lapply(zt_present, function(z) {
      idx <- which(zt_num == z)
      if (length(idx) == 0L) return(numeric(0))
      as.numeric(auc_ct[i, idx])
    })
    estimate_phaseR(Xg_zts, zt_present, period12, "Ftest")
  })

  amp_v   <- vapply(res, `[[`, numeric(1), "amp")
  mesor_v <- vapply(res, `[[`, numeric(1), "mesor")
  acro_v  <- vapply(res, `[[`, numeric(1), "acrophase")
  pval_v  <- vapply(res, `[[`, numeric(1), "p_value")
  rho_v   <- vapply(res, `[[`, numeric(1), "rho")
  pmac_v  <- vapply(res, `[[`, numeric(1), "p_value_macro")

  acro24  <- acro_v %% 24
  padj_v  <- stats::p.adjust(pval_v,  method = "BH")
  padj_mac<- stats::p.adjust(pmac_v,  method = "BH")

  # Per-ZT mean AUC
  R0_mat <- matrix(NA_real_, nrow = npaths, ncol = nzts)
  for (j in seq_along(zt_present)) {
    idx <- which(zt_num == zt_present[j])
    if (length(idx) > 0L)
      R0_mat[, j] <- rowMeans(auc_ct[, idx, drop = FALSE], na.rm = TRUE)
  }

  stats_df <- data.frame(
    Pathway         = path_names,
    Amp             = amp_v,
    Abs_Amp         = abs(amp_v),
    Mesor           = mesor_v,
    Acrophase       = acro_v,
    Acrophase_24    = acro24,
    Period          = rep(period, npaths),
    pvalue          = pval_v,
    pvalue_adj      = padj_v,
    Sine_corr       = rho_v,
    pvalue_corr     = pmac_v,
    pvalue_adj_corr = padj_mac,
    stringsAsFactors = FALSE
  )

  # Remove failed fits
  valid     <- !is.na(stats_df$pvalue) & !is.na(stats_df$pvalue_corr)
  stats_df  <- stats_df[valid, ]
  R0_mat    <- R0_mat[valid, , drop = FALSE]

  # Sort by adjusted correlation p-value
  ord       <- order(stats_df$pvalue_adj_corr, stats_df$pvalue_adj,
                     stats_df$Acrophase_24, -stats_df$Abs_Amp)
  stats_df  <- stats_df[ord, ]
  R0_mat    <- R0_mat[ord, , drop = FALSE]
  rownames(stats_df) <- NULL

  zt_df     <- data.frame(Pathway = stats_df$Pathway, stringsAsFactors = FALSE)
  for (j in seq_along(zt_labels)) zt_df[[zt_labels[j]]] <- R0_mat[, j]

  n_conf <- sum(stats_df$pvalue < 0.05 & stats_df$pvalue_corr < 0.05)
  message(sprintf("  Pathways tested: %d | Confident (p<0.05 both): %d",
                  nrow(stats_df), n_conf))

  list(stats = stats_df, zt_means = zt_df)
}


# ── 4. Write pathway results to Excel ─────────────────────────────────────────

#' Write pathway cosinor results to a formatted Excel workbook
#'
#' Creates a two-sheet workbook: "All_Pathways" (all tested) and
#' "Confident_Pathways" (p<0.05 for both F-test and Pearson correlation).
#' Columns are auto-width; confident rows are highlighted in the All sheet.
#'
#' @param results   List returned by \code{pathway_cosinor()} — must have \code{$stats}.
#' @param outpath   Full path for the output .xlsx file.
#' @param celltype  Cell type label used in the sheet title (default "").
#'
#' @return Invisibly: the \code{outpath}.
#' @export
write_pathway_results <- function(results, outpath, celltype = "") {
  for (pkg in c("openxlsx")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Install openxlsx first:  install.packages('openxlsx')")
  }

  stats_df <- results$stats
  conf_mask <- stats_df$pvalue < 0.05 & stats_df$pvalue_corr < 0.05
  conf_df   <- stats_df[conf_mask, ]

  wb <- openxlsx::createWorkbook()

  # ── Sheet 1: All pathways ──────────────────────────────────────────────────
  openxlsx::addWorksheet(wb, "All_Pathways")
  openxlsx::writeData(wb, "All_Pathways", stats_df, rowNames = FALSE)
  openxlsx::setColWidths(wb, "All_Pathways", cols = seq_len(ncol(stats_df)),
                         widths = "auto")

  # Highlight confident rows in yellow
  conf_rows <- which(conf_mask) + 1L  # +1 for header
  if (length(conf_rows) > 0L) {
    hi_style <- openxlsx::createStyle(fgFill = "#FFF2CC")
    openxlsx::addStyle(wb, "All_Pathways", style = hi_style,
                       rows = conf_rows, cols = seq_len(ncol(stats_df)),
                       gridExpand = TRUE, stack = TRUE)
  }

  # Header style
  hdr_style <- openxlsx::createStyle(
    fontColour = "#FFFFFF", fgFill = "#2F4F8F",
    halign = "center", textDecoration = "bold"
  )
  openxlsx::addStyle(wb, "All_Pathways", style = hdr_style,
                     rows = 1L, cols = seq_len(ncol(stats_df)),
                     gridExpand = TRUE, stack = TRUE)

  # ── Sheet 2: Confident pathways ──────────────────────────────────────────────
  openxlsx::addWorksheet(wb, "Confident_Pathways")
  openxlsx::writeData(wb, "Confident_Pathways", conf_df, rowNames = FALSE)
  openxlsx::setColWidths(wb, "Confident_Pathways", cols = seq_len(ncol(conf_df)),
                         widths = "auto")
  openxlsx::addStyle(wb, "Confident_Pathways", style = hdr_style,
                     rows = 1L, cols = seq_len(ncol(conf_df)),
                     gridExpand = TRUE, stack = TRUE)

  openxlsx::saveWorkbook(wb, outpath, overwrite = TRUE)
  message("Pathway results saved → ", outpath)
  invisible(outpath)
}


# ── 5. Gene regulatory network at a single ZT ─────────────────────────────────

#' Build a correlation-based GRN for a given set of genes at one ZT time point
#'
#' Computes pairwise Pearson correlations among cells at the specified ZT
#' using a supplied gene × cells dense sub-matrix.  Edges are retained
#' when |r| ≥ cor_thresh and the correlation p-value < p_thresh.
#'
#' @param expr_mat   genes × cells dense (or sparse) matrix for the full cell type.
#' @param gene_set   Character vector of gene names to include.
#' @param zt_num_vec Numeric ZT time vector aligned to expr_mat columns.
#' @param zt_target  Numeric ZT time to subset (e.g. 6 for ZT06).
#' @param cor_thresh Minimum |Pearson r| for an edge (default 0.2).
#' @param p_thresh   Maximum Pearson p-value for an edge (default 0.05).
#'
#' @return An \code{igraph} graph object (undirected, weighted by r).
#' @keywords internal
.build_zt_network <- function(
    expr_mat,
    gene_set,
    zt_num_vec,
    zt_target,
    cor_thresh = 0.2,
    p_thresh   = 0.05
) {
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Install igraph:  install.packages('igraph')")

  # Subset cells at this ZT
  cell_idx  <- which(zt_num_vec == zt_target)
  if (length(cell_idx) < 3L) {
    message(sprintf("  ZT%02d: only %d cells — skipping network", zt_target, length(cell_idx)))
    return(igraph::make_empty_graph(n = length(gene_set), directed = FALSE))
  }

  # Subset genes present in the matrix
  g_present <- intersect(gene_set, rownames(expr_mat))
  if (length(g_present) < 2L)
    return(igraph::make_empty_graph(n = 0L, directed = FALSE))

  sub_mat <- as.matrix(expr_mat[g_present, cell_idx, drop = FALSE])  # dense, small

  # Compute pairwise Pearson correlations
  ng  <- nrow(sub_mat)
  nc  <- ncol(sub_mat)
  df  <- nc - 2L

  # Correlation matrix
  rmat <- tryCatch(stats::cor(t(sub_mat), method = "pearson"), error = function(e) NULL)
  if (is.null(rmat)) return(igraph::make_empty_graph(n = ng, directed = FALSE))

  # t-statistic → p-value for each pair
  # t = r * sqrt((n-2) / (1-r^2)), df = n-2
  tmat <- rmat * sqrt(df / pmax(1 - rmat^2, 1e-10))
  pmat <- 2 * stats::pt(-abs(tmat), df = df)
  diag(rmat) <- 0; diag(pmat) <- 1

  # Edge list: keep upper triangle with |r| >= thresh AND p < thresh
  keep_mat <- (abs(rmat) >= cor_thresh) & (pmat < p_thresh)
  keep_mat[lower.tri(keep_mat, diag = TRUE)] <- FALSE

  edges <- which(keep_mat, arr.ind = TRUE)
  if (nrow(edges) == 0L) {
    g <- igraph::make_empty_graph(n = ng, directed = FALSE)
    igraph::V(g)$name <- g_present
    return(g)
  }

  edge_df <- data.frame(
    from   = g_present[edges[, 1]],
    to     = g_present[edges[, 2]],
    weight = rmat[edges],
    stringsAsFactors = FALSE
  )

  g <- igraph::graph_from_data_frame(edge_df, directed = FALSE,
                                      vertices = data.frame(name = g_present))
  g
}


# ── 6. GRN time-series plot ──────────────────────────────────────────────────

#' Plot GRN across ZT time points as a multi-panel figure
#'
#' Builds one network per ZT time point using \code{.build_zt_network()}, then
#' arranges them in a grid.  A fixed spring layout is computed from the pooled
#' (all-ZT) network so node positions are consistent across panels.
#'
#' @param obj          Seurat or SingleCellExperiment object.
#' @param circ_genes   Character vector of confident circadian genes (from T1).
#' @param pathway_genes Character vector of pathway member genes to overlay.
#' @param meta         data.frame of metadata (.get_meta(obj)).
#' @param celltype_col Metadata column for cell type.
#' @param zt_col       Metadata column for ZT strings.
#' @param tmeta        data.frame from build_tmeta().
#' @param target_ct    Cell type to analyse.
#' @param cor_thresh   Edge |r| threshold (default 0.2).
#' @param p_thresh     Edge p-value threshold (default 0.05).
#' @param use_norm     Use normalised slot (default TRUE).
#' @param outfile      Path for output PNG (default NULL = not saved, only returned).
#' @param ncol         Number of columns in the grid (default: auto).
#' @param node_size    Node size in ggraph (default 4).
#'
#' @return A ggplot2 object (the arranged grid).
#' @export
plot_grn_timeseries <- function(
    obj,
    circ_genes,
    pathway_genes,
    meta,
    celltype_col,
    zt_col,
    tmeta,
    target_ct,
    cor_thresh      = 0.2,
    p_thresh        = 0.05,
    use_norm        = TRUE,
    outfile         = NULL,
    ncol            = NULL,
    node_size       = 4,
    label_size      = 4,       # gene label font size (ggplot pts)
    edge_width_max  = 2.5,     # maximum edge line width (scaled by |r|)
    zt_title_size   = 18,      # ZT panel title font size
    zt_title_hjust  = 0.5      # 0 = left | 0.5 = centred | 1 = right
) {
  for (pkg in c("igraph", "ggraph", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Install ", pkg, ":  install.packages('", pkg, "')")
  }

  message("Building GRN time series for '", target_ct, "' ...")

  # All genes to include in the GRN (union of circadian + pathway genes)
  all_genes <- unique(c(circ_genes, pathway_genes))

  # Expression matrix for this cell type only (sparse → stays sparse until block)
  ct_mask  <- as.character(meta[[celltype_col]]) == target_ct
  zt_strs  <- as.character(meta[[zt_col]])[ct_mask]
  zt_num   <- tmeta$ZT_times[match(zt_strs, tmeta$zt_str)]

  expr_full <- .get_matrix(obj, use_normalized = use_norm)
  g_present <- intersect(all_genes, rownames(expr_full))
  if (length(g_present) < 2L)
    stop("Fewer than 2 target genes found in expression matrix.")

  # Dense sub-matrix for this cell type and gene set only (small — manageable)
  expr_ct   <- as.matrix(expr_full[g_present, ct_mask, drop = FALSE])
  rm(expr_full); gc(verbose = FALSE)

  # ZT time points
  zt_points <- sort(unique(zt_num[!is.na(zt_num)]))
  nzts      <- length(zt_points)
  zt_labels <- tmeta$zt_str[match(zt_points, tmeta$ZT_times)]

  message(sprintf("  %d genes × %d ZT time points", length(g_present), nzts))

  # ── Build per-ZT networks ───────────────────────────────────────────────────
  nets <- lapply(zt_points, function(zt) {
    .build_zt_network(expr_ct, g_present, zt_num, zt, cor_thresh, p_thresh)
  })
  names(nets) <- zt_labels

  # ── Fixed layout from pooled network (all ZT combined) ─────────────────────
  # Merge all edges to get a stable layout
  all_edges <- do.call(rbind, lapply(seq_along(nets), function(i) {
    g  <- nets[[i]]
    el <- igraph::as_data_frame(g, what = "edges")
    if (nrow(el) == 0L) return(NULL)
    el
  }))

  if (is.null(all_edges) || nrow(all_edges) == 0L) {
    # No edges anywhere — make empty grid
    message("  No edges found in any ZT network with current thresholds.")
    empty_p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
                        label = sprintf("No edges (|r|≥%.2f, p<%.2f)", cor_thresh, p_thresh),
                        size = 5) +
      ggplot2::theme_void()
    return(empty_p)
  }

  pooled_g <- igraph::graph_from_data_frame(all_edges, directed = FALSE,
                                             vertices = data.frame(name = g_present))
  set.seed(42)
  # layout_with_fr requires positive weights; use abs(r) so anti-correlations
  # still contribute to node spacing by magnitude without causing an error.
  layout_coords <- igraph::layout_with_fr(
    pooled_g,
    weights = abs(igraph::E(pooled_g)$weight)
  )
  rownames(layout_coords) <- igraph::V(pooled_g)$name

  # ── Classify nodes ─────────────────────────────────────────────────────────
  node_class <- setNames(
    ifelse(g_present %in% circ_genes & g_present %in% pathway_genes, "both",
    ifelse(g_present %in% circ_genes,  "circadian",
    ifelse(g_present %in% pathway_genes, "pathway", "other"))),
    g_present
  )
  class_colors <- c(
    circadian = "#E05C2F",   # orange-red
    pathway   = "#2F6DB5",   # blue
    both      = "#8B1A8B",   # purple
    other     = "#AAAAAA"    # grey
  )

  # ── One ggraph panel per ZT ─────────────────────────────────────────────────
  panels <- lapply(seq_along(nets), function(i) {
    g   <- nets[[i]]
    zt_lbl <- zt_labels[i]

    # Align node order to fixed layout
    vnames <- igraph::V(g)$name
    lmat   <- layout_coords[vnames, , drop = FALSE]

    # Positive / negative correlation → edge color
    el <- igraph::as_data_frame(g, what = "edges")

    p <- ggraph::ggraph(g, layout = lmat) +
      {if (nrow(el) > 0)
        ggraph::geom_edge_link(
          ggplot2::aes(color = weight,
                       alpha = abs(weight),
                       width = abs(weight)),
          show.legend = FALSE
        )
      else
        ggplot2::geom_blank()
      } +
      ggraph::scale_edge_color_gradient2(low = "#2166AC", mid = "#DDDDDD",
                                          high = "#B2182B", midpoint = 0,
                                          guide = "none") +
      ggraph::scale_edge_width(range = c(0.4, edge_width_max), guide = "none") +
      ggraph::scale_edge_alpha(range = c(0.35, 0.95),          guide = "none") +
      ggraph::geom_node_point(
        ggplot2::aes(color = node_class[vnames]), size = node_size
      ) +
      ggraph::geom_node_text(
        ggplot2::aes(label = vnames), size = label_size, repel = TRUE,
        max.overlaps = 30, fontface = "bold"
      ) +
      ggplot2::scale_color_manual(values = class_colors, name = "Gene type") +
      ggplot2::labs(title = zt_lbl,
                    subtitle = sprintf("%d edges", igraph::ecount(g))) +
      ggraph::theme_graph(base_family = "") +
      ggplot2::theme(
        legend.position = "none",
        plot.title      = ggplot2::element_text(size  = zt_title_size,
                                                face  = "bold",
                                                hjust = zt_title_hjust),
        plot.subtitle   = ggplot2::element_text(size = 11, color = "grey50",
                                                hjust = zt_title_hjust)
      )

    p
  })

  # ── Shared legend: nodes (left) + edge colour (right) ───────────────────────
  if (!requireNamespace("cowplot", quietly = TRUE))
    stop("Install cowplot:  install.packages('cowplot')")

  # Node legend
  node_leg_data <- data.frame(
    x = 1:4, y = 1,
    Gene_type = factor(c("circadian","pathway","both","other"),
                       levels = c("circadian","pathway","both","other"))
  )
  node_leg_p <- ggplot2::ggplot(node_leg_data,
                                 ggplot2::aes(x, y, color = Gene_type)) +
    ggplot2::geom_point(size = 5) +
    ggplot2::scale_color_manual(
      values = class_colors,
      labels = c("Circadian gene", "Pathway gene", "Both", "Other"),
      name   = "Nodes"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position  = "right",
      legend.title     = ggplot2::element_text(size = 12, face = "bold"),
      legend.text      = ggplot2::element_text(size = 12),
      legend.direction = "horizontal"
    )
  node_legend <- cowplot::get_legend(node_leg_p)

  # Edge colour legend — built from two coloured segments
  edge_leg_data <- data.frame(
    x    = c(0, 0),  xend = c(1, 1),
    y    = c(2, 1),  yend = c(2, 1),
    corr = factor(c("Positive (r > 0)", "Negative (r < 0)"),
                  levels = c("Positive (r > 0)", "Negative (r < 0)"))
  )
  edge_leg_p <- ggplot2::ggplot(edge_leg_data,
    ggplot2::aes(x=x, y=y, xend=xend, yend=yend, color=corr, linewidth=corr)) +
    ggplot2::geom_segment(lineend = "round") +
    ggplot2::scale_color_manual(
      values = c("Positive (r > 0)" = "#B2182B",
                 "Negative (r < 0)" = "#2166AC"),
      name   = "Edges"
    ) +
    ggplot2::scale_linewidth_manual(
      values = c("Positive (r > 0)" = 2.5,
                 "Negative (r < 0)" = 2.5),
      guide  = "none"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position  = "right",
      legend.title     = ggplot2::element_text(size = 12, face = "bold"),
      legend.text      = ggplot2::element_text(size = 12),
      legend.direction = "horizontal"
    )
  edge_legend <- cowplot::get_legend(edge_leg_p)

  # Combine both legends into one horizontal strip
  combined_legend <- cowplot::plot_grid(
    node_legend, edge_legend,
    ncol = 2, rel_widths = c(1, 0.6)
  )

  # ── Arrange into grid ───────────────────────────────────────────────────────
  nc   <- if (!is.null(ncol)) ncol else min(4L, nzts)
  grid <- cowplot::plot_grid(plotlist = panels, ncol = nc)
  final_plot <- cowplot::plot_grid(
    grid, combined_legend,
    ncol = 1, rel_heights = c(1, 0.1)
  )

  # ── Optional save ───────────────────────────────────────────────────────────
  if (!is.null(outfile)) {
    ggplot2::ggsave(
      filename = outfile,
      plot     = final_plot,
      width    = nc * 6,
      height   = ceiling(nzts / nc) * 6 + 1.2,
      dpi      = 180,
      bg       = "white"
    )
    message("GRN time-series saved → ", outfile)
  }

  invisible(final_plot)
}


# ── 7. Pathway circadian plot (violin / dots + cosine fit) ────────────────────

#' Plot a single pathway's AUC scores across ZT time points with cosine fit
#'
#' Mirrors \code{plot_gene_single()} but for pathway-level AUCell scores.
#' Shows per-cell AUC values (jittered dots or violin) per ZT time point,
#' per-ZT mean as orange points, and the fitted cosine curve in blue.
#'
#' @param auc_mat       Pathways × cells matrix from \code{auc_score_cells()}.
#' @param path_results  List returned by \code{pathway_cosinor()} with \code{$stats} and \code{$zt_means}.
#' @param meta          data.frame of cell-level metadata (\code{.get_meta(obj)}).
#' @param celltype_col  Metadata column for cell type.
#' @param zt_col        Metadata column for ZT strings.
#' @param tmeta         data.frame from \code{build_tmeta()}.
#' @param target_ct     Cell type name (must match celltype_col values).
#' @param target_pathway Exact pathway name (rowname in \code{auc_mat}).
#' @param period12      Logical. 12-hr (TRUE) or 24-hr (FALSE, default).
#' @param use_violin    Logical. Violin density (TRUE) or jittered dots (FALSE, default).
#' @param custom_zt     Optional character vector to filter ZT time points.
#'
#' @return A ggplot2 object.
#' @export
plot_pathway_single <- function(
    auc_mat,
    path_results,
    meta,
    celltype_col,
    zt_col,
    tmeta,
    target_ct,
    target_pathway,
    period12    = FALSE,
    use_violin  = FALSE,
    custom_zt   = NULL
) {
  # ── Find pathway in results ─────────────────────────────────────────────────
  stats_df <- path_results$stats
  pi_row   <- which(stats_df$Pathway == target_pathway)
  if (length(pi_row) == 0)
    stop(sprintf("Pathway '%s' not found in path_results$stats.", target_pathway))
  pi_row <- pi_row[1]

  amp_g   <- stats_df$Amp[pi_row]
  acro_g  <- stats_df$Acrophase[pi_row]
  mesor_g <- stats_df$Mesor[pi_row]
  acro_24 <- stats_df$Acrophase_24[pi_row]
  pval    <- stats_df$pvalue[pi_row]
  pcorr   <- stats_df$pvalue_corr[pi_row]
  period_val <- if (period12) 12 else 24

  # ── Subset cells ────────────────────────────────────────────────────────────
  ct_mask  <- as.character(meta[[celltype_col]]) == target_ct
  zt_strs  <- as.character(meta[[zt_col]])[ct_mask]
  if (!is.null(custom_zt)) {
    keep    <- zt_strs %in% custom_zt
    ct_mask[ct_mask][!keep] <- FALSE
    zt_strs <- zt_strs[keep]
  }
  zt_num   <- tmeta$ZT_times[match(zt_strs, tmeta$zt_str)]
  zt_present <- sort(unique(zt_num[!is.na(zt_num)]))
  n_cells  <- sum(ct_mask)

  # ── AUC scores for this pathway ─────────────────────────────────────────────
  pi_mat   <- which(rownames(auc_mat) == target_pathway)
  if (length(pi_mat) == 0)
    stop(sprintf("Pathway '%s' not found in auc_mat rownames.", target_pathway))
  auc_vals <- as.numeric(auc_mat[pi_mat[1], ct_mask])

  sc_data  <- data.frame(ZT = zt_num, AUC = auc_vals)
  sc_data  <- sc_data[!is.na(sc_data$ZT), ]

  # ── Per-ZT means from path_results$zt_means ─────────────────────────────────
  zt_df    <- path_results$zt_means
  pm_row   <- which(zt_df$Pathway == target_pathway)
  zt_cols  <- setdiff(colnames(zt_df), "Pathway")
  zt_numeric_m <- tmeta$ZT_times[match(zt_cols, tmeta$zt_str)]
  if (length(pm_row) > 0) {
    Rzts <- as.numeric(zt_df[pm_row[1], zt_cols])
  } else {
    Rzts <- tapply(sc_data$AUC, sc_data$ZT, mean, na.rm = TRUE)[as.character(zt_numeric_m)]
  }

  # ── Cosine curve ─────────────────────────────────────────────────────────────
  t0   <- min(zt_present); tf <- max(zt_present)
  tval <- seq(t0, tf, length.out = 300)
  fval <- amp_g * cos(2*pi*(tval - acro_g) / period_val) + mesor_g

  cos_df  <- data.frame(ZT = tval, AUC = fval)
  mean_df <- data.frame(ZT = zt_numeric_m, AUC = Rzts)

  # ── Colours (matches gene plot style) ───────────────────────────────────────
  col_cos  <- "#2278B5"; col_mean <- "#D65A0C"
  col_dots <- "#4D7FD1"; col_acro <- "#7A0000"; col_viol <- "#99CCEE"

  # ── Theme ────────────────────────────────────────────────────────────────────
  wt <- ggplot2::theme_bw() +
    ggplot2::theme(
      panel.background  = ggplot2::element_rect(fill="white", colour=NA),
      plot.background   = ggplot2::element_rect(fill="white", colour=NA),
      panel.grid.major  = ggplot2::element_line(colour="#d4d4d4"),
      panel.grid.minor  = ggplot2::element_blank(),
      axis.text         = ggplot2::element_text(colour="black", size=11),
      axis.title        = ggplot2::element_text(colour="black", size=11),
      plot.title        = ggplot2::element_text(colour="black", size=9, face="bold"),
      legend.background = ggplot2::element_rect(fill="white", colour="#b0b0b0"),
      legend.text       = ggplot2::element_text(colour="black", size=10),
      legend.title      = ggplot2::element_blank(),
      legend.position   = "top", legend.direction = "horizontal"
    )

  # Pretty pathway label: replace underscores, drop GOBP_ prefix
  path_label <- gsub("^GOBP_|^KEGG_|^REACTOME_|^HALLMARK_", "", target_pathway)
  path_label <- gsub("_", " ", path_label)
  path_label <- tools::toTitleCase(tolower(path_label))

  p <- ggplot2::ggplot() + wt +
    ggplot2::scale_x_continuous(breaks = zt_present) +
    ggplot2::labs(
      x     = "Zeitgeber Time (hrs)",
      y     = "AUCell score",
      title = sprintf("%s\nF-test p=%.3g | Corr p=%.3g | Phase ZT%.1f | N=%d",
                      path_label, pval, pcorr, acro_24, n_cells)
    )

  if (use_violin && nrow(sc_data) > 0) {
    sc_data$ZT_f <- factor(sc_data$ZT)
    y_max <- max(stats::quantile(sc_data$AUC, 0.99, na.rm=TRUE) * 1.15,
                 max(Rzts, na.rm=TRUE) * 3, na.rm=TRUE)
    p <- p +
      ggplot2::geom_violin(
        data = sc_data,
        ggplot2::aes(x=ZT, y=AUC, group=ZT_f, fill="Violin"),
        colour=col_cos, alpha=0.50, linewidth=0.4
      ) +
      ggplot2::scale_fill_manual(values=c("Violin"=col_viol))
  } else if (nrow(sc_data) > 0) {
    set.seed(42)
    p <- p +
      ggplot2::geom_jitter(
        data=sc_data,
        ggplot2::aes(x=ZT, y=AUC, colour="Cells"),
        width=0.35, alpha=0.3, size=0.8, inherit.aes=FALSE
      ) +
      ggplot2::scale_colour_manual(
        values=c("Cells"=col_dots, "Cosine fit"=col_cos, "Mean per ZT"=col_mean),
        breaks=c("Cosine fit","Mean per ZT")
      )
  }

  p <- p +
    ggplot2::geom_line(data=cos_df,
                       ggplot2::aes(x=ZT, y=AUC, colour="Cosine fit"),
                       linewidth=2.2, inherit.aes=FALSE) +
    ggplot2::geom_line(data=mean_df,
                       ggplot2::aes(x=ZT, y=AUC, colour="Mean per ZT"),
                       linewidth=1.8, inherit.aes=FALSE) +
    ggplot2::geom_point(data=mean_df,
                        ggplot2::aes(x=ZT, y=AUC, colour="Mean per ZT"),
                        size=3.5, inherit.aes=FALSE) +
    ggplot2::geom_vline(xintercept=acro_24, linetype="dashed",
                        colour=col_acro, linewidth=1.0) +
    ggplot2::annotate("text", x=acro_24, y=-Inf,
                      label=sprintf("ZT%.1f", acro_24),
                      vjust=-0.4, hjust=-0.1, colour=col_acro, size=3.5)

  if (!use_violin) {
    p <- p +
      ggplot2::scale_colour_manual(
        values = c("Cells"=col_dots, "Cosine fit"=col_cos, "Mean per ZT"=col_mean),
        breaks = c("Cosine fit","Mean per ZT")
      )
  } else {
    p <- p +
      ggplot2::scale_colour_manual(
        values = c("Cosine fit"=col_cos, "Mean per ZT"=col_mean),
        breaks = c("Cosine fit","Mean per ZT")
      )
  }

  p
}


#' Batch-save pathway circadian plots for top N confident pathways
#'
#' @param auc_mat      Pathways × cells matrix from \code{auc_score_cells()}.
#' @param path_results List from \code{pathway_cosinor()}.
#' @param meta         data.frame of cell-level metadata.
#' @param celltype_col Metadata column for cell type.
#' @param zt_col       Metadata column for ZT strings.
#' @param tmeta        data.frame from \code{build_tmeta()}.
#' @param target_ct    Cell type name.
#' @param n_top        Number of top confident pathways to plot (default 20).
#' @param period12     Logical. 12-hr (TRUE) or 24-hr (FALSE, default).
#' @param use_violin   Logical. Violin (TRUE) or jittered dots (FALSE, default).
#' @param outdir       Directory to save PNGs.
#' @param custom_zt    Optional character vector to filter ZT time points.
#'
#' @return Invisibly: character vector of saved file paths.
#' @export
save_batch_pathway_plots <- function(
    auc_mat,
    path_results,
    meta,
    celltype_col,
    zt_col,
    tmeta,
    target_ct,
    n_top      = 20L,
    period12   = FALSE,
    use_violin = FALSE,
    outdir,
    custom_zt  = NULL
) {
  stats_df <- path_results$stats
  conf_mask <- stats_df$pvalue < 0.05 & stats_df$pvalue_corr < 0.05
  conf_df   <- stats_df[conf_mask, ]

  if (nrow(conf_df) == 0) {
    message("No confident pathways to plot.")
    return(invisible(character(0)))
  }

  top_paths <- head(conf_df$Pathway, n_top)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  ct_safe   <- gsub("[^[:alnum:]_]", "_", trimws(target_ct))
  saved     <- character(length(top_paths))

  for (i in seq_along(top_paths)) {
    pw  <- top_paths[i]
    # Trim to 60 chars for filename only — plot title keeps the full name
    pw_safe <- substr(gsub("[^[:alnum:]_]", "_", pw), 1, 60)
    out_png <- file.path(outdir, sprintf("%s_%02d_%s.png", ct_safe, i, pw_safe))

    p <- tryCatch(
      plot_pathway_single(
        auc_mat, path_results, meta, celltype_col, zt_col, tmeta,
        target_ct, pw, period12, use_violin, custom_zt
      ),
      error = function(e) { message("  Skip ", pw, ": ", e$message); NULL }
    )
    if (!is.null(p)) {
      ggplot2::ggsave(out_png, plot=p, width=8, height=5, dpi=150, bg="white")
      saved[i] <- out_png
    }
  }

  message(sprintf("Saved %d pathway plots → %s", sum(nzchar(saved)), outdir))
  invisible(saved[nzchar(saved)])
}
