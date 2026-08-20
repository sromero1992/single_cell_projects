# =============================================================================
# de.R — differential expression / variability / enrichment helpers
# =============================================================================
# Ported from Script 07. Every value that used to be a script-level global is
# now a function argument whose DEFAULT equals the value Script 07 used, so
# behaviour is unchanged. Script 07 passes its config at the call sites that it
# actively tunes (contrast/sample columns, Enrichr DBs, ADD_ENSEMBL).
#
# NSE note: these use bare column names inside dplyr/ggplot2 (as the original
# script did); functionally correct, and the package imports ggplot2 + the pipe.
# =============================================================================

# -----------------------------------------------------------------------------
# Comparison column + Ensembl finalisation
# -----------------------------------------------------------------------------

#' Ensure the comparison column exists (use, or auto-build a compound "A_B")
#' @param obj A Seurat object.
#' @param col Column name; if a compound "ColA_ColB" whose parts both exist, it
#'   is built by pasting them with "_".
#' @return `obj`, guaranteed to carry `col` in its metadata.
#' @export
resolve_comparison_column <- function(obj, col) {
  if (col %in% colnames(obj@meta.data)) {
    message(paste0("  Comparison column '", col, "' found in metadata."))
    return(obj)
  }
  parts <- strsplit(col, "_", fixed = TRUE)[[1]]
  if (length(parts) >= 2 && all(parts %in% colnames(obj@meta.data))) {
    obj@meta.data[[col]] <- apply(obj@meta.data[, parts, drop = FALSE], 1,
                                  paste, collapse = "_")
    message(paste0("  Comparison column '", col, "' AUTO-BUILT from: ",
                   paste(parts, collapse = " + "), "."))
    return(obj)
  }
  stop(paste0("[ERROR] Comparison column '", col, "' is neither present in the ",
              "metadata nor buildable from existing columns.\n",
              "  Name it 'ColA_ColB' where both are existing columns.\n",
              "  Available metadata columns: ",
              paste(colnames(obj@meta.data), collapse = ", ")), call. = FALSE)
}

#' Map gene symbols to Ensembl IDs (mouse/human), NA-safe
#' @param symbols Character vector of gene symbols.
#' @param species "mouse" (default) or "human".
#' @return Character vector of Ensembl IDs, same length/order (NA where unmapped).
#' @export
symbol_to_ensembl <- function(symbols, species = "mouse") {
  symbols <- as.character(symbols)
  out <- rep(NA_character_, length(symbols))
  org_db <- if (tolower(species) %in% c("human", "hs")) "org.Hs.eg.db" else "org.Mm.eg.db"
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace(org_db, quietly = TRUE)) return(out)
  valid <- which(!is.na(symbols) & symbols != "")
  if (length(valid) == 0) return(out)
  keys <- unique(symbols[valid])
  db   <- getExportedValue(org_db, org_db)
  ens <- tryCatch(
    AnnotationDbi::mapIds(db, keys = keys, keytype = "SYMBOL",
                          column = "ENSEMBL", multiVals = "first"),
    error = function(e) stats::setNames(rep(NA_character_, length(keys)), keys))
  miss <- names(ens)[is.na(ens)]
  if (length(miss) > 0) {
    alias <- tryCatch(
      AnnotationDbi::mapIds(db, keys = miss, keytype = "ALIAS",
                            column = "ENSEMBL", multiVals = "first"),
      error = function(e) stats::setNames(rep(NA_character_, length(miss)), miss))
    ens[miss] <- alias
  }
  out[valid] <- unname(ens[symbols[valid]])
  out
}

#' Put 'gene' first and 'Ensembl' second in an export sheet
#' @param df A data.frame (a DE/DV/overlap sheet).
#' @param add_ensembl Add the Ensembl column (default TRUE).
#' @param species Passed to [symbol_to_ensembl()].
#' @return The reordered data.frame.
#' @export
finalize_sheet <- function(df, add_ensembl = TRUE, species = "mouse") {
  df <- as.data.frame(df)
  if (!"gene" %in% colnames(df)) return(df)          # e.g. placeholder sheets
  if (add_ensembl && !"Ensembl" %in% colnames(df))
    df$Ensembl <- symbol_to_ensembl(df$gene, species = species)
  front <- intersect(c("gene", "Ensembl"), colnames(df))
  df[, c(front, setdiff(colnames(df), front)), drop = FALSE]
}

#' Apply [finalize_sheet()] to every sheet in a list
#' @param sheet_list Named list of data.frames.
#' @param add_ensembl,species Passed to [finalize_sheet()].
#' @return The list, each sheet finalised.
#' @export
finalize_all <- function(sheet_list, add_ensembl = TRUE, species = "mouse")
  lapply(sheet_list, finalize_sheet, add_ensembl = add_ensembl, species = species)

# -----------------------------------------------------------------------------
# DE (MAST) and DV (SplineDV)
# -----------------------------------------------------------------------------

#' MAST differential expression between two groups within a Seurat object
#' @param so A (cell-type-subset) Seurat object.
#' @param ident.1,ident.2 The two levels of `group_by` to compare.
#' @param group_by Metadata column holding the comparison groups.
#' @param latent_vars Latent variable(s) for MAST (e.g. the sample column).
#' @param logfc_thresh,min_pct FindMarkers thresholds (defaults 0.5 / 0.10).
#' @param min_cells Skip if either group has fewer cells (default 20).
#' @return A tidy DE data.frame (gene, avg_log2FC, p_val_adj, confidence, ...) or NULL.
#' @export
run_mast <- function(so, ident.1, ident.2, group_by, latent_vars,
                     logfc_thresh = 0.5, min_pct = 0.10, min_cells = 20) {
  Seurat::Idents(so) <- so[[group_by, drop = TRUE]]
  n1 <- sum(Seurat::Idents(so) == ident.1)
  n2 <- sum(Seurat::Idents(so) == ident.2)
  if (n1 < min_cells || n2 < min_cells) {
    message(paste("  [SKIP] Too few cells:", ident.1, "=", n1, "|", ident.2, "=", n2,
                  "(min =", min_cells, ")"))
    return(NULL)
  }
  markers <- tryCatch(
    Seurat::FindMarkers(object = so, ident.1 = ident.1, ident.2 = ident.2,
                        test.use = "MAST", latent.vars = latent_vars,
                        logfc.threshold = logfc_thresh, min.pct = min_pct,
                        verbose = FALSE),
    error = function(e) { message(paste("  [ERROR] MAST failed:", e$message)); NULL })
  if (is.null(markers) || nrow(markers) == 0) return(NULL)

  markers %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(
      contrast   = paste0(ident.1, "_vs_", ident.2),
      direction  = ifelse(avg_log2FC > 0, paste0("up_in_", ident.1),
                          paste0("up_in_", ident.2)),
      confidence = dplyr::case_when(
        abs(avg_log2FC) >= 1.0 ~ "high",
        abs(avg_log2FC) >= 0.5 ~ "moderate",
        TRUE                   ~ "low")) %>%
    dplyr::arrange(p_val_adj, dplyr::desc(abs(avg_log2FC)))
}

#' SplineDV differential variability between two groups within one cell type
#' @param seurat_obj A Seurat object.
#' @param cell_type Cell type to test.
#' @param celltype_col,condition_col Metadata columns for cell type / condition.
#' @param group1,group2 The two condition levels to compare.
#' @param min_cells Skip if either group has fewer cells (default 20).
#' @return A tidy DV data.frame (gene, pval, padj, ...) or NULL.
#' @export
run_splinedv <- function(seurat_obj, cell_type, celltype_col, condition_col,
                         group1, group2, min_cells = 20) {
  if (!requireNamespace("SplineDV", quietly = TRUE))
    stop("[ERROR] SplineDV not installed. Run 00_rlibs_installation.R to install it.")

  cells_ct <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[celltype_col]] == cell_type]
  obj_ct   <- subset(seurat_obj, cells = cells_ct)
  obj_ct   <- subset(obj_ct, cells = rownames(obj_ct@meta.data)[
    obj_ct@meta.data[[condition_col]] %in% c(group1, group2)])

  n1 <- sum(obj_ct@meta.data[[condition_col]] == group1)
  n2 <- sum(obj_ct@meta.data[[condition_col]] == group2)
  if (n1 < min_cells || n2 < min_cells) {
    message(paste("  [SKIP SplineDV]", cell_type, "| insufficient cells:",
                  group1, "=", n1, "|", group2, "=", n2, "(min =", min_cells, ")"))
    return(NULL)
  }
  message(paste("  [SplineDV]", cell_type, "|", group1, "(n=", n1, ") vs",
                group2, "(n=", n2, ")"))

  cells_g1 <- rownames(obj_ct@meta.data)[obj_ct@meta.data[[condition_col]] == group1]
  cells_g2 <- rownames(obj_ct@meta.data)[obj_ct@meta.data[[condition_col]] == group2]
  counts_g1 <- Seurat::GetAssayData(obj_ct, assay = "RNA", layer = "counts")[, cells_g1, drop = FALSE]
  counts_g2 <- Seurat::GetAssayData(obj_ct, assay = "RNA", layer = "counts")[, cells_g2, drop = FALSE]

  dv_result <- tryCatch(SplineDV::splineDV(counts_g1, counts_g2),
    error = function(e) { message("    [ERROR] SplineDV failed: ", e$message); NULL })
  if (is.null(dv_result) || nrow(dv_result) == 0) return(NULL)
  dv_result <- as.data.frame(dv_result)

  if ("genes" %in% colnames(dv_result)) {
    dv_result <- dplyr::rename(dv_result, gene = genes)
  } else if ("Gene" %in% colnames(dv_result)) {
    dv_result <- dplyr::rename(dv_result, gene = Gene)
  } else if (!"gene" %in% colnames(dv_result)) {
    dv_result$gene <- rownames(dv_result)
  }
  if ("Pval" %in% colnames(dv_result)) {
    dv_result <- dplyr::rename(dv_result, pval = Pval)
  } else if ("PValue" %in% colnames(dv_result)) {
    dv_result <- dplyr::rename(dv_result, pval = PValue)
  } else if ("pvalue" %in% colnames(dv_result)) {
    dv_result <- dplyr::rename(dv_result, pval = pvalue)
  } else if ("p.value" %in% colnames(dv_result)) {
    dv_result <- dplyr::rename(dv_result, pval = p.value)
  }
  if (!"pval" %in% colnames(dv_result)) {
    message("    [WARN] No p-value column found in SplineDV output. Columns: ",
            paste(colnames(dv_result), collapse = ", "))
    dv_result$pval <- NA_real_
  }
  dv_result$padj <- if (any(!is.na(dv_result$pval)))
    stats::p.adjust(dv_result$pval, method = "BH") else NA_real_

  dv_result$cell_type  <- cell_type
  dv_result$comparison <- paste0(group1, "_vs_", group2)
  dv_result$n_cells_g1 <- n1
  dv_result$n_cells_g2 <- n2
  if (any(!is.na(dv_result$pval)))
    dv_result <- dplyr::arrange(dv_result, pval)
  dv_result
}

# -----------------------------------------------------------------------------
# Enrichr ORA
# -----------------------------------------------------------------------------

#' Enrichr over-representation analysis, flattened across databases
#' @param gene_list Character vector of gene symbols.
#' @param label Label for messages.
#' @param dbs Enrichr databases to query.
#' @param top_n Keep the top N terms per database (default 100).
#' @param padj Adjusted-p cutoff (default 0.05).
#' @param min_genes Skip if fewer than this many genes (default 5).
#' @param sleep Seconds to pause (Enrichr rate limit; default 2).
#' @param enabled If FALSE, returns NULL immediately.
#' @return A flat data.frame with a 'database' column, or NULL.
#' @export
run_enrichr_ora <- function(gene_list, label,
                            dbs = c("GO_Biological_Process_2026",
                                    "GO_Molecular_Function_2026",
                                    "KEGG_2026",
                                    "WikiPathways_2024_Mouse"),
                            top_n = 100, padj = 0.05, min_genes = 5,
                            sleep = 2, enabled = TRUE) {
  if (!enabled || length(gene_list) == 0) return(NULL)
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    message("  [SKIP Enrichr] enrichR not installed."); return(NULL)
  }
  gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
  if (length(gene_list) < min_genes) {
    message(paste("  [SKIP Enrichr] Too few genes for:", label)); return(NULL)
  }
  Sys.sleep(sleep)
  tryCatch({
    enrich_res <- enrichR::enrichr(gene_list, dbs)
    flat_rows  <- list()
    for (db in names(enrich_res)) {
      d <- enrich_res[[db]]
      if (!is.null(d) && nrow(d) > 0) {
        d <- d %>%
          dplyr::filter(Adjusted.P.value < padj) %>%
          dplyr::arrange(Adjusted.P.value) %>%
          utils::head(top_n) %>%
          dplyr::mutate(database = db, .before = 1)
        flat_rows[[db]] <- d
      }
    }
    if (length(flat_rows) > 0) {
      message(paste("  Enrichr done:", label))
      return(dplyr::bind_rows(flat_rows))
    }
    NULL
  }, error = function(e) { message(paste("  [ERROR Enrichr]", e$message)); NULL })
}

# -----------------------------------------------------------------------------
# Pathway heatmap / barplot / gene-hub heatmap
# -----------------------------------------------------------------------------

#' Cross-cell-type pathway heatmap (top GOBP terms, gene-count labels)
#' @param enr_list Named list of Enrichr flat data.frames.
#' @param title,out_file Plot title and output PNG path.
#' @param gobp_db GO-BP database name to filter to.
#' @param top_n,min_gene_overlap,redundancy_ratio,min_cell_types Selection knobs.
#' @param dpi Output resolution.
#' @return Invisibly, the ggplot (or NULL if skipped).
#' @export
make_pathway_heatmap <- function(enr_list, title, out_file,
                                 gobp_db = "GO_Biological_Process_2026",
                                 top_n = 20, min_gene_overlap = 3,
                                 redundancy_ratio = 0.5, min_cell_types = 2,
                                 dpi = 300) {
  if (length(enr_list) == 0) {
    message(paste("  [SKIP heatmap] No enrichr data for:", title)); return(invisible(NULL))
  }
  all_enr <- dplyr::bind_rows(enr_list)
  if (nrow(all_enr) == 0) return(invisible(NULL))
  all_enr <- all_enr %>% dplyr::filter(database == gobp_db)
  if (nrow(all_enr) == 0) {
    message(paste("  [SKIP heatmap] No GOBP results for:", title)); return(invisible(NULL))
  }
  all_enr <- all_enr %>%
    dplyr::mutate(gene_count = as.integer(sub("/.*", "", Overlap)),
                  neg_log10_padj = -log10(pmax(Adjusted.P.value, 1e-10))) %>%
    dplyr::filter(Adjusted.P.value < 0.05, gene_count >= min_gene_overlap)
  if (nrow(all_enr) == 0) {
    message(paste("  [SKIP heatmap] No significant pathways after filtering for:", title))
    return(invisible(NULL))
  }
  if ("direction" %in% colnames(all_enr)) {
    all_enr <- all_enr %>% dplyr::mutate(group_key = paste0(cell_type, "\n", contrast, "\n", direction))
  } else {
    all_enr <- all_enr %>% dplyr::mutate(group_key = paste0(cell_type, "\n", contrast))
  }

  deduplicate_pathways <- function(df) {
    if (nrow(df) <= 1) return(df)
    df <- df %>% dplyr::mutate(score = gene_count * neg_log10_padj) %>%
      dplyr::arrange(dplyr::desc(score))
    gene_lists <- lapply(df$Genes, function(g) unlist(strsplit(g, ";\\s*|,\\s*")))
    keep <- rep(TRUE, nrow(df)); n <- nrow(df)
    for (i in seq_len(n)) {
      if (!keep[i]) next
      if (i == n) break
      for (j in seq(i + 1, n)) {
        if (!keep[j]) next
        gi <- gene_lists[[i]]; gj <- gene_lists[[j]]
        shared <- length(intersect(gi, gj))
        pct_i  <- if (length(gi) > 0) shared / length(gi) else 0
        pct_j  <- if (length(gj) > 0) shared / length(gj) else 0
        if (pct_j > redundancy_ratio || pct_i > redundancy_ratio) keep[j] <- FALSE
      }
    }
    df[keep, ] %>% dplyr::select(-score)
  }

  all_enr_dedup <- all_enr %>% dplyr::group_by(group_key) %>%
    dplyr::group_modify(~ deduplicate_pathways(.x) %>% utils::head(top_n)) %>%
    dplyr::ungroup()

  median_gene_count <- stats::median(all_enr_dedup$gene_count, na.rm = TRUE)
  pathway_global <- all_enr_dedup %>% dplyr::group_by(Term) %>%
    dplyr::summarise(n_cell_types = dplyr::n_distinct(cell_type),
                     mean_gene_count = mean(gene_count, na.rm = TRUE),
                     mean_score = mean(gene_count * neg_log10_padj, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::filter(n_cell_types >= min_cell_types | mean_gene_count >= median_gene_count) %>%
    dplyr::arrange(dplyr::desc(n_cell_types), dplyr::desc(mean_score)) %>%
    utils::head(top_n)

  selected_terms <- pathway_global$Term
  if (length(selected_terms) == 0) {
    message(paste("  [SKIP heatmap] No cross-representative pathways for:", title))
    return(invisible(NULL))
  }
  message(paste("  Pathways selected:", length(selected_terms)))

  heatmap_df <- all_enr_dedup %>% dplyr::filter(Term %in% selected_terms) %>%
    dplyr::group_by(Term, group_key) %>%
    dplyr::summarise(neg_log10_padj = max(neg_log10_padj, na.rm = TRUE),
                     gene_count = max(gene_count, na.rm = TRUE), .groups = "drop") %>%
    tidyr::complete(Term, group_key, fill = list(neg_log10_padj = 0, gene_count = 0))

  term_order_wrapped <- stringr::str_wrap(pathway_global$Term, width = 35)
  heatmap_df$Term    <- stringr::str_wrap(heatmap_df$Term, width = 35)
  heatmap_df$Term    <- factor(heatmap_df$Term, levels = rev(term_order_wrapped))
  heatmap_df$group_key <- factor(heatmap_df$group_key)
  heatmap_df$label   <- ifelse(heatmap_df$gene_count > 0, as.character(heatmap_df$gene_count), "")
  n_terms  <- length(unique(heatmap_df$Term))
  n_groups <- length(unique(heatmap_df$group_key))

  p_heat <- ggplot(heatmap_df, aes(x = group_key, y = Term, fill = neg_log10_padj)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = label), size = 4.5, color = "black", fontface = "bold") +
    scale_fill_gradient(low = "white", high = "#d73027", name = "-log10\n(adj.p)", na.value = "grey95") +
    labs(title = title,
         subtitle = paste0("Top ", n_terms, " pathways | Gene count in cells | GOBP only"),
         x = NULL, y = NULL) +
    theme_bw(base_size = 16) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 13),
          panel.grid = element_blank(), plot.margin = margin(10, 10, 10, 20))

  ggsave(out_file, p_heat, width = max(10, n_groups * 1.6 + 4),
         height = max(8, n_terms * 0.7 + 5), dpi = dpi, bg = "white", limitsize = FALSE)
  message(paste("  Heatmap saved:", basename(out_file)))
  invisible(p_heat)
}

#' Top-N enriched GOBP pathways as a barplot (bar = odds ratio, colour = adj p)
#' @param enr_list Named list of Enrichr flat data.frames.
#' @param title,out_file Plot title and output PNG path.
#' @param top_n Number of pathways to show (default 10).
#' @param gobp_db GO-BP database name to filter to.
#' @param dpi Output resolution.
#' @return Invisibly, the ggplot (or NULL if skipped).
#' @export
make_pathway_barplot <- function(enr_list, title, out_file, top_n = 10,
                                 gobp_db = "GO_Biological_Process_2026", dpi = 300) {
  if (length(enr_list) == 0) {
    message(paste("  [SKIP barplot] No enrichr data for:", title)); return(invisible(NULL))
  }
  df <- dplyr::bind_rows(enr_list)
  if ("database" %in% colnames(df)) df <- dplyr::filter(df, database == gobp_db)
  df <- dplyr::filter(df, Adjusted.P.value < 0.05)
  if (nrow(df) == 0) {
    message(paste("  [SKIP barplot] No significant GOBP pathways for:", title)); return(invisible(NULL))
  }
  df$gene_count <- as.integer(sub("/.*", "", df$Overlap))
  xcol <- if ("Odds.Ratio" %in% colnames(df)) "Odds.Ratio" else
          if ("Combined.Score" %in% colnames(df)) "Combined.Score" else "gene_count"
  xlab <- switch(xcol, Odds.Ratio = "Odds Ratio", Combined.Score = "Combined Score", "Gene overlap")

  df <- df %>% dplyr::group_by(Term) %>%
    dplyr::slice_min(Adjusted.P.value, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>% dplyr::arrange(Adjusted.P.value) %>% utils::head(top_n)

  df$Term_wrapped <- stringr::str_wrap(df$Term, width = 45)
  df$Term_wrapped <- factor(df$Term_wrapped, levels = df$Term_wrapped[order(df[[xcol]])])

  p_bar <- ggplot(df, aes(x = .data[[xcol]], y = Term_wrapped, fill = Adjusted.P.value)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = gene_count), hjust = -0.3, size = 4, fontface = "bold") +
    scale_fill_gradient(low = "#B2182B", high = "#4393C3", name = "adj. p") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = title,
         subtitle = paste0("Top ", nrow(df), " GOBP pathways | bar = ", xlab,
                           " | colour = adj. p | number = gene count"),
         x = xlab, y = NULL) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0.5),
          panel.grid.major.y = element_blank())

  ggsave(out_file, p_bar, width = 11, height = max(3, nrow(df) * 0.5 + 2),
         dpi = dpi, bg = "white")
  message(paste("  Barplot saved:", basename(out_file)))
  invisible(p_bar)
}

#' Gene x pathway hub heatmap (hub genes across top GOBP pathways)
#' @param enr_list Named list of Enrichr flat data.frames.
#' @param de_list Named list of DE result data.frames.
#' @param title,out_file Plot title and output PNG path.
#' @param gobp_db GO-BP database name.
#' @param min_gene_overlap,top_n_pathways,top_n_genes,min_pathways,min_logfc,wrap_width Knobs.
#' @param dpi Output resolution.
#' @return Invisibly, the ggplot (or NULL if skipped).
#' @export
make_gene_hub_heatmap <- function(enr_list, de_list, title, out_file,
                                  gobp_db = "GO_Biological_Process_2026",
                                  min_gene_overlap = 3, top_n_pathways = 20,
                                  top_n_genes = 30, min_pathways = 2,
                                  min_logfc = 0.5, wrap_width = 30, dpi = 300) {
  if (length(enr_list) == 0 || length(de_list) == 0) {
    message(paste("  [SKIP gene hub] Insufficient data for:", title)); return(invisible(NULL))
  }
  all_enr <- dplyr::bind_rows(enr_list) %>%
    dplyr::filter(database == gobp_db) %>%
    dplyr::mutate(gene_count = as.integer(sub("/.*", "", Overlap)),
                  neg_log10_padj = -log10(pmax(Adjusted.P.value, 1e-10))) %>%
    dplyr::filter(Adjusted.P.value < 0.05, gene_count >= min_gene_overlap)
  if (nrow(all_enr) == 0) {
    message(paste("  [SKIP gene hub] No significant pathways for:", title)); return(invisible(NULL))
  }
  top_pathways <- all_enr %>% dplyr::group_by(Term) %>%
    dplyr::summarise(mean_gene_count = mean(gene_count, na.rm = TRUE),
                     mean_score = mean(gene_count * neg_log10_padj, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_score)) %>% utils::head(top_n_pathways)
  selected_pathways <- top_pathways$Term

  all_de <- dplyr::bind_rows(de_list) %>%
    dplyr::filter(p_val_adj <= 0.05, abs(avg_log2FC) >= min_logfc)
  if (nrow(all_de) == 0) {
    message(paste("  [SKIP gene hub] No qualifying DE genes for:", title)); return(invisible(NULL))
  }
  pathway_gene_lists <- all_enr %>% dplyr::filter(Term %in% selected_pathways) %>%
    dplyr::select(Term, Genes) %>% dplyr::distinct() %>% dplyr::rowwise() %>%
    dplyr::mutate(gene_vec = list(trimws(unlist(strsplit(Genes, ";|,"))))) %>%
    dplyr::ungroup()
  de_genes <- unique(all_de$gene)

  membership <- tidyr::expand_grid(gene = de_genes, Term = selected_pathways) %>%
    dplyr::left_join(pathway_gene_lists %>% dplyr::select(Term, gene_vec), by = "Term") %>%
    dplyr::rowwise() %>% dplyr::mutate(member = gene %in% gene_vec[[1]]) %>%
    dplyr::ungroup() %>% dplyr::select(gene, Term, member)

  gene_hub_scores <- membership %>% dplyr::filter(member) %>%
    dplyr::group_by(gene) %>% dplyr::summarise(n_pathways = dplyr::n_distinct(Term), .groups = "drop") %>%
    dplyr::filter(n_pathways >= min_pathways) %>%
    dplyr::left_join(all_de %>% dplyr::group_by(gene) %>%
                       dplyr::summarise(avg_log2FC = mean(avg_log2FC, na.rm = TRUE), .groups = "drop"),
                     by = "gene") %>%
    dplyr::mutate(hub_score = n_pathways * abs(avg_log2FC)) %>%
    dplyr::arrange(dplyr::desc(hub_score)) %>% utils::head(top_n_genes)
  if (nrow(gene_hub_scores) == 0) {
    message(paste("  [SKIP gene hub] No hub genes found (need >=", min_pathways,
                  "pathway memberships) for:", title)); return(invisible(NULL))
  }

  top_genes <- gene_hub_scores$gene
  selected_pathways_final <- selected_pathways[
    selected_pathways %in% (membership %>%
      dplyr::filter(gene %in% top_genes, member) %>% dplyr::pull(Term) %>% unique())]

  plot_df <- membership %>%
    dplyr::filter(gene %in% top_genes, Term %in% selected_pathways_final) %>%
    dplyr::left_join(all_de %>% dplyr::group_by(gene) %>%
                       dplyr::summarise(avg_log2FC = mean(avg_log2FC, na.rm = TRUE), .groups = "drop"),
                     by = "gene") %>%
    dplyr::mutate(fill_val = ifelse(member, avg_log2FC, NA_real_),
                  label = ifelse(member, sprintf("%.1f", avg_log2FC), ""))

  gene_order    <- gene_hub_scores$gene
  pathway_order <- stringr::str_wrap(selected_pathways_final, width = wrap_width)
  plot_df$gene <- factor(plot_df$gene, levels = rev(gene_order))
  plot_df$Term <- factor(stringr::str_wrap(plot_df$Term, width = wrap_width), levels = pathway_order)

  gene_labels <- gene_hub_scores %>%
    dplyr::mutate(label_full = paste0(gene, " (", n_pathways, ")")) %>%
    dplyr::select(gene, label_full)
  levels(plot_df$gene) <- gene_labels$label_full[match(levels(plot_df$gene), gene_labels$gene)]
  n_genes    <- length(unique(plot_df$gene))
  n_pathways <- length(unique(plot_df$Term))

  p_hub <- ggplot(plot_df, aes(x = Term, y = gene)) +
    geom_tile(aes(fill = fill_val), color = "white", linewidth = 0.4) +
    geom_text(aes(label = label), size = 3.8, color = "black", fontface = "bold") +
    scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                         midpoint = 0, name = "avg\nlog2FC", na.value = "grey92") +
    labs(title = title,
         subtitle = paste0("Genes: padj <= 0.05 | logFC >= ", min_logfc,
                           " | Red = up | Blue = down | Number = logFC | GOBP only"),
         x = NULL, y = NULL) +
    theme_bw(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5, size = 17, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 11),
          panel.grid = element_blank(), plot.margin = margin(10, 10, 10, 20))

  ggsave(out_file, p_hub, width = max(10, n_pathways * 1.4 + 4),
         height = max(8, n_genes * 0.5 + 5), dpi = dpi, bg = "white", limitsize = FALSE)
  message(paste("  Gene hub heatmap saved:", basename(out_file),
                "| Hub genes:", n_genes, "| Pathways:", n_pathways))
  invisible(p_hub)
}
