# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 12e: AUCell-SCORE THE Nr4a1 HIERARCHY ON THE UMAP
# =============================================================================
# Takes the TF -> FOCUS -> targets hierarchy defined by 12d and AUCell-scores its
# gene set(s) across the WHOLE object, then overlays the scores on the existing
# (HVG/Harmony) UMAP. This is a program-activity read of the regulatory
# neighborhood: an un-pruned, expression-driven cousin of a SCENIC regulon (no
# cisTarget motif validation), so it is orthogonal to the motif-based regulon and
# reads as "where in transcriptional space this program is coordinately active".
#
# To keep it DIRECTED/interpretable, each piece is scored separately:
#   <FOCUS>_upstream     the TFs that regulate FOCUS
#   <FOCUS>_downstream   the genes FOCUS regulates
#   <FOCUS>_all          upstream + FOCUS + downstream
#   <FOCUS>_down_mod<k>  each downstream co-expression module (from 12d, if present)
#
# SCOPE: by default the gene set (e.g. from the STEM-CELL GRN) is AUCell-scored
# across the whole COLONOCYTE lineage (NR4A1_HIER_SCOPE="Colonocytes") so you can
# follow a stem-origin program DOWNSTREAM through differentiation (stem -> TA ->
# absorptive) and across WT/Polyp/KO -- the relevant readout for an Lgr5-Cre KO of
# Nr4a1. Set NR4A1_HIER_SCOPE="whole" to score every cell, or "label" for just the
# label's own cells. The ranking is built WITHIN the scoped population and every
# plot (UMAP, violins, lineage trace) is drawn on it.
# CAVEAT: the set was inferred partly from the label's cells, so it trivially lights
# up there; the informative test is how it behaves across the rest of the lineage.
#
# INPUT : SCENIC/<label>/network/<FOCUS>_network_nodes[_augmented].csv  (from 12d)
#         <OUTPUT_DIR>/<PROJECT>_with_cell_scores.rds                    (expression + UMAP)
# OUTPUT (per label, under SCENIC/<label>/hierarchy_aucell/):
#   <FOCUS>_hierAUC_umap.png        each program on the RE-EMBEDDED colonocyte UMAP
#                                   (fresh HVG/PCA/Harmony/UMAP of the subset)
#   violin_AUC_<FOCUS>_<set>.png   faceted by colonocyte SUBTYPE, split by the 6
#   barplot_AUC_<FOCUS>_<set>.png  Genotype_sex conditions (Script 09 style + Wilcoxon)
#   <FOCUS>_hierAUC_percell.csv     per-cell AUC for each set + metadata
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
})

# --- Shared portable config (nr4a1 defaults; override via env vars, config.R). ----
# Robust config.R locator: works whether launched with Rscript from any dir,
# source()'d, or run interactively. Tries, in order: NR4A1_CONFIG env var; the
# dir of this script; the current working dir; the default install path.
.NR4A1_CFG <- local({
  cand <- Sys.getenv("NR4A1_CONFIG", "")
  a  <- commandArgs(FALSE); fa <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(fa)) cand <- c(cand, file.path(dirname(normalizePath(fa[1], mustWork = FALSE)), "config.R"))
  of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(of)) cand <- c(cand, file.path(dirname(normalizePath(of, mustWork = FALSE)), "config.R"))
  cand <- c(cand, "config.R",
            "/run/media/external/scdc_projects/2026_nr4a1_ack/r_process/config.R")
  hit <- cand[nzchar(cand) & file.exists(cand)]
  if (!length(hit)) stop("config.R not found. Set NR4A1_CONFIG=/full/path/to/config.R", call. = FALSE)
  hit[1]
})
if (!file.exists(.NR4A1_CFG)) stop("config.R not found at '", .NR4A1_CFG,
  "' - cd to the pipeline directory or set NR4A1_CONFIG.", call. = FALSE)
source(.NR4A1_CFG)

SCENIC_DIR  <- Sys.getenv("NR4A1_SCENIC_DIR", file.path(OUTPUT_DIR, "SCENIC"))
FOCUS_GENE  <- Sys.getenv("NR4A1_FOCUS_GENE", "Nr4a1-cust")  # must match 12d (reads its *_network_*.csv)
RDS_PATH    <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_with_cell_scores.rds"))
UMAP_REDUCTION   <- Sys.getenv("NR4A1_UMAP", "umap_harmony")
CELLTYPE_COLUMN  <- "CellType"
CONDITION_COLUMN <- "Genotype_sex"
# Scoring population for AUCell (ranks are built WITHIN this set, and all plots use
# it): "whole" = every cell (generalization test); "label" = only the label's SCENIC
# cells; ANY OTHER value = a compartment, matched in BROAD_COLUMN (e.g. "Colonocytes")
# or falling back to CELLTYPE_COLUMN -> score + plot across ALL cells of that lineage.
SCORE_SCOPE  <- Sys.getenv("NR4A1_HIER_SCOPE", "Colonocytes")
BROAD_COLUMN <- Sys.getenv("NR4A1_HIER_BROAD_COL", "CellType_broad")
# WHICH labels' GRN gene set to score. Default: the STEM-CELL and TUMOR-EPITHELIUM
# labels -- both are part of the colonocyte lineage; CD8/macrophages/etc. are NOT.
# Override with a regex (NR4A1_HIER_GENESET_PATTERN) or explicit list.
GENESET_PATTERN <- Sys.getenv("NR4A1_HIER_GENESET_PATTERN", "Stem|Tumor")
GENESET_LABELS  <- Filter(nzchar, trimws(strsplit(Sys.getenv("NR4A1_HIER_GENESET_LABELS", ""), ",")[[1]]))
# Lineage/potency trace: plot the program score vs a potency axis (from Script 09)
# and order cell types by mean potency (stem -> differentiated). Set to "" to skip.
POTENCY_COLUMN <- Sys.getenv("NR4A1_POTENCY_COL", "CytoTRACE2_Score")
AUC_MAXRANK_PCT <- as.numeric(Sys.getenv("NR4A1_HIER_AUCMAXRANK", "0.05"))
MIN_SET <- 4          # skip a gene set smaller than this
DPI_SETTING <- 300
LABELS <- NULL

if (!requireNamespace("AUCell", quietly = TRUE))
  stop("12e needs AUCell (already used by Script 09). Install via 00_rlibs_installation.R.",
       call. = FALSE)
have_patch <- requireNamespace("patchwork", quietly = TRUE)
if (!requireNamespace("ggpubr", quietly = TRUE))
  message("[note] ggpubr not installed -> condition violins without Wilcoxon brackets.")

# =============================================================================
# SHARED cell-score violin/barplot plotter -- IDENTICAL to Script 09's
# plot_score_comparison(): faceted by cell type, split by Genotype_sex with the
# defined levels + legend + within-sex Wilcoxon. Reused here so the hierarchy AUC
# plots match the rest of the pipeline exactly.
# =============================================================================
SPLIT_BY_COL    <- CONDITION_COLUMN
selected_groups <- c("WT_Female", "Polyp_Female", "Polyp_NR4a1_KO_Female",
                     "WT_Male",   "Polyp_Male",   "Polyp_NR4a1_KO_Male")
my_comparisons  <- list(
  c("Polyp_NR4a1_KO_Female", "Polyp_Female"),
  c("Polyp_Female",           "WT_Female"),
  c("Polyp_NR4a1_KO_Male",   "Polyp_Male"),
  c("Polyp_Male",            "WT_Male"))

plot_score_comparison <- function(seurat_obj, score_col, split_by = SPLIT_BY_COL,
                                  facet_by = "CellType_broad", groups = selected_groups,
                                  comparisons = my_comparisons, plot_type = "violin",
                                  save_dir = ".") {
  df <- FetchData(seurat_obj, vars = c(score_col, split_by, facet_by)) %>%
    tidyr::drop_na() %>% dplyr::filter(!!sym(split_by) %in% groups)
  if (nrow(df) == 0) return(NULL)
  df[[split_by]] <- factor(df[[split_by]], levels = groups)
  valid_comps <- Filter(function(x) all(x %in% unique(df[[split_by]])), comparisons)
  p <- ggplot(df, aes(x = !!sym(split_by), y = !!sym(score_col), fill = !!sym(split_by),
                      color = !!sym(split_by), group = !!sym(split_by)))
  if (plot_type == "barplot") {
    p <- p +
      geom_jitter(shape = 21, color = "black", stroke = 0.3, width = 0.2, size = 1.2,
                  alpha = 0.7, show.legend = FALSE) +
      stat_summary(fun = mean, geom = "bar", color = "black", linewidth = 0.6, alpha = 1, width = 0.7) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black", linewidth = 0.8)
  } else if (plot_type == "jitter") {
    p <- p +
      geom_jitter(shape = 21, color = "black", stroke = 0.3, width = 0.25, size = 1.2, alpha = 0.6) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, color = "black", linewidth = 0.8) +
      stat_summary(fun = mean, geom = "point", size = 3, color = "black")
  } else {
    p <- p +
      geom_violin(scale = "width", trim = TRUE, alpha = 0.7, color = "black", linewidth = 0.3) +
      geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.6, color = "black", linewidth = 0.3)
  }
  p <- p + facet_wrap(as.formula(paste("~", facet_by)), scales = "free_y", ncol = 4)
  if (requireNamespace("ggpubr", quietly = TRUE) && length(valid_comps))
    p <- p + ggpubr::stat_compare_means(
      comparisons = valid_comps, method = "wilcox.test", method.args = list(exact = FALSE),
      label = "p.signif",
      symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                         symbols = c("****", "***", "**", "*", "ns")),
      step.increase = 0.08, size = 4.5, bracket.size = 0.6, inherit.aes = TRUE)
  p <- p +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)), n.breaks = 4) +
    scale_fill_hue() + scale_color_hue() +
    labs(title = paste("Distribution of", score_col), x = NULL, y = score_col, fill = split_by) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          strip.text = element_text(face = "bold", size = 12),
          strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.8),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 13), axis.title.y = element_text(size = 13, face = "bold"),
          legend.position = "bottom", legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12), legend.key.size = unit(1.2, "cm"),
          panel.spacing = unit(1, "lines")) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE,
             override.aes = list(shape = 21, size = 4, color = "black", stroke = 0.5)), color = "none")
  n_facets <- length(unique(df[[facet_by]])); n_rows <- ceiling(n_facets / 4)
  plot_h <- max(5, n_rows * 3.5 + 1.5)
  out_path <- file.path(save_dir, paste0(plot_type, "_", score_col, ".png"))
  ggsave(out_path, p, width = 12, height = plot_h, dpi = 300, bg = "white")
  message("    saved: ", basename(out_path)); p
}

# --- load object + build AUCell rankings ONCE --------------------------------
if (!file.exists(RDS_PATH)) stop("RDS not found: ", RDS_PATH, call. = FALSE)
message("=== Loading object ===")
so <- readRDS(RDS_PATH)
DefaultAssay(so) <- "RNA"
if (inherits(so[["RNA"]], "Assay5") &&
    length(SeuratObject::Layers(so[["RNA"]], search = "counts")) > 1) so <- JoinLayers(so)
umap <- if (UMAP_REDUCTION %in% SeuratObject::Reductions(so)) UMAP_REDUCTION else
  { r <- grep("^umap", SeuratObject::Reductions(so), value = TRUE, ignore.case = TRUE)
    if (length(r)) r[1] else NA_character_ }

# --- re-embed a compartment subset with the SAME params as the colonocyte sub-
# annotation (Script 04): HVG=2000, 50 PCs, Harmony by SampleID, UMAP k=15,
# min.dist=0.3. Gives a fresh colonocyte UMAP that resolves the subtype regions,
# instead of the global broad-cell UMAP. No clustering/annotation needed here.
HARMONIZE_BY <- Sys.getenv("NR4A1_HIER_HARMONIZE_BY", "SampleID")
reembed <- function(obj, num_hvg = 2000, dims_pca = 50, min_dist = 0.3, kneigh = 15) {
  message("  re-embedding ", ncol(obj), " cells (HVG=", num_hvg, ", ", dims_pca,
          " PCs, Harmony by ", HARMONIZE_BY, ", UMAP k=", kneigh, ", min.dist=", min_dist, ")")
  DefaultAssay(obj) <- "RNA"
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = num_hvg, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = dims_pca, reduction.name = "pca", verbose = FALSE)
  if (HARMONIZE_BY %in% colnames(obj@meta.data) && requireNamespace("harmony", quietly = TRUE)) {
    obj <- harmony::RunHarmony(obj, group.by.vars = HARMONIZE_BY, reduction = "pca",
                               reduction.save = "harmony", verbose = FALSE); rd <- "harmony"
  } else { message("  [note] harmony/", HARMONIZE_BY, " unavailable -> PCA-only UMAP."); rd <- "pca" }
  obj <- RunUMAP(obj, dims = 1:dims_pca, reduction = rd, n.neighbors = kneigh,
                 min.dist = min_dist, n.epochs = 500, reduction.name = "umap_colon", verbose = FALSE)
  gc(); obj
}

# Build the scoring/plotting object ONCE from the chosen scope.
scope_l <- tolower(SCORE_SCOPE)
if (scope_l == "whole") {
  base_obj <- so; BASE_UMAP <- umap
  if (is.na(umap)) message("  [note] no UMAP reduction found; FeaturePlots will be skipped.")
} else {
  col <- if (BROAD_COLUMN %in% colnames(so@meta.data)) BROAD_COLUMN else
         if (CELLTYPE_COLUMN %in% colnames(so@meta.data)) CELLTYPE_COLUMN else NA_character_
  if (is.na(col)) stop("No ", BROAD_COLUMN, "/", CELLTYPE_COLUMN, " column to subset '",
                       SCORE_SCOPE, "'.", call. = FALSE)
  keep <- colnames(so)[as.character(so@meta.data[[col]]) == SCORE_SCOPE]
  if (length(keep) < 20) stop("Only ", length(keep), " cells in ", col, " == '", SCORE_SCOPE,
    "'. Available values: ", paste(utils::head(unique(as.character(so@meta.data[[col]])), 10),
    collapse = ", "), call. = FALSE)
  message("=== Subsetting ", length(keep), " '", SCORE_SCOPE,
          "' cells and RE-EMBEDDING (your colonocyte UMAP params) ===")
  base_obj <- subset(so, cells = keep)
  base_obj <- tryCatch(Seurat::DietSeurat(base_obj, assays = "RNA"), error = function(e) base_obj)
  base_obj <- reembed(base_obj)
  # lean it out: keep RNA + the new UMAP (drop scale.data/pca/harmony) so per-label
  # scoring copies are cheap.
  base_obj <- tryCatch(Seurat::DietSeurat(base_obj, assays = "RNA", dimreducs = "umap_colon"),
                       error = function(e) base_obj)
  BASE_UMAP <- "umap_colon"
  rm(so); gc()
}

label_cells <- function(ldir) {   # exact SCENIC cell set for a label
  f <- file.path(ldir, "aucell.csv")
  if (!file.exists(f)) return(NULL)
  if (requireNamespace("data.table", quietly = TRUE))
    as.character(data.table::fread(f, select = 1L, header = TRUE)[[1]])
  else rownames(utils::read.csv(f, row.names = 1, check.names = FALSE))
}

discover <- function() {
  d <- list.dirs(SCENIC_DIR, recursive = FALSE, full.names = FALSE)
  d[vapply(d, function(x) {
    nd <- file.path(SCENIC_DIR, x, "network")
    file.exists(file.path(nd, paste0(FOCUS_GENE, "_network_nodes_augmented.csv"))) ||
    file.exists(file.path(nd, paste0(FOCUS_GENE, "_network_nodes.csv")))
  }, logical(1))]
}

read_nodes <- function(ldir) {
  f_aug <- file.path(ldir, "network", paste0(FOCUS_GENE, "_network_nodes_augmented.csv"))
  f_bas <- file.path(ldir, "network", paste0(FOCUS_GENE, "_network_nodes.csv"))
  f <- if (file.exists(f_aug)) f_aug else f_bas
  utils::read.csv(f, stringsAsFactors = FALSE)
}

score_label <- function(label) {
  message("\n########## ", label, " ##########")
  ldir  <- file.path(SCENIC_DIR, label)
  nodes <- read_nodes(ldir)
  if (nrow(nodes) == 0) { message("  [skip] empty node table"); return(invisible()) }

  up   <- nodes$gene[nodes$role == "upstream"]
  down <- nodes$gene[nodes$role == "target"]
  gsets <- list()
  gsets[[paste0(FOCUS_GENE, "_upstream")]]   <- up
  gsets[[paste0(FOCUS_GENE, "_downstream")]] <- down
  gsets[[paste0(FOCUS_GENE, "_all")]]        <- unique(c(up, FOCUS_GENE, down))
  if ("module" %in% colnames(nodes)) {
    for (m in sort(unique(stats::na.omit(nodes$module)))) {
      g <- nodes$gene[nodes$role == "target" & !is.na(nodes$module) & nodes$module == m]
      gsets[[paste0(FOCUS_GENE, "_down_mod", m)]] <- g
    }
  }

  obj <- base_obj   # the colonocyte compartment, re-embedded once above

  # Expression for AUCell: counts if present, else the data (lognorm) layer -- the
  # saved object may not keep a counts layer; AUCell's ranking works on either.
  # Keep it SPARSE (dgCMatrix): densifying all colonocytes x ~20k genes is what
  # blows up RAM. AUCell_buildRankings accepts sparse input and ranks per block.
  exprMat <- tryCatch(SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts"),
                      error = function(e) NULL)
  if (is.null(exprMat) || !nrow(exprMat) || !ncol(exprMat))
    exprMat <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "data")
  gsets <- lapply(gsets, function(g) intersect(g, rownames(exprMat)))
  gsets <- gsets[vapply(gsets, length, 1L) >= MIN_SET]
  if (length(gsets) == 0) { message("  [skip] no gene set >= ", MIN_SET, " genes present"); return(invisible()) }
  message("  scoring ", length(gsets), " gene sets on ", ncol(exprMat),
          " ", SCORE_SCOPE, " cells (sparse)")

  nmax <- max(5L, as.integer(AUC_MAXRANK_PCT * nrow(exprMat)))
  # splitByBlocks=TRUE processes cells in chunks -> low peak memory. If an older
  # AUCell rejects the arg, retry without it.
  rk  <- tryCatch(
    AUCell::AUCell_buildRankings(exprMat, plotStats = FALSE, verbose = FALSE, splitByBlocks = TRUE),
    error = function(e) AUCell::AUCell_buildRankings(exprMat, plotStats = FALSE, verbose = FALSE))
  auc <- AUCell::getAUC(AUCell::AUCell_calcAUC(gsets, rk, aucMaxRank = nmax, verbose = FALSE))
  rm(rk, exprMat); gc()

  # attach scores to the SCORED object (every one of its cells has a value)
  for (sn in rownames(auc)) obj[[paste0("AUC_", sn)]] <- auc[sn, colnames(obj)]
  out <- file.path(ldir, "hierarchy_aucell"); if (!dir.exists(out)) dir.create(out, recursive = TRUE)
  cols <- paste0("AUC_", rownames(auc))

  # per-cell CSV (scored cells only)
  md <- intersect(c(CELLTYPE_COLUMN, CONDITION_COLUMN, "SampleID"), colnames(obj@meta.data))
  per <- data.frame(cell = colnames(obj), obj@meta.data[, c(md, cols), drop = FALSE],
                    check.names = FALSE, row.names = NULL)
  utils::write.csv(per, file.path(out, paste0(FOCUS_GENE, "_hierAUC_percell.csv")), row.names = FALSE)

  # FeaturePlot of each program on the RE-EMBEDDED colonocyte UMAP
  if (!is.na(BASE_UMAP) && BASE_UMAP %in% SeuratObject::Reductions(obj)) {
    fps <- lapply(cols, function(cn)
      FeaturePlot(obj, features = cn, reduction = BASE_UMAP, order = TRUE) +
        ggtitle(sub("^AUC_", "", cn)) + theme(plot.title = element_text(size = 10)))
    g <- if (have_patch) patchwork::wrap_plots(fps, ncol = 2) else fps[[1]]
    ggsave(file.path(out, paste0(FOCUS_GENE, "_hierAUC_umap.png")), g,
           width = 11, height = 3.2 * ceiling(length(cols) / 2) + 1,
           dpi = DPI_SETTING, bg = "white", limitsize = FALSE)
  }
  # Violins/barplots in the PIPELINE'S standard style (Script 09's plot_score_comparison):
  # faceted by colonocyte SUBTYPE (CELLTYPE_COLUMN), split by the 6 Genotype_sex
  # conditions with the defined levels + legend + within-sex Wilcoxon. One figure per
  # gene set (upstream TFs, downstream targets, all).
  v_sets <- paste0("AUC_", FOCUS_GENE, c("_upstream", "_downstream", "_all"))
  v_sets <- v_sets[v_sets %in% colnames(obj@meta.data)]
  for (vs in v_sets) {
    for (pt in c("violin", "barplot"))
      tryCatch(
        plot_score_comparison(obj, score_col = vs, facet_by = CELLTYPE_COLUMN,
                              split_by = CONDITION_COLUMN, plot_type = pt, save_dir = out),
        error = function(e) message("  [violin fail] ", vs, " (", pt, "): ", conditionMessage(e)))
  }
  # (lineage/potency-vs-celltype trace removed on request -- the comparison of
  #  interest is across CONDITIONS, handled by the faceted violins/barplots above.)
  message("  saved -> ", out)
}

labels <- discover()
if (length(GENESET_LABELS)) labels <- intersect(GENESET_LABELS, labels) else
  if (nzchar(GENESET_PATTERN)) labels <- grep(GENESET_PATTERN, labels, value = TRUE, ignore.case = TRUE)
if (!is.null(LABELS)) labels <- intersect(LABELS, labels)
if (length(labels) == 0)
  stop("No labels with 12d node CSVs under ", SCENIC_DIR,
       ". Run 12d_build_tf_network.py (and 12d_plot_tf_network.R for modules) first.",
       call. = FALSE)
message("Focus: ", FOCUS_GENE, " | scope: ", SCORE_SCOPE, " | labels: ",
        paste(labels, collapse = ", "))
for (lab in labels) tryCatch(score_label(lab),
  error = function(e) message("  [FAIL] ", lab, ": ", conditionMessage(e)))
message("\n=== 12e done === AUCell hierarchy scores + UMAP/violin under SCENIC/<label>/hierarchy_aucell/")
