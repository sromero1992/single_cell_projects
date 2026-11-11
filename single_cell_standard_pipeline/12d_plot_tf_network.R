# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 12d (R): PLOT TF -> FOCUS -> TARGETS  + downstream GRN
# =============================================================================
# Reads the node/edge CSVs from 12d_build_tf_network.py (the SCENIC backbone:
# upstream TFs -> focus -> targets) and, for the DOWNSTREAM target genes, computes
# a co-expression GRN from the SINGLE-CELL expression of that label's cells
# (Spearman correlation, thresholded), then draws the combined network:
#
#   upstream TFs --(directed, SCENIC)-->  FOCUS  --(directed, SCENIC)-->  targets
#                                                   \__ target<->target mesh __/
#                                                      (undirected, Spearman)
#
# WHY two sources: SCENIC's adjacency score is TF-centric regression (TF->target),
# so it is the wrong tool for target<->target wiring; that comes from expression.
# Backbone edges: solid = motif-confirmed regulon, dashed = adjacency-only.
# Co-expression edges: undirected arcs, width ~ |rho|; targets colored by Louvain
# co-expression module.
#
# INPUT : SCENIC/<label>/network/<FOCUS>_network_{nodes,edges}.csv   (from 12d .py)
#         SCENIC/<label>/aucell.csv                                  (the label's cells)
#         <OUTPUT_DIR>/<PROJECT>_with_cell_scores.rds                (expression)
# OUTPUT: SCENIC/<label>/network/<FOCUS>_network_ggraph.png
#         SCENIC/<label>/network/<FOCUS>_network_{nodes,edges}_augmented.csv
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2) })

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
FOCUS_GENE  <- Sys.getenv("NR4A1_FOCUS_GENE", "Nr4a1-cust")
RDS_PATH    <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_with_cell_scores.rds"))
ADD_COEXPR  <- toupper(Sys.getenv("NR4A1_NET_COEXPR", "TRUE")) %in% c("TRUE","1","YES")
COEXPR_R    <- as.numeric(Sys.getenv("NR4A1_NET_COEXPR_R", "0.3"))   # |assoc| threshold
# Downstream target<->target GRN method: "spearman" (default; rank correlation) or
# "pcnet" (scTenifoldNet's PC-regression network, same as Script 13). Both scaled
# to ~[-1,1] so COEXPR_R applies to either.
COEXPR_METHOD <- tolower(Sys.getenv("NR4A1_NET_COEXPR_METHOD", "spearman"))
DPI_SETTING <- 300
LABELS      <- NULL   # NULL = every label that has the network CSVs

# Genes to DROP from every network before plotting (nodes + their edges). The
# focus is the custom Nr4a1-cust feature, so the plain "Nr4a1" node is removed to
# avoid a duplicate hub. Comma-separated env override; "" disables.
DROP_GENES  <- Filter(nzchar, trimws(strsplit(
  Sys.getenv("NR4A1_NET_DROP_GENES", "Nr4a1"), ",")[[1]]))

# --- Per-condition GRN (one hub-and-ring PNG per Genotype_sex) -----------------
# For each condition: node FILL = mean expression in that condition's cells,
# node SIZE = % of cells expressing it, upstream TFs connect ONLY to the focus
# (static SCENIC backbone), target<->target ring edges = per-condition single-cell
# Spearman co-expression (rewires per condition).
PER_CONDITION_GRN <- toupper(Sys.getenv("NR4A1_NET_PER_CONDITION", "TRUE")) %in% c("TRUE","1","YES")
CONDITION_COLUMN  <- Sys.getenv("NR4A1_NET_CONDITION_COL", "Genotype_sex")
CONDITION_ORDER   <- Filter(nzchar, trimws(strsplit(Sys.getenv(
  "NR4A1_NET_CONDITION_ORDER",
  "WT_Female,Polyp_Female,Polyp_NR4a1_KO_Female,WT_Male,Polyp_Male,Polyp_NR4a1_KO_Male"),
  ",")[[1]]))
CELLTYPE_COLUMN   <- Sys.getenv("NR4A1_NET_CELLTYPE_COL", "CellType")
# Minimum cells to draw a per-condition panel. Node expression/% is meaningful at
# low N, so this is permissive; the co-expression RING separately needs >=20 cells
# (in .coexpr_edges) and is simply omitted when too few, rather than skipping the
# whole panel. Lowered so sparse groups (e.g. male tumour epithelium) still plot.
MIN_CELLS_LABEL <- as.integer(Sys.getenv("NR4A1_NET_MIN_CELLS_LABEL", "10"))  # total over the label
MIN_CELLS_COND  <- as.integer(Sys.getenv("NR4A1_NET_MIN_CELLS_COND",  "5"))   # per condition panel

for (pkg in c("ggraph","tidygraph","igraph"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("12d plot needs ggraph + tidygraph + igraph (added to 00_rlibs_installation.R).",
         call. = FALSE)

# --- load the Seurat object once (needed for co-expression or per-condition) --
so <- NULL
if (ADD_COEXPR || PER_CONDITION_GRN) {
  if (file.exists(RDS_PATH)) {
    suppressPackageStartupMessages(library(Seurat))
    message("Loading expression object: ", basename(RDS_PATH))
    so <- readRDS(RDS_PATH)
    SeuratObject::DefaultAssay(so) <- "RNA"
    if (inherits(so[["RNA"]], "Assay5") &&
        length(SeuratObject::Layers(so[["RNA"]], search = "data")) > 1)
      so <- SeuratObject::JoinLayers(so)
  } else {
    message("[note] RDS not found (", RDS_PATH,
            "); backbone only, no co-expression / per-condition plots.")
    ADD_COEXPR <- FALSE
    PER_CONDITION_GRN <- FALSE
  }
}

discover <- function() {
  if (!dir.exists(SCENIC_DIR)) stop("SCENIC_DIR not found: ", SCENIC_DIR, call. = FALSE)
  d <- list.dirs(SCENIC_DIR, recursive = FALSE, full.names = FALSE)
  d[vapply(d, function(x) file.exists(file.path(
    SCENIC_DIR, x, "network", paste0(FOCUS_GENE, "_network_nodes.csv"))), logical(1))]
}

# cell barcodes for a label = rows of its aucell.csv (exactly the SCENIC cell set)
label_cells <- function(ldir) {
  f <- file.path(ldir, "aucell.csv")
  if (!file.exists(f)) return(NULL)
  if (requireNamespace("data.table", quietly = TRUE))
    as.character(data.table::fread(f, select = 1L, header = TRUE)[[1]])
  else rownames(utils::read.csv(f, row.names = 1, check.names = FALSE))
}

# Spearman co-expression among `genes` on the label's cells -> edges + modules
downstream_coexpr <- function(label, genes) {
  if (is.null(so) || length(genes) < 3) return(NULL)
  cells <- label_cells(file.path(SCENIC_DIR, label))
  cells <- intersect(cells, colnames(so))
  if (length(cells) < 20) return(NULL)
  genes <- intersect(genes, rownames(so))
  if (length(genes) < 3) return(NULL)
  expr <- as.matrix(SeuratObject::GetAssayData(so, assay = "RNA", layer = "data")[genes, cells, drop = FALSE])
  keep <- apply(expr, 1, function(v) stats::sd(v) > 0)
  expr <- expr[keep, , drop = FALSE]
  if (nrow(expr) < 3) return(NULL)
  if (COEXPR_METHOD == "pcnet") {
    if (!requireNamespace("scTenifoldNet", quietly = TRUE))
      stop("NR4A1_NET_COEXPR_METHOD='pcnet' needs scTenifoldNet (used by Script 13).",
           call. = FALSE)
    rho <- as.matrix(scTenifoldNet::pcNet(
      expr, nComp = max(2, min(10, ncol(expr) - 1, nrow(expr) - 1)),
      scaleScores = TRUE, symmetric = TRUE, q = 0))
  } else {
    rho <- suppressWarnings(stats::cor(t(expr), method = "spearman"))
  }
  rho[!is.finite(rho)] <- 0
  diag(rho) <- 0
  ut <- which(upper.tri(rho) & abs(rho) >= COEXPR_R, arr.ind = TRUE)
  if (nrow(ut) == 0)
    return(list(edges = NULL,
                modules = setNames(rep(NA_integer_, nrow(expr)), rownames(expr))))
  ed <- data.frame(source = rownames(rho)[ut[, 1]],
                   target = rownames(rho)[ut[, 2]],
                   rho    = rho[ut], stringsAsFactors = FALSE)
  gmesh <- igraph::graph_from_data_frame(
    data.frame(ed$source, ed$target, weight = abs(ed$rho)),
    directed = FALSE, vertices = data.frame(name = rownames(expr)))
  mem <- igraph::membership(igraph::cluster_louvain(gmesh))
  list(edges = ed, modules = mem[rownames(expr)])
}

# ---- helpers for the per-condition hub-and-ring GRN --------------------------
# Co-expression edges (source,target,rho) among `genes` on an arbitrary cell set.
.coexpr_edges <- function(genes, cells) {
  genes <- intersect(genes, rownames(so)); cells <- intersect(cells, colnames(so))
  if (length(genes) < 3 || length(cells) < 20) return(NULL)
  expr <- as.matrix(SeuratObject::GetAssayData(so, assay = "RNA", layer = "data")[genes, cells, drop = FALSE])
  keep <- apply(expr, 1, function(v) stats::sd(v) > 0)
  expr <- expr[keep, , drop = FALSE]
  if (nrow(expr) < 3) return(NULL)
  if (COEXPR_METHOD == "pcnet") {
    rho <- as.matrix(scTenifoldNet::pcNet(
      expr, nComp = max(2, min(10, ncol(expr) - 1, nrow(expr) - 1)),
      scaleScores = TRUE, symmetric = TRUE, q = 0))
  } else {
    rho <- suppressWarnings(stats::cor(t(expr), method = "spearman"))
  }
  rho[!is.finite(rho)] <- 0; diag(rho) <- 0
  ut <- which(upper.tri(rho) & abs(rho) >= COEXPR_R, arr.ind = TRUE)
  if (nrow(ut) == 0) return(NULL)
  data.frame(source = rownames(rho)[ut[, 1]], target = rownames(rho)[ut[, 2]],
             rho = rho[ut], stringsAsFactors = FALSE)
}

# Cells belonging to a SCENIC label. Prefer its aucell.csv (exact SCENIC set); if
# that is absent (CLI run), reconstruct from "<CellType>__<split>" by matching the
# sanitized CellType and (if present) the sex/condition suffix in the label.
.label_cells_any <- function(label) {
  cells <- label_cells(file.path(SCENIC_DIR, label))
  if (!is.null(cells)) return(intersect(cells, colnames(so)))
  md <- so[[]]
  if (!CELLTYPE_COLUMN %in% colnames(md)) return(NULL)
  parts <- strsplit(label, "__", fixed = TRUE)[[1]]
  base <- parts[1]; split_val <- if (length(parts) > 1) parts[2] else NA_character_
  ct_vals <- as.character(md[[CELLTYPE_COLUMN]])
  uniq <- unique(ct_vals)
  san  <- setNames(gsub("[^A-Za-z0-9]+", "_", uniq), uniq)
  ct   <- names(san)[match(base, san)]
  if (is.na(ct)) return(NULL)
  bc <- rownames(md)[ct_vals == ct]
  if (!is.na(split_val) && CONDITION_COLUMN %in% colnames(md)) {
    gs  <- as.character(md[[CONDITION_COLUMN]])
    # EXACT match, never a loose suffix: "Male$" also matches "Female" (fe-MALE),
    # which was leaking Female cells into __Male labels. Match the full condition
    # label first, else the last underscore-delimited token (the sex) exactly.
    sel <- gs == split_val
    if (!any(sel)) sel <- tolower(sub("^.*_", "", gs)) == tolower(split_val)
    bc  <- intersect(bc, rownames(md)[sel])
  }
  bc
}

# One hub-and-ring PNG per condition: node fill = mean expr, size = % expressing,
# upstream TFs -> focus only (static SCENIC), target<->target ring = per-condition
# Spearman co-expression (rewires each condition).
plot_label_per_condition <- function(label) {
  ndir  <- file.path(SCENIC_DIR, label, "network")
  nodes <- utils::read.csv(file.path(ndir, paste0(FOCUS_GENE, "_network_nodes.csv")), stringsAsFactors = FALSE)
  edges <- utils::read.csv(file.path(ndir, paste0(FOCUS_GENE, "_network_edges.csv")), stringsAsFactors = FALSE)
  if (length(DROP_GENES)) {
    nodes <- nodes[!nodes$gene %in% DROP_GENES, , drop = FALSE]
    edges <- edges[!(edges$source %in% DROP_GENES | edges$target %in% DROP_GENES), , drop = FALSE]
  }
  if (nrow(nodes) == 0) { message("  [per-cond skip] ", label, ": empty after drop"); return(invisible()) }
  edges <- edges[edges$source %in% nodes$gene & edges$target %in% nodes$gene, , drop = FALSE]
  bb    <- edges[edges$source == FOCUS_GENE | edges$target == FOCUS_GENE, , drop = FALSE]  # static backbone only

  up <- nodes$gene[nodes$role == "upstream"]
  tg <- nodes$gene[nodes$role == "target"]
  foc <- FOCUS_GENE

  # --- hub-and-ring layout (shared across conditions) -------------------------
  nt <- length(tg)
  r  <- max(1.8, nt * 0.14)
  cx <- r + 1.6                       # ring sits entirely right of focus (x=0): never wraps into upstream
  coord <- data.frame(gene = foc, x = 0, y = 0, role = "focus", stringsAsFactors = FALSE)
  if (length(up)) {
    yy <- seq_along(up); yy <- (yy - mean(yy)) * 0.9
    coord <- rbind(coord, data.frame(gene = up, x = -3.0, y = yy, role = "upstream"))
  }
  if (nt) {
    th <- seq(0, 2 * pi, length.out = nt + 1)[-(nt + 1)]
    coord <- rbind(coord, data.frame(gene = tg, x = cx + r * cos(th), y = r * sin(th), role = "target"))
  }
  xy <- function(g) coord[match(g, coord$gene), c("x", "y")]

  # label anchors: upstream left, targets radially outward, focus above the hub
  coord$lx <- coord$x; coord$ly <- coord$y
  coord$hjust <- ifelse(coord$role == "upstream", 1.08, 0.5)
  isT <- coord$role == "target"
  if (any(isT)) {
    ang <- atan2(coord$y[isT], coord$x[isT] - cx)
    coord$lx[isT] <- coord$x[isT] + 0.30 * cos(ang)
    coord$ly[isT] <- coord$y[isT] + 0.30 * sin(ang)
    coord$hjust[isT] <- ifelse(cos(ang) >= 0, 0, 1)
  }
  coord$lx[coord$role == "upstream"] <- coord$x[coord$role == "upstream"] - 0.1
  coord$ly[coord$role == "focus"] <- 0.55

  # upstream -> focus (grey input) kept separate from focus -> target (orange GRN)
  bb_up <- bb[bb$target == FOCUS_GENE, , drop = FALSE]
  bb_dn <- bb[bb$source == FOCUS_GENE, , drop = FALSE]
  seg_up <- if (nrow(bb_up)) data.frame(
    x = xy(bb_up$source)$x, y = xy(bb_up$source)$y,
    xend = xy(bb_up$target)$x, yend = xy(bb_up$target)$y,
    confirmed = bb_up$confirmed, stringsAsFactors = FALSE) else NULL
  seg_dn <- if (nrow(bb_dn)) data.frame(
    x = xy(bb_dn$source)$x, y = xy(bb_dn$source)$y,
    xend = xy(bb_dn$target)$x, yend = xy(bb_dn$target)$y,
    stringsAsFactors = FALSE) else NULL

  md <- so[[]]
  if (!CONDITION_COLUMN %in% colnames(md)) { message("  [per-cond skip] ", label, ": no ", CONDITION_COLUMN); return(invisible()) }
  cells <- .label_cells_any(label)
  if (is.null(cells) || length(cells) < MIN_CELLS_LABEL) {
    message("  [per-cond skip] ", label, ": ", length(cells),
            " cells resolved (< MIN_CELLS_LABEL=", MIN_CELLS_LABEL, ")"); return(invisible()) }
  gs <- setNames(as.character(md[cells, CONDITION_COLUMN]), cells)
  conds <- intersect(CONDITION_ORDER, unique(gs))
  conds <- c(conds, setdiff(unique(gs[!is.na(gs)]), conds))
  message("  [per-cond] ", label, ": ", length(cells), " cells over conditions {",
          paste(sprintf("%s=%d", names(table(gs)), as.integer(table(gs))), collapse = ", "), "}")

  genes_in <- intersect(coord$gene, rownames(so))
  expr_all <- SeuratObject::GetAssayData(so, assay = "RNA", layer = "data")[genes_in, cells, drop = FALSE]

  for (cond in conds) {
    cc <- names(gs)[gs == cond]
    if (length(cc) < MIN_CELLS_COND) {
      message("  [per-cond] ", label, " / ", cond, ": ", length(cc),
              " cells (< MIN_CELLS_COND=", MIN_CELLS_COND, "), skip"); next }
    ex <- expr_all[, cc, drop = FALSE]
    mexpr <- Matrix::rowMeans(ex); pct <- Matrix::rowMeans(ex > 0) * 100
    nd <- coord
    nd$mean_expr <- unname(mexpr[match(nd$gene, names(mexpr))]); nd$mean_expr[is.na(nd$mean_expr)] <- 0
    nd$pct       <- unname(pct[match(nd$gene, names(pct))]);     nd$pct[is.na(nd$pct)] <- 0

    ring <- .coexpr_edges(tg, cc)
    seg_ring <- if (!is.null(ring) && nrow(ring)) data.frame(
      x = xy(ring$source)$x, y = xy(ring$source)$y,
      xend = xy(ring$target)$x, yend = xy(ring$target)$y,
      w = abs(ring$rho), stringsAsFactors = FALSE) else NULL

    p <- ggplot2::ggplot()
    if (!is.null(seg_ring))
      p <- p + ggplot2::geom_curve(data = seg_ring,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linewidth = w),
        colour = "#e6550d", alpha = 0.40, curvature = 0.18) +
        ggplot2::scale_linewidth_continuous(range = c(0.3, 1.6), guide = "none")
    if (!is.null(seg_dn))      # focus -> target: orange GRN hub spokes
      p <- p + ggplot2::geom_segment(data = seg_dn,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        colour = "#d94801", alpha = 0.75, linewidth = 0.4,
        arrow = grid::arrow(length = grid::unit(2.0, "mm"), type = "closed"))
    if (!is.null(seg_up))      # upstream TF -> focus: grey regulatory input
      p <- p + ggplot2::geom_segment(data = seg_up,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linetype = confirmed),
        colour = "grey45", alpha = 0.75,
        arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed"))
    p <- p +
      ggplot2::geom_point(data = nd,
        ggplot2::aes(x = x, y = y, fill = mean_expr, size = pct),
        shape = 21, colour = "grey20", stroke = 0.4) +
      ggplot2::geom_text(data = nd, ggplot2::aes(x = lx, y = ly, label = gene),
        hjust = nd$hjust, fontface = ifelse(nd$role == "focus", "bold", "plain"), size = 3.4) +
      ggplot2::scale_fill_gradientn(colours = c("#f7fbff", "#6baed6", "#08306b"), name = "mean expr") +
      ggplot2::scale_size_continuous(range = c(2, 11), limits = c(0, 100), name = "% expressing") +
      ggplot2::scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"), na.value = "22", name = "regulon") +
      ggplot2::coord_equal(clip = "off") +
      ggplot2::labs(title = paste0(FOCUS_GENE, " GRN — ", label, "  |  ", cond),
        subtitle = "node fill = mean expr, size = % expressing | grey = upstream TF input | orange = GRN (hub spokes + per-condition target co-expression)") +
      ggplot2::theme_void(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40", size = 8),
        legend.position = "right", plot.margin = ggplot2::margin(10, 60, 10, 40))

    out <- file.path(ndir, paste0(FOCUS_GENE, "_", gsub("[^A-Za-z0-9]+", "_", cond), "_grn.png"))
    ggplot2::ggsave(out, p, width = 11, height = 9, dpi = DPI_SETTING, bg = "white")
    message("  [per-cond] ", label, " / ", cond, ": ", length(cc), " cells, ",
            if (is.null(ring)) 0 else nrow(ring), " ring edges -> ", basename(out))
  }
}

plot_label <- function(label) {
  ndir  <- file.path(SCENIC_DIR, label, "network")
  nodes <- utils::read.csv(file.path(ndir, paste0(FOCUS_GENE, "_network_nodes.csv")),
                           stringsAsFactors = FALSE)
  edges <- utils::read.csv(file.path(ndir, paste0(FOCUS_GENE, "_network_edges.csv")),
                           stringsAsFactors = FALSE)
  if (length(DROP_GENES)) {
    nodes <- nodes[!nodes$gene %in% DROP_GENES, , drop = FALSE]
    edges <- edges[!(edges$source %in% DROP_GENES | edges$target %in% DROP_GENES), , drop = FALSE]
  }
  if (nrow(nodes) == 0 || nrow(edges) == 0) { message("  [skip] ", label, ": empty"); return(invisible()) }
  edges$directed <- TRUE
  nodes$module   <- NA_integer_

  # --- downstream co-expression GRN (Spearman) --------------------------------
  if (ADD_COEXPR) {
    tgt <- nodes$gene[nodes$role == "target"]
    ce  <- downstream_coexpr(label, tgt)
    if (!is.null(ce)) {
      if (!is.null(ce$modules))
        nodes$module[match(names(ce$modules), nodes$gene)] <- as.integer(ce$modules)
      if (!is.null(ce$edges) && nrow(ce$edges) > 0) {
        add <- data.frame(source = ce$edges$source, target = ce$edges$target,
                          importance = abs(ce$edges$rho), confirmed = NA,
                          edge_type = "coexpr", directed = FALSE,
                          stringsAsFactors = FALSE)
        edges <- rbind(edges[, names(add)], add)
      }
    }
  }

  # node colour group: role for TFs/focus, co-expression module for targets
  nodes$group <- ifelse(nodes$role == "target",
                        ifelse(is.na(nodes$module), "module: n/a",
                               paste0("module ", nodes$module)),
                        nodes$role)
  mods <- sort(unique(nodes$group[grepl("^module ", nodes$group)]))
  mod_cols <- if (length(mods)) setNames(scales::hue_pal()(length(mods)), mods) else character(0)
  col_vals <- c(upstream = "#4575b4", focus = "#d73027",
                `module: n/a` = "#cccccc", mod_cols)

  # layered coords: x = layer (0/1/2); spread y within each layer
  nodes$x <- nodes$layer
  nodes$y <- stats::ave(seq_len(nrow(nodes)), nodes$layer,
                        FUN = function(i) { r <- seq_along(i); r - mean(r) })
  edges <- edges[edges$source %in% nodes$gene & edges$target %in% nodes$gene, ]

  g   <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  lay <- ggraph::create_layout(g, layout = "manual", x = nodes$x, y = nodes$y)

  nlab <- union(nodes$gene[nodes$role %in% c("focus", "upstream")],
                utils::head(nodes$gene[order(-nodes$imp_from_focus)], 20))

  p <- ggraph::ggraph(lay) +
    # co-expression mesh among downstream targets (undirected arcs)
    ggraph::geom_edge_arc(aes(filter = !directed, edge_width = importance),
                          edge_colour = "#1b9e77", edge_alpha = 0.35, strength = 0.15) +
    # SCENIC backbone (directed), solid = regulon-confirmed, dashed = adjacency-only
    ggraph::geom_edge_link(aes(filter = directed, edge_width = importance,
                               edge_linetype = confirmed),
                           arrow = grid::arrow(length = grid::unit(2, "mm"), type = "closed"),
                           end_cap = ggraph::circle(2.5, "mm"),
                           edge_colour = "grey35", edge_alpha = 0.7) +
    ggraph::scale_edge_width(range = c(0.2, 2.0), guide = "none") +
    ggraph::scale_edge_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"),
                                        name = "regulon-confirmed", na.value = "22") +
    ggraph::geom_node_point(aes(color = group, size = role == "focus")) +
    scale_color_manual(values = col_vals, name = NULL) +
    scale_size_manual(values = c(`TRUE` = 7, `FALSE` = 3), guide = "none") +
    ggraph::geom_node_text(aes(label = ifelse(gene %in% nlab, gene, "")),
                           repel = TRUE, size = 3, max.overlaps = Inf) +
    labs(title = paste0(FOCUS_GENE, " regulatory neighborhood — ", label),
         subtitle = paste0("upstream TFs → ", FOCUS_GENE,
                           " → targets (SCENIC; solid=regulon, dashed=adjacency)  |  ",
                           "target↔target green arcs = Spearman co-expression (|rho|≥",
                           COEXPR_R, ")")) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
          legend.position = "bottom")

  ggsave(file.path(ndir, paste0(FOCUS_GENE, "_network_ggraph.png")), p,
         width = 13, height = 9, dpi = DPI_SETTING, bg = "white")
  # augmented tables (backbone + co-expression) for networkx / downstream use
  utils::write.csv(nodes[, c("gene","role","layer","is_tf","module")],
                   file.path(ndir, paste0(FOCUS_GENE, "_network_nodes_augmented.csv")),
                   row.names = FALSE)
  utils::write.csv(edges, file.path(ndir, paste0(FOCUS_GENE, "_network_edges_augmented.csv")),
                   row.names = FALSE)
  message("  ", label, ": ", sum(edges$directed), " backbone + ",
          sum(!edges$directed), " co-expression edges -> ggraph png + augmented CSVs")
}

labels <- if (is.null(LABELS)) discover() else intersect(LABELS, discover())
if (length(labels) == 0)
  stop("No labels with network CSVs under ", SCENIC_DIR,
       ". Run 12d_build_tf_network.py first.", call. = FALSE)
message("Focus gene: ", FOCUS_GENE, " | co-expression: ", ADD_COEXPR,
        " (|rho|>=", COEXPR_R, ") | labels: ", paste(labels, collapse = ", "))
for (lab in labels) {
  tryCatch(plot_label(lab),
    error = function(e) message("  [FAIL] ", lab, ": ", conditionMessage(e)))
  if (PER_CONDITION_GRN)
    tryCatch(plot_label_per_condition(lab),
      error = function(e) message("  [FAIL per-cond] ", lab, ": ", conditionMessage(e)))
}
message("\n=== 12d (R) done === combined GRN PNGs + augmented CSVs under SCENIC/<label>/network/")
if (PER_CONDITION_GRN)
  message("    per-condition GRNs: SCENIC/<label>/network/", FOCUS_GENE, "_<condition>_grn.png")
