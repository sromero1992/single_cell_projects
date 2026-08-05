# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 1b: CLUSTER-LEVEL DOUBLET REVIEW
# Version: 1.1 (UNIFIED build)
#
# Optional companion to 01_process_data.R. Runs on the processed object BEFORE
# annotation. Re-clusters at high resolution, averages the two doublet scores
# per cluster, flags clusters that cross a threshold, and (optionally) removes
# them and re-clusters the cleaned object so it is drop-in for Script 02.
#
# Reusable at the sub-type level too: point RDS_PATH at a subclustered object
# and set RECLUSTER = FALSE with CLUSTER_COL = the sub-cluster column; the
# DF_score / scDblFinder_score columns survive subsetting.
#
# REVIEW, DON'T AUTOPILOT: a high-scoring cluster can be proliferating cells or
# a real transitional state, not a doublet. Default ACTION = "flag". Look at
# the plots first; a cluster flagged by only ONE method is a review case.
# =============================================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(writexl)
set.seed(123)

# =============================================================================
# --- PART 1: USER CONFIGURATION (EDIT THIS SECTION) --------------------------
# =============================================================================
# --- 1.1: Paths --------------------------------------------------------------
PROJECT_NAME <- "Wu_Diet_project2"
ROOT_PATH    <- "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_wu_project2/r_process"
OUTPUT_DIR   <- file.path(ROOT_PATH, "seurat_output")
RDS_PATH   <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_processed_for_annotation.rds"))
REVIEW_DIR <- file.path(OUTPUT_DIR, "cluster_doublet_review")

# --- 1.2: Score columns (from Script 01, DOUBLET_METHOD = "both") ------------
DF_SCORE_COL    <- "DF_score"           # DoubletFinder pANN
SCDBL_SCORE_COL <- "scDblFinder_score"  # scDblFinder probability
DF_CLASS_COL    <- "DF_class"
SCDBL_CLASS_COL <- "scDblFinder_class"

# --- 1.3: Which clustering to score ------------------------------------------
# RECLUSTER = TRUE : fresh high-res clustering on harmony (isolates doublets).
# RECLUSTER = FALSE: use an existing column (set CLUSTER_COL to it).
RECLUSTER            <- TRUE
RECLUSTER_RESOLUTION <- 2.5            # 2.0-3.0 useful; try 2.5 then 3.0
RECLUSTER_REDUCTION  <- "harmony"
RECLUSTER_DIMS       <- 50            # match N_PCS_TO_USE from Script 01
CLUSTER_COL          <- "clusters_review"

# --- 1.4: Flagging rule ------------------------------------------------------
# All cutoffs are ADAPTIVE: for each method,
#     cutoff(k) = median + k * MAD  of the per-cluster means, floored.
# MAD (median absolute deviation) is not inflated by the doublet clusters
# themselves, so the cutoff tracks the clean bulk. Two k-levels are used:
#   WEAK_K   = a "moderately elevated" line (small k, lower line)
#   STRONG_K = an "extreme outlier" line     (large k, higher line)
#
# FLAG_RULE combines the two methods:
#   "hybrid"  (RECOMMENDED) flag a cluster if EITHER method is EXTREME
#             (>= STRONG_K line) OR BOTH methods are moderately elevated
#             (>= WEAK_K line). This catches lone extreme spikes (a cluster one
#             method is very confident about) AND concordant moderate clusters
#             (both methods mildly suspicious), while sparing a real cell type
#             that only one method mildly distrusts. Uses BOTH k's.
#   "either"  flag if EITHER method >= WEAK_K line. Aggressive. Uses WEAK_K.
#   "both"    flag only if BOTH methods >= WEAK_K line. Conservative. Uses WEAK_K.
FLAG_RULE <- "hybrid"

WEAK_K   <- 2.5    # moderate-outlier line (both-agree tier; also either/both)
STRONG_K <- 5.0    # extreme-outlier line (single-method spike; hybrid only)

# Absolute FLOOR: a cluster must clear this even if the adaptive line is lower.
# Protects a clean dataset (and real high-RNA cell types) from phantom flags.
DF_MEAN_THRESHOLD    <- 0.30   # DoubletFinder pANN floor
SCDBL_MEAN_THRESHOLD <- 0.30   # scDblFinder floor
FRACTION_THRESHOLD   <- 0.60   # scale-free backup: % of cluster called doublet

# --- 1.5: Action -------------------------------------------------------------
# "flag"   add cluster_dbl_flag column + plots + table; remove nothing.
# "remove" also drop flagged clusters, re-cluster the clean object at
#          ANNOTATION_RESOLUTION on the harmony embedding, and save it.
ACTION                <- "flag"
ANNOTATION_RESOLUTION <- 1.0
DPI_SETTING           <- 300

# =============================================================================
# --- PART 2: EXECUTION -------------------------------------------------------
# =============================================================================
if (!dir.exists(REVIEW_DIR)) dir.create(REVIEW_DIR, recursive = TRUE)

message("=== Script 1b: Cluster-level doublet review ===")
data <- readRDS(RDS_PATH)
DefaultAssay(data) <- "RNA"

# --- Clustering to score -----------------------------------------------------
if (RECLUSTER) {
  message(paste0("  Re-clustering at resolution ", RECLUSTER_RESOLUTION,
                 " on '", RECLUSTER_REDUCTION, "'..."))
  n_dims <- min(RECLUSTER_DIMS, ncol(Embeddings(data, RECLUSTER_REDUCTION)))
  data <- FindNeighbors(data, dims = 1:n_dims, reduction = RECLUSTER_REDUCTION,
                        graph.name = "review_nn", verbose = FALSE)
  data <- FindClusters(data, resolution = RECLUSTER_RESOLUTION,
                       graph.name = "review_nn", cluster.name = CLUSTER_COL,
                       verbose = FALSE)
}
data@meta.data[[CLUSTER_COL]] <- factor(as.character(data@meta.data[[CLUSTER_COL]]))
message(paste0("  Scoring ", nlevels(data@meta.data[[CLUSTER_COL]]),
               " clusters in '", CLUSTER_COL, "'."))

# --- Per-cluster stats --------------------------------------------------------
md <- data@meta.data
summ <- md %>%
  group_by(cluster = .data[[CLUSTER_COL]]) %>%
  summarise(
    N_Cells        = dplyr::n(),
    DF_mean        = mean(.data[[DF_SCORE_COL]],    na.rm = TRUE),
    scDbl_mean     = mean(.data[[SCDBL_SCORE_COL]], na.rm = TRUE),
    DF_frac_dbl    = mean(.data[[DF_CLASS_COL]]    == "Doublet", na.rm = TRUE),
    scDbl_frac_dbl = mean(.data[[SCDBL_CLASS_COL]] == "Doublet", na.rm = TRUE),
    .groups = "drop"
  )

# --- Adaptive cutoffs: median + k*MAD, floored (one weak + one strong each) --
madcut <- function(x, k, floor_val) max(median(x, na.rm = TRUE) + k * mad(x, na.rm = TRUE), floor_val)
DF_W <- madcut(summ$DF_mean,    WEAK_K,   DF_MEAN_THRESHOLD)
SC_W <- madcut(summ$scDbl_mean, WEAK_K,   SCDBL_MEAN_THRESHOLD)
DF_S <- madcut(summ$DF_mean,    STRONG_K, DF_MEAN_THRESHOLD)
SC_S <- madcut(summ$scDbl_mean, STRONG_K, SCDBL_MEAN_THRESHOLD)
message(sprintf("  Cutoffs: weak(DF>=%.2f, scDbl>=%.2f)  strong(DF>=%.2f, scDbl>=%.2f)  rule='%s'",
                DF_W, SC_W, DF_S, SC_S, FLAG_RULE))

# --- Flag (vectorised) -------------------------------------------------------
# frac backup reinforces the moderate (weak) tier for either/both/hybrid.
summ <- summ %>%
  mutate(
    df_weak    = DF_mean    >= DF_W | DF_frac_dbl    >= FRACTION_THRESHOLD,
    sc_weak    = scDbl_mean >= SC_W | scDbl_frac_dbl >= FRACTION_THRESHOLD,
    df_strong  = DF_mean    >= DF_S,
    sc_strong  = scDbl_mean >= SC_S,
    flagged    = case_when(
      FLAG_RULE == "both"   ~ df_weak & sc_weak,
      FLAG_RULE == "either" ~ df_weak | sc_weak,
      TRUE                  ~ (df_strong | sc_strong) | (df_weak & sc_weak)  # hybrid
    ),
    reason     = case_when(
      !flagged                 ~ "",
      df_strong & !sc_strong   ~ "DoubletFinder extreme",
      sc_strong & !df_strong   ~ "scDblFinder extreme",
      df_strong & sc_strong    ~ "both extreme",
      df_weak & sc_weak        ~ "both moderate",
      df_weak                  ~ "DoubletFinder only",
      TRUE                     ~ "scDblFinder only")
  )

flagged_clusters <- as.character(summ$cluster[summ$flagged])
data$cluster_dbl_flag <- ifelse(
  as.character(data@meta.data[[CLUSTER_COL]]) %in% flagged_clusters,
  "flagged_doublet_cluster", "keep")
n_flag_cells <- sum(data$cluster_dbl_flag == "flagged_doublet_cluster")

message(paste0("  FLAG_RULE='", FLAG_RULE, "' -> ", length(flagged_clusters),
               " cluster(s) flagged (", n_flag_cells, " cells, ",
               round(100 * n_flag_cells / ncol(data), 1), "%): ",
               paste(flagged_clusters, collapse = ", ")))

# --- Bar plot (your plot + threshold lines + flag marks) --------------------
lev <- levels(summ$cluster)
num <- suppressWarnings(as.numeric(lev))
ord <- if (!any(is.na(num))) lev[order(num)] else lev
bar_df <- summ %>%
  select(cluster, DoubletFinder = DF_mean, scDblFinder = scDbl_mean) %>%
  pivot_longer(-cluster, names_to = "Method", values_to = "Mean_Score") %>%
  mutate(cluster = factor(as.character(cluster), levels = ord))

p_bar <- ggplot(bar_df, aes(cluster, Mean_Score, fill = Method)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = DF_W, linetype = "dashed", color = "#F8766D") +
  geom_hline(yintercept = SC_W, linetype = "dashed", color = "#00BFC4") +
  geom_hline(yintercept = DF_S, linetype = "dotted", color = "#F8766D") +
  geom_hline(yintercept = SC_S, linetype = "dotted", color = "#00BFC4") +
  scale_fill_manual(values = c(DoubletFinder = "#F8766D", scDblFinder = "#00BFC4")) +
  theme_classic() +
  labs(title = "Mean doublet score per cluster",
       subtitle = paste0("Dashed = moderate (WEAK_K), dotted = extreme (STRONG_K). ",
                         "Flagged (", FLAG_RULE, "): ",
                         if (length(flagged_clusters)) paste(flagged_clusters, collapse = ", ") else "none"),
       x = CLUSTER_COL, y = "Mean score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(REVIEW_DIR, "01_mean_doublet_score_per_cluster.png"),
       p_bar, width = 13, height = 5, dpi = DPI_SETTING, bg = "white")

# --- Violin of raw scores ----------------------------------------------------
p_vln <- VlnPlot(data, features = c(SCDBL_SCORE_COL, DF_SCORE_COL),
                 group.by = CLUSTER_COL, pt.size = 0, ncol = 1) &
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
ggsave(file.path(REVIEW_DIR, "02_doublet_score_violin_per_cluster.png"),
       p_vln, width = 14, height = 8, dpi = DPI_SETTING, bg = "white")

# --- UMAP: clusters + flag ---------------------------------------------------
red <- if ("umap_harmony" %in% names(data@reductions)) "umap_harmony" else "umap"
p_um <- (DimPlot(data, group.by = CLUSTER_COL, reduction = red, label = TRUE) +
           NoLegend() + coord_fixed() + ggtitle(CLUSTER_COL)) +
        (DimPlot(data, group.by = "cluster_dbl_flag", reduction = red,
                 cols = c(keep = "grey85", flagged_doublet_cluster = "#B2182B")) +
           coord_fixed() + ggtitle("Flagged clusters"))
ggsave(file.path(REVIEW_DIR, "03_umap_flagged_clusters.png"),
       p_um, width = 15, height = 7, dpi = DPI_SETTING, bg = "white")

# --- Summary table -----------------------------------------------------------
out_tab <- summ %>%
  transmute(Cluster = as.character(cluster), N_Cells,
            DF_mean = round(DF_mean, 4), scDbl_mean = round(scDbl_mean, 4),
            DF_pct_doublet = round(100 * DF_frac_dbl, 1),
            scDbl_pct_doublet = round(100 * scDbl_frac_dbl, 1),
            Flagged = flagged, Reason = reason) %>%
  arrange(desc(Flagged), desc(pmax(DF_mean, scDbl_mean)))
write_xlsx(list(Per_Cluster = out_tab), file.path(REVIEW_DIR, "cluster_doublet_summary.xlsx"))

# --- Save flagged object -----------------------------------------------------
flag_path <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_processed_for_annotation_dbl_reviewed.rds"))
saveRDS(data, flag_path)
message(paste0("  Flagged object saved: ", basename(flag_path)))

# --- Remove + re-cluster (only when asked, and only if something flagged) ----
if (ACTION == "remove" && length(flagged_clusters) > 0) {
  data <- subset(data, subset = cluster_dbl_flag == "keep")
  message(paste0("  Removed flagged clusters -> ", ncol(data),
                 " cells. Re-clustering at resolution ", ANNOTATION_RESOLUTION, "..."))
  n_dims <- min(RECLUSTER_DIMS, ncol(Embeddings(data, "harmony")))
  data <- FindNeighbors(data, dims = 1:n_dims, reduction = "harmony",
                        graph.name = "harmony_nn", verbose = FALSE)
  data <- FindClusters(data, resolution = ANNOTATION_RESOLUTION,
                       graph.name = "harmony_nn", cluster.name = "clusters_harmony",
                       verbose = FALSE)
  data <- RunUMAP(data, dims = 1:n_dims, reduction = "harmony",
                  reduction.name = "umap_harmony", n.epochs = 500, verbose = FALSE)
  clean_path <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_processed_for_annotation_dbl_removed.rds"))
  saveRDS(data, clean_path)
  message(paste0("  Cleaned + re-clustered object saved: ", basename(clean_path)))
  message("  -> Point Script 02 RDS_PATH here.")
}

message("\n=== Script 1b complete. Review plots in: ", REVIEW_DIR, " ===")
message("  Open 01_mean_doublet_score_per_cluster.png and 03_umap_flagged_clusters.png,")
message("  confirm the flagged clusters are real artefacts, then set ACTION <- \"remove\".")
