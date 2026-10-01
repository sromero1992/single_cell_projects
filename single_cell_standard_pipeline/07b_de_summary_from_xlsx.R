# =============================================================================
# SCRIPT 07b: DE SUMMARY PLOT FROM SPREADSHEETS (post-DE, no MAST re-run)
# =============================================================================
# Rebuilds the "DE gene counts per cell type per contrast" summary bar plot by
# reading the per-cell-type MAST DE xlsx files already written by Script 07.
# Lets you re-threshold (padj / |log2FC|) WITHOUT re-running the DE, and orders
# the cell-type axis by lineage (subtypes grouped) for the subtype layer.
#
# INPUT : <DE_DIR>/<CellType>_MAST_DE_all_genes.xlsx  (or _MAST_DE.xlsx)
# OUTPUT: <DE_DIR>/SUMMARY_DE_gene_counts.(png|csv)   (overwrites/refreshes)
# =============================================================================
suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# =============================================================================
# --- CONFIG ------------------------------------------------------------------
# =============================================================================
# Shared portable config (nr4a1 defaults; override via env vars). Run from the
# pipeline directory, or set NR4A1_CONFIG=/full/path/to/config.R.
.NR4A1_CFG <- Sys.getenv("NR4A1_CONFIG", "config.R")
if (file.exists(.NR4A1_CFG)) source(.NR4A1_CFG) else
  message("[07b] config.R not found; using the literal DE_DIR below.")

# Which level folder to summarize ("subtypes" or "broad"). DE_DIR derives from the
# config OUTPUT_DIR unless you override NR4A1_DE_DIR or edit DE_DIR directly.
DE_LEVEL <- Sys.getenv("NR4A1_DE_LEVEL", "subtypes")
DE_DIR   <- Sys.getenv(
  "NR4A1_DE_DIR",
  if (exists("OUTPUT_DIR")) file.path(OUTPUT_DIR, "DE_results", DE_LEVEL)
  else "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_nr4a1_ack/r_process/seurat_output/DE_results/subtypes"
)

# Which files to read. "_MAST_DE_all_genes.xlsx" lets you re-threshold freely;
# "_MAST_DE.xlsx" is the already-significant subset (padj < 0.05).
FILE_SUFFIX <- "_MAST_DE_all_genes.xlsx"

# Significance thresholds for COUNTING up/down (re-tune here, no DE re-run).
PADJ_THRESH  <- 0.05
LOGFC_THRESH <- 0.5          # |avg_log2FC| >= this counts as up/down; 0 = any sig gene

MODE_LABEL   <- DE_LEVEL     # only used in the plot title (follows DE_LEVEL)

# Lineage ordering (subtype layer). Provide EITHER:
#   (a) LINEAGE_ORDER: an explicit vector of cell types in the order you want, or
#   (b) RDS_PATH + CELLTYPE_COLUMN: derive CellType -> CellType_broad from the object.
# Leave both NULL for alphabetical.
LINEAGE_ORDER   <- NULL
RDS_PATH        <- NULL      # e.g. ".../Nr4a1_s17_ack_with_cell_scores.rds"
CELLTYPE_COLUMN <- "CellType"

DPI_SETTING <- 300

# =============================================================================
# --- READ + COUNT ------------------------------------------------------------
# =============================================================================
files <- list.files(DE_DIR, pattern = paste0(gsub("\\.", "\\\\.", FILE_SUFFIX), "$"),
                     full.names = TRUE)
if (length(files) == 0)
  stop("No '", FILE_SUFFIX, "' files found in:\n  ", DE_DIR, call. = FALSE)
message("Found ", length(files), " DE files.")

rows <- list()
for (f in files) {
  ct <- sub(FILE_SUFFIX, "", basename(f), fixed = TRUE)
  for (sheet in openxlsx::getSheetNames(f)) {
    df <- tryCatch(openxlsx::read.xlsx(f, sheet = sheet), error = function(e) NULL)
    if (is.null(df) || !all(c("avg_log2FC", "p_val_adj") %in% colnames(df))) next
    sig <- df$p_val_adj < PADJ_THRESH & !is.na(df$p_val_adj)
    n_up <- sum(sig & df$avg_log2FC >=  LOGFC_THRESH, na.rm = TRUE)
    n_dn <- sum(sig & df$avg_log2FC <= -LOGFC_THRESH, na.rm = TRUE)
    rows[[paste(ct, sheet, sep = "|")]] <- data.frame(
      cell_type = ct, contrast = sheet,
      Up = n_up, Down = n_dn, stringsAsFactors = FALSE)
  }
}
counts <- dplyr::bind_rows(rows)
if (nrow(counts) == 0) stop("No usable sheets (need avg_log2FC + p_val_adj columns).", call. = FALSE)

de_summary <- counts %>%
  tidyr::pivot_longer(c("Up", "Down"), names_to = "direction", values_to = "n_genes")

# =============================================================================
# --- LINEAGE ORDER -----------------------------------------------------------
# =============================================================================
cts <- sort(unique(as.character(de_summary$cell_type)))
if (!is.null(LINEAGE_ORDER)) {
  cts <- c(intersect(LINEAGE_ORDER, cts), setdiff(cts, LINEAGE_ORDER))
} else if (!is.null(RDS_PATH) && file.exists(RDS_PATH)) {
  suppressPackageStartupMessages(library(Seurat))
  md <- readRDS(RDS_PATH)@meta.data
  if (all(c(CELLTYPE_COLUMN, "CellType_broad") %in% colnames(md))) {
    tab      <- table(as.character(md[[CELLTYPE_COLUMN]]), as.character(md$CellType_broad))
    broad_of <- colnames(tab)[max.col(tab, ties.method = "first")]
    names(broad_of) <- rownames(tab)
    key <- broad_of[cts]; key[is.na(key)] <- "zzz"
    cts <- cts[order(key, cts)]
  }
}
de_summary$cell_type <- factor(de_summary$cell_type, levels = rev(cts))
de_summary$direction <- factor(de_summary$direction, levels = c("Up", "Down"))
de_summary$n_signed  <- ifelse(de_summary$direction == "Down", -de_summary$n_genes, de_summary$n_genes)

# =============================================================================
# --- PLOT + SAVE -------------------------------------------------------------
# =============================================================================
p <- ggplot(de_summary, aes(x = n_signed, y = cell_type, fill = direction)) +
  geom_col() +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
  scale_fill_manual(values = c("Up" = "#d73027", "Down" = "#4575b4")) +
  scale_x_continuous(labels = abs) +
  facet_wrap(~ contrast, ncol = 2) +
  labs(title = paste0("DE Gene Counts per Cell Type — ", MODE_LABEL,
                      "  (padj<", PADJ_THRESH, ", |log2FC|>=", LOGFC_THRESH, ")"),
       x = "Number of Significant Genes", y = NULL, fill = "Direction") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        strip.text = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 11),
        legend.position = "top")

n_ct <- length(unique(de_summary$cell_type))
ggsave(file.path(DE_DIR, "SUMMARY_DE_gene_counts.png"), p,
       width = 14, height = max(6, n_ct * 0.4 + 3), dpi = DPI_SETTING)
write.csv(counts, file.path(DE_DIR, "SUMMARY_DE_gene_counts.csv"), row.names = FALSE)
message("Saved: SUMMARY_DE_gene_counts.png + .csv  in\n  ", DE_DIR)
