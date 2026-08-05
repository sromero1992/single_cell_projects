# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 08: CELLCHAT DIFFERENTIAL ANALYSIS
# Version: 1.0
# UNIFIED BUILD: part of unified_pipeline/. Consumes the object produced by
#   01_process_data.R v11.0. Doublet calls arrive standardised in the
#   'Doublet_Status' column regardless of which caller ran, so this script
#   requires no changes when DOUBLET_METHOD is switched.
#
# PURPOSE:
#   Sex-separated differential CellChat analysis across two comparisons:
#     1. Polyp vs WT      (disease effect)
#     2. KO vs Polyp      (NR4a1 rescue effect)
#   Run on BOTH broad cell types and sub-annotated cell types.
#   Followed by a focused re-run on the strongest differential interactors.
#   Includes immune checkpoint interaction panels (PD-L1/PD-1, CTLA4, etc.)
#
# COMPARISONS (sex-separated, Option B):
#   Female: Polyp_Female    vs WT_Female
#   Female: KO_Female       vs Polyp_Female
#   Male:   Polyp_Male      vs WT_Male
#   Male:   KO_Male         vs Polyp_Male
#
# OUTPUT STRUCTURE:
#   cellchat_results/
#     broad/
#       Female_Polyp_vs_WT/
#       Female_KO_vs_Polyp/
#       Male_Polyp_vs_WT/
#       Male_KO_vs_Polyp/
#     subtypes/
#       Female_Polyp_vs_WT/
#       ...
#     top_interactors/
#       broad/
#       subtypes/
#
# INPUT:  {PROJECT_NAME}_unified_annotated.rds  (output of 07_unifier.R)
# =============================================================================

library(Seurat)
library(CellChat)
library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(123)

# =============================================================================
# --- PART 1: CONFIGURATION ---------------------------------------------------
# =============================================================================
PROJECT_NAME <- "Wu_Diet_project2"
ROOT_PATH    <- "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_wu_project2/r_process"
OUTPUT_DIR   <- file.path(ROOT_PATH, "seurat_output")

# Input
UNIFIED_RDS  <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_unified_annotated.rds"))

# CellChat output root
CC_ROOT      <- file.path(OUTPUT_DIR, "cellchat_results")
CC_BROAD_DIR <- file.path(CC_ROOT, "broad")
CC_SUB_DIR   <- file.path(CC_ROOT, "subtypes")
CC_TOP_DIR   <- file.path(CC_ROOT, "top_interactors")

for (d in c(CC_ROOT, CC_BROAD_DIR, CC_SUB_DIR, CC_TOP_DIR,
            file.path(CC_TOP_DIR, "broad"),
            file.path(CC_TOP_DIR, "subtypes"))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# --- 1.1: Mode ---------------------------------------------------------------
# "broad"   → uses CellType_broad column
# "subtype" → uses CellType column (fine-grained)
# "both"    → runs both sequentially
MODE <- "both"   # <-- ACTION: change if needed

# --- 1.2: CellChat parameters ------------------------------------------------
MIN_CELLS       <- 10    # minimum cells per group for CellChat
DPI_SETTING     <- 300

# --- 1.3: Comparisons (sex-separated) ----------------------------------------
COMPARISONS <- list(
  Female_Polyp_vs_WT  = list(cond1 = "Polyp_Female",         cond2 = "WT_Female",
                             pos   = "Polyp_Female",         label = "Polyp_vs_WT_Female"),
  Female_KO_vs_Polyp  = list(cond1 = "Polyp_NR4a1_KO_Female", cond2 = "Polyp_Female",
                             pos   = "Polyp_NR4a1_KO_Female", label = "KO_vs_Polyp_Female"),
  Male_Polyp_vs_WT    = list(cond1 = "Polyp_Male",            cond2 = "WT_Male",
                             pos   = "Polyp_Male",            label = "Polyp_vs_WT_Male"),
  Male_KO_vs_Polyp    = list(cond1 = "Polyp_NR4a1_KO_Male",   cond2 = "Polyp_Male",
                             pos   = "Polyp_NR4a1_KO_Male",  label = "KO_vs_Polyp_Male")
)

# --- 1.4: Top interactors threshold ------------------------------------------
# Keep cell types in top X% of differential interaction strength for re-run
TOP_INTERACTOR_PCT <- 0.20   # top 20%

# --- 1.5: Immune checkpoint interactions to highlight ------------------------
# Names must match CellChatDB.mouse interaction_name EXACTLY (ALL CAPS, underscore)
CHECKPOINT_PAIRS <- c(
  "CD274_PDCD1",    # PD-L1 → PD-1: tumor/myeloid suppression of CD8+
  "PDCD1LG2_PDCD1", # PD-L2 → PD-1: alternative PD-1 ligand
  "CD80_CTLA4",     # CD80 → CTLA4: co-stimulation vs suppression
  "CD86_CTLA4",     # CD86 → CTLA4
  "ICOSL_CTLA4",    # ICOSL → CTLA4: additional co-inhibitory
  "CD80_CD274",     # CD80 → PD-L1: trans-endocytosis
  "PVR_TIGIT",      # CD155 → TIGIT: exhaustion checkpoint
  "NECTIN2_TIGIT",  # NECTIN2 → TIGIT: exhaustion checkpoint
  "NECTIN3_TIGIT"   # NECTIN3 → TIGIT: exhaustion checkpoint
)

# =============================================================================
# --- PART 2: LOAD DATA -------------------------------------------------------
# =============================================================================
message("=== Loading unified annotated object ===")
data <- readRDS(UNIFIED_RDS)
message(paste("  Loaded:", ncol(data), "cells"))
DefaultAssay(data) <- "RNA"

# =============================================================================
# --- PART 3: CELLCHAT PIPELINE FUNCTIONS -------------------------------------
# =============================================================================

# --- 3.1: Build one CellChat object for one condition -----------------------

build_cellchat <- function(so, condition, group_col, save_path = NULL) {
  
  # Skip if already computed and saved — validate before loading
  if (!is.null(save_path) && file.exists(save_path)) {
    cc_test <- tryCatch(readRDS(save_path), error = function(e) NULL)
    if (!is.null(cc_test)) {
      message(paste("  [SKIP] Already exists, loading:", basename(save_path)))
      return(cc_test)
    } else {
      message(paste("  [WARN] Corrupted RDS detected, recomputing:", basename(save_path)))
      file.remove(save_path)
    }
  }
  so_sub <- subset(so, Genotype_sex == condition)
  
  if (ncol(so_sub) < MIN_CELLS) {
    message(paste("  [SKIP] Too few cells for condition:", condition,
                  "(n =", ncol(so_sub), ")"))
    return(NULL)
  }
  
  so_sub[[group_col]] <- droplevels(factor(so_sub[[group_col, drop = TRUE]]))
  ct_counts <- table(so_sub[[group_col, drop = TRUE]])
  ct_keep   <- names(ct_counts)[ct_counts >= MIN_CELLS]
  
  if (length(ct_keep) < 2) {
    message(paste("  [SKIP] Fewer than 2 cell types with >=", MIN_CELLS,
                  "cells for:", condition))
    return(NULL)
  }
  
  so_sub <- subset(so_sub, !!sym(group_col) %in% ct_keep)
  so_sub[[group_col]] <- droplevels(factor(so_sub[[group_col, drop = TRUE]]))
  
  message(paste("  Building CellChat for:", condition,
                "| n =", ncol(so_sub),
                "| cell types =", length(ct_keep)))
  
  cc <- createCellChat(object   = so_sub,
                       meta     = so_sub@meta.data,
                       group.by = group_col,
                       assay    = "RNA")
  cc@DB <- CellChatDB.mouse
  cc    <- setIdent(cc, ident.use = group_col)
  cc    <- subsetData(cc)
  cc    <- identifyOverExpressedGenes(cc)
  cc    <- identifyOverExpressedInteractions(cc)
  cc    <- computeCommunProb(cc)
  cc    <- filterCommunication(cc, min.cells = MIN_CELLS)
  cc    <- computeCommunProbPathway(cc)
  cc    <- aggregateNet(cc)
  return(cc)
}

# --- 3.2: COMPUTATION — build and merge CellChat objects --------------------
run_comparison <- function(so, comparison, group_col, out_dir) {
  cond1 <- comparison$cond1
  cond2 <- comparison$cond2
  label <- comparison$label
  pos   <- comparison$pos
  
  comp_dir <- file.path(out_dir, label)
  if (!dir.exists(comp_dir)) dir.create(comp_dir, recursive = TRUE)
  
  rds_dir  <- file.path(out_dir, "rds")
  if (!dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE)
  cc1_path <- file.path(rds_dir, paste0("cellchat_", cond1, ".rds"))
  cc2_path <- file.path(rds_dir, paste0("cellchat_", cond2, ".rds"))
  
  merged_path <- file.path(comp_dir, "cellchat_merged.rds")
  if (file.exists(merged_path)) {
    cc_test <- tryCatch(readRDS(merged_path), error = function(e) NULL)
    if (!is.null(cc_test)) {
      message(paste("  [SKIP compute] Merged exists:", label))
      return(list(merged = cc_test, cond1 = cond1, cond2 = cond2,
                  label = label, pos = pos, comp_dir = comp_dir,
                  out_dir = out_dir))
    } else {
      message(paste("  [WARN] Corrupted merged RDS, recomputing:", label))
      file.remove(merged_path)
    }
  }
  
  message(paste("\n  Building CellChat:", cond1, "..."))
  cc1 <- build_cellchat(so, cond1, group_col, save_path = cc1_path)
  if (is.null(cc1)) return(NULL)
  
  message(paste("  Building CellChat:", cond2, "..."))
  cc2 <- build_cellchat(so, cond2, group_col, save_path = cc2_path)
  if (is.null(cc2)) return(NULL)
  
  cts1 <- levels(cc1@idents); cts2 <- levels(cc2@idents)
  shared_cts <- intersect(cts1, cts2)
  dropped    <- union(setdiff(cts1, shared_cts), setdiff(cts2, shared_cts))
  if (length(dropped) > 0) {
    message(paste("  [INFO] Dropping non-shared cell types:",
                  paste(dropped, collapse = ", ")))
    cc1 <- subsetCellChat(cc1, idents.use = shared_cts)
    cc2 <- subsetCellChat(cc2, idents.use = shared_cts)
  }
  
  saveRDS(cc1, cc1_path, compress = FALSE)
  saveRDS(cc2, cc2_path, compress = FALSE)
  
  cc_merged <- mergeCellChat(list(cond1 = cc1, cond2 = cc2),
                             add.names = c(cond1, cond2))
  rm(cc1, cc2); gc()
  saveRDS(cc_merged, merged_path, compress = FALSE)
  message(paste("  Merge saved:", label))
  
  return(list(merged = cc_merged, cond1 = cond1, cond2 = cond2,
              label = label, pos = pos, comp_dir = comp_dir,
              out_dir = out_dir))
}

# --- 3.3: PLOTTING — all plots from merged object ----------------------------
plot_comparison <- function(result) {
  if (is.null(result)) return(invisible(NULL))
  
  cc_merged  <- result$merged
  cond1      <- result$cond1
  cond2      <- result$cond2
  label      <- result$label
  pos        <- result$pos
  comp_dir   <- result$comp_dir
  out_dir    <- result$out_dir
  rds_dir    <- file.path(out_dir, "rds")
  
  message(paste("  Generating plots for:", label))
  
  # 01: Interaction count and weight
  tryCatch({
    p1 <- compareInteractions(cc_merged, show.legend = FALSE, group = c(1, 2))
    p2 <- compareInteractions(cc_merged, measure = "weight",
                              show.legend = FALSE, group = c(1, 2))
    png(file.path(comp_dir, "01_interaction_count_weight.png"),
        width = 10, height = 5, res = DPI_SETTING, units = "in")
    print(p1 + p2); dev.off()
  }, error = function(e) message(paste("  [SKIP 01]:", e$message)))
  
  # 02: Heatmap
  tryCatch({
    p_heat <- netVisual_heatmap(cc_merged, font.size = 13, font.size.title = 15)
    png(file.path(comp_dir, "02_heatmap_diff_interactions.png"),
        width = 8, height = 7, res = DPI_SETTING, units = "in")
    print(p_heat); dev.off()
  }, error = function(e) message(paste("  [SKIP 02]:", e$message)))
  
  # 03: Circle plots
  tryCatch({
    png(file.path(comp_dir, "03_circle_diff_interactions.png"),
        width = 10, height = 10, res = DPI_SETTING, units = "in")
    par(mfrow = c(1, 2), xpd = TRUE)
    netVisual_diffInteraction(cc_merged, comparison = c(1, 2),
                              weight.scale = TRUE, vertex.label.cex = 1.1,
                              vertex.size.max = 20, arrow.size = 0.8, margin = 0.1)
    netVisual_diffInteraction(cc_merged, comparison = c(1, 2), measure = "weight",
                              weight.scale = TRUE, vertex.label.cex = 1.1,
                              vertex.size.max = 20, arrow.size = 0.8, margin = 0.1)
    dev.off()
  }, error = function(e) message(paste("  [SKIP 03]:", e$message)))
  
  # 04-08: DEG mapping and bubble plots
  net.up   <- data.frame()
  net.down <- data.frame()
  
  net_up_path   <- file.path(comp_dir, "05_net_up.csv")
  net_down_path <- file.path(comp_dir, "06_net_down.csv")
  
  if (file.exists(net_up_path) && file.exists(net_down_path)) {
    message("  Loading existing net.up/net.down from CSV...")
    net.up   <- read.csv(net_up_path)
    net.down <- read.csv(net_down_path)
  } else {
    message("  Computing DEG mapping...")
    tryCatch({
      cc_merged <- identifyOverExpressedGenes(
        cc_merged, group.dataset = "datasets", pos.dataset = pos,
        features.name = "differential_genes", only.pos = FALSE,
        thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 0.05,
        group.DE.combined = FALSE)
      net <- netMappingDEG(cc_merged, features.name = "differential_genes",
                           variable.all = TRUE)
      write.csv(net, file.path(comp_dir, "04_net_mapping_DEG.csv"), row.names = FALSE)
      net.up   <- subsetCommunication(cc_merged, net = net, datasets = cond1,
                                      ligand.logFC = 0.1, thresh = 0.05,
                                      ligand.pvalues = 0.05)
      net.down <- subsetCommunication(cc_merged, net = net, datasets = cond2,
                                      ligand.logFC = -0.1, thresh = 0.05,
                                      ligand.pvalues = 0.05)
      net.up   <- net.up[net.up$prob > 0.01, ]
      net.down <- net.down[net.down$prob > 0.01, ]
      write.csv(net.up,   net_up_path,   row.names = FALSE)
      write.csv(net.down, net_down_path, row.names = FALSE)
    }, error = function(e) message(paste("  [SKIP DEG mapping]:", e$message)))
  }
  
  message(paste("  net.up:", nrow(net.up), "rows | net.down:", nrow(net.down), "rows"))
  
  # Bubble plots per source cell type
  if (nrow(net.up) > 0 || nrow(net.down) > 0) {
    filter_top_pairs <- function(net_df, top_n = 30, min_prob = 0.05) {
      if (nrow(net_df) == 0) return(data.frame(interaction_name = character(0)))
      net_df %>%
        dplyr::filter(prob >= min_prob) %>%
        dplyr::group_by(interaction_name) %>%
        dplyr::summarise(mean_prob = mean(prob, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(desc(mean_prob)) %>%
        head(top_n) %>%
        dplyr::select(interaction_name)
    }
    pairLR_up   <- filter_top_pairs(net.up)
    pairLR_down <- filter_top_pairs(net.down)
    
    source_names <- levels(cc_merged@idents[[cond1]])
    for (i in seq_along(source_names)) {
      src      <- source_names[i]
      safe_src <- gsub("[^A-Za-z0-9_]", "_", src)
      if (nrow(pairLR_up) > 0) {
        tryCatch({
          gg_up <- netVisual_bubble(cc_merged, pairLR.use = pairLR_up,
                                    sources.use = i, comparison = c(1, 2),
                                    remove.isolate = TRUE, color.heatmap = "viridis",
                                    title.name = paste0("Up in ", cond1, "\nSource: ", src),
                                    font.size = 14, font.size.title = 16, angle.x = 45) +
            theme(plot.title = element_text(hjust = 0.5),
                  legend.text = element_text(size = 12),
                  legend.title = element_text(size = 14, hjust = 0.5)) +
            guides(color = guide_colorbar(title = "Commun.\nProb."))
          h <- max(6, min(16, 4 + nrow(pairLR_up) * 0.25))
          ggsave(file.path(comp_dir, paste0("07_bubble_UP_source_", safe_src, ".png")),
                 gg_up, width = 9, height = h, dpi = DPI_SETTING)
        }, error = function(e) message(paste("  [SKIP bubble UP]", src, ":", e$message)))
      }
      if (nrow(pairLR_down) > 0) {
        tryCatch({
          gg_dn <- netVisual_bubble(cc_merged, pairLR.use = pairLR_down,
                                    sources.use = i, comparison = c(1, 2),
                                    remove.isolate = TRUE, color.heatmap = "viridis",
                                    title.name = paste0("Down in ", cond1, "\nSource: ", src),
                                    font.size = 14, font.size.title = 16, angle.x = 45) +
            theme(plot.title = element_text(hjust = 0.5),
                  legend.text = element_text(size = 12),
                  legend.title = element_text(size = 14, hjust = 0.5)) +
            guides(color = guide_colorbar(title = "Commun.\nProb."))
          h <- max(6, min(16, 4 + nrow(pairLR_down) * 0.25))
          ggsave(file.path(comp_dir, paste0("08_bubble_DOWN_source_", safe_src, ".png")),
                 gg_dn, width = 9, height = h, dpi = DPI_SETTING)
        }, error = function(e) message(paste("  [SKIP bubble DOWN]", src, ":", e$message)))
      }
    }
  }  # end if nrow(net.up) > 0 || nrow(net.down) > 0
  
  # 09-10: Immune checkpoint panel
  message("  Generating immune checkpoint plots...")
  checkpoint_dir <- file.path(comp_dir, "checkpoint_interactions")
  if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir)
  all_lr_pairs     <- rownames(cc_merged@DB$interaction)
  checkpoint_found <- CHECKPOINT_PAIRS[CHECKPOINT_PAIRS %in% all_lr_pairs]
  message(paste("  Checkpoint pairs found:", length(checkpoint_found),
                "of", length(CHECKPOINT_PAIRS)))
  if (length(checkpoint_found) > 0) {
    tryCatch({
      gg_cp <- netVisual_bubble(cc_merged,
                                pairLR.use = data.frame(interaction_name = checkpoint_found),
                                comparison = c(1, 2), remove.isolate = TRUE,
                                color.heatmap = "Spectral",
                                title.name = paste0("Immune Checkpoints\n",
                                                    cond1, " vs ", cond2),
                                font.size = 14, font.size.title = 16, angle.x = 45) +
        theme(plot.title = element_text(hjust = 0.5),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 14, hjust = 0.5)) +
        guides(color = guide_colorbar(title = "Commun.\nProb."))
      h <- max(7, min(18, 4 + length(checkpoint_found) * 0.4))
      ggsave(file.path(checkpoint_dir, "09_checkpoint_bubble_all.png"),
             gg_cp, width = 12, height = h, dpi = DPI_SETTING)
      message("  Checkpoint plot saved.")
    }, error = function(e) message(paste("  [SKIP checkpoint]:", e$message)))
    
    # PD-L1 chord — load individual RDS
    pdl1_pair <- checkpoint_found[grepl("CD274|PDCD1", checkpoint_found)]
    if (length(pdl1_pair) > 0) {
      tryCatch({
        cc1_single <- readRDS(file.path(rds_dir, paste0("cellchat_", cond1, ".rds")))
        png(file.path(checkpoint_dir, "10_PDL1_PD1_chord.png"),
            width = 8, height = 8, res = DPI_SETTING, units = "in")
        netVisual_chord_gene(cc1_single,
                             pairLR.use = data.frame(interaction_name = pdl1_pair),
                             title.name = paste0("PD-L1 -> PD-1\n", cond1))
        dev.off()
        rm(cc1_single); gc()
      }, error = function(e) {
        if (dev.cur() > 1) dev.off()
        message(paste("  [SKIP PDL1 chord]:", e$message))
      })
    }
  } else {
    message("  [WARN] No checkpoint pairs found in LRsig.")
  }
  
  # 11-12: Top pairs chord diagrams (prob x logFC scored)
  message("  Generating top pairs chord diagrams...")
  tryCatch({
    cc1_path <- file.path(rds_dir, paste0("cellchat_", cond1, ".rds"))
    cc2_path <- file.path(rds_dir, paste0("cellchat_", cond2, ".rds"))
    
    if (nrow(net.up) > 0 && file.exists(cc1_path)) {
      top_up <- net.up %>%
        dplyr::filter(prob >= 0.05) %>%
        dplyr::mutate(score = prob * abs(ligand.logFC)) %>%
        dplyr::arrange(desc(score)) %>% head(20) %>%
        dplyr::pull(interaction_name) %>% unique()
      if (length(top_up) > 0) {
        cc1_single <- readRDS(cc1_path)
        png(file.path(comp_dir, "11_chord_top_pairs_UP.png"),
            width = 10, height = 10, res = DPI_SETTING, units = "in")
        netVisual_chord_gene(cc1_single,
                             pairLR.use   = data.frame(interaction_name = top_up),
                             title.name   = paste0("Top UP pairs\n", cond1),
                             legend.pos.x = 8)
        dev.off()
        rm(cc1_single); gc()
        message(paste("  Chord UP saved:", length(top_up), "pairs"))
      }
    }
    if (nrow(net.down) > 0 && file.exists(cc2_path)) {
      top_dn <- net.down %>%
        dplyr::filter(prob >= 0.05) %>%
        dplyr::mutate(score = prob * abs(ligand.logFC)) %>%
        dplyr::arrange(desc(score)) %>% head(20) %>%
        dplyr::pull(interaction_name) %>% unique()
      if (length(top_dn) > 0) {
        cc2_single <- readRDS(cc2_path)
        png(file.path(comp_dir, "12_chord_top_pairs_DOWN.png"),
            width = 10, height = 10, res = DPI_SETTING, units = "in")
        netVisual_chord_gene(cc2_single,
                             pairLR.use   = data.frame(interaction_name = top_dn),
                             title.name   = paste0("Top DOWN pairs\n", cond2),
                             legend.pos.x = 8)
        dev.off()
        rm(cc2_single); gc()
        message(paste("  Chord DOWN saved:", length(top_dn), "pairs"))
      }
    }
  }, error = function(e) message(paste("  [SKIP chord]:", e$message)))
  
  # ==========================================================================
  # --- PART A: Signaling Role Scatter (2D sending/receiving) ----------------
  # ==========================================================================
  message("  Generating signaling role scatter plots...")
  scatter_dir <- file.path(comp_dir, "signaling_roles")
  if (!dir.exists(scatter_dir)) dir.create(scatter_dir)
  
  tryCatch({
    object.list <- list(cc_merged@net[[cond1]], cc_merged@net[[cond2]])
    # Use individual single objects loaded from RDS for proper scatter
    cc1_single <- readRDS(file.path(rds_dir, paste0("cellchat_", cond1, ".rds")))
    cc2_single <- readRDS(file.path(rds_dir, paste0("cellchat_", cond2, ".rds")))
    obj_list   <- list(cc1_single, cc2_single)
    names(obj_list) <- c(cond1, cond2)
    
    # Control dot size consistently across datasets
    num.link <- sapply(obj_list, function(x) {
      rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
    })
    weight.MinMax <- c(min(num.link), max(num.link))
    
    # (A) Overall signaling role scatter — all cell types
    gg_scatter <- lapply(seq_along(obj_list), function(i) {
      netAnalysis_signalingRole_scatter(obj_list[[i]],
                                        title        = names(obj_list)[i],
                                        weight.MinMax = weight.MinMax)
    })
    p_scatter <- patchwork::wrap_plots(plots = gg_scatter, ncol = 2)
    ggsave(file.path(scatter_dir, "A1_signaling_role_scatter.png"),
           p_scatter, width = 14, height = 7, dpi = DPI_SETTING)
    message("  Signaling role scatter saved.")
    
    # (B) Signaling changes per cell type — auto-detect top changed cell types
    ct_names <- levels(cc1_single@idents)
    for (ct in ct_names) {
      safe_ct <- gsub("[^A-Za-z0-9_]", "_", ct)
      tryCatch({
        gg_changes <- netAnalysis_signalingChanges_scatter(
          cc_merged, idents.use = ct)
        ggsave(file.path(scatter_dir,
                         paste0("B_signaling_changes_", safe_ct, ".png")),
               gg_changes, width = 8, height = 6, dpi = DPI_SETTING)
      }, error = function(e) NULL)  # skip silently if cell type has no changes
    }
    message("  Signaling changes scatter saved per cell type.")
    
    rm(cc1_single, cc2_single); gc()
  }, error = function(e) message(paste("  [SKIP scatter]:", e$message)))
  
  # ==========================================================================
  # --- PART B: Information Flow / Pathway Ranking ---------------------------
  # ==========================================================================
  message("  Generating pathway ranking plots...")
  ranking_dir <- file.path(comp_dir, "pathway_ranking")
  if (!dir.exists(ranking_dir)) dir.create(ranking_dir)
  
  tryCatch({
    gg_rank1 <- rankNet(cc_merged, mode = "comparison", measure = "weight",
                        sources.use = NULL, targets.use = NULL,
                        stacked = TRUE,  do.stat = TRUE)
    gg_rank2 <- rankNet(cc_merged, mode = "comparison", measure = "weight",
                        sources.use = NULL, targets.use = NULL,
                        stacked = FALSE, do.stat = TRUE)
    p_rank <- gg_rank1 + gg_rank2
    ggsave(file.path(ranking_dir, "C_pathway_information_flow.png"),
           p_rank, width = 16, height = 10, dpi = DPI_SETTING)
    message("  Pathway ranking saved.")
  }, error = function(e) message(paste("  [SKIP ranking]:", e$message)))
  
  # ==========================================================================
  # --- PART C: Per-pathway circle + chord + heatmap -------------------------
  # ==========================================================================
  message("  Generating per-pathway visualizations...")
  pathway_dir <- file.path(comp_dir, "pathway_plots")
  if (!dir.exists(pathway_dir)) dir.create(pathway_dir)
  
  tryCatch({
    cc1_single <- readRDS(file.path(rds_dir, paste0("cellchat_", cond1, ".rds")))
    cc2_single <- readRDS(file.path(rds_dir, paste0("cellchat_", cond2, ".rds")))
    obj_list   <- list(cc1_single, cc2_single)
    names(obj_list) <- c(cond1, cond2)
    
    # Get active pathways in either condition
    pathways_active <- union(cc1_single@netP$pathways,
                             cc2_single@netP$pathways)
    weight.max <- getMaxWeight(obj_list, slot.name = "netP",
                               attribute = pathways_active[1])  # for scaling
    
    message(paste("  Plotting", length(pathways_active), "pathways..."))
    
    for (pw in pathways_active) {
      safe_pw <- gsub("[^A-Za-z0-9_]", "_", pw)
      pw_dir  <- file.path(pathway_dir, safe_pw)
      if (!dir.exists(pw_dir)) dir.create(pw_dir)
      
      # Get max weight for this pathway for consistent scaling
      tryCatch({
        wmax <- getMaxWeight(obj_list, slot.name = "netP", attribute = pw)
      }, error = function(e) { wmax <- NULL })
      
      # Circle plot — side by side
      tryCatch({
        png(file.path(pw_dir, paste0(safe_pw, "_circle.png")),
            width = 14, height = 7, res = DPI_SETTING, units = "in")
        par(mfrow = c(1, 2), xpd = TRUE)
        for (i in seq_along(obj_list)) {
          netVisual_aggregate(obj_list[[i]], signaling = pw,
                              layout = "circle",
                              edge.weight.max = if (!is.null(wmax)) wmax[1] else NULL,
                              edge.width.max  = 10,
                              signaling.name  = paste(pw, names(obj_list)[i]))
        }
        dev.off()
      }, error = function(e) {
        if (dev.cur() > 1) dev.off()
      })
      
      # Chord diagram — side by side
      tryCatch({
        png(file.path(pw_dir, paste0(safe_pw, "_chord.png")),
            width = 14, height = 7, res = DPI_SETTING, units = "in")
        par(mfrow = c(1, 2), xpd = TRUE)
        for (i in seq_along(obj_list)) {
          netVisual_aggregate(obj_list[[i]], signaling = pw,
                              layout = "chord",
                              signaling.name = paste(pw, names(obj_list)[i]))
        }
        dev.off()
      }, error = function(e) {
        if (dev.cur() > 1) dev.off()
      })
      
      # Heatmap — side by side
      tryCatch({
        ht <- lapply(seq_along(obj_list), function(i) {
          netVisual_heatmap(obj_list[[i]], signaling = pw,
                            color.heatmap = "Spectral",
                            title.name = paste(pw, "signaling", names(obj_list)[i]))
        })
        png(file.path(pw_dir, paste0(safe_pw, "_heatmap.png")),
            width = 14, height = 6, res = DPI_SETTING, units = "in")
        ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
        dev.off()
      }, error = function(e) {
        if (dev.cur() > 1) dev.off()
      })
      
      # Gene expression violin — uses original Seurat object
      tryCatch({
        cc_merged@meta$datasets <- factor(
          cc_merged@meta$datasets, levels = c(cond2, cond1))
        gg_expr <- plotGeneExpression(cc_merged, signaling = pw,
                                      split.by = "datasets",
                                      colors.ggplot = TRUE,
                                      type = "violin")
        ggsave(file.path(pw_dir, paste0(safe_pw, "_gene_expression.png")),
               gg_expr, width = 12, height = 6, dpi = DPI_SETTING)
      }, error = function(e) NULL)
    }
    
    rm(cc1_single, cc2_single); gc()
    message(paste("  Per-pathway plots saved to:", pathway_dir))
  }, error = function(e) message(paste("  [SKIP pathway plots]:", e$message)))
  message("  Generating pathway activity dotplot...")
  tryCatch({
    pathways_use <- union(cc_merged@netP[[cond1]]$pathways,
                          cc_merged@netP[[cond2]]$pathways)
    if (length(pathways_use) > 0) {
      gg_path <- netVisual_bubble(cc_merged,
                                  sources.use = NULL, targets.use = NULL,
                                  signaling = pathways_use, comparison = c(1, 2),
                                  angle.x = 45, remove.isolate = TRUE,
                                  title.name = paste0("Pathway Activity\n",
                                                      cond1, " vs ", cond2),
                                  font.size = 12, font.size.title = 14)
      h <- max(8, min(24, 4 + length(pathways_use) * 0.35))
      ggsave(file.path(comp_dir, "13_pathway_activity_dotplot.png"),
             gg_path, width = 16, height = h, dpi = DPI_SETTING)
      message(paste("  Pathway dotplot saved:", length(pathways_use), "pathways"))
    }
  }, error = function(e) message(paste("  [SKIP pathway dotplot]:", e$message)))
  
  return(invisible(NULL))
}

run_top_interactors <- function(so, result, group_col, out_dir) {
  if (is.null(result)) return(invisible(NULL))
  
  cc_merged <- result$merged
  label     <- result$label
  
  # Get differential interaction matrix (count)
  diff_mat <- cc_merged@net[[result$cond1]]$count -
    cc_merged@net[[result$cond2]]$count
  
  # Absolute differential interaction strength per cell type
  ct_strength <- rowSums(abs(diff_mat)) + colSums(abs(diff_mat))
  threshold   <- quantile(ct_strength, 1 - TOP_INTERACTOR_PCT)
  top_cts     <- names(ct_strength)[ct_strength >= threshold]
  
  message(paste("  Top interactors (top", TOP_INTERACTOR_PCT * 100, "%):",
                paste(top_cts, collapse = ", ")))
  
  top_dir <- file.path(out_dir, paste0(label, "_top_interactors"))
  if (!dir.exists(top_dir)) dir.create(top_dir, recursive = TRUE)
  
  # Subset CellChat to top interactors only
  tryCatch({
    cc_top <- subsetCellChat(cc_merged, idents.use = top_cts)
    saveRDS(cc_top, file.path(top_dir, "cellchat_top_merged.rds"))
    
    # Re-run overview plots on subset
    p_heat_top <- netVisual_heatmap(cc_top, font.size = 13, font.size.title = 15)
    png(file.path(top_dir, "01_heatmap_top_interactors.png"),
        width = 8, height = 7, res = DPI_SETTING, units = "in")
    print(p_heat_top)
    dev.off()
    
    png(file.path(top_dir, "02_circle_top_interactors.png"),
        width = 10, height = 10, res = DPI_SETTING, units = "in")
    par(mfrow = c(1, 2), xpd = TRUE)
    netVisual_diffInteraction(cc_top, comparison = c(1, 2),
                              weight.scale = TRUE,
                              vertex.label.cex = 1.3,
                              vertex.size.max  = 25,
                              arrow.size       = 0.9,
                              margin           = 0.1)
    netVisual_diffInteraction(cc_top, comparison = c(1, 2),
                              measure      = "weight",
                              weight.scale = TRUE,
                              vertex.label.cex = 1.3,
                              vertex.size.max  = 25,
                              arrow.size       = 0.9,
                              margin           = 0.1)
    dev.off()
    
    # Checkpoint interactions on top interactors
    checkpoint_found <- CHECKPOINT_PAIRS[
      toupper(CHECKPOINT_PAIRS) %in% toupper(rownames(cc_top@LR$LRsig))
    ]
    if (length(checkpoint_found) > 0) {
      tryCatch({
        pairLR_cp <- data.frame(interaction_name = checkpoint_found,
                                stringsAsFactors = FALSE)
        gg_cp_top <- netVisual_bubble(cc_top,
                                      pairLR.use     = pairLR_cp,
                                      comparison     = c(1, 2),
                                      remove.isolate = TRUE,
                                      color.heatmap = "Spectral",
                                      title.name     = paste0(
                                        "Checkpoints — Top Interactors\n", label),
                                      font.size = 14, angle.x = 45)
        ggsave(file.path(top_dir, "03_checkpoint_top_interactors.png"),
               gg_cp_top, width = 12, height = 8, dpi = DPI_SETTING)
      }, error = function(e) {
        message(paste("  [SKIP top checkpoint]:", e$message))
      })
    }
    
    write.csv(data.frame(cell_type = top_cts, strength = ct_strength[top_cts]),
              file.path(top_dir, "top_interactors_list.csv"), row.names = FALSE)
    
  }, error = function(e) {
    message(paste("  [SKIP top interactors re-run]:", e$message))
  })
}

# =============================================================================
# --- PART 4: RUN BROAD -------------------------------------------------------
# =============================================================================
run_broad <- MODE %in% c("broad", "both")
run_sub   <- MODE %in% c("subtype", "both")

if (run_broad) {
  message("\n\n=== BROAD CELLCHAT ANALYSIS ===")
  data$group_col_broad <- data$CellType_broad
  
  broad_results <- list()
  for (comp_name in names(COMPARISONS)) {
    message(paste("\n--- Broad comparison:", comp_name, "---"))
    result <- run_comparison(data, COMPARISONS[[comp_name]],
                             group_col = "group_col_broad",
                             out_dir   = CC_BROAD_DIR)
    broad_results[[comp_name]] <- result
    
    # Plot — always runs, regenerates plots without recomputing
    plot_comparison(result)
    
    if (!is.null(result)) {
      message("  Running top interactors re-run...")
      tryCatch(
        run_top_interactors(data, result,
                            group_col = "group_col_broad",
                            out_dir   = file.path(CC_TOP_DIR, "broad")),
        error = function(e) message(paste("  [SKIP top interactors]:", e$message))
      )
      result$merged <- NULL
    }
    gc()
  }
  message("\n=== BROAD CELLCHAT COMPLETE ===")
}

if (run_sub) {
  message("\n\n=== SUBTYPE CELLCHAT ANALYSIS ===")
  
  sub_results <- list()
  for (comp_name in names(COMPARISONS)) {
    message(paste("\n--- Subtype comparison:", comp_name, "---"))
    result <- run_comparison(data, COMPARISONS[[comp_name]],
                             group_col = "CellType",
                             out_dir   = CC_SUB_DIR)
    sub_results[[comp_name]] <- result
    
    # Plot — always runs, regenerates plots without recomputing
    plot_comparison(result)
    
    if (!is.null(result)) {
      message("  Running top interactors re-run...")
      tryCatch(
        run_top_interactors(data, result,
                            group_col = "CellType",
                            out_dir   = file.path(CC_TOP_DIR, "subtypes")),
        error = function(e) message(paste("  [SKIP top interactors]:", e$message))
      )
      result$merged <- NULL
    }
    gc()
  }
  message("\n=== SUBTYPE CELLCHAT COMPLETE ===")
}

# =============================================================================
# --- PART 6: FINAL MESSAGE ---------------------------------------------------
# =============================================================================
message(paste0(
  "\n=== CELLCHAT ANALYSIS COMPLETE ===\n",
  "  Mode:    ", MODE, "\n",
  "  Output:  ", CC_ROOT, "\n",
  "  Comparisons:\n",
  paste0("    ", names(COMPARISONS), collapse = "\n"), "\n",
  "\n  Output structure:\n",
  "    cellchat_results/broad/{comparison}/\n",
  "    cellchat_results/subtypes/{comparison}/\n",
  "    cellchat_results/top_interactors/broad/\n",
  "    cellchat_results/top_interactors/subtypes/\n"
))