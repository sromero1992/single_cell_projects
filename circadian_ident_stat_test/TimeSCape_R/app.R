# ============================================================================
# TimeSCape v0.2 — Shiny GUI (Seurat-native)
# Run with: shiny::runApp("app.R")
# ============================================================================

library(shiny)
library(bslib)
library(ggplot2)
library(rhandsontable)
library(pheatmap)
library(grid)
library(future)
library(future.apply)

# Source pipeline functions
source("R/estimate_phaseR.R")
source("R/run_timescape.R")
source("R/generate_heatmap.R")
source("R/plot_gene.R")

# ── UI ───────────────────────────────────────────────────────────────────────

ui <- bslib::page_sidebar(
  title = tags$span(
    tags$b("TimeSCape"), " v0.2  –  Circadian Rhythm Detection"
  ),
  theme = bslib::bs_theme(bootswatch = "flatly", base_font = "Helvetica Neue"),

  # ── TOP BAR: load data ────────────────────────────────────────────────────
  tags$head(tags$style(HTML("
    .section-header { color: #1a6ebd; font-weight: bold; margin-top: 8px; }
    .status-ok  { color: #27ae60; font-size: 0.85em; }
    .status-err { color: #c0392b; font-size: 0.85em; }
    .status-info{ color: #2980b9; font-size: 0.85em; }
    hr { margin: 6px 0; }
  "))),

  # ── SIDEBAR ───────────────────────────────────────────────────────────────
  sidebar = bslib::sidebar(
    width = 370,
    open  = "desktop",

    # ① LOAD DATA
    p(class="section-header", "① Load Seurat Object"),
    fileInput("file_seurat", NULL, accept = ".rds",
              buttonLabel = "Browse .rds…", placeholder = "No file selected",
              width = "100%"),
    uiOutput("ui_load_status"),
    hr(),

    # ② DEFINE ZT METADATA
    p(class="section-header", "② Define ZT Time Metadata"),
    uiOutput("ui_col_selectors"),   # cell-type col + ZT col pickers (appear after load)
    actionButton("btn_build_tmeta", "Build / Preview ZT Table",
                 width = "100%", class = "btn-primary btn-sm"),
    uiOutput("ui_tmeta_status"),
    hr(),

    # ③ ANALYSIS SETTINGS
    p(class="section-header", "③ Analysis Settings"),
    uiOutput("ui_celltype_sel"),
    selectInput("sel_norm", "Normalization:",
                choices = c(
                  "Library size → log1p  (recommended)" = "lib_size",
                  "Seurat NormalizedData slot"           = "seurat",
                  "None  (raw counts or pre-normalized)" = "none"
                ),
                selected = "lib_size", width = "100%"),
    helpText(style="font-size:0.78em; color:#555;",
      "lib_size: per-cell CPM×10k+log1p (safe across stages/replicates). ",
      "seurat: use slot already computed by NormalizeData()/SCTransform. ",
      "none: pass counts as-is (use if already normalized externally)."),
    checkboxInput("chk_period12", "Use 12-hr period  (default: 24-hr)", value = FALSE),
    hr(),

    # ④ RUN ANALYSIS
    p(class="section-header", "④ Run Analysis  →  writes ALL + Confident"),
    fluidRow(
      column(6, actionButton("btn_run",     "▶  Run Analysis",  width="100%", class="btn-success btn-sm")),
      column(6, actionButton("btn_run_all", "▶▶  All Types",    width="100%", class="btn-success btn-sm"))
    ),
    uiOutput("ui_run_status"),
    hr(),

    # ⑤ GENERATE HEATMAP
    p(class="section-header", "⑤ Generate Heatmap"),
    checkboxInput("chk_strict_heat", "Show confident genes only  (recommended)", value = TRUE),
    checkboxInput("chk_circ_heat",   "Restrict to core circadian gene set",       value = FALSE),
    textInput("txt_custom_label", "Custom label:", value = "", width = "100%"),
    actionButton("btn_heatmap", "Generate Heatmap", width = "100%", class = "btn-info btn-sm"),
    hr(),

    # ⑥ GENE EXPLORER
    p(class="section-header", "⑥ Gene Explorer"),
    fluidRow(
      column(7, selectInput("sel_batch_type", "Batch plots:",
                            choices = c("Confident genes" = "1",
                                        "Non-confident"   = "2",
                                        "Classical circ." = "3"),
                            width = "100%")),
      column(5, br(), actionButton("btn_batch", "Save (Batch)",
                                   width = "100%", class = "btn-secondary btn-sm"))
    ),
    selectInput("sel_gene", "Single gene:", choices = NULL, width = "100%"),
    uiOutput("ui_gene_group_sel"),   # appears only when group_col is set
    textInput("txt_gene", "Or type gene name:", value = "", width = "100%"),
    checkboxInput("chk_scdata", "Overlay single-cell data", value = TRUE),
    radioButtons("rad_style", "SC style:",
                 choices = c("Violin" = "violin", "Dots" = "dots"),
                 selected = "dots", inline = TRUE),
    actionButton("btn_plot_gene", "Plot Single Gene",
                 width = "100%", class = "btn-warning btn-sm")
  ),

  # ── MAIN PANEL ────────────────────────────────────────────────────────────
  bslib::layout_columns(
    col_widths = 12,
    bslib::card(
      bslib::card_header(
        class = "d-flex justify-content-between align-items-center",
        span("Plot"),
        span(
          downloadButton("btn_save_png", "💾 PNG",  class = "btn-outline-secondary btn-sm"),
          downloadButton("btn_save_pdf", "💾 PDF",  class = "btn-outline-secondary btn-sm"),
          downloadButton("btn_save_svg", "💾 SVG",  class = "btn-outline-secondary btn-sm")
        )
      ),
      plotOutput("main_plot", height = "560px")
    )
  )
)


# ── SERVER ───────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  rv <- reactiveValues(
    seurat_obj   = NULL,
    seurat_path  = NULL,
    tmeta        = NULL,
    outdir       = getwd(),
    T1           = NULL,    # last analysis T1 (stats table)
    T2           = NULL,    # last analysis T2 (per-ZT means)
    all_results  = NULL,    # full list from run_timescape()
    current_plot = NULL,    # ggplot or "heatmap" tag for save
    current_celltype = NULL,
    status_run   = "Not yet run"
  )

  # ── ① LOAD SEURAT RDS ─────────────────────────────────────────────────────
  output$ui_load_status <- renderUI({ NULL })

  observeEvent(input$file_seurat, {
    req(input$file_seurat)
    withProgress(message = "Loading Seurat object…", value = 0.5, {
      tryCatch({
        obj <- readRDS(input$file_seurat$datapath)
        if (!inherits(obj, "Seurat"))
          stop("The .rds file does not contain a Seurat object.")
        rv$seurat_obj  <- obj
        rv$seurat_path <- input$file_seurat$datapath
        rv$outdir      <- dirname(input$file_seurat$datapath)

        ncells <- ncol(obj); ngenes <- nrow(obj)
        meta_cols <- colnames(obj@meta.data)
        output$ui_load_status <- renderUI(
          tags$p(class = "status-ok",
                 sprintf("✓ Loaded: %d cells × %d genes | %d metadata columns",
                         ncells, ngenes, length(meta_cols)))
        )

        # Set parallel workers (use half the cores, min 2)
        nw <- max(2, floor(parallel::detectCores() / 2))
        future::plan(future::multisession, workers = nw)

      }, error = function(e) {
        output$ui_load_status <- renderUI(
          tags$p(class = "status-err", paste("✗", e$message))
        )
      })
    })
  })

  # ── ② COLUMN SELECTORS (appear once Seurat object is loaded) ─────────────
  output$ui_col_selectors <- renderUI({
    req(rv$seurat_obj)
    meta_cols <- colnames(rv$seurat_obj@meta.data)
    # Guess defaults
    ct_default <- if ("cell_type"    %in% meta_cols) "cell_type"
                  else if ("celltype" %in% meta_cols) "celltype"
                  else meta_cols[1]
    zt_default <- if ("ZT_str"  %in% meta_cols) "ZT_str"
                  else if ("ZT"  %in% meta_cols) "ZT"
                  else if (any(grepl("(?i)^ZT", meta_cols, perl=TRUE)))
                       meta_cols[grepl("(?i)^ZT", meta_cols, perl=TRUE)][1]
                  else meta_cols[1]
    tagList(
      selectInput("sel_celltype_col", "Cell type column:",
                  choices = meta_cols, selected = ct_default, width = "100%"),
      selectInput("sel_zt_col", "ZT time column:",
                  choices = meta_cols, selected = zt_default, width = "100%"),
      selectInput("sel_group_col", "Group by column  (optional):",
                  choices = c("None", meta_cols), selected = "None", width = "100%"),
      helpText(style="font-size:0.78em; color:#555;",
        "Split analysis by a 2nd variable (e.g. cancer stage, replicate, treatment). ",
        "Each cell-type × group pair is analysed independently and saved to its own folder.")
    )
  })

  output$ui_tmeta_status <- renderUI({ NULL })

  # ── Helper: active group column (NULL when "None") ─────────────────────────
  active_group_col <- reactive({
    gc <- input$sel_group_col
    if (is.null(gc) || gc == "None" || !isTruthy(gc)) NULL else gc
  })

  # ── Group value picker in Gene Explorer (only visible when group_col set) ──
  output$ui_gene_group_sel <- renderUI({
    req(rv$seurat_obj)
    gc <- active_group_col()
    if (is.null(gc)) return(NULL)
    grp_vals <- sort(unique(as.character(rv$seurat_obj@meta.data[[gc]])))
    selectInput("sel_gene_group",
                paste0(gc, " (group value to plot):"),
                choices = grp_vals, width = "100%")
  })

  # ── Build / Preview ZT table ──────────────────────────────────────────────
  observeEvent(input$btn_build_tmeta, {
    req(rv$seurat_obj, input$sel_zt_col)
    tryCatch({
      tmeta <- build_tmeta_from_seurat(rv$seurat_obj, input$sel_zt_col)
      rv$tmeta <- tmeta

      # Show editable preview in a modal
      showModal(modalDialog(
        title   = "ZT Time Metadata  –  verify & adjust if needed",
        size    = "m",
        p("Numeric ZT hours are auto-parsed from the ZT string column. ",
          "Edit ZT_times if your format differs from 'ZTXX'."),
        rhandsontable::rHandsontableOutput("modal_tmeta_tbl"),
        footer = tagList(
          actionButton("btn_tmeta_save", "✔  Confirm & Close", class = "btn-success"),
          modalButton("Cancel")
        )
      ))
    }, error = function(e) {
      output$ui_tmeta_status <- renderUI(
        tags$p(class = "status-err", paste("✗", e$message))
      )
    })
  })

  output$modal_tmeta_tbl <- rhandsontable::renderRHandsontable({
    req(rv$tmeta)
    rhandsontable::rhandsontable(rv$tmeta, stretchH = "all", rowHeaders = NULL)
  })

  observeEvent(input$btn_tmeta_save, {
    updated <- rhandsontable::hot_to_r(input$modal_tmeta_tbl)
    if (!is.null(updated)) rv$tmeta <- updated
    removeModal()
    n <- nrow(rv$tmeta)
    output$ui_tmeta_status <- renderUI(
      tags$p(class = "status-ok",
             sprintf("✓ Tmeta set: %d time points (%s)", n,
                     paste(rv$tmeta$zt_str, collapse=", ")))
    )
    # Populate cell-type dropdown
    ct_vals <- sort(unique(as.character(
      rv$seurat_obj@meta.data[[input$sel_celltype_col]])))
    updateSelectInput(session, "sel_celltype_run", choices = ct_vals)
  })

  # ── ③ CELL TYPE DROPDOWN (populated after tmeta is confirmed) ─────────────
  output$ui_celltype_sel <- renderUI({
    selectInput("sel_celltype_run", "Cell Type:",
                choices = if (!is.null(rv$seurat_obj))
                  sort(unique(as.character(
                    rv$seurat_obj@meta.data[[
                      if (!is.null(input$sel_celltype_col)) input$sel_celltype_col
                      else colnames(rv$seurat_obj@meta.data)[1]
                    ]])))
                else "Load data first",
                width = "100%")
  })

  output$ui_run_status <- renderUI({
    tags$p(class = "status-info", rv$status_run)
  })

  # ── ④ RUN ANALYSIS (single cell type) ─────────────────────────────────────
  observeEvent(input$btn_run, {
    req(rv$seurat_obj, rv$tmeta, input$sel_celltype_run)
    ct <- input$sel_celltype_run
    rv$status_run <- paste("Running:", ct, "…")
    rv$current_celltype <- ct

    withProgress(message = paste("Analysing", ct, "…"), value = 0, {
      tryCatch({
        gc <- active_group_col()
        res <- run_timescape(
          seurat_obj     = rv$seurat_obj,
          celltype_col   = input$sel_celltype_col,
          zt_col         = input$sel_zt_col,
          tmeta          = rv$tmeta,
          rm_low_conf    = TRUE,
          period12       = input$chk_period12,
          custom_celltype = ct,
          plot_heat      = FALSE,
          norm_str       = input$sel_norm,
          outdir         = rv$outdir,
          group_col      = gc
        )
        rv$all_results <- res
        # When group_col active, results are keyed as "CellType_Group"
        # Use first matching key to populate gene dropdown
        matching_key <- names(res)[startsWith(names(res),
                           gsub("[^[:alnum:]_]", "_", trimws(ct)))]
        first_key <- if (length(matching_key) > 0) matching_key[1] else ct
        if (!is.null(res[[first_key]])) {
          rv$T1 <- res[[first_key]]$T1
          rv$T2 <- res[[first_key]]$T2
          n_conf <- sum(res[[first_key]]$T1$pvalue < 0.05 &
                        res[[first_key]]$T1$pvalue_corr < 0.05)
          rv$status_run <- sprintf("✓ Done: %d genes tested, %d confident → %s",
                                   nrow(rv$T1), n_conf,
                                   if (!is.null(gc)) paste0(ct, " (all groups)") else ct)
          # Populate gene dropdown with confident genes first
          conf_genes <- rv$T1$Genes[rv$T1$pvalue < 0.05 & rv$T1$pvalue_corr < 0.05]
          all_genes  <- rv$T1$Genes
          updateSelectInput(session, "sel_gene",
                            choices = c(conf_genes, setdiff(all_genes, conf_genes)))
        }
      }, error = function(e) {
        rv$status_run <- paste("✗ Error:", e$message)
      })
    })
  })

  # ── RUN ALL CELL TYPES ─────────────────────────────────────────────────────
  observeEvent(input$btn_run_all, {
    req(rv$seurat_obj, rv$tmeta)
    rv$status_run <- "Running all cell types…"
    withProgress(message = "Running all cell types…", value = 0, {
      tryCatch({
        res <- run_timescape(
          seurat_obj   = rv$seurat_obj,
          celltype_col = input$sel_celltype_col,
          zt_col       = input$sel_zt_col,
          tmeta        = rv$tmeta,
          rm_low_conf  = TRUE,
          period12     = input$chk_period12,
          plot_heat    = TRUE,
          norm_str     = input$sel_norm,
          outdir       = rv$outdir,
          group_col    = active_group_col()
        )
        rv$all_results <- res
        rv$status_run  <- sprintf("✓ All types done (%d combinations)", length(res))
      }, error = function(e) {
        rv$status_run <- paste("✗ Error:", e$message)
      })
    })
  })

  # ── ⑤ GENERATE HEATMAP ────────────────────────────────────────────────────
  observeEvent(input$btn_heatmap, {
    req(rv$T1, rv$current_celltype)
    ct      <- rv$current_celltype
    ct_safe <- gsub("[^[:alnum:]_]", "_", trimws(ct))
    gc      <- active_group_col()
    grp_safe <- if (!is.null(gc) && isTruthy(input$sel_gene_group))
                  gsub("[^[:alnum:]_]", "_", trimws(input$sel_gene_group))
                else NULL
    combo_name <- if (!is.null(grp_safe)) paste0(ct_safe, "_", grp_safe) else ct_safe
    ct_outdir <- file.path(rv$outdir, combo_name)
    per_label <- if (input$chk_period12) "_period_12_" else "_period_24_"

    tryCatch({
      ph <- generate_heatmap(
        celltype    = combo_name,
        strict      = input$chk_strict_heat,
        custom_name = input$txt_custom_label,
        circ        = input$chk_circ_heat,
        period12    = input$chk_period12,
        outdir      = ct_outdir,
        return_obj  = TRUE   # returns pheatmap object for display
      )
      output$main_plot <- renderPlot({
        grid::grid.newpage()
        grid::grid.draw(ph$gtable)
      }, bg = "white")
      rv$current_plot <- list(type = "heatmap", ph = ph)
    }, error = function(e) {
      showNotification(paste("Heatmap error:", e$message), type = "error", duration = 8)
    })
  })

  # ── ⑥ BATCH GENE PLOTS ────────────────────────────────────────────────────
  observeEvent(input$btn_batch, {
    req(rv$T1, rv$tmeta, rv$current_celltype)
    ct      <- rv$current_celltype
    ct_safe <- gsub("[^[:alnum:]_]", "_", trimws(ct))
    gc      <- active_group_col()
    grp_safe <- if (!is.null(gc) && isTruthy(input$sel_gene_group))
                  gsub("[^[:alnum:]_]", "_", trimws(input$sel_gene_group))
                else NULL
    combo_name <- if (!is.null(grp_safe)) paste0(ct_safe, "_", grp_safe) else ct_safe
    ct_outdir <- file.path(rv$outdir, combo_name)
    withProgress(message = "Saving batch plots…", value = 0, {
      tryCatch({
        save_batch_plots(
          tmeta      = rv$tmeta,
          cust_cells = combo_name,
          plot_type  = as.integer(input$sel_batch_type),
          period12   = input$chk_period12,
          outdir     = ct_outdir
        )
        showNotification("Batch plots saved!", type = "message")
      }, error = function(e) {
        showNotification(paste("Batch error:", e$message), type = "error", duration = 8)
      })
    })
  })

  # ── PLOT SINGLE GENE ──────────────────────────────────────────────────────
  observeEvent(input$btn_plot_gene, {
    req(rv$T1, rv$tmeta, rv$current_celltype)

    gene <- if (nchar(trimws(input$txt_gene)) > 0) trimws(input$txt_gene)
            else input$sel_gene
    if (is.null(gene) || gene == "") {
      showNotification("Select or type a gene name.", type = "warning"); return()
    }

    ct      <- rv$current_celltype
    ct_safe <- gsub("[^[:alnum:]_]", "_", trimws(ct))
    gc      <- active_group_col()
    grp_safe <- if (!is.null(gc) && isTruthy(input$sel_gene_group))
                  gsub("[^[:alnum:]_]", "_", trimws(input$sel_gene_group))
                else NULL
    combo_name <- if (!is.null(grp_safe)) paste0(ct_safe, "_", grp_safe) else ct_safe
    ct_outdir <- file.path(rv$outdir, combo_name)

    tryCatch({
      p <- plot_gene_single(
        tmeta         = rv$tmeta,
        cust_cells    = combo_name,
        period12      = input$chk_period12,
        cust_gene     = gene,
        print_scdata  = input$chk_scdata,
        sce           = if (input$chk_scdata) rv$seurat_obj else NULL,
        celltype_col  = input$sel_celltype_col,
        zt_col        = input$sel_zt_col,
        use_violin    = (input$rad_style == "violin"),
        outdir        = ct_outdir
      )
      output$main_plot <- renderPlot({ print(p) }, bg = "white")
      rv$current_plot  <- list(type = "ggplot", p = p)
    }, error = function(e) {
      showNotification(paste("Plot error:", e$message), type = "error", duration = 8)
    })
  })

  # ── SAVE FIGURE DOWNLOADS ─────────────────────────────────────────────────
  make_download <- function(fmt) {
    downloadHandler(
      filename = function() paste0("TimeSCape_", format(Sys.time(),"%Y%m%d_%H%M%S"), ".", fmt),
      content  = function(file) {
        req(rv$current_plot)
        cp <- rv$current_plot
        if (cp$type == "ggplot") {
          ggplot2::ggsave(file, cp$p, width=10, height=6.5, dpi=200,
                          device=fmt, bg="white")
        } else {
          # heatmap
          if (fmt == "pdf") {
            pdf(file, width=10, height=8, bg="white")
          } else if (fmt == "svg") {
            svg(file, width=10, height=8, bg="white")
          } else {
            png(file, width=1600, height=1300, res=160, bg="white")
          }
          grid::grid.newpage(); grid::grid.draw(cp$ph$gtable)
          dev.off()
        }
      }
    )
  }
  output$btn_save_png <- make_download("png")
  output$btn_save_pdf <- make_download("pdf")
  output$btn_save_svg <- make_download("svg")

  # ── WELCOME PLOT ──────────────────────────────────────────────────────────
  output$main_plot <- renderPlot({
    ggplot2::ggplot() +
      ggplot2::annotate("text", x=0.5, y=0.55,
                        label="TimeSCape v0.2",
                        size=9, fontface="bold", colour="#1a6ebd") +
      ggplot2::annotate("text", x=0.5, y=0.42,
                        label="Load a Seurat .rds  →  pick metadata columns  →  run analysis",
                        size=4.5, colour="#555555") +
      ggplot2::theme_void() +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill="white", colour=NA))
  }, bg = "white")
}

shinyApp(ui, server)
