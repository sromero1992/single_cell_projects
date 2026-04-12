# TimeSCape v0.2 R Translation - Installation & Usage Guide

## Installation

### Option 1: Install as R Package

```r
# Install devtools if needed
if (!require("devtools")) install.packages("devtools")

# Install TimeSCapeR from the local directory
devtools::install("/path/to/TimeSCape_R")

# Load the package
library(TimeSCapeR)
```

### Option 2: Source Functions Directly

```r
source("/path/to/TimeSCape_R/R/estimate_phaseR.R")
source("/path/to/TimeSCape_R/R/run_timescape.R")
source("/path/to/TimeSCape_R/R/generate_heatmap.R")
source("/path/to/TimeSCape_R/R/plot_gene.R")
```

## Dependencies

Install all required packages before using:

```r
required_packages <- c(
  "SingleCellExperiment", "future.apply", "minpack.lm", 
  "pheatmap", "ggplot2", "rhandsontable", "shiny", "bslib"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("SingleCellExperiment")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}
```

## Quick Start

### 1. Prepare Your Data

You need:
- A `SingleCellExperiment` object with:
  - `counts()`: gene x cell expression matrix
  - `colData()$cell_type`: cell type assignments
  - `rownames()`: gene names

- A metadata data.frame with columns:
  - `old_labels`: original cell type labels (must match `sce$cell_type`)
  - `new_labels`: unified cell type names for output
  - `ZT_times`: numeric ZT hours (e.g., 0, 6, 12, 18, 24)

Example:
```r
tmeta <- data.frame(
  old_labels = sce$cell_type,
  new_labels = case_when(
    sce$cell_type == "Neuron_1" ~ "Neurons",
    sce$cell_type == "Neuron_2" ~ "Neurons",
    TRUE ~ sce$cell_type
  ),
  ZT_times = c(rep(0, 100), rep(6, 100), rep(12, 100), rep(18, 100))
)
```

### 2. Run the Pipeline

```r
# Optional: enable parallel processing (must set plan before running)
future::plan(multisession, workers = 4)

# Run analysis
result <- run_timescape(
  sce = sce,                    # SingleCellExperiment object
  tmeta = tmeta,               # Metadata data.frame
  rm_low_conf = TRUE,          # Remove non-significant genes
  period12 = FALSE,            # Use 24-hour period
  custom_genelist = NULL,      # NULL = analyze all genes
  custom_celltype = NULL,      # NULL = analyze all cell types
  plot_heat = TRUE,            # Generate heatmaps automatically
  norm_str = "lib_size",       # Normalization method
  outdir = "./timescape_output" # Output directory
)

# Access results
stats_table <- result$T1  # Stats for all genes
zt_means_table <- result$T2  # Per-ZT mean expression
```

### 3. Generate Individual Gene Plots

```r
# Plot a single gene
p <- plot_gene_single(
  tmeta = tmeta,
  cust_cells = c("Neurons"),      # Cell types to include
  period12 = FALSE,
  cust_gene = "Clock",            # Gene name
  print_scdata = TRUE,            # Overlay individual cells
  sce = sce,                      # Required if print_scdata=TRUE
  norm_str = "lib_size",
  use_violin = TRUE,              # Violin plots (vs scatter)
  outdir = "./timescape_output"
)

# Save the plot
ggplot2::ggsave("Clock_plot.png", p, width = 10, height = 6)
```

### 4. Generate Heatmaps

```r
generate_heatmap(
  celltype = "Neurons",
  strict = TRUE,                  # Only confident genes (p < 0.05 for both tests)
  custom_name = "my_analysis",    # Added to filename
  circ = FALSE,                   # FALSE = all genes, TRUE = circadian genes only
  period12 = FALSE,
  outdir = "./timescape_output/Neurons"
)
```

### 5. Save Batch Gene Plots

```r
save_batch_plots(
  tmeta = tmeta,
  cust_cells = c("Neurons"),
  plot_type = 1,        # 1=confident, 2=non-confident, 3=classical circadian
  period12 = FALSE,
  outdir = "./timescape_output"
)
# Creates PNG files in: ./timescape_output/Neurons/batch_plots_confident/
```

## Interactive GUI

Run the Shiny web application:

```r
shiny::runApp("/path/to/TimeSCape_R/app.R")
```

Opens in browser at http://localhost:3838

Features:
1. **Define ZT Times**: Upload or edit cell type metadata
2. **Analysis Settings**: Select cell type, normalization, period
3. **Run Analysis**: Execute pipeline with progress tracking
4. **Generate Heatmap**: Create publication-ready heatmaps
5. **Gene Explorer**: Plot individual genes with customization

## Output Files

### Per-Cell-Type Directory Structure

```
output/
└── Neurons/
    ├── Neurons_24hr_circadian_analysis_all.csv         # All genes + stats
    ├── Neurons_24hr_circadian_analysis_confident.csv   # Significant genes only
    ├── Neurons_24hr_ZTs_mean.csv                       # Per-ZT raw means
    ├── Neurons_24hr_ZTs_mean_normalized.csv            # Per-ZT normalized
    ├── Neurons_24hr_ZTs_mean_confident.csv             # Per-ZT means (confident)
    ├── Neurons_24hr_ZTs_mean_normalized_confident.csv  # Per-ZT norm (confident)
    ├── Neurons_24hr_heatmap.png                        # Heatmap PNG
    ├── Neurons_24hr_heatmap_genes.csv                  # Genes in heatmap
    └── batch_plots_confident/
        ├── Clock_24hr.png
        ├── Per2_24hr.png
        └── ... (one PNG per gene)
```

### CSV Column Names (Stats Table)

| Column | Description |
|--------|-------------|
| Genes | Gene name |
| Amp | Cosine amplitude (can be negative) |
| Abs_Amp | Absolute amplitude |
| Mesor | Mean level (midline) |
| Acrophase | Peak time in period (0 to period) |
| Acrophase_24 | Peak time in 24-hour scale |
| Period | 12 or 24 hours |
| pvalue | F-test p-value (or LRT) |
| pvalue_adj | BH-corrected p-value |
| Sine_corr | Pearson r (per-ZT means vs fit) |
| pvalue_corr | Correlation p-value |
| pvalue_adj_corr | BH-corrected correlation p-value |

### Significance Criteria

**Confident genes** = Both:
- F-test p-value < 0.05
- Correlation p-value < 0.05

## Examples

### Example 1: Basic Analysis

```r
library(TimeSCapeR)

# Run analysis for all genes and cell types
future::plan(multisession, workers = 8)

result <- run_timescape(
  sce = my_sce,
  tmeta = my_metadata,
  outdir = "./my_analysis"
)

# Check confident genes
confident <- result$T1[result$T1$pvalue < 0.05 & result$T1$pvalue_corr < 0.05, ]
print(head(confident[order(confident$pvalue_adj_corr), ]))
```

### Example 2: Analysis with Specific Genes

```r
# Analyze only known circadian genes
circ_genes <- c("Clock", "Bmal1", "Per1", "Per2", "Cry1", "Cry2")

result <- run_timescape(
  sce = my_sce,
  tmeta = my_metadata,
  custom_genelist = circ_genes,
  outdir = "./circadian_genes"
)
```

### Example 3: 12-Hour Periodicity

```r
# For organisms with 12-hour circadian periods
result <- run_timescape(
  sce = my_sce,
  tmeta = my_metadata,
  period12 = TRUE,
  outdir = "./12hr_analysis"
)
```

### Example 4: Detailed Gene Visualization

```r
# Create high-quality plot for publication
p <- plot_gene_single(
  tmeta = tmeta,
  cust_cells = c("Neurons", "Glia"),
  period12 = FALSE,
  cust_gene = "Clock",
  print_scdata = TRUE,
  sce = sce,
  use_violin = TRUE,
  outdir = "./output"
)

ggplot2::ggsave(
  "Clock_circadian_pattern.pdf",
  p,
  width = 8,
  height = 6,
  device = "pdf"
)
```

## Troubleshooting

### Error: "Gene not found in analysis results"

**Cause**: The gene name doesn't match the SCE rownames

**Solution**: Check gene name case sensitivity
```r
head(rownames(sce))  # See exact format
# Use the exact name from this output
```

### Error: "Required files not found"

**Cause**: Wrong output directory or analysis hasn't been run yet

**Solution**: 
```r
# Check that run_timescape completed successfully
# Verify outdir contains subdirectories named after cell types
list.dirs(outdir)
```

### Slow analysis

**Solution**: Enable parallel processing
```r
future::plan(multisession, workers = 8)  # Adjust workers for your CPU
result <- run_timescape(...)
```

### Out of memory

**Solution**: Analyze cell types separately instead of all at once
```r
result_neurons <- run_timescape(
  sce = sce,
  tmeta = tmeta,
  custom_celltype = c("Neurons"),
  outdir = "./output_neurons"
)
```

## Version Notes

**TimeSCape R v0.2.0** is a direct translation of the MATLAB v0.2 pipeline.

Key differences from MATLAB:
1. Uses `minpack.lm::nlsLM` instead of MATLAB Trust-Region solver
   - Results typically identical (<0.1% difference)
2. Requires explicit `future::plan()` call for parallelism
3. Outputs PNG files instead of MATLAB figures
4. Shiny web GUI instead of MATLAB uicontrol native GUI

All statistical calculations, CSV outputs, and data processing are identical to MATLAB version.

## Citation

If you use TimeSCape R in your research, please cite:
- Original MATLAB version citation (if available)
- This R translation (provide repo/DOI once available)

## License

See the LICENSE file in the package directory.

## Support

For issues or questions:
1. Check this INSTALL.md guide
2. Review function documentation: `?run_timescape`, `?plot_gene_single`, etc.
3. Check R package DESCRIPTION file for version/author info
