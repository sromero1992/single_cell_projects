"""
TimeSCape: Circadian Rhythm Detection in Single-Cell RNA-seq Data
=================================================================
Python package — v0.2

Selim Romero · Texas A&M University · ssromerogon@tamu.edu

Mirrors the MATLAB v0.2 and R (TimeSCape_R) implementations.
Statistical results are numerically equivalent across all three platforms.
"""

from .pipeline   import run_timescape, build_tmeta
from .core       import estimate_phase_r
from .normalize  import normalize_adata
from .visualize  import generate_heatmap, plot_gene_single, save_batch_plots
from .utils      import bh_adjust, wrap_acrophase
from .pathway    import (pull_genesets, phase_bin_analysis,
                         auc_score_cells, pathway_cosinor,
                         write_pathway_results)
from .grn        import select_hub_genes, plot_grn_timeseries

__version__ = "0.2.0"
__author__  = "Selim Romero"
__email__   = "ssromerogon@tamu.edu"

__all__ = [
    # Core pipeline
    "run_timescape",
    "build_tmeta",
    "estimate_phase_r",
    "normalize_adata",
    # Visualization
    "generate_heatmap",
    "plot_gene_single",
    "save_batch_plots",
    # Pathway enrichment
    "pull_genesets",
    "phase_bin_analysis",
    "auc_score_cells",
    "pathway_cosinor",
    "write_pathway_results",
    # GRN
    "select_hub_genes",
    "plot_grn_timeseries",
    # Utilities
    "bh_adjust",
    "wrap_acrophase",
]
