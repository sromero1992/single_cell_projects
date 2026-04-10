# Author: Selim Romero, Texas A&M University
"""
QUBO-based Differential RNA Co-expression Analysis for Single-Cell Data

This package implements a quantum-inspired optimization pipeline for identifying
gene subnetworks with maximum differential co-expression between conditions
(e.g., Knockout vs. Wild Type) within biological pathways.

The algorithm casts the subnetwork selection as a Quadratic Unconstrained Binary
Optimization (QUBO) problem, solvable via simulated annealing or quantum annealers.

Key features:
- Pathway gene fetching from KEGG and Gene Ontology databases
- Differential co-expression computation via Gram matrices
- Mutual Nearest Neighbor (MNN) graph construction
- QUBO assembly and solving via simulated annealing
- Visualization of co-expression networks and heatmaps
"""

__version__ = "1.0.0"
__author__ = "Selim Romero, Texas A&M University"

from qubo_dr.pipeline import run_pipeline, run_pipeline_classical
from qubo_dr.plot import (
    plot_hub_network,
    plot_gene_network,
    plot_condition_heatmaps,
    plot_coexpr_heatmap,
)

__all__ = [
    "run_pipeline",
    "run_pipeline_classical",
    "plot_hub_network",
    "plot_gene_network",
    "plot_condition_heatmaps",
    "plot_coexpr_heatmap",
]
