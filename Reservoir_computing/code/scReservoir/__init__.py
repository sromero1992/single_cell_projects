"""
scReservoir: Reservoir Computing Framework for Single-Cell RNA-seq GRN Inference

A Python implementation of reservoir computing (Echo State Networks) adapted
for single-cell RNA-seq data to infer gene regulatory networks, detect
attractor states, and model developmental landscapes.
"""

__version__ = '0.1.0'
__author__ = 'scReservoir Contributors'
__license__ = 'MIT'

from .reservoir import ScReservoir, ScGraphReservoir
from .grn import ScReservoirGRN
from .landscape import ScReservoirLandscape
from .utils import preprocess, order_by_pseudotime, compute_gene_velocity
from .plotting import (
    plot_grn_heatmap,
    plot_energy_landscape,
    plot_fate_probabilities,
    plot_attractor_genes
)

__all__ = [
    'ScReservoir',
    'ScGraphReservoir',
    'ScReservoirGRN',
    'ScReservoirLandscape',
    'preprocess',
    'order_by_pseudotime',
    'compute_gene_velocity',
    'plot_grn_heatmap',
    'plot_energy_landscape',
    'plot_fate_probabilities',
    'plot_attractor_genes'
]
