#!/usr/bin/env python
# Author: Selim Romero, Texas A&M University
"""
Command-line interface for QUBO differential co-expression pipeline.

Provides a user-friendly entry point with argument parsing, pathway listing,
and full pipeline execution from the command line.

Example usage:
    # Run pipeline with KEGG Wnt signaling pathway (human)
    python run_pipeline.py \\
        --data counts.h5ad \\
        --pathway-source kegg \\
        --pathway-id 04310 \\
        --organism hsa \\
        --K 30 \\
        --num-reads 2000

    # Run with GO Biological Process
    python run_pipeline.py \\
        --data counts.h5ad \\
        --pathway-source gobp \\
        --pathway-id "Wnt signaling pathway" \\
        --organism human \\
        --K 30

    # List available KEGG pathways for human
    python run_pipeline.py --list-kegg-pathways

    # Search GO BP terms by keyword
    python run_pipeline.py --list-gobp-terms wnt
"""

import argparse
import sys
import numpy as np

try:
    import anndata
    HAS_ANNDATA = True
except ImportError:
    HAS_ANNDATA = False

from qubo_dr import run_pipeline
from qubo_dr.pathway import (
    fetch_kegg_pathway_list,
    search_kegg_pathways,
    list_gobp_terms,
)


def load_data(data_path):
    """
    Load count matrix from .h5ad or .npy file.

    Parameters
    ----------
    data_path : str
        Path to .h5ad (AnnData) or .npy (NumPy) file.

    Returns
    -------
    X : np.ndarray, shape (genes, cells)
        Count matrix.
    g : list of str
        Gene names.
    batch_id : np.ndarray, shape (cells,)
        Condition labels per cell.

    Raises
    ------
    ValueError
        If file format not recognized or data not found.
    """
    if data_path.endswith('.h5ad'):
        if not HAS_ANNDATA:
            raise ImportError("anndata required for .h5ad files: pip install anndata")

        print(f"Loading data from {data_path}")
        adata = anndata.read_h5ad(data_path)

        X = adata.X.T.toarray() if hasattr(adata.X, 'toarray') else adata.X.T
        g = adata.var_names.tolist()

        if 'batch' not in adata.obs:
            raise ValueError("Data must have 'batch' column in .obs")

        batch_id = adata.obs['batch'].values

        print(f"  Shape: {X.shape} (genes x cells)")
        print(f"  Conditions: {np.unique(batch_id)}")

        return X, g, batch_id

    elif data_path.endswith('.npy'):
        print(f"Loading data from {data_path}")
        X = np.load(data_path)
        print(f"  Shape: {X.shape} (genes x cells)")

        # User must provide gene names and batch_id separately
        raise ValueError(
            ".npy files require manual data loading. "
            "Use the Python API directly: from qubo_dr import run_pipeline"
        )

    else:
        raise ValueError(f"Unknown file format: {data_path}. Use .h5ad or .npy")


def main():
    """Parse arguments and execute pipeline or utility functions."""
    parser = argparse.ArgumentParser(
        description="QUBO differential co-expression analysis for scRNA-seq data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # KEGG pathway (Wnt signaling, human)
  python run_pipeline.py --data counts.h5ad --pathway-source kegg --pathway-id 04310

  # GO Biological Process pathway
  python run_pipeline.py --data counts.h5ad --pathway-source gobp --pathway-id "Wnt signaling pathway"

  # List KEGG pathways
  python run_pipeline.py --list-kegg-pathways

  # Search GO terms
  python run_pipeline.py --list-gobp-terms wnt
        """
    )

    # Pipeline arguments
    parser.add_argument(
        '--data', type=str, default=None,
        help='Path to count matrix (.h5ad or .npy file)'
    )
    parser.add_argument(
        '--pathway-source', type=str, choices=['kegg', 'gobp'], default=None,
        help='Pathway database: kegg or gobp'
    )
    parser.add_argument(
        '--pathway-id', type=str, default=None,
        help='KEGG pathway ID (e.g., 04310) or GO term name/ID'
    )
    parser.add_argument(
        '--organism', type=str, default='human',
        help='Organism for pathway mapping. KEGG: hsa, mmu, dme, etc. GOBP: human, mouse'
    )
    parser.add_argument(
        '--K', type=int, default=30,
        help='Target subnetwork size (default 30)'
    )
    parser.add_argument(
        '--ko-label', type=str, default='KO',
        help='KO condition label in batch column (default KO)'
    )
    parser.add_argument(
        '--wt-label', type=str, default='WT',
        help='WT condition label in batch column (default WT)'
    )
    parser.add_argument(
        '--method', type=str, choices=['knn', 'mnn'], default='mnn',
        help='Graph construction method: knn or mnn (default mnn)'
    )
    parser.add_argument(
        '--neighbors', type=int, default=30,
        help='Number of neighbors for KNN/MNN (default 30)'
    )
    parser.add_argument(
        '--num-reads', type=int, default=1000,
        help='SA iterations (default 1000)'
    )
    parser.add_argument(
        '--no-plot', action='store_true',
        help='Suppress plot generation'
    )
    parser.add_argument(
        '--outfile', type=str, default='qubo_genes_solution.txt',
        help='Output filename for selected genes (default qubo_genes_solution.txt)'
    )

    # Utility arguments
    parser.add_argument(
        '--list-kegg-pathways', action='store_true',
        help='List available KEGG pathways for organism and exit'
    )
    parser.add_argument(
        '--list-gobp-terms', type=str, default=None, metavar='KEYWORD',
        help='Search GO BP terms by keyword and exit'
    )

    args = parser.parse_args()

    # Handle utility commands
    if args.list_kegg_pathways:
        organism = _infer_kegg_organism(args.organism)
        print(f"KEGG pathways for {organism}:")
        try:
            pathways = fetch_kegg_pathway_list(organism=organism)
            for path_id, path_name in sorted(pathways.items())[:50]:
                print(f"  {path_id}: {path_name}")
            if len(pathways) > 50:
                print(f"  ... and {len(pathways) - 50} more")
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
        sys.exit(0)

    if args.list_gobp_terms:
        print(f"GO Biological Process terms matching '{args.list_gobp_terms}':")
        try:
            terms = list_gobp_terms(args.list_gobp_terms, organism=args.organism)
            for term_name in sorted(terms.keys())[:50]:
                print(f"  {term_name}")
            if len(terms) > 50:
                print(f"  ... and {len(terms) - 50} more")
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
        sys.exit(0)

    # Run pipeline
    if args.data is None:
        print("Error: --data required for pipeline execution", file=sys.stderr)
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    if args.pathway_source is None or args.pathway_id is None:
        print(
            "Error: --pathway-source and --pathway-id required",
            file=sys.stderr
        )
        sys.exit(1)

    try:
        X, g, batch_id = load_data(args.data)
    except Exception as e:
        print(f"Error loading data: {e}", file=sys.stderr)
        sys.exit(1)

    # Infer KEGG organism code if needed
    organism = args.organism
    if args.pathway_source == 'kegg':
        organism = _infer_kegg_organism(args.organism)

    print("\n" + "=" * 60)
    print("QUBO Differential Co-expression Pipeline")
    print("=" * 60 + "\n")

    try:
        results = run_pipeline(
            X=X,
            g=g,
            batch_id=batch_id,
            genelist=None,
            K=args.K,
            source=args.pathway_source,
            pathway_id=args.pathway_id,
            organism=organism,
            ko_label=args.ko_label,
            wt_label=args.wt_label,
            method=args.method,
            n_neighbors=args.neighbors,
            num_reads=args.num_reads,
            plotit=not args.no_plot,
            outfile=args.outfile,
        )

        print("\n" + "=" * 60)
        print(f"Results saved to: {args.outfile}")
        print(f"Selected genes: {results['n_selected']}/{results['n_pathway_genes']}")
        print("=" * 60)

        if not args.no_plot:
            print("\nTo view plots, save them with:")
            print("  import matplotlib.pyplot as plt")
            print("  for fig, ax in results['figures']: fig.savefig(...)")

    except Exception as e:
        print(f"Pipeline error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


def _infer_kegg_organism(organism):
    """
    Convert human-readable organism names to KEGG codes.

    Parameters
    ----------
    organism : str
        Organism name ('human', 'mouse', etc.) or KEGG code ('hsa', 'mmu', etc.).

    Returns
    -------
    kegg_code : str
        KEGG organism code.
    """
    mapping = {
        'human': 'hsa',
        'mouse': 'mmu',
        'fly': 'dme',
        'worm': 'cel',
        'yeast': 'sce',
        'rat': 'rno',
    }

    organism_lower = organism.lower()
    return mapping.get(organism_lower, organism_lower)


if __name__ == '__main__':
    main()
