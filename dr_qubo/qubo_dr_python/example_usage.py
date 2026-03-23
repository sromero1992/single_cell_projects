#!/usr/bin/env python
# Author: Selim Romero, Texas A&M University
"""
Example usage of the qubo_dr package.

Demonstrates the Python API for running the QUBO differential co-expression
pipeline on synthetic scRNA-seq data.
"""

import numpy as np
import matplotlib.pyplot as plt
from qubo_dr import run_pipeline
from qubo_dr.pathway import get_pathway_genes, list_gobp_terms


def generate_synthetic_data(n_genes=500, n_cells_ko=150, n_cells_wt=150):
    """
    Generate synthetic scRNA-seq data for testing.

    Parameters
    ----------
    n_genes : int
        Number of genes
    n_cells_ko : int
        Number of KO condition cells
    n_cells_wt : int
        Number of WT condition cells

    Returns
    -------
    X : np.ndarray, shape (genes, cells)
        Count matrix
    g : list of str
        Gene names
    batch_id : np.ndarray
        Condition labels
    """
    np.random.seed(42)

    n_cells = n_cells_ko + n_cells_wt

    # Generate synthetic counts (Poisson-like)
    X = np.random.poisson(lam=5, size=(n_genes, n_cells)).astype(float)

    # Add condition-specific signal to first 50 genes
    X[:50, :n_cells_ko] *= 1.5  # Elevated in KO
    X[:50, n_cells_ko:] *= 0.8  # Reduced in WT

    # Add noise
    X += np.random.normal(0, 0.5, X.shape)
    X = np.clip(X, 0, None)

    # Generate gene names (Wnt pathway + random)
    wnt_genes = [
        'WNT1', 'WNT2', 'WNT3', 'WNT4', 'WNT5A', 'WNT6', 'WNT7A', 'WNT7B',
        'WNT8A', 'WNT8B', 'WNT9A', 'WNT9B', 'WNT10A', 'WNT10B', 'WNT11',
        'GSK3B', 'CTNNB1', 'TCF7L2', 'LEF1', 'APC', 'AXIN1', 'AXIN2', 'CKB',
    ]
    other_genes = [f'GENE_{i:04d}' for i in range(n_genes - len(wnt_genes))]
    g = wnt_genes + other_genes

    # Generate batch labels
    batch_id = np.array(['KO'] * n_cells_ko + ['WT'] * n_cells_wt, dtype=str)

    return X, g, batch_id


def example_1_synthetic_data():
    """Example 1: Run pipeline on synthetic data with custom gene list."""
    print("\n" + "=" * 70)
    print("Example 1: Synthetic Data with Custom Gene List")
    print("=" * 70 + "\n")

    # Generate synthetic data
    X, g, batch_id = generate_synthetic_data(n_genes=500, n_cells_ko=150, n_cells_wt=150)

    print(f"Generated synthetic data:")
    print(f"  Shape: {X.shape} (genes × cells)")
    print(f"  Genes: {len(g)}")
    print(f"  Conditions: {np.unique(batch_id)}")

    # Define custom gene list (Wnt pathway subset)
    genelist = ['WNT1', 'WNT2', 'WNT3', 'WNT5A', 'GSK3B', 'CTNNB1', 'TCF7L2', 'LEF1']

    print(f"\nRunning pipeline with custom gene list: {genelist}")

    try:
        results = run_pipeline(
            X=X,
            g=g,
            batch_id=batch_id,
            genelist=genelist,
            K=5,
            ko_label='KO',
            wt_label='WT',
            method='mnn',
            n_neighbors=20,
            num_reads=500,
            plotit=True,
            outfile='/tmp/qubo_results_example1.txt',
        )

        print(f"\nResults:")
        print(f"  Selected genes: {results['sub_g_net']}")
        print(f"  Subnetwork size: {results['n_selected']}")
        print(f"  Best energy: {results['best_energy']:.4f}")
        print(f"  Generated {len(results['figures'])} figures")

        # Save figures
        for i, (fig, ax) in enumerate(results['figures']):
            if fig is not None:
                fig.savefig(f'/tmp/figure_example1_{i}.png', dpi=100, bbox_inches='tight')
                print(f"  Saved figure {i} to /tmp/figure_example1_{i}.png")

        return results

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return None


def example_2_pathway_lookup():
    """Example 2: Fetch and explore pathways from databases."""
    print("\n" + "=" * 70)
    print("Example 2: Pathway Lookup (KEGG and GO)")
    print("=" * 70 + "\n")

    # Try KEGG (if online)
    print("Fetching KEGG pathways...")
    try:
        # Note: This requires internet access
        genes = get_pathway_genes('kegg', '04310', organism='hsa')
        print(f"  KEGG pathway 04310 (Wnt signaling): {len(genes)} genes")
        print(f"  Sample genes: {genes[:10]}")
    except Exception as e:
        print(f"  Could not fetch KEGG: {e}")

    # Try GOBP (if gseapy available)
    print("\nSearching GO Biological Process terms...")
    try:
        terms = list_gobp_terms('wnt', organism='human')
        if terms:
            term_list = list(terms.keys())[:5]
            print(f"  Found {len(terms)} terms matching 'wnt'")
            print(f"  Sample terms: {term_list}")
        else:
            print("  No terms found")
    except Exception as e:
        print(f"  Could not search GO BP: {e}")


def example_3_advanced_pipeline():
    """Example 3: Advanced pipeline with cell state integration."""
    print("\n" + "=" * 70)
    print("Example 3: Advanced Pipeline with Cell State")
    print("=" * 70 + "\n")

    # Generate synthetic data with pseudotime
    X, g, batch_id = generate_synthetic_data(n_genes=500, n_cells_ko=150, n_cells_wt=150)
    n_cells = X.shape[1]

    # Simulate pseudotime trajectory
    pseudotime = np.linspace(0, 1, n_cells)
    # Add condition-specific bias
    pseudotime[batch_id == 'KO'] *= 1.2
    pseudotime[batch_id == 'WT'] *= 0.9

    print(f"Generated synthetic data with pseudotime")
    print(f"  Pseudotime range: {pseudotime.min():.3f} - {pseudotime.max():.3f}")

    # Define gene list
    genelist = ['WNT1', 'WNT2', 'WNT3', 'WNT5A', 'GSK3B', 'CTNNB1', 'TCF7L2', 'LEF1', 'APC']

    print(f"\nRunning pipeline with cell state integration...")

    try:
        results = run_pipeline(
            X=X,
            g=g,
            batch_id=batch_id,
            genelist=genelist,
            K=6,
            ko_label='KO',
            wt_label='WT',
            method='mnn',
            n_neighbors=20,
            num_reads=500,
            use_cell_state=True,
            cell_state=pseudotime,
            plotit=True,
            outfile='/tmp/qubo_results_example3.txt',
        )

        print(f"\nResults:")
        print(f"  Selected genes: {results['sub_g_net']}")
        print(f"  Subnetwork size: {results['n_selected']}")
        print(f"  Best energy: {results['best_energy']:.4f}")

        return results

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return None


def example_4_network_analysis():
    """Example 4: Post-processing network analysis."""
    print("\n" + "=" * 70)
    print("Example 4: Network Analysis from Results")
    print("=" * 70 + "\n")

    # Run pipeline to get graph
    X, g, batch_id = generate_synthetic_data()
    genelist = ['WNT1', 'WNT2', 'WNT3', 'WNT5A', 'GSK3B', 'CTNNB1']

    try:
        results = run_pipeline(
            X=X, g=g, batch_id=batch_id,
            genelist=genelist, K=5,
            num_reads=300, plotit=False,
            outfile='/tmp/qubo_results_example4.txt',
        )

        print(f"Selected genes: {results['sub_g_net']}")

        # Analyze network
        if results['G_graph'] is not None:
            import networkx as nx

            G = results['G_graph']
            print(f"\nNetwork analysis:")
            print(f"  Nodes: {G.number_of_nodes()}")
            print(f"  Edges: {G.number_of_edges()}")

            # Degree centrality
            degree_cent = nx.degree_centrality(G)
            print(f"\n  Degree centrality:")
            for node_idx, score in sorted(degree_cent.items(), key=lambda x: x[1], reverse=True):
                gene = results['sub_g_net'][node_idx]
                print(f"    {gene}: {score:.3f}")

            # Betweenness centrality
            between_cent = nx.betweenness_centrality(G)
            print(f"\n  Betweenness centrality:")
            for node_idx, score in sorted(between_cent.items(), key=lambda x: x[1], reverse=True)[:3]:
                gene = results['sub_g_net'][node_idx]
                print(f"    {gene}: {score:.3f}")

        return results

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    """Run all examples."""
    print("\n" + "=" * 70)
    print("QUBO Differential Co-expression Package - Examples")
    print("=" * 70)

    # Example 1: Synthetic data with custom gene list
    results1 = example_1_synthetic_data()

    # Example 2: Pathway lookup
    example_2_pathway_lookup()

    # Example 3: Advanced pipeline with cell state
    results3 = example_3_advanced_pipeline()

    # Example 4: Network analysis
    results4 = example_4_network_analysis()

    print("\n" + "=" * 70)
    print("Examples Complete!")
    print("=" * 70)
    print("\nCheck /tmp/figure_example*.png and /tmp/qubo_results_*.txt for outputs")


if __name__ == '__main__':
    main()
