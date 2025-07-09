#!/usr/bin/env python
"""
Optimized MEDALT runner script
Implements all improvements based on original algorithm
"""

import sys
import os
import numpy as np
import pandas as pd
from medalt_implementation import (
    MEDALT, LineageSpeciationAnalysis, load_scRNA_data, logger
)
from medalt_visualizations_fixed import (
    MEDALTTreeVisualizer, MEDALTHeatmapVisualizer, MEDALTManhattanPlot
)

def main():
    """Run optimized MEDALT analysis"""
    
    # Input and output files
    input_file = 'scRNA.CNV.txt'
    tree_output = 'CNV.tree.txt'
    lsa_output = 'gene.LSA.txt'
    
    logger.info("Starting optimized MEDALT analysis...")
    
    # Check if input file exists
    if not os.path.exists(input_file):
        logger.error(f"Input file {input_file} not found!")
        sys.exit(1)
    
    # Load data (no windowing to preserve individual genes)
    logger.info("Loading CNV data...")
    cnv_matrix, cell_names, gene_names = load_scRNA_data(input_file, window_size=1)
    
    logger.info(f"Loaded {len(cell_names)} cells and {len(gene_names)} genes")
    
    # Build MEDALT tree with optimized algorithm
    logger.info("Building MEDALT tree...")
    medalt = MEDALT(cnv_matrix)
    tree = medalt.build_tree()
    
    # Save tree in correct format
    logger.info("Saving tree...")
    with open(tree_output, 'w') as f:
        f.write("from\tto\tdist\n")
        for parent, child, data in tree.edges(data=True):
            parent_name = parent if parent == 'root' else cell_names[parent]
            child_name = cell_names[child]
            f.write(f"{parent_name}\t{child_name}\t{int(data['weight'])}\n")
    
    logger.info(f"Tree saved to {tree_output}")
    
    # Run LSA with optimized permutation tests
    logger.info("Running optimized Lineage Speciation Analysis...")
    lsa = LineageSpeciationAnalysis(tree, cnv_matrix, 
                                   n_permutations=100,  # Reasonable for demo
                                   data_type='RNA')
    results = lsa.run_analysis(min_lineage_size=5, n_jobs=1, gene_names=gene_names)
    
    # Save LSA results
    if not results.empty:
        results.to_csv(lsa_output, sep='\t', index=False)
        logger.info(f"LSA results saved to {lsa_output}")
        logger.info(f"Found {len(results)} significant associations")
        
        # Show top results
        print("\nTop 10 LSA results:")
        print(results.head(10)[['region', 'Score', 'pvalue', 'adjustp', 'cell', 'CNA']].to_string(index=False))
        
        # Show summary statistics
        n_significant = len(results[results['adjustp'] < 0.05])
        n_amp = len(results[results['CNA'] == 'AMP'])
        n_del = len(results[results['CNA'] == 'DEL'])
        print(f"\nSummary:")
        print(f"  Total associations: {len(results)}")
        print(f"  Significant (FDR < 0.05): {n_significant}")
        print(f"  Amplifications: {n_amp}")
        print(f"  Deletions: {n_del}")
        
    else:
        logger.warning("No significant associations found")
        # Create empty file with headers
        empty_df = pd.DataFrame(columns=['region', 'Score', 'pvalue', 'adjustp', 'cell', 'depth', 'subtreesize', 'CNA'])
        empty_df.to_csv(lsa_output, sep='\t', index=False)
        results = empty_df
    
    # Generate visualizations
    logger.info("Generating visualizations...")
    
    # Initialize visualizer
    visualizer = MEDALTTreeVisualizer(tree, cnv_matrix, cell_names, gene_names, results)
    
    # Create tree visualization
    logger.info("Creating tree visualization...")
    visualizer.plot_tree('singlecell.tree.pdf', figsize=(15, 12))
    
    # Create LSA tree visualization
    if not results.empty:
        logger.info("Creating LSA tree visualization...")
        visualizer.plot_lsa_tree('LSA.tree.pdf', figsize=(15, 12))
    
    # Create CNV heatmap
    logger.info("Creating CNV heatmap...")
    heatmap_viz = MEDALTHeatmapVisualizer(cnv_matrix, cell_names, gene_names, tree)
    heatmap_viz.plot_cnv_heatmap('cnv_heatmap.pdf', figsize=(20, 12))
    
    # Create Manhattan plot
    if not results.empty:
        logger.info("Creating Manhattan plot...")
        manhattan_viz = MEDALTManhattanPlot(results, gene_names)
        manhattan_viz.plot_manhattan('lsa_manhattan.pdf', figsize=(15, 8))
    
    logger.info("Analysis complete!")
    
    # Print comparison with desired outputs
    print("\n" + "="*60)
    print("COMPARISON WITH DESIRED OUTPUTS")
    print("="*60)
    
    # Check if outputs exist
    if os.path.exists('desired_output/CNV.tree.txt'):
        print("\nTree structure comparison:")
        print(f"  Generated tree edges: {tree.number_of_edges()}")
        
        # Count desired edges
        with open('desired_output/CNV.tree.txt', 'r') as f:
            desired_edges = sum(1 for line in f) - 1  # subtract header
        print(f"  Desired tree edges: {desired_edges}")
    
    if os.path.exists('desired_output/gene.LSA.txt'):
        print("\nLSA results comparison:")
        desired_lsa = pd.read_csv('desired_output/gene.LSA.txt', sep='\t')
        print(f"  Generated LSA results: {len(results)}")
        print(f"  Desired LSA results: {len(desired_lsa)}")
        
        if not results.empty and not desired_lsa.empty:
            # Compare significant results
            gen_sig = len(results[results['adjustp'] < 0.05])
            des_sig = len(desired_lsa[desired_lsa['adjustp'] < 0.05])
            print(f"  Generated significant: {gen_sig}")
            print(f"  Desired significant: {des_sig}")
            
            # Compare genes
            gen_genes = set(results['region'].unique())
            des_genes = set(desired_lsa['region'].unique())
            overlap = len(gen_genes.intersection(des_genes))
            print(f"  Gene overlap: {overlap}/{len(des_genes)} ({overlap/len(des_genes)*100:.1f}%)")
    
    print("\nFiles generated:")
    print(f"  - {tree_output}")
    print(f"  - {lsa_output}")
    print(f"  - singlecell.tree.pdf")
    if not results.empty:
        print(f"  - LSA.tree.pdf")
        print(f"  - lsa_manhattan.pdf")
    print(f"  - cnv_heatmap.pdf")

if __name__ == '__main__':
    main()