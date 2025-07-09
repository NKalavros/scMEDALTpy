#!/usr/bin/env python
"""
Simple MEDALT runner script
Processes scRNA.CNV.txt and generates outputs matching the expected format
"""

import sys
import os
import numpy as np
import pandas as pd
from medalt_implementation import (
    MEDALT, LineageSpeciationAnalysis, load_scRNA_data, logger
)

def main():
    """Run MEDALT analysis on scRNA.CNV.txt"""
    
    # Input and output files
    input_file = 'scRNA.CNV.txt'
    tree_output = 'CNV.tree.txt'
    lsa_output = 'gene.LSA.txt'
    
    logger.info("Starting MEDALT analysis...")
    
    # Check if input file exists
    if not os.path.exists(input_file):
        logger.error(f"Input file {input_file} not found!")
        sys.exit(1)
    
    # Load data (no windowing to preserve individual genes)
    logger.info("Loading CNV data...")
    cnv_matrix, cell_names, gene_names = load_scRNA_data(input_file, window_size=1)
    
    logger.info(f"Loaded {len(cell_names)} cells and {len(gene_names)} genes")
    
    # Build MEDALT tree
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
            f.write(f"{parent_name}\\t{child_name}\\t{int(data['weight'])}\\n")
    
    logger.info(f"Tree saved to {tree_output}")
    
    # Run LSA
    logger.info("Running Lineage Speciation Analysis...")
    lsa = LineageSpeciationAnalysis(tree, cnv_matrix, n_permutations=50)  # Reduced for speed
    results = lsa.run_analysis(min_lineage_size=5, n_jobs=1, gene_names=gene_names)
    
    # Save results
    if not results.empty:
        results.to_csv(lsa_output, sep='\t', index=False)
        logger.info(f"LSA results saved to {lsa_output}")
        logger.info(f"Found {len(results)} significant associations")
        
        print("\\nTop 10 results:")
        print(results.head(10).to_string(index=False))
    else:
        logger.warning("No significant associations found")
        # Create empty file with headers
        empty_df = pd.DataFrame(columns=['region', 'Score', 'pvalue', 'adjustp', 'cell', 'depth', 'subtreesize', 'CNA'])
        empty_df.to_csv(lsa_output, sep='\t', index=False)
    
    logger.info("Analysis complete!")

if __name__ == '__main__':
    main()