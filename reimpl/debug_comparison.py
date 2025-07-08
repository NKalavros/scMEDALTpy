#!/usr/bin/env python3
"""
Debug script to compare R and Python MEDALT outputs
"""

import pandas as pd
import numpy as np
import os
import sys
from collections import defaultdict

def compare_trees(r_tree_path, py_tree_path):
    """Compare tree structures between R and Python outputs."""
    print("\n=== Tree Comparison ===")
    
    # Read R tree
    r_tree = pd.read_csv(r_tree_path, sep='\t')
    print(f"R tree: {len(r_tree)} edges")
    
    # Read Python tree  
    py_tree = pd.read_csv(py_tree_path, sep='\t')
    # Rename columns to match R format
    py_tree.columns = ['from', 'to', 'dist']
    print(f"Python tree: {len(py_tree)} edges")
    
    # Compare edge weights distribution
    print("\nEdge weight statistics:")
    print(f"R - min: {r_tree['dist'].min()}, max: {r_tree['dist'].max()}, "
          f"mean: {r_tree['dist'].mean():.2f}, zeros: {(r_tree['dist'] == 0).sum()}")
    print(f"Py - min: {py_tree['dist'].min()}, max: {py_tree['dist'].max()}, "
          f"mean: {py_tree['dist'].mean():.2f}, zeros: {(py_tree['dist'] == 0).sum()}")
    
    # Find common edges
    r_edges = set(zip(r_tree['from'], r_tree['to']))
    py_edges = set(zip(py_tree['from'], py_tree['to']))
    
    common = r_edges & py_edges
    r_only = r_edges - py_edges
    py_only = py_edges - r_edges
    
    print(f"\nCommon edges: {len(common)}")
    print(f"R-only edges: {len(r_only)}")
    print(f"Python-only edges: {len(py_only)}")
    
    # Check root nodes
    r_roots = set(r_tree['from']) - set(r_tree['to'])
    py_roots = set(py_tree['from']) - set(py_tree['to'])
    
    print(f"\nR root nodes: {r_roots}")
    print(f"Python root nodes: {py_roots}")
    
    return r_tree, py_tree

def compare_binned_cnv(r_cnv_path, py_cnv_path):
    """Compare binned CNV profiles."""
    print("\n=== Binned CNV Comparison ===")
    
    # Read CNV files
    r_cnv = pd.read_csv(r_cnv_path, sep='\t', index_col=0)
    py_cnv = pd.read_csv(py_cnv_path, sep='\t', index_col=0)
    
    print(f"R CNV shape: {r_cnv.shape}")
    print(f"Python CNV shape: {py_cnv.shape}")
    
    # Compare cells
    r_cells = set(r_cnv.index)
    py_cells = set(py_cnv.index)
    common_cells = r_cells & py_cells
    
    print(f"\nCommon cells: {len(common_cells)}/{len(r_cells)}")
    
    # Compare regions
    print(f"\nR regions (first 5): {list(r_cnv.columns[:5])}")
    print(f"Python regions (first 5): {list(py_cnv.columns[:5])}")
    
    # Compare CNV values for common cells
    if common_cells:
        sample_cell = list(common_cells)[0]
        print(f"\nCNV values for cell {sample_cell}:")
        print(f"R: min={r_cnv.loc[sample_cell].min()}, max={r_cnv.loc[sample_cell].max()}, "
              f"mean={r_cnv.loc[sample_cell].mean():.2f}")
        print(f"Py: min={py_cnv.loc[sample_cell].min()}, max={py_cnv.loc[sample_cell].max()}, "
              f"mean={py_cnv.loc[sample_cell].mean():.2f}")
    
    return r_cnv, py_cnv

def compare_lsa_results(r_lsa_path, py_lsa_path):
    """Compare LSA results."""
    print("\n=== LSA Results Comparison ===")
    
    # Read LSA results
    r_lsa = pd.read_csv(r_lsa_path, sep='\t')
    py_lsa = pd.read_csv(py_lsa_path, sep='\t')
    
    print(f"R LSA: {len(r_lsa)} significant CNAs")
    print(f"Python LSA: {len(py_lsa)} significant CNAs")
    
    # Compare cells with significant CNAs
    r_cells = r_lsa['cell'].unique()
    py_cells = py_lsa['cell'].unique()
    
    print(f"\nR cells with CNAs: {len(r_cells)}")
    print(f"Python cells with CNAs: {len(py_cells)}")
    print(f"Common cells: {len(set(r_cells) & set(py_cells))}")
    
    # Compare p-value distributions
    print(f"\nR p-values: min={r_lsa['pvalue'].min():.6f}, max={r_lsa['pvalue'].max():.6f}")
    print(f"Python p-values: min={py_lsa['pvalue'].min():.6f}, max={py_lsa['pvalue'].max():.6f}")
    
    # Compare CNA types
    r_amp = (r_lsa['CNA'] == 'AMP').sum()
    r_del = (r_lsa['CNA'] == 'DEL').sum()
    py_amp = (py_lsa['CNA'] == 'AMP').sum()
    py_del = (py_lsa['CNA'] == 'DEL').sum()
    
    print(f"\nR CNAs: {r_amp} AMP, {r_del} DEL")
    print(f"Python CNAs: {py_amp} AMP, {py_del} DEL")
    
    # Compare regions
    r_regions = set(r_lsa['region'])
    py_regions = set(py_lsa['region'])
    
    print(f"\nR regions: {len(r_regions)}")
    print(f"Python regions: {len(py_regions)}")
    print(f"Common regions: {len(r_regions & py_regions)}")
    
    # Show top results
    print("\nTop 5 R results:")
    print(r_lsa[['region', 'cell', 'Score', 'pvalue', 'CNA']].head())
    
    print("\nTop 5 Python results:")
    print(py_lsa[['region', 'cell', 'Score', 'pvalue', 'CNA']].head())
    
    # Calculate consistency metrics
    print("\n=== Consistency Analysis ===")
    
    # Find matching cell-region pairs (since region names might differ)
    r_cell_regions = set(zip(r_lsa['cell'], r_lsa['region']))
    py_cell_regions = set(zip(py_lsa['cell'], py_lsa['region']))
    
    # Also check by cell only (more lenient)
    common_cells_with_cnas = set(r_cells) & set(py_cells)
    
    if common_cells_with_cnas:
        consistency_pct = (len(common_cells_with_cnas) / len(set(r_cells) | set(py_cells))) * 100
        print(f"Cell-level consistency: {consistency_pct:.1f}%")
        print(f"Cells with CNAs in both: {sorted(common_cells_with_cnas)}")
    
    # Check for matching chromosomes/CNA types in common cells
    if common_cells_with_cnas:
        chr_matches = 0
        cna_type_matches = 0
        
        for cell in common_cells_with_cnas:
            r_cell_data = r_lsa[r_lsa['cell'] == cell]
            py_cell_data = py_lsa[py_lsa['cell'] == cell]
            
            # Extract chromosomes (simplified)
            r_chrs = set([region.split(':')[0] if ':' in region else region.split('_')[0] 
                         for region in r_cell_data['region']])
            py_chrs = set([region.split(':')[0] if ':' in region else region.split('_')[0] 
                          for region in py_cell_data['region']])
            
            if r_chrs & py_chrs:  # Any common chromosomes
                chr_matches += 1
                
            # Check CNA types
            r_cna_types = set(r_cell_data['CNA'])
            py_cna_types = set(py_cell_data['CNA'])
            
            if r_cna_types & py_cna_types:  # Any common CNA types
                cna_type_matches += 1
        
        if common_cells_with_cnas:
            chr_consistency = (chr_matches / len(common_cells_with_cnas)) * 100
            cna_consistency = (cna_type_matches / len(common_cells_with_cnas)) * 100
            
            print(f"Chromosome-level consistency: {chr_consistency:.1f}%")
            print(f"CNA type consistency: {cna_consistency:.1f}%")
    
    return r_lsa, py_lsa

def check_permutations(perm_dir, n_expected=500):
    """Check if permutation files were generated correctly."""
    print("\n=== Permutation Check ===")
    
    if not os.path.exists(perm_dir):
        print(f"Permutation directory {perm_dir} does not exist!")
        return False
    
    # Count permutation files
    cnv_files = [f for f in os.listdir(perm_dir) if f.endswith('.CNV.txt') and 'permute' in f]
    gene_files = [f for f in os.listdir(perm_dir) if f.endswith('.gene.CNV.txt') and 'permute' in f]
    
    print(f"Found {len(cnv_files)} binned CNV files")
    print(f"Found {len(gene_files)} gene CNV files")
    print(f"Expected: {n_expected} of each")
    
    # Check a sample permutation file
    if cnv_files:
        sample_file = os.path.join(perm_dir, 'permute.1.CNV.txt')
        if os.path.exists(sample_file):
            perm_cnv = pd.read_csv(sample_file, sep='\t', index_col=0)
            print(f"\nSample permutation shape: {perm_cnv.shape}")
            print(f"Value range: [{perm_cnv.min().min()}, {perm_cnv.max().max()}]")
    
    return len(cnv_files) >= n_expected

def check_visualizations(r_dir, py_dir):
    """Check visualization files in both directories."""
    
    # Expected visualization files
    r_viz_files = [
        'singlecell.tree.pdf',  # R naming convention
        'LSA.tree.pdf'          # R naming convention
    ]
    
    py_viz_files = [
        'singlecell_tree.pdf',           # Python naming convention
        'LSA_tree.pdf',                  # Python naming convention
        'MEDALT_visualization_report.pdf' # Python comprehensive report
    ]
    
    print("R visualization files:")
    for viz_file in r_viz_files:
        path = os.path.join(r_dir, viz_file)
        exists = os.path.exists(path)
        print(f"  {viz_file}: {'✓' if exists else '✗'}")
        if exists:
            size = os.path.getsize(path)
            print(f"    Size: {size:,} bytes")
    
    print("\nPython visualization files:")
    for viz_file in py_viz_files:
        path = os.path.join(py_dir, viz_file)
        exists = os.path.exists(path)
        print(f"  {viz_file}: {'✓' if exists else '✗'}")
        if exists:
            size = os.path.getsize(path)
            print(f"    Size: {size:,} bytes")
    
    # Check for any additional PDF files
    try:
        r_pdfs = [f for f in os.listdir(r_dir) if f.endswith('.pdf')]
        py_pdfs = [f for f in os.listdir(py_dir) if f.endswith('.pdf')]
        
        if r_pdfs:
            print(f"\nAll R PDF files: {r_pdfs}")
        if py_pdfs:
            print(f"All Python PDF files: {py_pdfs}")
            
    except OSError:
        print("Could not list directory contents for PDF check")

def main():
    """Main comparison function."""
    if len(sys.argv) < 3:
        print("Usage: python debug_comparison.py <R_output_dir> <Python_output_dir>")
        sys.exit(1)
    
    r_dir = sys.argv[1]
    py_dir = sys.argv[2]
    
    print(f"Comparing R output: {r_dir}")
    print(f"With Python output: {py_dir}")
    
    # Compare trees
    r_tree_path = os.path.join(r_dir, 'CNV.tree.txt')
    py_tree_path = os.path.join(py_dir, '3_CNV.tree.txt')
    
    if os.path.exists(r_tree_path) and os.path.exists(py_tree_path):
        r_tree, py_tree = compare_trees(r_tree_path, py_tree_path)
    else:
        print("Tree files not found!")
    
    # Compare binned CNV
    # Locate binned CNV files in both R and Python outputs
    py_cnv_files = [f for f in os.listdir(py_dir) if f.startswith('2_') and ('_bin_' in f or f.endswith('.CNV.txt'))]
    r_cnv_files = [f for f in os.listdir(r_dir) if f.startswith('2_') and ('_bin_' in f or f.endswith('.CNV.txt'))]

    if py_cnv_files and r_cnv_files:
        py_cnv_path = os.path.join(py_dir, py_cnv_files[0])
        r_cnv_path = os.path.join(r_dir, r_cnv_files[0])
        compare_binned_cnv(r_cnv_path, py_cnv_path)
    else:
        print("Binned CNV files not found in one or both directories – skipping CNV comparison.")
    
    # Compare LSA results
    r_lsa_path = os.path.join(r_dir, 'segmental.LSA.txt')
    py_lsa_path = os.path.join(py_dir, 'segmental.LSA.txt')
    
    if os.path.exists(r_lsa_path) and os.path.exists(py_lsa_path):
        r_lsa, py_lsa = compare_lsa_results(r_lsa_path, py_lsa_path)
    else:
        print("LSA files not found!")
    
    # Check permutations
    py_perm_dir = os.path.join(py_dir, 'permutation')
    check_permutations(py_perm_dir)
    
    # Check visualizations
    print("\n=== Visualization Files Check ===")
    check_visualizations(r_dir, py_dir)
    
    print("\n=== Summary ===")
    print("Comparison complete! Key areas analyzed:")
    print("1. Tree structure (edges, weights, topology)")
    print("2. Binning strategy and CNV profiles")
    print("3. LSA results consistency (cells, chromosomes, CNA types)")
    print("4. Permutation file generation")
    print("5. Visualization file availability")
    print("\nFor detailed biological validation, compare:")
    print("- Common cells with significant CNAs")
    print("- Chromosome-level patterns (AMP/DEL)")
    print("- Tree topology and phylogenetic relationships")
    print("- Generated visualization files")
    print("\nThe Python implementation now includes comprehensive visualizations!")
    print("Check the generated PDF files for graphical results comparison.")

if __name__ == '__main__':
    main()