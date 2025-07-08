#!/usr/bin/env python
"""
Example usage of MEDALT for scRNA-seq CNV analysis
"""

import numpy as np
import pandas as pd
import os
from medalt_implementation import (
    MED, MEDALT, LineageSpeciationAnalysis,
    load_scRNA_data
)
from medalt_utils import MEDALTPipeline, MEDALTVisualizer, ChromosomeMapper


def example_basic_usage():
    """Basic example using the core components directly"""
    print("=== Basic MEDALT Usage Example ===\n")
    
    # 1. Demonstrate MED calculation
    print("1. Minimal Event Distance Calculation:")
    profile_a = np.array([2, 2, 2, 3, 3, 1, 1])
    profile_b = np.array([2, 3, 3, 3, 3, 2, 2])
    
    distance = MED.compute(profile_a, profile_b)
    print(f"   Profile A: {profile_a}")
    print(f"   Profile B: {profile_b}")
    print(f"   MED = {distance}")
    
    # Example with homozygous deletion
    profile_c = np.array([0, 0, 2, 2, 2])
    profile_d = np.array([1, 1, 2, 2, 2])
    
    distance2 = MED.compute(profile_c, profile_d)
    print(f"\n   Profile C (with CN=0): {profile_c}")
    print(f"   Profile D: {profile_d}")
    print(f"   MED = {distance2} (infinity due to homozygous deletion)")
    
    # 2. Small example of tree construction
    print("\n2. Tree Construction Example:")
    # Create a small synthetic dataset
    cnv_matrix = np.array([
        [2, 2, 2, 2, 2],  # Cell 0: Normal
        [2, 3, 3, 2, 2],  # Cell 1: Small amplification
        [2, 3, 3, 3, 3],  # Cell 2: Larger amplification
        [1, 1, 2, 2, 2],  # Cell 3: Deletion
        [1, 1, 1, 2, 2],  # Cell 4: Larger deletion
    ])
    
    # Build tree
    medalt = MEDALT(cnv_matrix)
    tree = medalt.build_tree()
    
    print("   Tree edges:")
    for parent, child, data in tree.edges(data=True):
        print(f"   {parent} -> {child} (distance: {data['weight']:.1f})")


def example_full_pipeline():
    """Example using the full pipeline with example data"""
    print("\n\n=== Full Pipeline Example ===\n")
    
    # Check if example data exists
    example_file = 'example/scRNA.CNV.txt'
    
    if not os.path.exists(example_file):
        print("Creating synthetic example data...")
        create_example_data(example_file)
    
    # Run the complete pipeline
    print("Running MEDALT pipeline on example data...")
    pipeline = MEDALTPipeline(output_dir='example_output')
    
    # Run with fewer permutations for quick example
    tree, results = pipeline.run_analysis(
        input_file=example_file,
        data_type='R',
        genome_version='hg19',
        window_size=30,
        n_permutations=50,  # Reduced for example
        min_lineage_size=3,  # Smaller for example
        n_jobs=2  # Limit parallelism for example
    )
    
    print("\nPipeline completed!")
    print(f"Output files saved to: example_output/")
    print("\nGenerated files:")
    for file in os.listdir('example_output'):
        print(f"  - {file}")
    
    if not results.empty:
        print(f"\nFound {len(results)} significant CNAs")
        print("\nTop results:")
        print(results[['region', 'CNA', 'pvalue', 'adjustp', 'cell']].head())


def create_example_data(filepath: str):
    """Create synthetic scRNA-seq CNV data for testing"""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    
    np.random.seed(42)
    n_cells = 50
    n_genes = 1000
    
    # Create base profiles (relative copy numbers)
    # Normal cells (around 1.0)
    normal_profile = np.random.normal(1.0, 0.1, (20, n_genes))
    normal_profile = np.clip(normal_profile, 0.5, 1.5)
    
    # Create lineage with amplification
    amp_profile = normal_profile[:15].copy()
    # Add amplification in region 200-400
    amp_profile[:, 200:400] = np.random.normal(1.5, 0.1, (15, 200))
    
    # Create lineage with deletion
    del_profile = normal_profile[:15].copy()
    # Add deletion in region 600-800
    del_profile[:, 600:800] = np.random.normal(0.5, 0.1, (15, 200))
    
    # Combine all profiles
    all_profiles = np.vstack([normal_profile, amp_profile, del_profile])
    
    # Add some noise
    all_profiles += np.random.normal(0, 0.05, all_profiles.shape)
    all_profiles = np.clip(all_profiles, 0.1, 2.0)
    
    # Create DataFrame
    cell_names = [f'cell_{i}' for i in range(n_cells)]
    gene_names = [f'gene_{i}' for i in range(n_genes)]
    
    df = pd.DataFrame(all_profiles.T, index=gene_names, columns=cell_names)
    
    # Save to file
    df.to_csv(filepath, sep='\t')
    print(f"Created example data with {n_cells} cells and {n_genes} genes")
    print(f"- 20 normal cells")
    print(f"- 15 cells with amplification (genes 200-400)")
    print(f"- 15 cells with deletion (genes 600-800)")


def example_custom_analysis():
    """Example of custom analysis using individual components"""
    print("\n\n=== Custom Analysis Example ===\n")
    
    # Load data
    example_file = 'example/scRNA.CNV.txt'
    if not os.path.exists(example_file):
        create_example_data(example_file)
    
    cnv_matrix, cell_names, bin_names = load_scRNA_data(example_file, window_size=50)
    
    print(f"Loaded data: {cnv_matrix.shape[0]} cells, {cnv_matrix.shape[1]} bins")
    
    # Build tree with custom parameters
    print("\nBuilding tree...")
    medalt = MEDALT(cnv_matrix)
    tree = medalt.build_tree()
    
    # Analyze tree structure
    print(f"\nTree statistics:")
    print(f"  Total nodes: {tree.number_of_nodes()}")
    print(f"  Total edges: {tree.number_of_edges()}")
    
    # Find major lineages
    large_subtrees = []
    for node in tree.nodes():
        if node == 'root':
            continue
        descendants = list(nx.descendants(tree, node))
        if len(descendants) >= 5:  # Lineages with at least 5 cells
            large_subtrees.append((node, len(descendants)))
    
    print(f"\nMajor lineages (≥5 cells):")
    for node, size in sorted(large_subtrees, key=lambda x: x[1], reverse=True)[:5]:
        depth = nx.shortest_path_length(tree, 'root', node)
        print(f"  Cell {node}: {size} descendants, depth {depth}")
    
    # Run targeted LSA on specific lineages
    if large_subtrees:
        print("\nRunning LSA on largest lineage...")
        
        # Get chromosome boundaries
        chr_mapper = ChromosomeMapper('hg19')
        chr_boundaries = chr_mapper.get_chromosome_boundaries(cnv_matrix.shape[1])
        
        # Run LSA with custom parameters
        lsa = LineageSpeciationAnalysis(
            tree, cnv_matrix, chr_boundaries, 
            n_permutations=100  # Quick analysis
        )
        
        # Get results
        results = lsa.run_analysis(min_lineage_size=5, n_jobs=2)
        
        if not results.empty:
            print(f"\nFound {len(results)} CNAs in lineages")
            print("\nMost significant CNAs:")
            print(results[['region', 'CNA', 'pvalue', 'cell', 'subtreesize']].head())
        
        # Visualize specific lineage
        largest_lineage_root = large_subtrees[0][0]
        subtree = tree.subgraph(
            nx.descendants(tree, largest_lineage_root) | {largest_lineage_root}
        )
        
        print(f"\nVisualizing largest lineage (rooted at cell {largest_lineage_root})...")
        visualizer = MEDALTVisualizer()
        visualizer.plot_tree(subtree, 'largest_lineage.png')


def example_parallel_evolution():
    """Example showing how to detect parallel evolution"""
    print("\n\n=== Parallel Evolution Detection Example ===\n")
    
    # Create synthetic data with parallel evolution
    n_cells = 60
    n_bins = 100
    
    # Base diploid state
    cnv_matrix = np.full((n_cells, n_bins), 2, dtype=int)
    
    # Lineage 1: cells 10-25 with amplification at bins 20-30
    cnv_matrix[10:25, 20:30] = 3
    
    # Lineage 2: cells 35-50 with same amplification (parallel evolution)
    cnv_matrix[35:50, 20:30] = 3
    
    # Different alterations to make lineages distinct
    cnv_matrix[10:25, 60:70] = 1  # Deletion in lineage 1
    cnv_matrix[35:50, 80:90] = 1  # Different deletion in lineage 2
    
    print("Created synthetic data with parallel amplification in two lineages")
    
    # Build tree and run analysis
    medalt = MEDALT(cnv_matrix)
    tree = medalt.build_tree()
    
    # This would be extended with PLSA implementation
    print("\nParallel evolution analysis would identify:")
    print("  - Amplification at bins 20-30 in multiple independent lineages")
    print("  - Statistical significance of this convergent evolution")
    
    # Save example
    np.save('parallel_evolution_example.npy', cnv_matrix)
    print("\nSaved parallel evolution example to 'parallel_evolution_example.npy'")


if __name__ == '__main__':
    # Run all examples
    print("MEDALT Example Scripts")
    print("=" * 50)
    
    # Basic usage
    example_basic_usage()
    
    # Full pipeline
    example_full_pipeline()
    
    # Custom analysis
    example_custom_analysis()
    
    # Parallel evolution
    example_parallel_evolution()
    
    print("\n\nAll examples completed!")
    print("\nTo run MEDALT on your own data:")
    print("  python medalt_utils.py -I your_data.txt -D R -G hg19 -O output_dir")