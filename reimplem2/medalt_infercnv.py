#!/usr/bin/env python
"""
MEDALT Guide for inferCNV Output Analysis
Practical examples and utilities for analyzing inferCNV results
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Optional

from medalt_optimized import (
    InferCNVReader, OptimizedMEDALT, OptimizedLSA, 
    InferCNVPipeline, compute_med_fast
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class InferCNVAnalysisGuide:
    """
    Complete guide for analyzing inferCNV output with MEDALT
    """
    
    def __init__(self, infercnv_dir: str):
        """
        Initialize with inferCNV output directory
        
        Args:
            infercnv_dir: Path to inferCNV output directory
        """
        self.infercnv_dir = Path(infercnv_dir)
        self.validate_infercnv_output()
    
    def validate_infercnv_output(self):
        """Check that required inferCNV files exist"""
        required_files = [
            'infercnv.observations.txt',
            'infercnv.references.txt'
        ]
        
        optional_files = [
            'infercnv.observation_groupings.txt',
            'infercnv.observations_dendrogram.txt'
        ]
        
        logger.info(f"Checking inferCNV output in {self.infercnv_dir}")
        
        for file in required_files:
            if not (self.infercnv_dir / file).exists():
                raise FileNotFoundError(f"Required file {file} not found!")
            logger.info(f"✓ Found {file}")
        
        for file in optional_files:
            if (self.infercnv_dir / file).exists():
                logger.info(f"✓ Found optional file {file}")
    
    def explore_data_characteristics(self) -> Dict:
        """
        Explore the characteristics of the inferCNV data
        """
        obs_file = self.infercnv_dir / 'infercnv.observations.txt'
        
        # Get dimensions without loading full data
        with open(obs_file, 'r') as f:
            header = f.readline().strip().split('\t')
            n_cells = len(header) - 1  # First column is gene names
            
            # Count genes
            n_genes = sum(1 for _ in f)
        
        logger.info(f"Data dimensions: {n_cells} cells × {n_genes} genes")
        
        # Sample data to check value ranges
        sample_df = pd.read_csv(obs_file, sep='\t', index_col=0, nrows=100)
        
        value_stats = {
            'min': sample_df.values.min(),
            'max': sample_df.values.max(),
            'mean': sample_df.values.mean(),
            'std': sample_df.values.std()
        }
        
        logger.info(f"Value statistics (sample): {value_stats}")
        
        # Check for existing groupings
        groupings = {}
        groupings_file = self.infercnv_dir / 'infercnv.observation_groupings.txt'
        if groupings_file.exists():
            group_df = pd.read_csv(groupings_file, sep='\t', header=None)
            groupings = dict(zip(group_df[0], group_df[1]))
            unique_groups = len(set(groupings.values()))
            logger.info(f"Found {unique_groups} cell groups in inferCNV groupings")
        
        return {
            'n_cells': n_cells,
            'n_genes': n_genes,
            'value_stats': value_stats,
            'has_groupings': len(groupings) > 0,
            'n_groups': len(set(groupings.values())) if groupings else 0
        }
    
    def select_optimal_parameters(self, data_stats: Dict) -> Dict:
        """
        Suggest optimal parameters based on data characteristics
        """
        n_cells = data_stats['n_cells']
        n_genes = data_stats['n_genes']
        
        # Gene selection
        if n_genes > 10000:
            top_k_genes = 1000
            logger.info(f"Large gene set ({n_genes}), selecting top 1000 variance genes")
        elif n_genes > 5000:
            top_k_genes = min(1000, n_genes // 5)
            logger.info(f"Medium gene set, selecting top {top_k_genes} genes")
        else:
            top_k_genes = n_genes // 3
            logger.info(f"Small gene set, selecting top {top_k_genes} genes")
        
        # Tree building strategy
        if n_cells > 5000:
            subsample_size = 2000
            logger.info(f"Large dataset ({n_cells} cells), using subsampling (n={subsample_size})")
        else:
            subsample_size = None
            logger.info(f"Small dataset, building complete tree")
        
        # Permutations - balance accuracy and speed
        if n_cells > 10000:
            n_permutations = 300
            logger.info("Very large dataset, using 300 permutations for speed")
        else:
            n_permutations = 500
            logger.info("Using 500 permutations for better statistical power")
        
        # Lineage size
        min_lineage_size = max(20, int(n_cells * 0.01))
        logger.info(f"Minimum lineage size: {min_lineage_size}")
        
        return {
            'top_k_genes': top_k_genes,
            'subsample_size': subsample_size,
            'n_permutations': n_permutations,
            'min_lineage_size': min_lineage_size
        }
    
    def run_basic_analysis(self, output_dir: str = 'medalt_analysis'):
        """
        Run basic MEDALT analysis with automatic parameter selection
        """
        logger.info("Starting basic MEDALT analysis...")
        
        # Explore data
        data_stats = self.explore_data_characteristics()
        
        # Get optimal parameters
        params = self.select_optimal_parameters(data_stats)
        
        # Run pipeline
        pipeline = InferCNVPipeline(output_dir)
        
        tree, results = pipeline.run_infercnv_analysis(
            observations_file=str(self.infercnv_dir / 'infercnv.observations.txt'),
            reference_file=str(self.infercnv_dir / 'infercnv.references.txt'),
            **params
        )
        
        return tree, results
    
    def advanced_analysis_with_groups(self, output_dir: str = 'medalt_grouped'):
        """
        Advanced analysis incorporating inferCNV cell groupings
        """
        logger.info("Running advanced analysis with cell groupings...")
        
        # Load groupings if available
        groupings_file = self.infercnv_dir / 'infercnv.observation_groupings.txt'
        if not groupings_file.exists():
            logger.warning("No groupings file found, running standard analysis")
            return self.run_basic_analysis(output_dir)
        
        # Load data with groupings
        obs_file = self.infercnv_dir / 'infercnv.observations.txt'
        ref_file = self.infercnv_dir / 'infercnv.references.txt'
        
        # Load groupings
        group_df = pd.read_csv(groupings_file, sep='\t', header=None)
        cell_groups = dict(zip(group_df[0], group_df[1]))
        
        # Load CNV data
        cnv_matrix, cell_names, gene_names, gene_indices = InferCNVReader.load_observations(
            str(obs_file), top_k_genes=1000
        )
        
        # Create group labels array
        group_labels = np.array([cell_groups.get(cell, 'Unknown') for cell in cell_names])
        unique_groups = np.unique(group_labels)
        
        logger.info(f"Found {len(unique_groups)} groups: {unique_groups}")
        
        # Build tree
        medalt = OptimizedMEDALT(cnv_matrix)
        tree = medalt.build_tree_fast(subsample_size=2000)
        
        # Run group-aware LSA
        results = self._group_aware_lsa(tree, cnv_matrix, group_labels, gene_names)
        
        # Save results
        Path(output_dir).mkdir(exist_ok=True)
        results.to_csv(Path(output_dir) / 'group_aware_results.txt', sep='\t', index=False)
        
        # Visualize
        self.visualize_groups_on_tree(tree, group_labels, cell_names, output_dir)
        
        return tree, results
    
    def _group_aware_lsa(self, tree, cnv_matrix, group_labels, gene_names):
        """LSA that considers pre-defined groups"""
        # Standard LSA first
        lsa = OptimizedLSA(tree, cnv_matrix, n_permutations=500)
        standard_results = lsa.run_analysis_fast()
        
        # Add group enrichment analysis
        if not standard_results.empty:
            standard_results['gene_name'] = standard_results['gene_idx'].apply(
                lambda x: gene_names[x]
            )
            
            # For each significant lineage, check group enrichment
            group_enrichments = []
            
            for _, row in standard_results.iterrows():
                lineage_root = row['lineage_root']
                
                # Get cells in lineage
                descendants = set(nx.descendants(tree, lineage_root))
                descendants.add(lineage_root)
                lineage_cells = [n for n in descendants if n != 'root']
                
                # Calculate group enrichment
                lineage_groups = group_labels[lineage_cells]
                group_counts = pd.Series(lineage_groups).value_counts()
                
                # Find most enriched group
                total_cells = len(lineage_cells)
                max_group = group_counts.index[0]
                max_fraction = group_counts.iloc[0] / total_cells
                
                group_enrichments.append({
                    'enriched_group': max_group,
                    'group_fraction': max_fraction
                })
            
            # Add to results
            enrichment_df = pd.DataFrame(group_enrichments)
            standard_results = pd.concat([standard_results, enrichment_df], axis=1)
        
        return standard_results
    
    def visualize_groups_on_tree(self, tree, group_labels, cell_names, output_dir):
        """Visualize cell groups on the lineage tree"""
        import networkx as nx
        from matplotlib import cm
        
        # Create color map for groups
        unique_groups = np.unique(group_labels)
        colors = cm.get_cmap('tab10')(np.linspace(0, 1, len(unique_groups)))
        group_colors = dict(zip(unique_groups, colors))
        
        # Assign colors to nodes
        node_colors = []
        for node in tree.nodes():
            if node == 'root':
                node_colors.append('white')
            else:
                group = group_labels[node]
                node_colors.append(group_colors[group])
        
        # Plot
        plt.figure(figsize=(12, 10))
        pos = nx.spring_layout(tree, k=2, iterations=50)
        
        nx.draw_networkx_nodes(tree, pos, node_color=node_colors, 
                              node_size=50, alpha=0.8)
        nx.draw_networkx_edges(tree, pos, alpha=0.3, arrows=False)
        
        # Add legend
        for group, color in group_colors.items():
            plt.scatter([], [], c=[color], label=group, s=100)
        
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.title("Cell Groups on Lineage Tree")
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(Path(output_dir) / 'groups_on_tree.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def compare_with_infercnv_clusters(self, medalt_tree, output_dir):
        """
        Compare MEDALT lineages with inferCNV's hierarchical clustering
        """
        dendrogram_file = self.infercnv_dir / 'infercnv.observations_dendrogram.txt'
        
        if not dendrogram_file.exists():
            logger.warning("No inferCNV dendrogram found for comparison")
            return
        
        logger.info("Comparing MEDALT lineages with inferCNV clustering...")
        
        # This would require parsing the dendrogram and comparing
        # For now, we'll create a comparison summary
        
        # Get MEDALT lineages
        medalt_lineages = self._extract_major_lineages(medalt_tree)
        
        summary = {
            'n_medalt_lineages': len(medalt_lineages),
            'lineage_sizes': [len(lin['cells']) for lin in medalt_lineages]
        }
        
        logger.info(f"MEDALT found {summary['n_medalt_lineages']} major lineages")
        logger.info(f"Lineage sizes: {summary['lineage_sizes']}")
        
        return summary
    
    def _extract_major_lineages(self, tree, min_size=50):
        """Extract major lineages from tree"""
        lineages = []
        
        for node in tree.nodes():
            if node == 'root':
                continue
            
            descendants = set(nx.descendants(tree, node))
            descendants.add(node)
            
            if len(descendants) >= min_size:
                # Check if this is not a subset of existing lineage
                is_subset = False
                for existing in lineages:
                    if descendants.issubset(existing['cells']):
                        is_subset = True
                        break
                
                if not is_subset:
                    lineages.append({
                        'root': node,
                        'cells': descendants,
                        'size': len(descendants)
                    })
        
        return lineages
    
    def export_for_downstream(self, tree, results, output_dir):
        """
        Export results in formats suitable for downstream analysis
        """
        output_path = Path(output_dir)
        
        # 1. Cell lineage assignments
        lineage_assignments = self._get_lineage_assignments(tree)
        lineage_df = pd.DataFrame(lineage_assignments.items(), 
                                 columns=['cell_id', 'lineage'])
        lineage_df.to_csv(output_path / 'cell_lineage_assignments.txt', 
                         sep='\t', index=False)
        
        # 2. Significant genes per lineage
        if not results.empty:
            sig_genes = results[results['qvalue'] < 0.05]
            sig_genes.to_csv(output_path / 'significant_genes_by_lineage.txt',
                           sep='\t', index=False)
        
        # 3. Tree structure for Cytoscape
        cytoscape_edges = []
        for parent, child, data in tree.edges(data=True):
            cytoscape_edges.append({
                'source': parent,
                'target': child,
                'weight': data.get('weight', 0),
                'interaction': 'lineage'
            })
        
        cyto_df = pd.DataFrame(cytoscape_edges)
        cyto_df.to_csv(output_path / 'tree_for_cytoscape.txt', 
                      sep='\t', index=False)
        
        logger.info(f"Exported results to {output_path}")
    
    def _get_lineage_assignments(self, tree):
        """Assign each cell to its major lineage"""
        # Find major branching points
        major_branches = []
        
        for node in tree.nodes():
            if node == 'root':
                continue
            
            # Major branch if has many descendants
            descendants = list(nx.descendants(tree, node))
            if len(descendants) >= 50:
                major_branches.append((node, len(descendants)))
        
        # Sort by size
        major_branches.sort(key=lambda x: x[1], reverse=True)
        
        # Assign cells
        assignments = {}
        for i, (branch_root, _) in enumerate(major_branches[:10]):  # Top 10 lineages
            descendants = nx.descendants(tree, branch_root)
            descendants.add(branch_root)
            
            for cell in descendants:
                if cell not in assignments and cell != 'root':
                    assignments[cell] = f'Lineage_{i+1}'
        
        # Assign remaining cells
        for node in tree.nodes():
            if node != 'root' and node not in assignments:
                assignments[node] = 'Other'
        
        return assignments


def interactive_analysis_example():
    """
    Example of interactive analysis workflow
    """
    print("=== MEDALT Interactive Analysis for inferCNV ===\n")
    
    # Get inferCNV directory
    infercnv_dir = input("Enter path to inferCNV output directory: ").strip()
    
    if not Path(infercnv_dir).exists():
        print(f"Error: Directory {infercnv_dir} not found!")
        return
    
    # Initialize guide
    guide = InferCNVAnalysisGuide(infercnv_dir)
    
    # Explore data
    print("\nExploring data characteristics...")
    data_stats = guide.explore_data_characteristics()
    
    # Ask for analysis type
    print("\nSelect analysis type:")
    print("1. Basic analysis (automatic parameters)")
    print("2. Advanced analysis with cell groups")
    print("3. Custom analysis")
    
    choice = input("Enter choice (1-3): ").strip()
    
    if choice == '1':
        tree, results = guide.run_basic_analysis()
        print(f"\nAnalysis complete! Found {len(results)} significant associations")
        
    elif choice == '2':
        tree, results = guide.advanced_analysis_with_groups()
        print(f"\nGroup-aware analysis complete!")
        
    elif choice == '3':
        # Custom parameters
        print("\nEnter custom parameters:")
        top_k = int(input(f"Number of top variance genes (default 1000): ") or 1000)
        n_perm = int(input(f"Number of permutations (default 500): ") or 500)
        min_lin = int(input(f"Minimum lineage size (default 20): ") or 20)
        
        # Run with custom parameters
        pipeline = InferCNVPipeline('custom_analysis')
        tree, results = pipeline.run_infercnv_analysis(
            observations_file=str(Path(infercnv_dir) / 'infercnv.observations.txt'),
            top_k_genes=top_k,
            n_permutations=n_perm,
            min_lineage_size=min_lin
        )
    
    # Export results
    if input("\nExport results for downstream analysis? (y/n): ").lower() == 'y':
        guide.export_for_downstream(tree, results, 'exported_results')
        print("Results exported to 'exported_results' directory")


def batch_analysis_example():
    """
    Example for batch processing multiple inferCNV outputs
    """
    import glob
    
    # Find all inferCNV output directories
    infercnv_dirs = glob.glob('*/infercnv.observations.txt')
    infercnv_dirs = [Path(p).parent for p in infercnv_dirs]
    
    logger.info(f"Found {len(infercnv_dirs)} inferCNV outputs")
    
    all_results = []
    
    for dir_path in infercnv_dirs:
        sample_name = dir_path.name
        logger.info(f"\nProcessing {sample_name}...")
        
        try:
            guide = InferCNVAnalysisGuide(str(dir_path))
            tree, results = guide.run_basic_analysis(
                output_dir=f'batch_results/{sample_name}'
            )
            
            # Add sample name to results
            if not results.empty:
                results['sample'] = sample_name
                all_results.append(results)
                
        except Exception as e:
            logger.error(f"Failed to process {sample_name}: {e}")
    
    # Combine results
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined.to_csv('batch_results/all_samples_combined.txt', sep='\t', index=False)
        logger.info(f"\nCombined results saved for {len(all_results)} samples")


def main():
    """Main entry point with examples"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='MEDALT Guide for inferCNV Output Analysis'
    )
    parser.add_argument('infercnv_dir', help='Path to inferCNV output directory')
    parser.add_argument('--mode', choices=['basic', 'advanced', 'export'],
                       default='basic', help='Analysis mode')
    parser.add_argument('--output', default='medalt_results', 
                       help='Output directory')
    parser.add_argument('--interactive', action='store_true',
                       help='Run interactive analysis')
    
    args = parser.parse_args()
    
    if args.interactive:
        interactive_analysis_example()
        return
    
    # Initialize guide
    guide = InferCNVAnalysisGuide(args.infercnv_dir)
    
    if args.mode == 'basic':
        tree, results = guide.run_basic_analysis(args.output)
        
    elif args.mode == 'advanced':
        tree, results = guide.advanced_analysis_with_groups(args.output)
        
    elif args.mode == 'export':
        # Run analysis first
        tree, results = guide.run_basic_analysis(args.output)
        # Then export
        guide.export_for_downstream(tree, results, args.output)
    
    print(f"\nAnalysis complete! Results saved to {args.output}")
    
    if not results.empty:
        print(f"\nTop 5 significant associations:")
        print(results[['gene_name', 'direction', 'enrichment_score', 
                      'pvalue', 'qvalue']].head())


if __name__ == '__main__':
    main()