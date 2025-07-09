#!/usr/bin/env python
"""
MEDALT Pipeline Runner and Utilities
Provides high-level interface for running MEDALT analysis
"""

import argparse
import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional
import networkx as nx
from datetime import datetime
import logging

# Import core MEDALT components
from medalt_implementation import (
    MED, MEDALT, LineageSpeciationAnalysis,
    load_scRNA_data, logger
)


class ChromosomeMapper:
    """Map genomic bins to chromosome coordinates"""
    
    def __init__(self, genome_version: str = 'hg19'):
        self.genome_version = genome_version
        self.chr_info = self._load_chromosome_info()
    
    def _load_chromosome_info(self) -> Dict:
        """Load chromosome information for the genome"""
        # Simplified chromosome sizes for hg19/hg38
        # In real implementation, load from reference files
        chr_sizes = {
            'chr1': 249250621, 'chr2': 242193529, 'chr3': 198295559,
            'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
            'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
            'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
            'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
            'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
            'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
            'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
        }
        return chr_sizes
    
    def get_chromosome_boundaries(self, n_bins: int) -> List[Tuple[int, int]]:
        """
        Get chromosome boundaries for given number of bins
        
        Returns:
            List of (start_bin, end_bin) tuples for each chromosome
        """
        total_genome_size = sum(self.chr_info.values())
        bin_size = total_genome_size // n_bins
        
        boundaries = []
        current_bin = 0
        
        for chr_name in ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']:
            if chr_name not in self.chr_info:
                continue
            
            chr_bins = self.chr_info[chr_name] // bin_size
            if chr_bins > 0:
                boundaries.append((current_bin, current_bin + chr_bins))
                current_bin += chr_bins
        
        return boundaries
    
    def bin_to_genomic_location(self, bin_idx: int, n_total_bins: int) -> str:
        """Convert bin index to genomic location string"""
        boundaries = self.get_chromosome_boundaries(n_total_bins)
        
        for i, (start, end) in enumerate(boundaries):
            if start <= bin_idx < end:
                chr_idx = i + 1 if i < 22 else ('X' if i == 22 else 'Y')
                chr_name = f'chr{chr_idx}'
                
                # Determine band (simplified)
                position_in_chr = (bin_idx - start) / (end - start)
                if position_in_chr < 0.5:
                    arm = 'p'
                    band = int(position_in_chr * 30)  # Simplified banding
                else:
                    arm = 'q'
                    band = int((position_in_chr - 0.5) * 30)
                
                return f"{chr_name}:{arm}{band}"
        
        return f"bin_{bin_idx}"


class MEDALTVisualizer:
    """Visualization utilities for MEDALT results"""
    
    @staticmethod
    def plot_tree(tree: nx.DiGraph, output_file: str = 'tree_visualization.png',
                  highlight_nodes: List = None, node_colors: Dict = None):
        """
        Visualize the lineage tree
        
        Args:
            tree: NetworkX directed graph
            output_file: Output filename
            highlight_nodes: List of nodes to highlight
            node_colors: Dict mapping nodes to colors
        """
        plt.figure(figsize=(12, 10))
        
        # Calculate layout
        pos = nx.nx_agraph.graphviz_layout(tree, prog='dot')
        
        # Default node colors
        if node_colors is None:
            node_colors = {}
        
        default_color = 'lightblue'
        colors = [node_colors.get(node, default_color) for node in tree.nodes()]
        
        # Highlight specific nodes
        if highlight_nodes:
            node_sizes = [300 if node in highlight_nodes else 100 for node in tree.nodes()]
        else:
            node_sizes = 100
        
        # Draw the tree
        nx.draw(tree, pos, node_color=colors, node_size=node_sizes,
                with_labels=True, arrows=True, arrowsize=20,
                font_size=8, font_weight='bold')
        
        # Add edge labels (distances)
        edge_labels = nx.get_edge_attributes(tree, 'weight')
        edge_labels = {k: f'{v:.1f}' for k, v in edge_labels.items()}
        nx.draw_networkx_edge_labels(tree, pos, edge_labels, font_size=6)
        
        plt.title("MEDALT Lineage Tree", fontsize=16)
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Tree visualization saved to {output_file}")
    
    @staticmethod
    def plot_cnv_heatmap(cnv_matrix: np.ndarray, tree: nx.DiGraph,
                        output_file: str = 'cnv_heatmap.png',
                        cell_names: List[str] = None):
        """
        Plot CNV heatmap ordered by tree structure
        
        Args:
            cnv_matrix: Copy number matrix
            tree: NetworkX directed graph
            output_file: Output filename
            cell_names: List of cell names
        """
        # Order cells by tree traversal
        ordered_cells = []
        
        def traverse(node):
            if node != 'root' and node not in ordered_cells:
                ordered_cells.append(node)
            for child in tree.successors(node):
                traverse(child)
        
        traverse('root')
        
        # Reorder matrix
        ordered_matrix = cnv_matrix[ordered_cells]
        
        # Create heatmap
        plt.figure(figsize=(15, 10))
        
        # Define color map (0=white, 1=blue, 2=gray, 3+=red gradient)
        cmap = plt.cm.colors.ListedColormap(['white', 'blue', 'gray', 'salmon', 'red', 'darkred'])
        bounds = [0, 0.5, 1.5, 2.5, 3.5, 4.5, 6]
        norm = plt.cm.colors.BoundaryNorm(bounds, cmap.N)
        
        # Plot heatmap
        ax = sns.heatmap(ordered_matrix, cmap=cmap, norm=norm,
                        cbar_kws={'label': 'Copy Number'},
                        yticklabels=False, xticklabels=False)
        
        # Add cell labels if provided
        if cell_names:
            ordered_names = [cell_names[i] for i in ordered_cells]
            ax.set_yticklabels(ordered_names, fontsize=6)
        
        plt.title("Copy Number Heatmap (Cells ordered by lineage)", fontsize=14)
        plt.xlabel("Genomic Bins", fontsize=12)
        plt.ylabel("Cells", fontsize=12)
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"CNV heatmap saved to {output_file}")
    
    @staticmethod
    def plot_lsa_results(lsa_results: pd.DataFrame, output_file: str = 'lsa_manhattan.png',
                        significance_threshold: float = 0.05):
        """
        Create Manhattan-like plot for LSA results
        
        Args:
            lsa_results: DataFrame with LSA results
            output_file: Output filename
            significance_threshold: P-value threshold for significance
        """
        if lsa_results.empty:
            logger.warning("No LSA results to plot")
            return
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10), 
                                       gridspec_kw={'height_ratios': [1, 1]})
        
        # Separate amplifications and deletions
        amp_results = lsa_results[lsa_results['CNA'] == 'AMP'].copy()
        del_results = lsa_results[lsa_results['CNA'] == 'DEL'].copy()
        
        # Extract bin numbers from region names
        amp_results['bin'] = amp_results['region'].str.extract('(\d+)').astype(int)
        del_results['bin'] = del_results['region'].str.extract('(\d+)').astype(int)
        
        # Plot amplifications
        if not amp_results.empty:
            ax1.scatter(amp_results['bin'], -np.log10(amp_results['pvalue']),
                       c='red', alpha=0.6, s=30)
            ax1.axhline(y=-np.log10(significance_threshold), color='black',
                       linestyle='--', alpha=0.5, label=f'p={significance_threshold}')
            ax1.set_title('Amplifications', fontsize=12)
            ax1.set_ylabel('-log10(p-value)', fontsize=10)
            ax1.legend()
        
        # Plot deletions
        if not del_results.empty:
            ax2.scatter(del_results['bin'], -np.log10(del_results['pvalue']),
                       c='blue', alpha=0.6, s=30)
            ax2.axhline(y=-np.log10(significance_threshold), color='black',
                       linestyle='--', alpha=0.5, label=f'p={significance_threshold}')
            ax2.set_title('Deletions', fontsize=12)
            ax2.set_ylabel('-log10(p-value)', fontsize=10)
            ax2.set_xlabel('Genomic Bin', fontsize=10)
            ax2.legend()
        
        plt.suptitle('LSA Results: Significant CNAs', fontsize=16)
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"LSA Manhattan plot saved to {output_file}")


class MEDALTPipeline:
    """Complete MEDALT analysis pipeline"""
    
    def __init__(self, output_dir: str = 'medalt_output'):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # Set up logging to file
        log_file = os.path.join(output_dir, 'medalt_analysis.log')
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        ))
        logger.addHandler(file_handler)
        
    def run_analysis(self, input_file: str, data_type: str = 'R',
                    genome_version: str = 'hg19', window_size: int = 30,
                    n_permutations: int = 500, min_lineage_size: int = 5,
                    n_jobs: int = -1):
        """
        Run complete MEDALT analysis pipeline
        
        Args:
            input_file: Path to input CNV file
            data_type: 'R' for scRNA-seq, 'D' for scDNA-seq
            genome_version: Genome version (hg19 or hg38)
            window_size: Window size for scRNA-seq data
            n_permutations: Number of permutations for LSA
            min_lineage_size: Minimum lineage size for LSA
            n_jobs: Number of parallel jobs
        """
        logger.info(f"Starting MEDALT analysis on {input_file}")
        logger.info(f"Parameters: data_type={data_type}, genome={genome_version}, "
                   f"window_size={window_size}, permutations={n_permutations}")
        
        # Load data
        if data_type == 'R':
            cnv_matrix, cell_names, gene_names = load_scRNA_data(input_file, window_size)
        else:
            # For scDNA-seq, implement appropriate loader
            raise NotImplementedError("scDNA-seq data loading not yet implemented")
        
        # Get chromosome boundaries
        chr_mapper = ChromosomeMapper(genome_version)
        chr_boundaries = chr_mapper.get_chromosome_boundaries(cnv_matrix.shape[1])
        
        # Build MEDALT tree
        logger.info("Building MEDALT tree...")
        medalt = MEDALT(cnv_matrix, chr_boundaries)
        tree = medalt.build_tree()
        
        # Save tree
        tree_file = os.path.join(self.output_dir, 'CNV.tree.txt')
        self._save_tree(tree, tree_file, cell_names)
        
        # Run LSA
        logger.info("Running Lineage Speciation Analysis...")
        lsa = LineageSpeciationAnalysis(tree, cnv_matrix, chr_boundaries, n_permutations)
        lsa_results = lsa.run_analysis(min_lineage_size, n_jobs, gene_names)
        
        # Map bins to genomic locations (only if using windowed data)
        if not lsa_results.empty and window_size > 1:
            lsa_results['genomic_location'] = lsa_results['region'].apply(
                lambda x: chr_mapper.bin_to_genomic_location(
                    gene_names.index(x) if x in gene_names else 0, cnv_matrix.shape[1]
                )
            )
        
        # Save LSA results
        if not lsa_results.empty:
            # Gene-level results
            gene_file = os.path.join(self.output_dir, 'gene.LSA.txt')
            lsa_results.to_csv(gene_file, sep='\t', index=False)
            
            # Segmental results (aggregate by chromosome arm)
            seg_results = self._aggregate_to_segments(lsa_results)
            seg_file = os.path.join(self.output_dir, 'segmental.LSA.txt')
            seg_results.to_csv(seg_file, sep='\t', index=False)
        
        # Generate visualizations
        logger.info("Generating visualizations...")
        visualizer = MEDALTVisualizer()
        
        # Tree visualization
        tree_plot = os.path.join(self.output_dir, 'singlecell.tree.pdf')
        visualizer.plot_tree(tree, tree_plot)
        
        # CNV heatmap
        heatmap_file = os.path.join(self.output_dir, 'cnv_heatmap.png')
        visualizer.plot_cnv_heatmap(cnv_matrix, tree, heatmap_file, cell_names)
        
        # LSA results plot
        if not lsa_results.empty:
            lsa_plot = os.path.join(self.output_dir, 'LSA.tree.pdf')
            visualizer.plot_lsa_results(lsa_results, lsa_plot)
        
        # Save summary statistics
        self._save_summary(tree, lsa_results)
        
        logger.info(f"Analysis complete! Results saved to {self.output_dir}")
        
        return tree, lsa_results
    
    def _save_tree(self, tree: nx.DiGraph, output_file: str, cell_names: List[str] = None):
        """Save tree in MEDALT format"""
        with open(output_file, 'w') as f:
            f.write("from\tto\tdist\n")
            for parent, child, data in tree.edges(data=True):
                parent_name = parent if parent == 'root' else f'cell_{parent}'
                child_name = f'cell_{child}'
                
                if cell_names and parent != 'root':
                    parent_name = cell_names[parent]
                if cell_names and isinstance(child, int):
                    child_name = cell_names[child]
                
                f.write(f"{parent_name}\t{child_name}\t{int(data['weight'])}\n")
    
    def _aggregate_to_segments(self, lsa_results: pd.DataFrame) -> pd.DataFrame:
        """Aggregate bin-level results to chromosome segments"""
        if 'genomic_location' not in lsa_results.columns:
            return lsa_results
        
        # Extract chromosome and arm
        lsa_results['chromosome'] = lsa_results['genomic_location'].str.split(':').str[0]
        lsa_results['arm'] = lsa_results['genomic_location'].str.extract(':([pq])')
        
        # Group by chromosome arm
        grouped = lsa_results.groupby(['chromosome', 'arm', 'CNA', 'cell']).agg({
            'score': 'mean',
            'pvalue': 'min',  # Most significant p-value
            'adjustp': 'min',
            'depth': 'first',
            'subtreesize': 'first'
        }).reset_index()
        
        # Create segment names
        grouped['region'] = grouped['chromosome'] + ':' + grouped['arm']
        
        return grouped
    
    def _save_summary(self, tree: nx.DiGraph, lsa_results: pd.DataFrame):
        """Save analysis summary"""
        summary = {
            'timestamp': datetime.now().isoformat(),
            'tree_stats': {
                'n_nodes': tree.number_of_nodes(),
                'n_edges': tree.number_of_edges(),
                'max_depth': max(nx.shortest_path_length(tree, 'root').values()),
                'n_leaves': len([n for n in tree.nodes() if tree.out_degree(n) == 0])
            },
            'lsa_stats': {
                'n_significant': len(lsa_results[lsa_results['adjustp'] < 0.05]) if not lsa_results.empty else 0,
                'n_amplifications': len(lsa_results[lsa_results['CNA'] == 'AMP']) if not lsa_results.empty else 0,
                'n_deletions': len(lsa_results[lsa_results['CNA'] == 'DEL']) if not lsa_results.empty else 0
            }
        }
        
        summary_file = os.path.join(self.output_dir, 'analysis_summary.json')
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)


def main():
    """Command line interface for MEDALT"""
    parser = argparse.ArgumentParser(description='MEDALT: Minimal Event Distance Aneuploidy Lineage Tree')
    
    # Required arguments
    parser.add_argument('-I', '--input', required=True,
                       help='Input CNV file')
    parser.add_argument('-D', '--datatype', required=True, choices=['D', 'R'],
                       help='Data type: D for scDNA-seq, R for scRNA-seq')
    parser.add_argument('-G', '--genome', required=True, choices=['hg19', 'hg38'],
                       help='Genome version')
    
    # Optional arguments
    parser.add_argument('-O', '--outpath', default='medalt_output',
                       help='Output directory (default: medalt_output)')
    parser.add_argument('-W', '--windows', type=int, default=30,
                       help='Window size for scRNA-seq (default: 30)')
    parser.add_argument('-R', '--permutation', type=int, default=500,
                       help='Number of permutations (default: 500)')
    parser.add_argument('-M', '--minsize', type=int, default=5,
                       help='Minimum lineage size (default: 5)')
    parser.add_argument('-J', '--jobs', type=int, default=-1,
                       help='Number of parallel jobs (default: all CPUs)')
    
    args = parser.parse_args()
    
    # Run pipeline
    pipeline = MEDALTPipeline(args.outpath)
    tree, results = pipeline.run_analysis(
        input_file=args.input,
        data_type=args.datatype,
        genome_version=args.genome,
        window_size=args.windows,
        n_permutations=args.permutation,
        min_lineage_size=args.minsize,
        n_jobs=args.jobs
    )
    
    # Print summary
    if not results.empty:
        print(f"\nAnalysis complete!")
        print(f"Tree nodes: {tree.number_of_nodes()}")
        print(f"Significant CNAs: {len(results[results['adjustp'] < 0.05])}")
        print(f"\nTop 5 results:")
        print(results.head())
    else:
        print("\nAnalysis complete, but no significant CNAs found.")


if __name__ == '__main__':
    main()