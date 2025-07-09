#!/usr/bin/env python
"""
MEDALT Visualization Module - Fixed version
Implements visualizations matching the original MEDALT output
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap
import seaborn as sns
import networkx as nx
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

class MEDALTTreeVisualizer:
    """Tree visualization matching original MEDALT style"""
    
    def __init__(self, tree: nx.DiGraph, cnv_matrix: np.ndarray, 
                 cell_names: List[str], gene_names: List[str],
                 lsa_results: pd.DataFrame = None):
        """
        Initialize visualizer
        
        Args:
            tree: NetworkX directed graph from MEDALT
            cnv_matrix: Copy number matrix
            cell_names: List of cell names
            gene_names: List of gene names
            lsa_results: LSA results DataFrame
        """
        self.tree = tree
        self.cnv_matrix = cnv_matrix
        self.cell_names = cell_names
        self.gene_names = gene_names
        self.lsa_results = lsa_results
        
        # Create cell name mapping
        self.cell_name_map = {i: name for i, name in enumerate(cell_names)}
        self.cell_name_map['root'] = 'root'
    
    def plot_tree(self, output_file: str = 'singlecell.tree.pdf', 
                  figsize: Tuple[int, int] = (15, 12),
                  show_distances: bool = True,
                  node_size: int = 200) -> None:
        """
        Plot tree visualization matching original MEDALT style
        """
        fig, ax = plt.subplots(figsize=figsize)
        
        # Create layout
        pos = self._create_layout()
        
        # Get node colors
        node_colors = self._get_node_colors()
        
        # Draw edges
        nx.draw_networkx_edges(self.tree, pos, 
                             edge_color='#666666',
                             arrows=True, 
                             arrowsize=20,
                             arrowstyle='->', 
                             width=1.5,
                             alpha=0.7,
                             ax=ax)
        
        # Draw nodes
        nx.draw_networkx_nodes(self.tree, pos,
                             node_color=node_colors,
                             node_size=node_size,
                             alpha=0.8,
                             linewidths=1,
                             edgecolors='black',
                             ax=ax)
        
        # Add node labels
        node_labels = {node: self._get_node_label(node) for node in self.tree.nodes()}
        nx.draw_networkx_labels(self.tree, pos, 
                              labels=node_labels,
                              font_size=8,
                              font_weight='bold',
                              ax=ax)
        
        # Add edge labels if requested
        if show_distances:
            edge_labels = nx.get_edge_attributes(self.tree, 'weight')
            edge_labels = {k: f'{int(v)}' for k, v in edge_labels.items()}
            nx.draw_networkx_edge_labels(self.tree, pos, 
                                        edge_labels=edge_labels,
                                        font_size=6,
                                        ax=ax)
        
        # Set title and styling
        ax.set_title('MEDALT Single-Cell Lineage Tree', fontsize=16, fontweight='bold')
        ax.set_aspect('equal')
        ax.axis('off')
        
        # Add legend
        self._add_legend(ax)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight', format='pdf')
        plt.close()
        
        print(f"Tree visualization saved to {output_file}")
    
    def plot_lsa_tree(self, output_file: str = 'LSA.tree.pdf',
                      figsize: Tuple[int, int] = (15, 12),
                      significance_threshold: float = 0.05) -> None:
        """
        Plot LSA tree with significant CNAs highlighted
        """
        if self.lsa_results is None or self.lsa_results.empty:
            print("No LSA results available for visualization")
            return
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Create layout
        pos = self._create_layout()
        
        # Get node colors based on LSA significance
        node_colors = self._get_lsa_node_colors(significance_threshold)
        
        # Draw edges
        nx.draw_networkx_edges(self.tree, pos, 
                             edge_color='#666666',
                             arrows=True, 
                             arrowsize=20,
                             arrowstyle='->', 
                             width=1.5,
                             alpha=0.7,
                             ax=ax)
        
        # Draw nodes
        nx.draw_networkx_nodes(self.tree, pos,
                             node_color=node_colors,
                             node_size=300,
                             alpha=0.8,
                             linewidths=2,
                             edgecolors='black',
                             ax=ax)
        
        # Add node labels
        node_labels = {node: self._get_node_label(node) for node in self.tree.nodes()}
        nx.draw_networkx_labels(self.tree, pos, 
                              labels=node_labels,
                              font_size=8,
                              font_weight='bold',
                              ax=ax)
        
        # Add LSA annotations
        self._add_lsa_annotations(ax, pos, significance_threshold)
        
        # Set title and styling
        ax.set_title('MEDALT LSA Tree - Significant CNAs', fontsize=16, fontweight='bold')
        ax.set_aspect('equal')
        ax.axis('off')
        
        # Add LSA legend
        self._add_lsa_legend(ax, significance_threshold)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight', format='pdf')
        plt.close()
        
        print(f"LSA tree visualization saved to {output_file}")
    
    def _create_layout(self) -> Dict:
        """Create hierarchical layout for tree visualization"""
        # Try graphviz layout first
        try:
            pos = nx.nx_agraph.graphviz_layout(self.tree, prog='dot')
        except:
            # Fallback to spring layout
            pos = nx.spring_layout(self.tree, k=3, iterations=50, seed=42)
            
            # Make it more hierarchical
            if 'root' in pos:
                depths = nx.single_source_shortest_path_length(self.tree, 'root')
                max_depth = max(depths.values()) if depths else 1
                
                for node in pos:
                    if node in depths:
                        pos[node] = (pos[node][0], max_depth - depths[node])
        
        return pos
    
    def _get_node_colors(self) -> List[str]:
        """Get node colors for basic tree visualization"""
        colors = []
        for node in self.tree.nodes():
            if node == 'root':
                colors.append('#333333')
            else:
                colors.append('#87CEEB')
        return colors
    
    def _get_lsa_node_colors(self, significance_threshold: float) -> List[str]:
        """Get node colors based on LSA significance"""
        colors = []
        
        # Find significant nodes
        significant_nodes = set()
        if self.lsa_results is not None and not self.lsa_results.empty:
            sig_results = self.lsa_results[self.lsa_results['adjustp'] < significance_threshold]
            significant_nodes = set(sig_results['cell'].unique())
        
        for node in self.tree.nodes():
            if node == 'root':
                colors.append('#333333')
            elif self.cell_name_map.get(node, '') in significant_nodes:
                colors.append('#FF6B6B')
            else:
                colors.append('#87CEEB')
        
        return colors
    
    def _get_node_label(self, node) -> str:
        """Get label for a node"""
        if node == 'root':
            return 'Root'
        elif isinstance(node, int) and node < len(self.cell_names):
            name = self.cell_names[node]
            if len(name) > 12:
                return name[:9] + '...'
            return name
        else:
            return str(node)
    
    def _add_lsa_annotations(self, ax, pos, significance_threshold: float) -> None:
        """Add LSA-specific annotations"""
        if self.lsa_results is None or self.lsa_results.empty:
            return
        
        sig_results = self.lsa_results[self.lsa_results['adjustp'] < significance_threshold]
        
        for node in self.tree.nodes():
            if node == 'root':
                continue
            
            cell_name = self.cell_name_map.get(node, '')
            if not cell_name:
                continue
            
            # Get significant CNAs for this cell
            cell_cnas = sig_results[sig_results['cell'] == cell_name]
            
            if not cell_cnas.empty:
                # Count AMPs and DELs
                amp_count = len(cell_cnas[cell_cnas['CNA'] == 'AMP'])
                del_count = len(cell_cnas[cell_cnas['CNA'] == 'DEL'])
                
                # Create summary text
                summary_parts = []
                if amp_count > 0:
                    summary_parts.append(f"AMP: {amp_count}")
                if del_count > 0:
                    summary_parts.append(f"DEL: {del_count}")
                
                if summary_parts:
                    summary_text = '\\n'.join(summary_parts)
                    x, y = pos[node]
                    ax.annotate(summary_text, 
                              xy=(x, y), 
                              xytext=(x + 0.2, y + 0.2),
                              fontsize=7,
                              fontweight='bold',
                              bbox=dict(boxstyle='round,pad=0.3', 
                                       facecolor='orange', 
                                       alpha=0.8),
                              arrowprops=dict(arrowstyle='->', 
                                            connectionstyle='arc3,rad=0.2'))
    
    def _add_legend(self, ax) -> None:
        """Add legend for basic tree visualization"""
        legend_elements = [
            plt.Line2D([0], [0], marker='o', color='w', 
                      markerfacecolor='#333333', markersize=10, 
                      label='Root'),
            plt.Line2D([0], [0], marker='o', color='w', 
                      markerfacecolor='#87CEEB', markersize=10, 
                      label='Cell')
        ]
        
        ax.legend(handles=legend_elements, loc='upper right', 
                 bbox_to_anchor=(1, 1))
    
    def _add_lsa_legend(self, ax, significance_threshold: float) -> None:
        """Add legend for LSA tree visualization"""
        legend_elements = [
            plt.Line2D([0], [0], marker='o', color='w', 
                      markerfacecolor='#333333', markersize=10, 
                      label='Root'),
            plt.Line2D([0], [0], marker='o', color='w', 
                      markerfacecolor='#87CEEB', markersize=10, 
                      label='Normal Cell'),
            plt.Line2D([0], [0], marker='o', color='w', 
                      markerfacecolor='#FF6B6B', markersize=10, 
                      label=f'Significant CNAs (p < {significance_threshold})')
        ]
        
        ax.legend(handles=legend_elements, loc='upper right', 
                 bbox_to_anchor=(1, 1))


class MEDALTHeatmapVisualizer:
    """CNV heatmap visualization"""
    
    def __init__(self, cnv_matrix: np.ndarray, cell_names: List[str], 
                 gene_names: List[str], tree: nx.DiGraph = None):
        self.cnv_matrix = cnv_matrix
        self.cell_names = cell_names
        self.gene_names = gene_names
        self.tree = tree
    
    def plot_cnv_heatmap(self, output_file: str = 'cnv_heatmap.pdf',
                        figsize: Tuple[int, int] = (20, 12),
                        order_by_tree: bool = True) -> None:
        """
        Plot CNV heatmap ordered by tree structure
        """
        fig, ax = plt.subplots(figsize=figsize)
        
        # Order cells by tree if available
        if order_by_tree and self.tree is not None:
            ordered_indices = self._get_tree_ordered_indices()
            ordered_matrix = self.cnv_matrix[ordered_indices]
            ordered_names = [self.cell_names[i] for i in ordered_indices]
        else:
            ordered_matrix = self.cnv_matrix
            ordered_names = self.cell_names
        
        # Create custom colormap
        colors = ['white', '#4575b4', '#999999', '#f46d43', '#d73027', '#a50026']
        cmap = ListedColormap(colors)
        
        # Create heatmap
        im = ax.imshow(ordered_matrix, cmap=cmap, aspect='auto', 
                      vmin=0, vmax=6, interpolation='nearest')
        
        # Set labels
        ax.set_xlabel('Genes', fontsize=12)
        ax.set_ylabel('Cells', fontsize=12)
        ax.set_title('Copy Number Variation Heatmap', fontsize=16, fontweight='bold')
        
        # Set cell labels on y-axis
        if len(ordered_names) <= 50:
            ax.set_yticks(range(len(ordered_names)))
            ax.set_yticklabels(ordered_names, fontsize=8)
        else:
            ax.set_yticks([])
        
        # Set gene labels on x-axis
        if len(self.gene_names) <= 50:
            ax.set_xticks(range(len(self.gene_names)))
            ax.set_xticklabels(self.gene_names, rotation=90, fontsize=8)
        else:
            step = max(1, len(self.gene_names) // 20)
            tick_positions = range(0, len(self.gene_names), step)
            ax.set_xticks(tick_positions)
            ax.set_xticklabels([self.gene_names[i] for i in tick_positions], 
                              rotation=90, fontsize=8)
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.6)
        cbar.set_label('Copy Number', rotation=270, labelpad=15)
        cbar.set_ticks([0, 1, 2, 3, 4, 5, 6])
        cbar.set_ticklabels(['0', '1', '2', '3', '4', '5', '6+'])
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight', format='pdf')
        plt.close()
        
        print(f"CNV heatmap saved to {output_file}")
    
    def _get_tree_ordered_indices(self) -> List[int]:
        """Get cell indices ordered by tree traversal"""
        if self.tree is None:
            return list(range(len(self.cell_names)))
        
        ordered_indices = []
        
        def dfs_traverse(node):
            if node != 'root' and isinstance(node, int):
                ordered_indices.append(node)
            
            children = sorted(self.tree.successors(node))
            for child in children:
                dfs_traverse(child)
        
        dfs_traverse('root')
        
        # Add any missing indices
        for i in range(len(self.cell_names)):
            if i not in ordered_indices:
                ordered_indices.append(i)
        
        return ordered_indices


class MEDALTManhattanPlot:
    """Manhattan plot for LSA results"""
    
    def __init__(self, lsa_results: pd.DataFrame, gene_names: List[str] = None):
        self.lsa_results = lsa_results
        self.gene_names = gene_names
    
    def plot_manhattan(self, output_file: str = 'lsa_manhattan.pdf',
                      figsize: Tuple[int, int] = (15, 8),
                      significance_threshold: float = 0.05) -> None:
        """
        Create Manhattan plot for LSA results
        """
        if self.lsa_results.empty:
            print("No LSA results available for Manhattan plot")
            return
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, 
                                       gridspec_kw={'height_ratios': [1, 1]})
        
        # Separate amplifications and deletions
        amp_results = self.lsa_results[self.lsa_results['CNA'] == 'AMP'].copy()
        del_results = self.lsa_results[self.lsa_results['CNA'] == 'DEL'].copy()
        
        # Add x-axis positions
        if self.gene_names:
            gene_to_pos = {gene: i for i, gene in enumerate(self.gene_names)}
            
            if not amp_results.empty:
                amp_results['x_pos'] = amp_results['region'].map(gene_to_pos)
                amp_results = amp_results.dropna(subset=['x_pos'])
            
            if not del_results.empty:
                del_results['x_pos'] = del_results['region'].map(gene_to_pos)
                del_results = del_results.dropna(subset=['x_pos'])
        else:
            if not amp_results.empty:
                amp_results['x_pos'] = range(len(amp_results))
            if not del_results.empty:
                del_results['x_pos'] = range(len(del_results))
        
        # Plot amplifications
        if not amp_results.empty:
            ax1.scatter(amp_results['x_pos'], -np.log10(amp_results['adjustp']),
                       c='red', alpha=0.7, s=30, label='Amplifications')
            ax1.axhline(y=-np.log10(significance_threshold), color='black',
                       linestyle='--', alpha=0.5, 
                       label=f'Significance (p = {significance_threshold})')
        
        ax1.set_ylabel('-log10(adjusted p-value)', fontsize=12)
        ax1.set_title('Amplifications', fontsize=14, fontweight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot deletions
        if not del_results.empty:
            ax2.scatter(del_results['x_pos'], -np.log10(del_results['adjustp']),
                       c='blue', alpha=0.7, s=30, label='Deletions')
            ax2.axhline(y=-np.log10(significance_threshold), color='black',
                       linestyle='--', alpha=0.5, 
                       label=f'Significance (p = {significance_threshold})')
        
        ax2.set_xlabel('Genomic Position', fontsize=12)
        ax2.set_ylabel('-log10(adjusted p-value)', fontsize=12)
        ax2.set_title('Deletions', fontsize=14, fontweight='bold')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.suptitle('LSA Results: Manhattan Plot', fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight', format='pdf')
        plt.close()
        
        print(f"Manhattan plot saved to {output_file}")