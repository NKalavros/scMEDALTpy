#!/usr/bin/env python
"""
MEDALT Visualization Module
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
        
        Args:
            output_file: Output filename
            figsize: Figure size
            show_distances: Whether to show edge distances
            node_size: Size of nodes
        """
        fig, ax = plt.subplots(figsize=figsize)
        
        # Create hierarchical layout
        pos = self._create_hierarchical_layout()
        
        # Get node colors based on significance
        node_colors = self._get_node_colors()
        
        # Draw edges first
        edge_colors = ['#666666'] * len(self.tree.edges())
        nx.draw_networkx_edges(self.tree, pos, 
                             edge_color=edge_colors,
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
        
        # Add edge labels (distances) if requested
        if show_distances:
            edge_labels = nx.get_edge_attributes(self.tree, 'weight')
            edge_labels = {k: f'{int(v)}' for k, v in edge_labels.items()}
            nx.draw_networkx_edge_labels(self.tree, pos, 
                                        edge_labels=edge_labels,
                                        font_size=6,
                                        ax=ax)
        
        # Add CNA annotations
        self._add_cna_annotations(ax, pos)
        
        # Set title and remove axes
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
        
        Args:
            output_file: Output filename
            figsize: Figure size
            significance_threshold: Significance threshold for highlighting
        """
        if self.lsa_results is None or self.lsa_results.empty:
            print("No LSA results available for visualization")
            return
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Create hierarchical layout
        pos = self._create_hierarchical_layout()
        
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
        
        # Draw nodes with LSA-based coloring
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
        
        # Set title and remove axes
        ax.set_title('MEDALT LSA Tree - Significant CNAs', fontsize=16, fontweight='bold')
        ax.set_aspect('equal')
        ax.axis('off')
        
        # Add LSA legend
        self._add_lsa_legend(ax, significance_threshold)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight', format='pdf')
        plt.close()
        
        print(f"LSA tree visualization saved to {output_file}")
    
    def _create_hierarchical_layout(self) -> Dict:
        """Create hierarchical layout for tree visualization"""
        # Use graphviz layout if available, otherwise fallback to spring layout
        try:
            pos = nx.nx_agraph.graphviz_layout(self.tree, prog='dot')
        except:
            # Fallback to hierarchical spring layout
            pos = nx.spring_layout(self.tree, k=3, iterations=50)
            
            # Adjust positions to be more hierarchical
            if 'root' in pos:
                # Get depths
                depths = nx.single_source_shortest_path_length(self.tree, 'root')
                
                # Adjust y-coordinates based on depth
                max_depth = max(depths.values()) if depths else 1
                for node in pos:
                    if node in depths:
                        pos[node] = (pos[node][0], max_depth - depths[node])
        
        return pos\n    \n    def _get_node_colors(self) -> List[str]:\n        \"\"\"Get node colors for basic tree visualization\"\"\"\n        colors = []\n        \n        for node in self.tree.nodes():\n            if node == 'root':\n                colors.append('#333333')  # Dark gray for root\n            else:\n                colors.append('#87CEEB')  # Light blue for cells\n        \n        return colors\n    \n    def _get_lsa_node_colors(self, significance_threshold: float) -> List[str]:\n        \"\"\"Get node colors based on LSA significance\"\"\"\n        colors = []\n        \n        # Create color map for significant nodes\n        significant_nodes = set()\n        if self.lsa_results is not None and not self.lsa_results.empty:\n            sig_results = self.lsa_results[self.lsa_results['adjustp'] < significance_threshold]\n            significant_nodes = set(sig_results['cell'].unique())\n        \n        for node in self.tree.nodes():\n            if node == 'root':\n                colors.append('#333333')  # Dark gray for root\n            elif self.cell_name_map.get(node, '') in significant_nodes:\n                colors.append('#FF6B6B')  # Red for significant nodes\n            else:\n                colors.append('#87CEEB')  # Light blue for normal cells\n        \n        return colors\n    \n    def _get_node_label(self, node) -> str:\n        \"\"\"Get label for a node\"\"\"\n        if node == 'root':\n            return 'Root'\n        elif isinstance(node, int) and node < len(self.cell_names):\n            # Shorten long cell names\n            name = self.cell_names[node]\n            if len(name) > 12:\n                return name[:9] + '...'\n            return name\n        else:\n            return str(node)\n    \n    def _add_cna_annotations(self, ax, pos) -> None:\n        \"\"\"Add CNA annotations to nodes\"\"\"\n        if self.lsa_results is None or self.lsa_results.empty:\n            return\n        \n        # Get top CNAs per cell\n        for node in self.tree.nodes():\n            if node == 'root':\n                continue\n            \n            cell_name = self.cell_name_map.get(node, '')\n            if not cell_name:\n                continue\n            \n            # Get significant CNAs for this cell\n            cell_cnas = self.lsa_results[\n                (self.lsa_results['cell'] == cell_name) & \n                (self.lsa_results['adjustp'] < 0.05)\n            ].head(3)  # Top 3 CNAs\n            \n            if not cell_cnas.empty:\n                # Create annotation text\n                annotations = []\n                for _, cna in cell_cnas.iterrows():\n                    cna_text = f\"{cna['region']}({cna['CNA']})\"\n                    annotations.append(cna_text)\n                \n                if annotations:\n                    ann_text = '\\n'.join(annotations)\n                    x, y = pos[node]\n                    ax.annotate(ann_text, \n                              xy=(x, y), \n                              xytext=(x + 0.15, y + 0.15),\n                              fontsize=6,\n                              bbox=dict(boxstyle='round,pad=0.3', \n                                       facecolor='yellow', \n                                       alpha=0.7),\n                              arrowprops=dict(arrowstyle='->', \n                                            connectionstyle='arc3,rad=0.2'))\n    \n    def _add_lsa_annotations(self, ax, pos, significance_threshold: float) -> None:\n        \"\"\"Add LSA-specific annotations\"\"\"\n        if self.lsa_results is None or self.lsa_results.empty:\n            return\n        \n        sig_results = self.lsa_results[self.lsa_results['adjustp'] < significance_threshold]\n        \n        for node in self.tree.nodes():\n            if node == 'root':\n                continue\n            \n            cell_name = self.cell_name_map.get(node, '')\n            if not cell_name:\n                continue\n            \n            # Get significant CNAs for this cell\n            cell_cnas = sig_results[sig_results['cell'] == cell_name]\n            \n            if not cell_cnas.empty:\n                # Count AMPs and DELs\n                amp_count = len(cell_cnas[cell_cnas['CNA'] == 'AMP'])\n                del_count = len(cell_cnas[cell_cnas['CNA'] == 'DEL'])\n                \n                # Create summary text\n                summary_parts = []\n                if amp_count > 0:\n                    summary_parts.append(f\"AMP: {amp_count}\")\n                if del_count > 0:\n                    summary_parts.append(f\"DEL: {del_count}\")\n                \n                if summary_parts:\n                    summary_text = '\\n'.join(summary_parts)\n                    x, y = pos[node]\n                    ax.annotate(summary_text, \n                              xy=(x, y), \n                              xytext=(x + 0.2, y + 0.2),\n                              fontsize=7,\n                              fontweight='bold',\n                              bbox=dict(boxstyle='round,pad=0.3', \n                                       facecolor='orange', \n                                       alpha=0.8),\n                              arrowprops=dict(arrowstyle='->', \n                                            connectionstyle='arc3,rad=0.2'))\n    \n    def _add_legend(self, ax) -> None:\n        \"\"\"Add legend for basic tree visualization\"\"\"\n        legend_elements = [\n            plt.Line2D([0], [0], marker='o', color='w', \n                      markerfacecolor='#333333', markersize=10, \n                      label='Root'),\n            plt.Line2D([0], [0], marker='o', color='w', \n                      markerfacecolor='#87CEEB', markersize=10, \n                      label='Cell')\n        ]\n        \n        ax.legend(handles=legend_elements, loc='upper right', \n                 bbox_to_anchor=(1, 1))\n    \n    def _add_lsa_legend(self, ax, significance_threshold: float) -> None:\n        \"\"\"Add legend for LSA tree visualization\"\"\"\n        legend_elements = [\n            plt.Line2D([0], [0], marker='o', color='w', \n                      markerfacecolor='#333333', markersize=10, \n                      label='Root'),\n            plt.Line2D([0], [0], marker='o', color='w', \n                      markerfacecolor='#87CEEB', markersize=10, \n                      label='Normal Cell'),\n            plt.Line2D([0], [0], marker='o', color='w', \n                      markerfacecolor='#FF6B6B', markersize=10, \n                      label=f'Significant CNAs (p < {significance_threshold})')\n        ]\n        \n        ax.legend(handles=legend_elements, loc='upper right', \n                 bbox_to_anchor=(1, 1))


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
        \"\"\"
        Plot CNV heatmap ordered by tree structure
        
        Args:
            output_file: Output filename
            figsize: Figure size
            order_by_tree: Whether to order cells by tree structure
        \"\"\"
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
        # 0=white (homozygous deletion), 1=blue (loss), 2=gray (diploid), 3+=red (amplification)
        colors = ['white', '#4575b4', '#999999', '#f46d43', '#d73027', '#a50026']
        bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 6.5]
        cmap = ListedColormap(colors)
        
        # Create heatmap
        im = ax.imshow(ordered_matrix, cmap=cmap, aspect='auto', 
                      vmin=0, vmax=6, interpolation='nearest')
        
        # Set labels
        ax.set_xlabel('Genes', fontsize=12)
        ax.set_ylabel('Cells', fontsize=12)
        ax.set_title('Copy Number Variation Heatmap', fontsize=16, fontweight='bold')
        
        # Set cell labels on y-axis
        if len(ordered_names) <= 50:  # Only show labels if not too many cells
            ax.set_yticks(range(len(ordered_names)))
            ax.set_yticklabels(ordered_names, fontsize=8)
        else:
            ax.set_yticks([])\n        \n        # Set gene labels on x-axis (only show subset)\n        if len(self.gene_names) <= 50:\n            ax.set_xticks(range(len(self.gene_names)))\n            ax.set_xticklabels(self.gene_names, rotation=90, fontsize=8)\n        else:\n            # Show every 10th gene\n            step = max(1, len(self.gene_names) // 20)\n            tick_positions = range(0, len(self.gene_names), step)\n            ax.set_xticks(tick_positions)\n            ax.set_xticklabels([self.gene_names[i] for i in tick_positions], \n                              rotation=90, fontsize=8)\n        \n        # Add colorbar\n        cbar = plt.colorbar(im, ax=ax, shrink=0.6)\n        cbar.set_label('Copy Number', rotation=270, labelpad=15)\n        cbar.set_ticks([0, 1, 2, 3, 4, 5, 6])\n        cbar.set_ticklabels(['0', '1', '2', '3', '4', '5', '6+'])\n        \n        plt.tight_layout()\n        plt.savefig(output_file, dpi=300, bbox_inches='tight', format='pdf')\n        plt.close()\n        \n        print(f\"CNV heatmap saved to {output_file}\")\n    \n    def _get_tree_ordered_indices(self) -> List[int]:\n        \"\"\"Get cell indices ordered by tree traversal\"\"\"\n        if self.tree is None:\n            return list(range(len(self.cell_names)))\n        \n        ordered_indices = []\n        \n        def dfs_traverse(node):\n            if node != 'root' and isinstance(node, int):\n                ordered_indices.append(node)\n            \n            # Get children and sort them for consistent ordering\n            children = sorted(self.tree.successors(node))\n            for child in children:\n                dfs_traverse(child)\n        \n        dfs_traverse('root')\n        \n        # Add any missing indices\n        for i in range(len(self.cell_names)):\n            if i not in ordered_indices:\n                ordered_indices.append(i)\n        \n        return ordered_indices


class MEDALTManhattanPlot:
    \"\"\"Manhattan plot for LSA results\"\"\"\n    \n    def __init__(self, lsa_results: pd.DataFrame, gene_names: List[str] = None):\n        self.lsa_results = lsa_results\n        self.gene_names = gene_names\n    \n    def plot_manhattan(self, output_file: str = 'lsa_manhattan.pdf',\n                      figsize: Tuple[int, int] = (15, 8),\n                      significance_threshold: float = 0.05) -> None:\n        \"\"\"\n        Create Manhattan plot for LSA results\n        \n        Args:\n            output_file: Output filename\n            figsize: Figure size\n            significance_threshold: Significance threshold line\n        \"\"\"\n        if self.lsa_results.empty:\n            print(\"No LSA results available for Manhattan plot\")\n            return\n        \n        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, \n                                       gridspec_kw={'height_ratios': [1, 1]})\n        \n        # Separate amplifications and deletions\n        amp_results = self.lsa_results[self.lsa_results['CNA'] == 'AMP'].copy()\n        del_results = self.lsa_results[self.lsa_results['CNA'] == 'DEL'].copy()\n        \n        # Add x-axis positions\n        if self.gene_names:\n            # Map gene names to positions\n            gene_to_pos = {gene: i for i, gene in enumerate(self.gene_names)}\n            \n            if not amp_results.empty:\n                amp_results['x_pos'] = amp_results['region'].map(gene_to_pos)\n                amp_results = amp_results.dropna(subset=['x_pos'])\n            \n            if not del_results.empty:\n                del_results['x_pos'] = del_results['region'].map(gene_to_pos)\n                del_results = del_results.dropna(subset=['x_pos'])\n        else:\n            # Use index positions\n            if not amp_results.empty:\n                amp_results['x_pos'] = range(len(amp_results))\n            if not del_results.empty:\n                del_results['x_pos'] = range(len(del_results))\n        \n        # Plot amplifications\n        if not amp_results.empty:\n            ax1.scatter(amp_results['x_pos'], -np.log10(amp_results['adjustp']),\n                       c='red', alpha=0.7, s=30, label='Amplifications')\n            ax1.axhline(y=-np.log10(significance_threshold), color='black',\n                       linestyle='--', alpha=0.5, \n                       label=f'Significance (p = {significance_threshold})')\n        \n        ax1.set_ylabel('-log10(adjusted p-value)', fontsize=12)\n        ax1.set_title('Amplifications', fontsize=14, fontweight='bold')\n        ax1.legend()\n        ax1.grid(True, alpha=0.3)\n        \n        # Plot deletions\n        if not del_results.empty:\n            ax2.scatter(del_results['x_pos'], -np.log10(del_results['adjustp']),\n                       c='blue', alpha=0.7, s=30, label='Deletions')\n            ax2.axhline(y=-np.log10(significance_threshold), color='black',\n                       linestyle='--', alpha=0.5, \n                       label=f'Significance (p = {significance_threshold})')\n        \n        ax2.set_xlabel('Genomic Position', fontsize=12)\n        ax2.set_ylabel('-log10(adjusted p-value)', fontsize=12)\n        ax2.set_title('Deletions', fontsize=14, fontweight='bold')\n        ax2.legend()\n        ax2.grid(True, alpha=0.3)\n        \n        plt.suptitle('LSA Results: Manhattan Plot', fontsize=16, fontweight='bold')\n        plt.tight_layout()\n        plt.savefig(output_file, dpi=300, bbox_inches='tight', format='pdf')\n        plt.close()\n        \n        print(f\"Manhattan plot saved to {output_file}\")