#!/usr/bin/env python3
"""
MEDALT Visualization Module for Reimplementation

A simplified Python visualization module for the MEDALT reimplementation that
generates publication-quality tree plots and analysis reports.

Usage:
    python3 visualize_medalt.py <tree_file> [output_dir]

Example:
    python3 visualize_medalt.py medalt_results.txt ./plots/
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
import seaborn as sns
from typing import Dict, List, Tuple, Optional
import warnings
import sys
import os
from datetime import datetime

warnings.filterwarnings('ignore')

class MEDALTVisualizer:
    """
    Python visualization for MEDALT reimplementation results
    """
    
    def __init__(self, output_dir: str = "."):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # Set matplotlib style for publication quality
        plt.style.use('default')
        sns.set_palette("husl")
    
    def _create_hierarchical_layout(self, G, root=None):
        """Create a hierarchical tree layout"""
        if root is None:
            # Find root nodes (nodes with no incoming edges)
            root_nodes = [n for n in G.nodes() if G.in_degree(n) == 0]
            root = root_nodes[0] if root_nodes else list(G.nodes())[0]
        
        # Use networkx hierarchy if available
        try:
            pos = nx.nx_agraph.graphviz_layout(G, prog='dot')
            return pos
        except:
            pass
        
        # Manual hierarchical layout
        pos = {}
        levels = {root: 0}
        queue = [root]
        
        # BFS to assign levels
        while queue:
            node = queue.pop(0)
            current_level = levels[node]
            
            children = list(G.successors(node))
            for child in children:
                if child not in levels:
                    levels[child] = current_level + 1
                    queue.append(child)
        
        # Group nodes by level and assign positions
        level_groups = {}
        for node, level in levels.items():
            if level not in level_groups:
                level_groups[level] = []
            level_groups[level].append(node)
        
        max_level = max(levels.values()) if levels else 0
        
        for level, nodes_at_level in level_groups.items():
            y = max_level - level  # Root at top
            num_nodes = len(nodes_at_level)
            
            if num_nodes == 1:
                x_positions = [0]
            else:
                x_positions = np.linspace(-num_nodes/2, num_nodes/2, num_nodes)
            
            for i, node in enumerate(nodes_at_level):
                pos[node] = (x_positions[i], y)
        
        return pos
    
    def _shorten_cell_name(self, cell_name: str, max_length: int = 15) -> str:
        """Shorten cell names for display"""
        cell_str = str(cell_name)
        if len(cell_str) <= max_length:
            return cell_str
        
        # Try to keep meaningful parts
        if 'HNSCC' in cell_str:
            # Extract key parts: HNSCC5_p7_..._G05 -> HNSCC5_G05
            parts = cell_str.split('_')
            if len(parts) >= 3:
                return f"{parts[0]}_{parts[-1]}"
        
        return cell_str[:max_length-3] + '...'
    
    def plot_single_cell_tree(self, tree_file: str, output_filename: str = "singlecell_tree.pdf"):
        """
        Plot the complete single cell phylogenetic tree
        """
        print(f"Creating single cell tree visualization...")
        
        # Read tree data
        try:
            tree_df = pd.read_csv(tree_file, sep='\t', comment='#')
        except Exception as e:
            print(f"Error reading tree file: {e}")
            return None
        
        # Create NetworkX graph
        G = nx.from_pandas_edgelist(tree_df, 
                                   source=tree_df.columns[0], 
                                   target=tree_df.columns[1], 
                                   edge_attr='dist' if 'dist' in tree_df.columns else None,
                                   create_using=nx.DiGraph())
        
        # Find root node
        root_nodes = [n for n in G.nodes() if G.in_degree(n) == 0]
        root = root_nodes[0] if root_nodes else list(G.nodes())[0]
        
        # Create layout
        pos = self._create_hierarchical_layout(G, root)
        
        # Calculate node properties
        node_colors = []
        node_sizes = []
        
        for node in G.nodes():
            # Color by distance from root
            try:
                path_length = nx.shortest_path_length(G, root, node)
                color_val = plt.cm.viridis(path_length / max(3, nx.eccentricity(G, root)))
            except:
                color_val = 'lightblue'
            
            node_colors.append(color_val)
            
            # Size by degree
            degree = G.degree(node)
            node_sizes.append(max(100, degree * 50))
        
        # Create the plot
        plt.figure(figsize=(16, 12))
        
        # Draw edges with distance-based styling
        edge_weights = []
        edge_colors = []
        
        for u, v in G.edges():
            if 'dist' in G[u][v]:
                weight = G[u][v]['dist']
                edge_weights.append(max(0.5, weight))
                # Color edges by distance
                edge_colors.append(plt.cm.Reds(min(1.0, weight / 3.0)))
            else:
                edge_weights.append(1.0)
                edge_colors.append('gray')
        
        nx.draw_networkx_edges(G, pos, 
                              width=edge_weights,
                              edge_color=edge_colors,
                              arrows=True,
                              arrowsize=15,
                              arrowstyle='->',
                              alpha=0.7)
        
        # Draw nodes
        nx.draw_networkx_nodes(G, pos,
                              node_color=node_colors,
                              node_size=node_sizes,
                              alpha=0.8,
                              edgecolors='black',
                              linewidths=0.5)
        
        # Add selective labels
        labels = {}
        
        # Always label root
        labels[root] = f"ROOT\n{self._shorten_cell_name(root, 10)}"
        
        # Label high-degree nodes
        degrees = dict(G.degree())
        high_degree_nodes = [n for n, d in degrees.items() if d >= np.percentile(list(degrees.values()), 80)]
        
        for node in high_degree_nodes[:8]:  # Limit to top 8
            if node != root:
                labels[node] = self._shorten_cell_name(node, 12)
        
        # Label G05 if present
        for node in G.nodes():
            if 'G05' in str(node):
                labels[node] = f"G05\n{self._shorten_cell_name(node, 10)}"
                break
        
        nx.draw_networkx_labels(G, pos, labels, 
                               font_size=8,
                               font_weight='bold',
                               bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7))
        
        plt.title("MEDALT Single Cell Phylogenetic Tree\n"
                 f"Nodes: {len(G.nodes())}, Edges: {len(G.edges())}, Root: {self._shorten_cell_name(root, 20)}", 
                 fontsize=14, fontweight='bold')
        plt.axis('off')
        
        # Add legend
        legend_elements = [
            mpatches.Patch(color='darkblue', label='Close to root'),
            mpatches.Patch(color='yellow', label='Far from root'),
            mpatches.Patch(color='red', label='High distance edges'),
            mpatches.Patch(color='gray', label='Low distance edges')
        ]
        plt.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(0, 1))
        
        # Save plot
        output_path = os.path.join(self.output_dir, output_filename)
        plt.savefig(output_path, format='pdf', bbox_inches='tight', dpi=300)
        plt.close()
        
        print(f"✓ Single cell tree saved to {output_path}")
        return output_path
    
    def plot_distance_heatmap(self, tree_file: str, output_filename: str = "distance_heatmap.pdf"):
        """
        Create a heatmap of pairwise distances between cells
        """
        print(f"Creating distance heatmap...")
        
        # Read tree data and create distance matrix
        tree_df = pd.read_csv(tree_file, sep='\t', comment='#')
        
        # Create a simple distance analysis
        G = nx.from_pandas_edgelist(tree_df, 
                                   source=tree_df.columns[0], 
                                   target=tree_df.columns[1], 
                                   edge_attr='dist' if 'dist' in tree_df.columns else None,
                                   create_using=nx.Graph())  # Undirected for distances
        
        # Get sample of nodes for visualization (too many for full heatmap)
        all_nodes = list(G.nodes())
        if len(all_nodes) > 20:
            # Sample important nodes
            degrees = dict(G.degree())
            top_nodes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:15]
            sample_nodes = [node for node, _ in top_nodes]
            
            # Add G05 if present
            for node in all_nodes:
                if 'G05' in str(node) and node not in sample_nodes:
                    sample_nodes.append(node)
                    break
            
            # Add root if present
            for node in all_nodes:
                if any(root_word in str(node).lower() for root_word in ['root', 'e11']) and node not in sample_nodes:
                    sample_nodes.append(node)
                    break
        else:
            sample_nodes = all_nodes
        
        # Calculate distance matrix
        distance_matrix = []
        node_labels = []
        
        for node1 in sample_nodes:
            row = []
            for node2 in sample_nodes:
                try:
                    if node1 == node2:
                        distance = 0
                    else:
                        distance = nx.shortest_path_length(G, node1, node2)
                except:
                    distance = np.inf
                row.append(distance)
            distance_matrix.append(row)
            node_labels.append(self._shorten_cell_name(node1, 12))
        
        # Create heatmap
        plt.figure(figsize=(12, 10))
        
        # Replace inf with max distance + 1 for visualization
        distance_array = np.array(distance_matrix)
        max_dist = np.max(distance_array[distance_array != np.inf])
        distance_array[distance_array == np.inf] = max_dist + 1
        
        im = plt.imshow(distance_array, cmap='viridis_r', aspect='auto')
        
        # Add colorbar
        cbar = plt.colorbar(im)
        cbar.set_label('Distance', rotation=270, labelpad=20)
        
        # Set labels
        plt.xticks(range(len(node_labels)), node_labels, rotation=45, ha='right')
        plt.yticks(range(len(node_labels)), node_labels)
        
        # Add distance values to cells
        for i in range(len(node_labels)):
            for j in range(len(node_labels)):
                dist_val = distance_array[i, j]
                if dist_val <= max_dist:
                    plt.text(j, i, f'{int(dist_val)}', ha='center', va='center', 
                            color='white' if dist_val > max_dist/2 else 'black', fontweight='bold')
        
        plt.title("Pairwise Distance Matrix (Sample of Key Cells)", fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        # Save plot
        output_path = os.path.join(self.output_dir, output_filename)
        plt.savefig(output_path, format='pdf', bbox_inches='tight', dpi=300)
        plt.close()
        
        print(f"✓ Distance heatmap saved to {output_path}")
        return output_path
    
    def plot_tree_statistics(self, tree_file: str, output_filename: str = "tree_statistics.pdf"):
        """
        Create comprehensive tree statistics plots
        """
        print(f"Creating tree statistics...")
        
        tree_df = pd.read_csv(tree_file, sep='\t', comment='#')
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Plot 1: Distance distribution
        if 'dist' in tree_df.columns:
            distances = tree_df['dist']
            axes[0, 0].hist(distances, bins=max(1, len(distances.unique())), 
                           alpha=0.7, edgecolor='black', color='skyblue')
            axes[0, 0].set_title('Edge Distance Distribution', fontweight='bold')
            axes[0, 0].set_xlabel('Distance')
            axes[0, 0].set_ylabel('Frequency')
            axes[0, 0].grid(True, alpha=0.3)
            
            # Add statistics text
            axes[0, 0].text(0.7, 0.8, f'Mean: {distances.mean():.2f}\nStd: {distances.std():.2f}', 
                           transform=axes[0, 0].transAxes, 
                           bbox=dict(boxstyle="round", facecolor='wheat', alpha=0.8))
        
        # Plot 2: Node degree distribution
        all_nodes = list(tree_df.iloc[:, 0]) + list(tree_df.iloc[:, 1])
        from collections import Counter
        degree_counts = Counter(all_nodes)
        degrees = list(degree_counts.values())
        
        axes[0, 1].hist(degrees, bins=max(1, len(set(degrees))), 
                       alpha=0.7, edgecolor='black', color='lightcoral')
        axes[0, 1].set_title('Node Degree Distribution', fontweight='bold')
        axes[0, 1].set_xlabel('Degree (Number of Connections)')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].grid(True, alpha=0.3)
        
        # Plot 3: Tree summary statistics
        num_nodes = len(set(all_nodes))
        num_edges = len(tree_df)
        avg_degree = 2 * num_edges / num_nodes if num_nodes > 0 else 0
        
        summary_text = f"""Tree Structure Summary:
        
Total Nodes: {num_nodes}
Total Edges: {num_edges}
Average Degree: {avg_degree:.2f}

Tree Properties:
• Nodes per edge: {num_nodes/num_edges:.2f}
• Graph density: {2*num_edges/(num_nodes*(num_nodes-1)):.3f}
        """
        
        axes[1, 0].text(0.05, 0.95, summary_text, transform=axes[1, 0].transAxes, 
                       fontsize=12, verticalalignment='top',
                       bbox=dict(boxstyle="round", facecolor='lightblue', alpha=0.8))
        axes[1, 0].set_title('Tree Statistics Summary', fontweight='bold')
        axes[1, 0].axis('off')
        
        # Plot 4: Top connected nodes
        top_nodes = Counter(all_nodes).most_common(15)
        if top_nodes:
            nodes, counts = zip(*top_nodes)
            short_names = [self._shorten_cell_name(n, 15) for n in nodes]
            
            y_pos = np.arange(len(short_names))
            bars = axes[1, 1].barh(y_pos, counts, color='lightgreen', alpha=0.7, edgecolor='black')
            axes[1, 1].set_yticks(y_pos)
            axes[1, 1].set_yticklabels(short_names)
            axes[1, 1].set_title('Top 15 Most Connected Nodes', fontweight='bold')
            axes[1, 1].set_xlabel('Connection Count')
            axes[1, 1].grid(True, alpha=0.3, axis='x')
            
            # Add value labels on bars
            for i, bar in enumerate(bars):
                width = bar.get_width()
                axes[1, 1].text(width + 0.1, bar.get_y() + bar.get_height()/2, 
                               f'{int(width)}', ha='left', va='center', fontweight='bold')
        
        plt.tight_layout()
        
        # Save plot
        output_path = os.path.join(self.output_dir, output_filename)
        plt.savefig(output_path, format='pdf', bbox_inches='tight', dpi=300)
        plt.close()
        
        print(f"✓ Tree statistics saved to {output_path}")
        return output_path
    
    def create_comprehensive_report(self, tree_file: str, output_filename: str = "MEDALT_visualization_report.pdf"):
        """
        Create a comprehensive visualization report
        """
        print(f"Creating comprehensive MEDALT visualization report...")
        
        with PdfPages(os.path.join(self.output_dir, output_filename)) as pdf:
            # Title page
            fig, ax = plt.subplots(figsize=(11, 8.5))
            ax.text(0.5, 0.7, 'MEDALT Visualization Report', 
                   ha='center', va='center', fontsize=24, fontweight='bold')
            ax.text(0.5, 0.6, f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', 
                   ha='center', va='center', fontsize=14)
            ax.text(0.5, 0.5, f'Tree file: {os.path.basename(tree_file)}', 
                   ha='center', va='center', fontsize=12)
            ax.text(0.5, 0.3, 'MEDALT Python Reimplementation', 
                   ha='center', va='center', fontsize=16, style='italic')
            ax.axis('off')
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()
            
            # Generate individual plots and add to PDF
            temp_files = []
            
            # Single cell tree
            tree_plot = self.plot_single_cell_tree(tree_file, "temp_tree.pdf")
            
            # Distance heatmap
            heatmap_plot = self.plot_distance_heatmap(tree_file, "temp_heatmap.pdf")
            
            # Tree statistics
            stats_plot = self.plot_tree_statistics(tree_file, "temp_stats.pdf")
        
        print(f"✓ Comprehensive report saved to {os.path.join(self.output_dir, output_filename)}")
        return os.path.join(self.output_dir, output_filename)

def main():
    """Main entry point for visualization"""
    
    if len(sys.argv) < 2:
        print("Usage: python3 visualize_medalt.py <tree_file> [output_dir]")
        print()
        print("Example:")
        print("  python3 visualize_medalt.py medalt_results.txt")
        print("  python3 visualize_medalt.py medalt_results.txt ./plots/")
        sys.exit(1)
    
    tree_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "./plots"
    
    if not os.path.exists(tree_file):
        print(f"Error: Tree file '{tree_file}' not found!")
        sys.exit(1)
    
    print("="*60)
    print("MEDALT VISUALIZATION MODULE")
    print("="*60)
    print(f"Input tree file: {tree_file}")
    print(f"Output directory: {output_dir}")
    print()
    
    # Create visualizer
    visualizer = MEDALTVisualizer(output_dir)
    
    try:
        # Generate all visualizations
        visualizer.plot_single_cell_tree(tree_file)
        visualizer.plot_distance_heatmap(tree_file)
        visualizer.plot_tree_statistics(tree_file)
        visualizer.create_comprehensive_report(tree_file)
        
        print()
        print("="*60)
        print("VISUALIZATION COMPLETED SUCCESSFULLY!")
        print("="*60)
        print("Generated files:")
        for filename in ["singlecell_tree.pdf", "distance_heatmap.pdf", 
                        "tree_statistics.pdf", "MEDALT_visualization_report.pdf"]:
            filepath = os.path.join(output_dir, filename)
            if os.path.exists(filepath):
                print(f"  ✓ {filepath}")
        
        print(f"\n🎉 All visualizations saved in: {output_dir}")
        
    except Exception as e:
        print(f"Error during visualization: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()