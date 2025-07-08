#!/usr/bin/env python3
"""
LSA Visualization for MEDALT

Creates segmental LSA visualizations equivalent to the R implementation.
Generates publication-quality plots showing lineage-associated copy number alterations.

Usage:
    python3 lsa_visualization.py <lsa_file> <tree_file> [output_dir]

Example:
    python3 lsa_visualization.py segmental.LSA.txt tree.txt ./plots/
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import networkx as nx
import sys
import os
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

class LSAVisualizer:
    """
    LSA visualization implementation
    """
    
    def __init__(self, output_dir: str = "./lsa_plots"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # Set plotting style
        plt.style.use('default')
        sns.set_palette("husl")
    
    def _create_lsa_network(self, lsa_df: pd.DataFrame, tree_df: pd.DataFrame) -> pd.DataFrame:
        """
        Create network connections for LSA visualization
        """
        # Get cells with significant LSA events
        lsa_cells = set(lsa_df['cell'].unique())
        
        # Filter tree to include LSA-relevant connections
        lsa_connections = []
        
        for _, row in tree_df.iterrows():
            from_cell = str(row.iloc[0])
            to_cell = str(row.iloc[1])
            distance = row.iloc[2] if len(row) > 2 else 1
            
            # Include connection if either cell has LSA events
            if from_cell in lsa_cells or to_cell in lsa_cells:
                lsa_connections.append({
                    'from': from_cell,
                    'to': to_cell,
                    'weight': max(1, distance)
                })
        
        return pd.DataFrame(lsa_connections)
    
    def _assign_node_properties(self, G: nx.Graph, lsa_df: pd.DataFrame) -> Tuple[List, List, Dict]:
        """
        Assign colors, sizes, and labels to nodes based on LSA events
        """
        node_colors = []
        node_sizes = []
        node_labels = {}
        
        # Count LSA events per cell
        lsa_counts = lsa_df['cell'].value_counts()
        
        # Get unique cells and assign colors
        lsa_cells = lsa_df['cell'].unique()
        colors = plt.cm.Set3(np.linspace(0, 1, len(lsa_cells)))
        cell_color_map = {cell: colors[i] for i, cell in enumerate(lsa_cells)}
        
        for node in G.nodes():
            if node in lsa_counts:
                # LSA cell - color by cell type, size by event count
                count = lsa_counts[node]
                color = cell_color_map[node]
                size = max(200, count * 100)  # Scale size by number of events
                
                # Create label with top LSA events
                cell_events = lsa_df[lsa_df['cell'] == node].sort_values('pvalue')
                top_events = cell_events.head(2)
                
                label_parts = [str(node)[:15]]  # Shortened cell name
                for _, event in top_events.iterrows():
                    region = str(event['region'])[:10]  # Shortened region
                    cna_type = event['CNA']
                    label_parts.append(f"{region}:{cna_type}")
                
                node_labels[node] = '\n'.join(label_parts)
                
            else:
                # Non-LSA cell - gray color, smaller size
                color = 'lightgray'
                size = 100
                # No label for non-LSA cells
            
            node_colors.append(color)
            node_sizes.append(size)
        
        return node_colors, node_sizes, node_labels
    
    def plot_lsa_tree(self, lsa_file: str, tree_file: str, 
                     output_filename: str = "segmental_LSA_tree.pdf") -> Optional[str]:
        """
        Create LSA tree visualization showing cells with significant events
        """
        print(f"Creating LSA tree visualization...")
        
        try:
            # Load data
            lsa_df = pd.read_csv(lsa_file, sep='\t')
            tree_df = pd.read_csv(tree_file, sep='\t', comment='#')
            
            print(f"  LSA events: {len(lsa_df)}")
            print(f"  LSA cells: {len(lsa_df['cell'].unique())}")
            
        except Exception as e:
            print(f"Error loading LSA data: {e}")
            return None
        
        if len(lsa_df) == 0:
            print("No LSA events to visualize")
            return None
        
        # Create LSA network
        lsa_network = self._create_lsa_network(lsa_df, tree_df)
        
        if len(lsa_network) == 0:
            print("No LSA network connections found")
            return None
        
        # Create NetworkX graph
        G = nx.from_pandas_edgelist(lsa_network, 
                                   source='from', 
                                   target='to', 
                                   edge_attr='weight', 
                                   create_using=nx.DiGraph())
        
        # Assign node properties
        node_colors, node_sizes, node_labels = self._assign_node_properties(G, lsa_df)
        
        # Create layout
        try:
            # Try hierarchical layout
            pos = nx.nx_agraph.graphviz_layout(G, prog='dot')
        except:
            try:
                # Fallback to spring layout with good parameters
                pos = nx.spring_layout(G, k=3, iterations=100, seed=42)
            except:
                # Final fallback
                pos = nx.circular_layout(G)
        
        # Create the plot
        plt.figure(figsize=(14, 10))
        
        # Draw edges
        edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
        nx.draw_networkx_edges(G, pos, 
                              width=[w*0.5 for w in edge_weights],
                              edge_color='gray',
                              arrows=True,
                              arrowsize=15,
                              arrowstyle='->',
                              alpha=0.6)
        
        # Draw nodes
        nx.draw_networkx_nodes(G, pos,
                              node_color=node_colors,
                              node_size=node_sizes,
                              alpha=0.8,
                              edgecolors='black',
                              linewidths=1)
        
        # Add labels for LSA cells only
        if node_labels:
            # Position labels slightly offset from nodes
            label_pos = {node: (x, y+0.1) for node, (x, y) in pos.items() if node in node_labels}
            
            nx.draw_networkx_labels(G, label_pos, node_labels,
                                   font_size=8,
                                   font_weight='bold',
                                   bbox=dict(boxstyle="round,pad=0.3", 
                                           facecolor="white", 
                                           alpha=0.8,
                                           edgecolor='black'))
        
        plt.title("Lineage Speciation Analysis (LSA) Tree\n"
                 f"Cells with Significant Copy Number Alterations", 
                 fontsize=14, fontweight='bold')
        
        # Create legend
        legend_elements = []
        
        # CNA type legend
        legend_elements.extend([
            mpatches.Patch(color='red', label='Amplification (AMP)'),
            mpatches.Patch(color='blue', label='Deletion (DEL)'),
            mpatches.Patch(color='lightgray', label='No significant LSA')
        ])
        
        # Size legend
        legend_elements.extend([
            plt.scatter([], [], s=100, c='black', label='1-2 LSA events'),
            plt.scatter([], [], s=200, c='black', label='3-4 LSA events'),
            plt.scatter([], [], s=400, c='black', label='5+ LSA events')
        ])
        
        plt.legend(handles=legend_elements, 
                  loc='upper left', 
                  bbox_to_anchor=(0, 1),
                  frameon=True,
                  fancybox=True,
                  shadow=True)
        
        plt.axis('off')
        
        # Save plot
        output_path = os.path.join(self.output_dir, output_filename)
        plt.savefig(output_path, format='pdf', bbox_inches='tight', dpi=300)
        plt.close()
        
        print(f"✓ LSA tree saved to {output_path}")
        return output_path
    
    def plot_lsa_summary(self, lsa_file: str, output_filename: str = "LSA_summary.pdf") -> Optional[str]:
        """
        Create summary plots for LSA results
        """
        print(f"Creating LSA summary plots...")
        
        try:
            lsa_df = pd.read_csv(lsa_file, sep='\t')
        except Exception as e:
            print(f"Error loading LSA data: {e}")
            return None
        
        if len(lsa_df) == 0:
            print("No LSA events to summarize")
            return None
        
        # Create subplot figure
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('Lineage Speciation Analysis (LSA) Summary', fontsize=16, fontweight='bold')
        
        # Plot 1: CNA type distribution
        cna_counts = lsa_df['CNA'].value_counts()
        colors = ['red' if cna == 'AMP' else 'blue' for cna in cna_counts.index]
        
        axes[0, 0].pie(cna_counts.values, labels=cna_counts.index, autopct='%1.1f%%', 
                      colors=colors, startangle=90)
        axes[0, 0].set_title('CNA Type Distribution', fontweight='bold')
        
        # Plot 2: Chromosome distribution
        lsa_df['chromosome'] = lsa_df['region'].str.extract(r'(chr\d+|chrX|chrY)')
        chr_counts = lsa_df['chromosome'].value_counts().head(10)
        
        bars = axes[0, 1].bar(range(len(chr_counts)), chr_counts.values, 
                             color='skyblue', edgecolor='black', alpha=0.7)
        axes[0, 1].set_xticks(range(len(chr_counts)))
        axes[0, 1].set_xticklabels(chr_counts.index, rotation=45)
        axes[0, 1].set_title('LSA Events by Chromosome', fontweight='bold')
        axes[0, 1].set_ylabel('Number of Events')
        
        # Add value labels on bars
        for bar, value in zip(bars, chr_counts.values):
            axes[0, 1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                           str(value), ha='center', va='bottom', fontweight='bold')
        
        # Plot 3: P-value distribution
        axes[0, 2].hist(lsa_df['pvalue'], bins=20, alpha=0.7, color='green', 
                       edgecolor='black')
        axes[0, 2].axvline(0.05, color='red', linestyle='--', linewidth=2, 
                          label='p=0.05 threshold')
        axes[0, 2].set_title('P-value Distribution', fontweight='bold')
        axes[0, 2].set_xlabel('P-value')
        axes[0, 2].set_ylabel('Frequency')
        axes[0, 2].legend()
        
        # Plot 4: Score distribution by CNA type
        for cna_type in lsa_df['CNA'].unique():
            data = lsa_df[lsa_df['CNA'] == cna_type]['Score']
            color = 'red' if cna_type == 'AMP' else 'blue'
            axes[1, 0].hist(data, alpha=0.6, label=cna_type, bins=15, 
                           color=color, edgecolor='black')
        
        axes[1, 0].set_title('Score Distribution by CNA Type', fontweight='bold')
        axes[1, 0].set_xlabel('LSA Score')
        axes[1, 0].set_ylabel('Frequency')
        axes[1, 0].legend()
        
        # Plot 5: Subtree size vs Score
        colors = ['red' if cna == 'AMP' else 'blue' for cna in lsa_df['CNA']]
        scatter = axes[1, 1].scatter(lsa_df['subtreesize'], lsa_df['Score'], 
                                    c=colors, alpha=0.6, s=50, edgecolors='black')
        axes[1, 1].set_title('Subtree Size vs LSA Score', fontweight='bold')
        axes[1, 1].set_xlabel('Subtree Size')
        axes[1, 1].set_ylabel('LSA Score')
        axes[1, 1].grid(True, alpha=0.3)
        
        # Add trend line
        z = np.polyfit(lsa_df['subtreesize'], lsa_df['Score'], 1)
        p = np.poly1d(z)
        axes[1, 1].plot(lsa_df['subtreesize'], p(lsa_df['subtreesize']), 
                       "r--", alpha=0.8, linewidth=2)
        
        # Plot 6: Top cells with LSA events
        cell_counts = lsa_df['cell'].value_counts().head(10)
        
        if len(cell_counts) > 0:
            y_pos = np.arange(len(cell_counts))
            bars = axes[1, 2].barh(y_pos, cell_counts.values, color='orange', 
                                  alpha=0.7, edgecolor='black')
            axes[1, 2].set_yticks(y_pos)
            
            # Shorten cell names for display
            short_names = [name[:20] + '...' if len(name) > 20 else name 
                          for name in cell_counts.index]
            axes[1, 2].set_yticklabels(short_names)
            
            axes[1, 2].set_title('Top 10 Cells with LSA Events', fontweight='bold')
            axes[1, 2].set_xlabel('Number of LSA Events')
            
            # Add value labels
            for bar, value in zip(bars, cell_counts.values):
                axes[1, 2].text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2,
                               str(value), ha='left', va='center', fontweight='bold')
        
        plt.tight_layout()
        
        # Save plot
        output_path = os.path.join(self.output_dir, output_filename)
        plt.savefig(output_path, format='pdf', bbox_inches='tight', dpi=300)
        plt.close()
        
        print(f"✓ LSA summary saved to {output_path}")
        return output_path
    
    def create_lsa_visualizations(self, lsa_file: str, tree_file: str) -> List[str]:
        """
        Create all LSA visualizations
        """
        output_files = []
        
        # Check if LSA file exists and has data
        try:
            lsa_df = pd.read_csv(lsa_file, sep='\t')
            if len(lsa_df) == 0:
                print("No LSA events found - skipping LSA visualizations")
                return []
        except:
            print("Could not load LSA file - skipping LSA visualizations")
            return []
        
        # Create LSA tree
        tree_plot = self.plot_lsa_tree(lsa_file, tree_file)
        if tree_plot:
            output_files.append(tree_plot)
        
        # Create LSA summary
        summary_plot = self.plot_lsa_summary(lsa_file)
        if summary_plot:
            output_files.append(summary_plot)
        
        return output_files

def main():
    """Main entry point"""
    
    if len(sys.argv) < 3:
        print("Usage: python3 lsa_visualization.py <lsa_file> <tree_file> [output_dir]")
        print()
        print("Example:")
        print("  python3 lsa_visualization.py segmental.LSA.txt tree.txt")
        print("  python3 lsa_visualization.py segmental.LSA.txt tree.txt ./plots/")
        sys.exit(1)
    
    lsa_file = sys.argv[1]
    tree_file = sys.argv[2]
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "./lsa_plots"
    
    # Check if files exist
    for file_path in [lsa_file, tree_file]:
        if not os.path.exists(file_path):
            print(f"Error: File '{file_path}' not found!")
            sys.exit(1)
    
    print("="*60)
    print("LSA VISUALIZATION MODULE")
    print("="*60)
    print(f"LSA file: {lsa_file}")
    print(f"Tree file: {tree_file}")
    print(f"Output directory: {output_dir}")
    print()
    
    # Create visualizer
    visualizer = LSAVisualizer(output_dir)
    
    # Generate visualizations
    output_files = visualizer.create_lsa_visualizations(lsa_file, tree_file)
    
    if output_files:
        print()
        print("="*60)
        print("LSA VISUALIZATION COMPLETED!")
        print("="*60)
        print("Generated files:")
        for file_path in output_files:
            print(f"  ✓ {file_path}")
        print(f"\n🎉 LSA visualizations saved in: {output_dir}")
    else:
        print("\n⚠️  No LSA visualizations generated")

if __name__ == "__main__":
    main()