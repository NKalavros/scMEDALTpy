import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
import seaborn as sns
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

class MEDALTVisualizer:
    """
    Python implementation of MEDALT visualization functions from R
    """
    
    def __init__(self, output_path: str = "."):
        self.output_path = output_path
    
    def _create_tree_layout(self, G, root):
        """
        Create a hierarchical tree layout manually
        """
        pos = {}
        
        # BFS to assign levels
        levels = {root: 0}
        queue = [root]
        
        while queue:
            node = queue.pop(0)
            current_level = levels[node]
            
            # Get children (successors in directed graph)
            children = list(G.successors(node))
            for child in children:
                if child not in levels:
                    levels[child] = current_level + 1
                    queue.append(child)
        
        # Group nodes by level
        level_groups = {}
        for node, level in levels.items():
            if level not in level_groups:
                level_groups[level] = []
            level_groups[level].append(node)
        
        # Assign positions
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
        
    def plot_single_cell_tree(self, tree_file: str, output_filename: str = "singlecell_tree.pdf"):
        """
        Plot cell tree figure (equivalent to R lines 42-51)
        """
        # Read tree data
        celltree = pd.read_csv(tree_file, sep='\t')
        
        # Create nodes dataframe
        all_nodes = set(celltree.iloc[:, 0].astype(str)) | set(celltree.iloc[:, 1].astype(str))
        nodes = pd.DataFrame({'id': list(all_nodes), 'size': [5] * len(all_nodes)})
        
        # Set node colors - root node is black, others are lightblue
        nodes['color'] = 'lightblue'
        root_nodes = set(celltree.iloc[:, 0].astype(str)) - set(celltree.iloc[:, 1].astype(str))
        nodes.loc[nodes['id'].isin(root_nodes), 'color'] = 'black'
        
        # Create networkx graph
        G = nx.from_pandas_edgelist(celltree, source=celltree.columns[0], 
                                   target=celltree.columns[1], create_using=nx.DiGraph())
        
        # Plot
        plt.figure(figsize=(15, 10))
        
        # Use hierarchical tree layout
        try:
            # Try graphviz hierarchical layout first
            pos = nx.nx_agraph.graphviz_layout(G, prog='dot')
        except:
            try:
                # Fallback to manual hierarchical layout
                pos = self._create_tree_layout(G, 'root' if 'root' in G.nodes() else list(G.nodes())[0])
            except:
                # Final fallback to spring layout but with better parameters
                pos = nx.spring_layout(G, k=3, iterations=100)
        
        # Draw nodes
        node_colors = [nodes[nodes['id'] == node]['color'].iloc[0] for node in G.nodes()]
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=300, alpha=0.8)
        
        # Draw edges
        nx.draw_networkx_edges(G, pos, edge_color='gray', arrows=True, 
                              arrowsize=10, arrowstyle='->', alpha=0.6)
        
        # Add labels for important nodes (root and nodes with many connections)
        important_nodes = {}
        degree_centrality = nx.degree_centrality(G)
        sorted_nodes = sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)
        
        # Label root and top connected nodes
        for node, centrality in sorted_nodes[:min(10, len(sorted_nodes))]:
            if node == 'root' or centrality > 0.1:
                # Shorten cell names for display
                if len(str(node)) > 15:
                    label = str(node)[:12] + '...'
                else:
                    label = str(node)
                important_nodes[node] = label
        
        if important_nodes:
            nx.draw_networkx_labels(G, pos, important_nodes, font_size=6, font_weight='bold')
        
        plt.title("Single Cell Phylogenetic Tree")
        plt.axis('off')
        
        # Save to PDF
        output_file = f"{self.output_path}/{output_filename}"
        plt.savefig(output_file, format='pdf', bbox_inches='tight', dpi=300)
        plt.close()
        
        print(f"Single cell tree saved to {output_file}")
        return output_file
    
    def plot_lsa_tree(self, lsa_file: str, tree_file: str, output_filename: str = "LSA_tree.pdf"):
        """
        Plot LSA Tree with genomic annotations (equivalent to R lines 287-319)
        Only shows cells with significant CNAs and their connections
        """
        # Read LSA results
        try:
            allsig = pd.read_csv(lsa_file, sep='\t')
        except FileNotFoundError:
            print(f"LSA file {lsa_file} not found. Skipping LSA tree visualization.")
            return None
            
        if allsig.empty:
            print("No significant LSA results found. Skipping LSA tree visualization.")
            return None
            
        # Read tree data
        celltree = pd.read_csv(tree_file, sep='\t')
        
        # Create LSA network connections
        lsa_network = self._create_lsa_network(allsig, celltree)
        
        if lsa_network.empty:
            print("No LSA network connections found. Skipping visualization.")
            return None
        
        # Create nodes dataframe
        all_nodes = set(lsa_network['from']) | set(lsa_network['to'])
        nodes = pd.DataFrame({'id': list(all_nodes), 'size': [5] * len(all_nodes)})
        
        # Calculate node sizes based on number of CNAs
        cna_counts = allsig['cell'].value_counts()
        for i, node_id in enumerate(nodes['id']):
            if node_id in cna_counts:
                nodes.loc[i, 'size'] = max(5, cna_counts[node_id] * 5)
        
        # Normalize node sizes
        max_size = nodes['size'].max()
        nodes['size'] = (nodes['size'] / max_size) * 1000  # Scale for visualization
        
        # Set node colors
        nodes['color'] = 'gray'
        unique_cells = allsig['cell'].unique()
        colors = plt.cm.Set3(np.linspace(0, 1, len(unique_cells)))
        
        for i, cell in enumerate(unique_cells):
            mask = nodes['id'] == cell
            if mask.any():
                # Convert color to string format
                color_val = colors[i]
                nodes.loc[mask, 'color'] = f'#{int(color_val[0]*255):02x}{int(color_val[1]*255):02x}{int(color_val[2]*255):02x}'
        
        # Create annotations for nodes
        annotations = {}
        for node_id in nodes['id']:
            if node_id in allsig['cell'].values:
                cell_cnas = allsig[allsig['cell'] == node_id].copy()
                cell_cnas = cell_cnas.sort_values('pvalue')
                
                # Create annotation string (up to 3 CNAs)
                cna_strs = []
                for _, row in cell_cnas.head(3).iterrows():
                    cna_strs.append(f"{row['region']}:{row['CNA']}")
                
                annotations[node_id] = ';'.join(cna_strs)
            else:
                annotations[node_id] = ''
        
        # Create networkx graph
        G = nx.from_pandas_edgelist(lsa_network, source='from', target='to', 
                                   edge_attr='weight', create_using=nx.DiGraph())
        
        # Plot
        plt.figure(figsize=(12, 12))
        
        # Use hierarchical layout
        try:
            # Try graphviz hierarchical layout first
            pos = nx.nx_agraph.graphviz_layout(G, prog='dot')
        except:
            try:
                # Fallback to manual hierarchical layout
                pos = self._create_tree_layout(G, 'root' if 'root' in G.nodes() else list(G.nodes())[0])
            except:
                # Final fallback to spring layout but with better parameters
                pos = nx.spring_layout(G, k=3, iterations=100)
        
        # Draw nodes
        node_colors = [nodes[nodes['id'] == node]['color'].iloc[0] for node in G.nodes()]
        node_sizes = [nodes[nodes['id'] == node]['size'].iloc[0] for node in G.nodes()]
        
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.8)
        
        # Draw edges
        edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
        nx.draw_networkx_edges(G, pos, edge_color='gray', arrows=True, 
                              arrowsize=15, arrowstyle='->', alpha=0.6, width=edge_weights)
        
        # Add labels - only for nodes with annotations or root
        labels = {}
        for node in G.nodes():
            if node == 'root':
                labels[node] = 'ROOT'
            elif annotations.get(node, ''):
                # Show cell name (shortened) and top CNA
                cell_short = str(node)[:10] + '...' if len(str(node)) > 13 else str(node)
                top_cna = annotations[node].split(';')[0]  # Just the first CNA
                labels[node] = f"{cell_short}\n{top_cna}"
        
        if labels:
            nx.draw_networkx_labels(G, pos, labels, font_size=6, font_weight='bold')
        
        plt.title("LSA Tree with Genomic Annotations")
        plt.axis('off')
        
        # Save to PDF
        output_file = f"{self.output_path}/{output_filename}"
        plt.savefig(output_file, format='pdf', bbox_inches='tight', dpi=300)
        plt.close()
        
        print(f"LSA tree saved to {output_file}")
        return output_file
    
    def _create_lsa_network(self, allsig: pd.DataFrame, celltree: pd.DataFrame) -> pd.DataFrame:
        """
        Create LSA network connections (equivalent to R CNAconnect function)
        """
        # This is a simplified version - the actual R function is complex
        # We'll create connections based on the cell tree structure
        
        significant_cells = set(allsig['cell'].unique())
        
        # Filter tree to only include significant cells
        lsa_connections = []
        
        for _, row in celltree.iterrows():
            from_cell = str(row.iloc[0])
            to_cell = str(row.iloc[1])
            distance = row.iloc[2] if len(row) > 2 else 1
            
            # Include connection if either cell has significant LSA
            if from_cell in significant_cells or to_cell in significant_cells:
                lsa_connections.append({
                    'from': from_cell,
                    'to': to_cell,
                    'weight': max(1, distance)
                })
        
        return pd.DataFrame(lsa_connections)
    
    def create_comprehensive_report(self, tree_file: str, lsa_file: str = None, 
                                  output_filename: str = "MEDALT_visualization_report.pdf"):
        """
        Create a comprehensive visualization report with multiple plots
        """
        with PdfPages(f"{self.output_path}/{output_filename}") as pdf:
            # Plot 1: Single cell tree
            self.plot_single_cell_tree(tree_file, "temp_singlecell.pdf")
            
            # Plot 2: LSA tree (if available)
            if lsa_file and pd.read_csv(lsa_file, sep='\t').shape[0] > 0:
                self.plot_lsa_tree(lsa_file, tree_file, "temp_lsa.pdf")
                
                # Plot 3: LSA summary statistics
                self._plot_lsa_summary(lsa_file, pdf)
            
            # Plot 4: Tree statistics
            self._plot_tree_statistics(tree_file, pdf)
        
        print(f"Comprehensive report saved to {self.output_path}/{output_filename}")
        return f"{self.output_path}/{output_filename}"
    
    def _plot_lsa_summary(self, lsa_file: str, pdf: PdfPages):
        """
        Create summary plots for LSA results
        """
        allsig = pd.read_csv(lsa_file, sep='\t')
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Plot 1: CNA type distribution
        cna_counts = allsig['CNA'].value_counts()
        axes[0, 0].pie(cna_counts.values, labels=cna_counts.index, autopct='%1.1f%%')
        axes[0, 0].set_title('CNA Type Distribution')
        
        # Plot 2: Chromosome distribution
        allsig['chromosome'] = allsig['region'].str.extract(r'(chr\d+|chrX|chrY)')
        chr_counts = allsig['chromosome'].value_counts().head(10)
        axes[0, 1].bar(chr_counts.index, chr_counts.values)
        axes[0, 1].set_title('Top 10 Chromosomes with LSA Events')
        axes[0, 1].tick_params(axis='x', rotation=45)
        
        # Plot 3: P-value distribution
        axes[1, 0].hist(allsig['pvalue'], bins=20, alpha=0.7, edgecolor='black')
        axes[1, 0].set_title('P-value Distribution')
        axes[1, 0].set_xlabel('P-value')
        axes[1, 0].set_ylabel('Frequency')
        
        # Plot 4: Score distribution by CNA type
        for cna_type in allsig['CNA'].unique():
            data = allsig[allsig['CNA'] == cna_type]['Score']
            axes[1, 1].hist(data, alpha=0.5, label=cna_type, bins=15)
        axes[1, 1].set_title('Score Distribution by CNA Type')
        axes[1, 1].set_xlabel('Score')
        axes[1, 1].set_ylabel('Frequency')
        axes[1, 1].legend()
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _plot_tree_statistics(self, tree_file: str, pdf: PdfPages):
        """
        Create tree statistics plots
        """
        celltree = pd.read_csv(tree_file, sep='\t')
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Plot 1: Distance distribution
        if 'dist' in celltree.columns:
            distances = celltree['dist']
            axes[0, 0].hist(distances, bins=20, alpha=0.7, edgecolor='black')
            axes[0, 0].set_title('Edge Distance Distribution')
            axes[0, 0].set_xlabel('Distance')
            axes[0, 0].set_ylabel('Frequency')
        
        # Plot 2: Node degree distribution
        all_nodes = list(celltree.iloc[:, 0]) + list(celltree.iloc[:, 1])
        from collections import Counter
        degree_counts = Counter(all_nodes)
        degrees = list(degree_counts.values())
        
        axes[0, 1].hist(degrees, bins=max(1, len(set(degrees))), alpha=0.7, edgecolor='black')
        axes[0, 1].set_title('Node Degree Distribution')
        axes[0, 1].set_xlabel('Degree')
        axes[0, 1].set_ylabel('Frequency')
        
        # Plot 3: Tree structure summary
        num_nodes = len(set(all_nodes))
        num_edges = len(celltree)
        
        axes[1, 0].text(0.1, 0.7, f'Total Nodes: {num_nodes}', fontsize=14, transform=axes[1, 0].transAxes)
        axes[1, 0].text(0.1, 0.5, f'Total Edges: {num_edges}', fontsize=14, transform=axes[1, 0].transAxes)
        axes[1, 0].text(0.1, 0.3, f'Average Degree: {2*num_edges/num_nodes:.2f}', fontsize=14, transform=axes[1, 0].transAxes)
        axes[1, 0].set_title('Tree Statistics')
        axes[1, 0].axis('off')
        
        # Plot 4: Top connected nodes
        top_nodes = Counter(all_nodes).most_common(10)
        if top_nodes:
            nodes, counts = zip(*top_nodes)
            axes[1, 1].barh(range(len(nodes)), counts)
            axes[1, 1].set_yticks(range(len(nodes)))
            axes[1, 1].set_yticklabels([str(n)[:20] + '...' if len(str(n)) > 20 else str(n) for n in nodes])
            axes[1, 1].set_title('Top 10 Connected Nodes')
            axes[1, 1].set_xlabel('Connection Count')
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

def create_medalt_visualizations(tree_file: str, lsa_file: str = None, output_path: str = "."):
    """
    Main function to create MEDALT visualizations
    """
    visualizer = MEDALTVisualizer(output_path)
    
    # Create individual plots
    visualizer.plot_single_cell_tree(tree_file)
    
    if lsa_file:
        visualizer.plot_lsa_tree(lsa_file, tree_file)
    
    # Create comprehensive report
    visualizer.create_comprehensive_report(tree_file, lsa_file)
    
    print("MEDALT visualizations completed!")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python visualization.py <tree_file> [lsa_file] [output_path]")
        sys.exit(1)
    
    tree_file = sys.argv[1]
    lsa_file = sys.argv[2] if len(sys.argv) > 2 else None
    output_path = sys.argv[3] if len(sys.argv) > 3 else "."
    
    create_medalt_visualizations(tree_file, lsa_file, output_path)