import os
import copy
import numpy as np
import pandas as pd
from collections import deque
from datetime import datetime as dt
import networkx as nx

def read_cnv(in_seg_path):
    """
    Reads a segmental CNV file and returns:
    - node_dic: dict mapping cell names to lists of lists (copy number per chromosome segment)
    - root: cell name for the diploid (or near-diploid) root
    """
    df = pd.read_csv(in_seg_path, sep="\t", index_col=0)
    df.index = df.index.map(str)  # Ensure index is always string cell names
    
    # Ensure copy numbers are integers
    df = df.round().astype(int)
    
    # Organize by chromosome
    chr_scan = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chr_data = {}
    
    for col in df.columns:
        # Extract chromosome from column name (chr7:p22.3 -> chr7 or chr7_1 -> chr7)
        if ':' in col:  # Cytoband format
            chr_match = col.split(':')[0]
        else:  # Simple format
            chr_match = col.split('_')[0]
        
        if chr_match in chr_scan:
            if chr_match not in chr_data:
                chr_data[chr_match] = []
            chr_data[chr_match].append(col)
    
    # Build node dictionary
    cell_lst = list(df.index)
    node_dic = {cell: [] for cell in cell_lst}
    
    for chr_i in sorted(chr_data.keys(), key=lambda x: (len(x), x)):
        chr_cols = sorted(chr_data[chr_i])  # Sort segments within chromosome
        for cell in cell_lst:
            node_dic[cell].append(list(df.loc[cell, chr_cols]))
    
    # Find root (most diploid cell)
    deviations = (df - 2).abs().sum(axis=1)
    diploids = deviations[deviations == 0].index.tolist()
    
    if diploids:
        root = diploids[0]
        print(f"Found diploid root: {root}")
    else:
        # Add artificial diploid root
        root = "root"
        diploid_profile = []
        for chr_i in sorted(chr_data.keys(), key=lambda x: (len(x), x)):
            diploid_profile.append([2] * len(chr_data[chr_i]))
        node_dic[root] = diploid_profile
        print(f"No diploid cell found. Added artificial root.")
    
    return node_dic, str(root)

def dist(node1, node2):
    """Calculate total distance between two nodes."""
    d = 0
    for i in range(len(node1)):
        d += disthelper(node1[i], node2[i])
    return d

def disthelper(node1, node2):
    """Helper function to handle chromosome-level distance."""
    if 0 in node1 or 0 in node2:
        return zerodisthelper(node1, node2)
    return distcalc(node1, node2)

def distcalc(node1, node2):
    """Calculate distance between two chromosome segments."""
    assert len(node1) == len(node2)
    
    if len(node1) == 1:
        return abs(node1[0] - node2[0])
    
    d = 0
    newlist = copy.deepcopy(node1)
    
    for i in range(len(node2)):
        newlist[i] -= node2[i]
    
    while newlist:
        if newlist[0] == 0:
            newlist.pop(0)
        elif newlist[0] > 0:
            k = 0
            for i in range(len(newlist)):
                if newlist[i] > 0:
                    k = i
                else:
                    break
            for i in range(k + 1):
                newlist[i] -= 1
            d += 1
        elif newlist[0] < 0:
            k = 0
            for i in range(len(newlist)):
                if newlist[i] < 0:
                    k = i
                else:
                    break
            for i in range(k + 1):
                newlist[i] += 1
            d += 1
    
    return abs(d)

def zerodisthelper(node1, node2):
    """Handle distance calculation when zeros are present."""
    n1 = copy.deepcopy(node1)
    n2 = copy.deepcopy(node2)
    temp1 = []
    temp2 = []
    
    while n1:
        x1 = n1.pop()
        x2 = n2.pop()
        if x1 == 0:
            if x2 == 0:
                temp1.append(x1)
                temp2.append(x2)
            else:
                return 1000000
        else:
            temp1.append(x1)
            temp2.append(x2)
    
    return distcalc(temp1, temp2)

def create_tree(nodes, node_list, root, proximity=True, len_threshold=30, df_cor=None, debug=False):
    """
    Build a directed graph using NetworkX for all pairwise distances.
    """
    G = nx.DiGraph()
    
    # Calculate all pairwise distances
    edge_buffer = []
    
    if debug:
        print(f"DEBUG: Computing distances for {len(node_list)} nodes...")
    
    dist_matrix = {}
    for i, A_node in enumerate(node_list):
        if debug and i % 10 == 0:
            print(f"DEBUG: Processing node {i}/{len(node_list)}")
        
        for B_node in node_list:
            if A_node == B_node:
                continue
            
            # Check spatial constraints if coordinates provided
            if df_cor is not None:
                physical_dist = (
                    (df_cor.loc[A_node, "coor_x"] - df_cor.loc[B_node, "coor_x"]) ** 2 +
                    (df_cor.loc[A_node, "coor_y"] - df_cor.loc[B_node, "coor_y"]) ** 2
                )
                if physical_dist >= len_threshold ** 2:
                    continue
            
            # Calculate genetic distance
            w = dist(nodes[A_node], nodes[B_node])
            dist_matrix[(A_node, B_node)] = w
            edge_buffer.append((w, str(A_node), str(B_node)))
    
    if debug:
        weights = [e[0] for e in edge_buffer]
        print(f"DEBUG: Distance statistics - min: {min(weights)}, max: {max(weights)}, "
              f"mean: {np.mean(weights):.2f}")
        print(f"DEBUG: Zero-weight edges: {sum(1 for w in weights if w == 0)}")
    
    # Sort edges deterministically
    edge_buffer.sort(key=lambda x: (x[0], x[1], x[2]))
    
    # Add edges to graph
    for w, src, dst in edge_buffer:
        G.add_edge(src, dst, weight=w)
    
    if debug:
        print(f"DEBUG: Created graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
    
    return G

def compute_rdmst(G, root, debug=False):
    """
    Compute minimum spanning arborescence from the given root.
    """
    G2 = G.copy()
    
    # Ensure root has no incoming edges
    in_edges = list(G2.in_edges(root))
    if debug:
        print(f"DEBUG: Removing {len(in_edges)} incoming edges to root {root}")
    G2.remove_edges_from(in_edges)
    
    # Ensure all nodes are reachable from root
    if root not in G2:
        G2.add_node(root)
    
    # Add edges from root to unreachable nodes with high weight
    reachable = nx.ancestors(G2, root) | {root}
    unreachable = set(G2.nodes()) - reachable
    
    if debug and unreachable:
        print(f"DEBUG: {len(unreachable)} nodes unreachable from root, adding high-weight edges")
    
    for node in unreachable:
        if node != root:
            # Add high-weight edge from root
            G2.add_edge(root, node, weight=1000000)
    
    # Compute minimum spanning arborescence
    try:
        mst = nx.minimum_spanning_arborescence(G2)
    except nx.NetworkXError as e:
        if debug:
            print(f"DEBUG: Error computing MST: {e}")
        # Fallback: create star topology from root
        mst = nx.DiGraph()
        mst.add_node(root)
        for node in G2.nodes():
            if node != root:
                weight = G2[root][node]['weight'] if G2.has_edge(root, node) else 1000000
                mst.add_edge(root, node, weight=weight)
    
    if debug:
        print(f"DEBUG: MST has {mst.number_of_nodes()} nodes and {mst.number_of_edges()} edges")
        weights = [d['weight'] for _, _, d in mst.edges(data=True)]
        if weights:
            print(f"DEBUG: MST edge weights - min: {min(weights)}, max: {max(weights)}, "
                  f"mean: {np.mean(weights):.2f}")
    
    # Convert to nested dict format
    tree_dict = {}
    for node in mst.nodes():
        tree_dict[node] = {}
        for neighbor in mst.successors(node):
            tree_dict[node][neighbor] = mst[node][neighbor]['weight']
    
    return tree_dict, mst

def graph_rank_dist(g, startnode):
    """Calculate graph distance from start node."""
    dist = {}
    for node in g:
        dist[node] = float('inf')
    dist[startnode] = 0
    queue = deque([startnode])
    
    while queue:
        node = queue.popleft()
        for nbr in g.get(node, {}):
            if dist[nbr] == float('inf'):
                dist[nbr] = dist[node] + 1
                queue.append(nbr)
    
    return dist

def reverse__graph(g):
    """Reverse all edges in the graph."""
    rev_graph = {}
    for node in g:
        rev_graph[node] = {}
    
    for a_node in g:
        for b_node in g[a_node]:
            rev_graph[b_node][a_node] = g[a_node][b_node]
    
    return rev_graph