import os
import copy
import numpy as np
from collections import deque
from datetime import datetime as dt
import networkx as nx


def read_cnv(in_seg_path):
    """
    Reads a segmental CNV file and returns:
    - node_dic: dict mapping cell names to lists of lists (copy number per chromosome segment)
    - root: cell name for the diploid (or near-diploid) root. If no perfectly diploid
      cell exists, choose the cell with the *smallest total deviation* from copy-number 2.
    """
    import pandas as pd
    df = pd.read_csv(in_seg_path, sep="\t", index_col=0)
    df.index = df.index.map(str)  # Ensure index is always string cell names
    # Ensure copy numbers are integers.
    df = df.round().astype(int)

    # Build node_dic: list of segment lists (each segment as length-1 list) to
    # preserve compatibility with downstream `dist` routine without relying on
    # a particular column-naming convention (e.g., "chr1_1" vs "chr1:p11").
    node_dic = {cell: [[int(cn)] for cn in df.loc[cell].tolist()] for cell in df.index}

    # Identify a diploid (or the most diploid-like) cell to serve as root
    deviations = (df - 2).abs().sum(axis=1)
    diploids = deviations[deviations == 0].index.tolist()
    if diploids:
        root = diploids[0]
    else:
        root = deviations.idxmin()
        print(f"No perfectly diploid cell found. Using {root} (least deviated) as root.")
    return node_dic, str(root)

def dist(node1, node2):
    d = 0
    for i in range(0, len(node1)):
        d = d + disthelper(node1[i], node2[i])
    return d

def disthelper(node1, node2):
    if 0 in node1 or 0 in node2:
        return zerodisthelper(node1, node2)
    return distcalc(node1, node2)

def distcalc(node1, node2):
    assert len(node1) == len(node2)
    if len(node1) == 1:
        return abs(node1[0] - node2[0])
    else:
        d = 0
        newlist = copy.deepcopy(node1)
        for i in range(0, len(node2)):
            newlist[i] -= node2[i]
        while newlist:
            if newlist[0] == 0:
                newlist.pop(0)
            elif newlist[0] > 0:
                k = 0
                for i in range(0, len(newlist)):
                    if newlist[i] > 0:
                        k = i
                    else:
                        break
                for i in range(0, k + 1):
                    newlist[i] -= 1
                d += 1
            elif newlist[0] < 0:
                k = 0
                for i in range(0, len(newlist)):
                    if newlist[i] < 0:
                        k = i
                    else:
                        break
                for i in range(0, k + 1):
                    newlist[i] += 1
                d += 1
        return abs(d)

def zerodisthelper(node1, node2):
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

def create_tree(nodes, node_list, root, proximity=True, len_threshold=30, df_cor=None):
    """
    Build a directed graph using NetworkX for all pairwise distances.
    Edges are inserted in a deterministic, lexicographically-sorted order so that
    NetworkX's minimum_spanning_arborescence produces reproducible results even
    when multiple edges share the same weight.
    Returns a NetworkX DiGraph.
    """
    G = nx.DiGraph()

    # Collect candidate edges first so we can sort them deterministically.
    edge_buffer = []  # (weight, src, dst)
    for A_node in node_list:
        for B_node in node_list:
            if A_node == B_node:
                continue
            if df_cor is not None:
                physical_dist = (
                    (df_cor.loc[A_node, "coor_x"] - df_cor.loc[B_node, "coor_x"]) ** 2
                    + (df_cor.loc[A_node, "coor_y"] - df_cor.loc[B_node, "coor_y"])
                )
                if physical_dist >= len_threshold ** 2:
                    continue
            w = dist(nodes[A_node], nodes[B_node])
            edge_buffer.append((w, str(A_node), str(B_node)))  # ensure str for lexicographic sort

    # Sort by weight first, then lexicographically by (src, dst).
    edge_buffer.sort(key=lambda x: (x[0], x[1], x[2]))

    # Insert edges in the sorted order.
    for w, src, dst in edge_buffer:
        G.add_edge(src, dst, weight=w)

    return G

def compute_rdmst(G, root):
    """Return a minimum spanning arborescence oriented away from the given root.

    This wraps NetworkX's minimum_spanning_arborescence but first guarantees that
    the provided `root` is the (unique) root of the arborescence by removing any
    incoming edges to `root`.  As a result, the returned DiGraph already has
    edges oriented outwards from `root`, so we no longer need the extra BFS
    re-orientation step used previously.
    """
    # Work on a copy so we don't mutate callers' graph.
    G2 = G.copy()

    # Remove all edges that point **into** the designated root so that the root
    # is forced to have in-degree 0 in the arborescence.
    in_edges_to_root = list(G2.in_edges(root, data=True))
    G2.remove_edges_from([(u, v) for u, v, *_ in in_edges_to_root])

    # Compute the minimum spanning arborescence.
    mst = nx.minimum_spanning_arborescence(G2)

    # Convert to nested-dict representation expected by downstream code.
    tree_dict = {u: {v: d["weight"] for v, d in mst[u].items()} for u in mst.nodes()}

    return tree_dict, mst

def graph_rank_dist(g, startnode):
    dist = {}
    for node in g:
        dist[node] = float('inf')
    dist[startnode] = 0
    queue = deque([startnode])
    while queue:
        node = queue.popleft()
        for nbr in g[node]:
            if dist[nbr] == float('inf'):
                dist[nbr] = dist[node] + 1
                queue.append(nbr)
    return dist

def reverse__graph(g):
    rev_graph = {}
    for node in g:
        rev_graph[node] = {}
    for a_node in g:
        for b_node in g[a_node]:
            rev_graph[b_node][a_node] = g[a_node][b_node]
    return rev_graph

def contract_graph(in_graph, root):
    rg = in_graph.copy()
    for node in rg:
        if not node == root:
            minimum = min(rg[node].values())
            for in_nbr in rg[node]:
                rg[node][in_nbr] -= minimum
    return rg

def get_rdst_graph(rg, root):
    # For each node, keep only zero-weight incoming edges
    g = {k: {} for k in rg}
    for node in rg:
        for nbr in rg[node]:
            if rg[node][nbr] == 0:
                g[node][nbr] = 0
    return g

def find_graphloop(rdst_candidate):
    # Find a cycle in the graph, if any
    visited = set()
    stack = []
    parent = {}
    for node in rdst_candidate:
        if node in visited:
            continue
        stack.append(node)
        while stack:
            curr = stack.pop()
            if curr in visited:
                continue
            visited.add(curr)
            for nbr in rdst_candidate[curr]:
                if nbr not in visited:
                    parent[nbr] = curr
                    stack.append(nbr)
                elif parent.get(curr, None) != nbr:
                    # Found a cycle
                    cycle = [nbr, curr]
                    p = parent.get(curr, None)
                    while p and p != nbr:
                        cycle.append(p)
                        p = parent.get(p, None)
                    cycle.reverse()
                    return cycle
    return None

def prune____graph(g, cycle):
    # Remove the cycle from the graph
    g2 = copy.deepcopy(g)
    for node in cycle:
        g2.pop(node, None)
    return g2, cycle

def expand___graph(g_contract, g_new_rdst, loop_node, loop_repr):
    # Re-expand the contracted cycle
    g_expanded = copy.deepcopy(g_contract)
    for node in loop_node:
        g_expanded[node] = {}
    for node in g_new_rdst:
        for nbr in g_new_rdst[node]:
            g_expanded[node][nbr] = g_new_rdst[node][nbr]
    return g_expanded 