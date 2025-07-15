import numpy as np
from numba import jit
from collections import *
from copy import *
import sys

INF = float('inf')

def reverse_g(g):
    rev_graph = {}
    for node in g:
        rev_graph[node] = {}
    for a_node in g:
        for b_node in g[a_node]:
            rev_graph[b_node][a_node] = g[a_node][b_node]
    return rev_graph

def update_dege(rg, root):
    for node, nbrs in rg.items():
        if node != root and nbrs:
            # find minimum incoming weight in C
            min_val = min(nbrs.values())
            # subtract it from every incoming edge
            for in_nbr in nbrs:
                nbrs[in_nbr] -= min_val


def compute_rdst_candidate(rg, root):
    candidate = {node: {} for node in rg}
    for node, nbrs in rg.items():
        if node != root and nbrs:
            # get the minimum incoming weight
            min_val = min(nbrs.values())
            # pick the first edge that attains it
            for in_nbr, w in nbrs.items():
                if w == min_val:
                    candidate[node][in_nbr] = w
                    break
    return candidate


def get_cycle(rdst_candidate):
    node_unvisited = []
    for node in rdst_candidate:  
        node_unvisited.append(node)
    while not node_unvisited == []:
        start_node = node_unvisited.pop()  
        stack = []  
        trail = []  
        stack.append(start_node)
        while not len(stack) == 0:
            node = stack.pop(-1)
            for nbr in rdst_candidate[node]:
                if nbr in trail:
                    return tuple(trail[trail.index(nbr):]) 
                else:
                    stack.append(nbr)
                    trail.append(nbr)
                    if nbr in node_unvisited:
                        node_unvisited.remove(nbr)
    return False


def contract_cycle(g, cycle):
    cstar = max(g.keys()) + "1"
    contracted_graph = {}
    contracted_graph[cstar] = {}
    for node in g:
        if not node in cycle:
            contracted_graph[node] = {}
    for node in g:
        for nbr in g[node]:
            if node in cycle:
                if nbr in cycle:
                    pass
                else: 
                    if nbr in contracted_graph[cstar]:
                        contracted_graph[cstar][nbr] = min(contracted_graph[cstar][nbr], g[node][nbr])
                    else:
                        contracted_graph[cstar][nbr] = g[node][nbr]
            else:
                if nbr in cycle:
                    if cstar in contracted_graph[node]:
                        contracted_graph[node][cstar] = min(contracted_graph[node][cstar], g[node][nbr])
                    else:
                        contracted_graph[node][cstar] = g[node][nbr]
                else:
                    contracted_graph[node][nbr] = g[node][nbr]
    return contracted_graph, cstar


def expand_graph(g, rdst_candidate, cycle, cstar):
    restored_graph = {}
    for node in g:
        restored_graph[node] = {}
    for node in rdst_candidate:  
        for nbr in rdst_candidate[node]:
            if node == cstar:  
                min = float('inf')
                for orig in cycle:
                    if nbr in g[orig]:
                        if g[orig][nbr] < min:
                            min = g[orig][nbr]
                            point = orig
                restored_graph[point][nbr] = min
            else:
                if nbr == cstar:  
                    min = float('inf')
                    for orig_nbr in g[node]:
                        if orig_nbr in cycle:
                            if g[node][orig_nbr] < min:
                                min = g[node][orig_nbr]
                                start_pt = orig_nbr  
                    restored_graph[node][start_pt] = min
                else: 
                    restored_graph[node][nbr] = g[node][nbr]
    for index in range(len(cycle) - 1): 
        restored_graph[cycle[index + 1]][cycle[index]] = g[cycle[index + 1]][cycle[index]]
    restored_graph[cycle[0]][cycle[-1]] = g[cycle[0]][cycle[-1]]
    index = cycle.index(start_pt)
    if index == len(cycle) - 1:
        restored_graph[cycle[0]].pop(cycle[index])
    else:
        restored_graph[cycle[index + 1]].pop(cycle[index])
    return restored_graph


def bfs(g, startnode):
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


def compute_rdmst(g, root):
    """Compute the rooted directed minimum spanning tree (RDMST) using Edmonds' algorithm"""
    from datetime import datetime as dt_
    import os
    try:
        import psutil
        def pc_mem():
            process = psutil.Process(os.getpid())
            mem_info = process.memory_info()
            return round(mem_info.rss/1024/1024/1024, 2)
    except ImportError:
        def pc_mem():
            return 0.0
    
    if root not in g:
        print "The root node does not exist"
        return

    distances = bfs(g, root)
    for node in g:
        if distances[node] == float('inf'):
            print "The root does not reach every other node in the graph"
            return

    t0 = dt_.now()
    print "Inferring single cell evolution tree..."
    rgraph = reverse_g(g)
    print 'Reversed graph created'
    rdmst = compute_rdmst_helper(g, rgraph, root)

    # Calculate weight
    rdmst_weight = 0
    for node in rdmst:
        for nbr in rdmst[node]:
            rdmst[node][nbr] = g[node][nbr]
            rdmst_weight += rdmst[node][nbr]
    
    elapsed = dt_.now() - t0
    print "RDMST computed in %.2f seconds" % elapsed.total_seconds()
    
    return (rdmst, rdmst_weight)


def compute_rdmst_helper(g, rgraph, root, depth=0):
    # progress tracking
    sys.stderr.write("compute_rdmst_helper depth %d: processing %d nodes\n" % (depth, len(g)))
    # subtract min incoming edge on rgraph
    update_dege(rgraph, root)
    # pick one candidate incoming edge per node
    rdst_candidate = compute_rdst_candidate(rgraph, root)
    cycle = get_cycle(rdst_candidate)

    if not cycle:
        # simple case: reverse the candidate graph back
        return reverse_g(rdst_candidate)
    else:
        # prepare contracted graph and its reverse
        g_copy = deepcopy(rgraph)
        g_copy = reverse_g(g_copy)
        (contracted_g, cstar) = contract_cycle(g_copy, cycle)
        contracted_rgraph = reverse_g(contracted_g)

        # recurse with increased depth
        new_rdmst = compute_rdmst_helper(contracted_g, contracted_rgraph, root, depth+1)
        # expand into original rgraph context
        return expand_graph(rgraph, new_rdmst, cycle, cstar)


@jit(nopython=True)
def rdmst_edmonds_optimized(adj_matrix, root_idx):
    n = adj_matrix.shape[0]
    min_edges = np.full(n, np.inf, dtype=np.float32)
    parent    = np.full(n, -1,   dtype=np.int32)
    for i in range(n):
        if i != root_idx:
            for j in range(n):
                w = adj_matrix[j, i]
                if w < min_edges[i]:
                    min_edges[i] = w
                    parent[i]    = j
    visited  = np.zeros(n, dtype=np.bool_)
    in_cycle = np.zeros(n, dtype=np.bool_)
    cycles   = []
    for i in range(n):
        if i != root_idx and not visited[i]:
            cycle = find_cycle_from_node(parent, i, visited, in_cycle)
            if cycle is not None:
                cycles.append(cycle)
    if len(cycles) == 0:
        result = np.full((n, n), np.inf, dtype=np.float32)
        for i in range(n):
            pi = parent[i]
            if pi >= 0:
                result[pi, i] = min_edges[i]
        return result
    else:
        return contract_and_solve_optimized(adj_matrix, root_idx, cycles, min_edges, parent)

@jit(nopython=True)
def find_cycle_from_node(parent, start, visited, in_cycle):
    path    = []
    current = start
    while current != -1 and not visited[current]:
        if current in path:
            idx   = path.index(current)
            cyc   = path[idx:]
            for node in cyc:
                in_cycle[node] = True
            return cyc
        path.append(current)
        visited[current] = True
        current = parent[current]
    return None

@jit(nopython=True)
def contract_and_solve_optimized(adj_matrix, root_idx, cycles, min_edges, parent):
    # only handle single cycle in Numba; fallback for >1
    if len(cycles) == 1:
        return contract_single_cycle_optimized(adj_matrix, root_idx, cycles[0], min_edges, parent)
    else:
        return rdmst_edmonds_fallback(adj_matrix, root_idx)

@jit(nopython=True)
def contract_single_cycle_optimized(adj_matrix, root_idx, cycle, min_edges, parent):
    n           = adj_matrix.shape[0]
    in_cycle    = set(cycle)
    non_cycle   = [i for i in range(n) if i not in in_cycle]
    new_size    = len(non_cycle) + 1
    old2new     = {uc:i for i,uc in enumerate(non_cycle)}
    super_node  = new_size - 1
    contracted  = np.full((new_size, new_size), np.inf, dtype=np.float32)
    # copy non-cycle edges
    for i in range(len(non_cycle)):
        oi = non_cycle[i]
        for j in range(len(non_cycle)):
            oj = non_cycle[j]
            contracted[i,j] = adj_matrix[oi,oj]
    # edges between cycle & non-cycle
    for i in range(len(non_cycle)):
        oi = non_cycle[i]
        # to cycle
        m2c = np.inf
        for c in cycle:
            w = adj_matrix[oi,c]
            if w < m2c: m2c = w
        if m2c < np.inf: contracted[i,super_node] = m2c
        # from cycle
        mfc = np.inf
        for c in cycle:
            w = adj_matrix[c,oi]
            if w < mfc: mfc = w
        if mfc < np.inf: contracted[super_node,i] = mfc
    new_root = old2new[root_idx] if root_idx in old2new else super_node
    sol      = rdmst_edmonds_optimized(contracted, new_root)
    # expand back
    result = np.full((n,n), np.inf, dtype=np.float32)
    for i in range(len(non_cycle)):
        for j in range(len(non_cycle)):
            if sol[i,j] < np.inf:
                old_i = new_to_old[i]
                old_j = new_to_old[j]
                result[old_i, old_j] = contracted_result[i, j]
    
    # Handle cycle edges
    for i in range(len(cycle) - 1):
        result[cycle[i+1], cycle[i]] = adj_matrix[cycle[i+1], cycle[i]]
    result[cycle[0], cycle[-1]] = adj_matrix[cycle[0], cycle[-1]]
    
    return result


@jit(nopython=True)
def rdmst_edmonds_fallback(adj_matrix, root_idx):
    n       = adj_matrix.shape[0]
    result  = np.full((n,n), np.inf, dtype=np.float32)
    visited = np.zeros(n, dtype=np.bool_)
    visited[root_idx] = True
    for _ in range(n-1):
        me = np.inf; mi = -1; mj = -1
        for i in range(n):
            if visited[i]:
                for j in range(n):
                    w = adj_matrix[i,j]
                    if not visited[j] and w < me:
                        me, mi, mj = w, i, j
        if mi >= 0:
            result[mi,mj] = me
            visited[mj]  = True
    return result
