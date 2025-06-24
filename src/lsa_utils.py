import numpy as np
import pandas as pd
from collections import defaultdict, deque
from statsmodels.stats.multitest import multipletests

def get_subtree(tree, root):
    """
    Given a tree (dict of dicts), return all nodes in the subtree rooted at 'root'.
    """
    nodes = set()
    queue = deque([root])
    while queue:
        node = queue.popleft()
        nodes.add(node)
        children = tree.get(node, {})
        for child in children:
            if child not in nodes:
                queue.append(child)
    return nodes

def lineage_cfl(cnv_df, lineage_cells):
    """
    Calculate cumulative fold level (CFL) for a set of cells (lineage) for each feature (bin/gene).
    Returns a pandas Series (index: feature, value: mean CNV across lineage).
    """
    return cnv_df.loc[list(lineage_cells)].mean(axis=0)

def lineage_score(cnv_df, lineage_cells, amp_cutoff=0.5, del_cutoff=-0.5):
    """
    Calculate AMP/DEL scores for a lineage (mean CNV for each feature, split by cutoff).
    Returns a DataFrame with AMP and DEL scores for each feature.
    """
    cfl = cnv_df.loc[lineage_cells].mean(axis=0)
    amp = cfl.where(cfl >= amp_cutoff, 0)
    dele = cfl.where(cfl <= del_cutoff, 0)
    return pd.DataFrame({'AMP': amp, 'DEL': dele})

def get_depths(tree, root):
    """
    Compute depth of each node in the tree (distance from root).
    Returns a dict: node -> depth
    """
    depths = {root: 0}
    queue = deque([(root, 0)])
    while queue:
        node, d = queue.popleft()
        for child in tree.get(node, {}):
            if child not in depths:
                depths[child] = d + 1
                queue.append((child, d + 1))
    return depths

def get_children_count(tree):
    """
    Compute number of children for each node in the tree.
    Returns a dict: node -> number of children
    """
    children_count = {node: len(children) for node, children in tree.items()}
    return children_count

def split_tree(tree, node):
    """
    Return all descendants (subtree) of a node, including itself.
    """
    return get_subtree(tree, node)

def find_nas(x):
    """
    Utility to check for NAs in a vector (returns 0 if any NA, else 1).
    """
    return 0 if pd.isnull(x).any() else 1

def refine_cna(res_df, tree, realcell_df):
    """
    Remove redundant lineages: if two lineages are associated with the same CNA and one lineage's cell-set is a subset of another's, drop the smaller.
    Keep all maximal (non-subset) lineages per region.
    """
    keep_idx = []
    for region in res_df['region'].unique():
        sub = res_df[res_df['region'] == region]
        # Build lineage sets
        lineage_sets = {row['cell']: split_tree(tree, row['cell']) for _, row in sub.iterrows()}
        cells = list(lineage_sets.keys())
        redundant = set()
        for i, c1 in enumerate(cells):
            for j, c2 in enumerate(cells):
                if i == j or c2 in redundant:
                    continue
                if lineage_sets[c1].issubset(lineage_sets[c2]):
                    redundant.add(c1)
                    break
        keep_cells = [c for c in cells if c not in redundant]
        keep_idx.extend(sub[sub['cell'].isin(keep_cells)].index.tolist())
    return res_df.loc[keep_idx].reset_index(drop=True)

def merge_lsa(res_df, tree, realcell_df):
    """
    Merge overlapping significant CNAs (basic version).
    res_df: DataFrame with columns ['cell', 'region', ...]
    tree: dict of dicts
    realcell_df: DataFrame with cell info (must include 'cell' and 'depth')
    Returns: merged DataFrame
    """
    # This is a placeholder for more advanced merging logic
    return res_df.drop_duplicates(['cell', 'region'])

def compute_lsa(tree, cnv_df, root):
    """
    For each non-root node, extract the lineage (subtree), calculate CFL for each feature.
    Returns a DataFrame: columns=[cell, depth, subtreesize, feature, cfl]
    """
    results = []
    for node in tree:
        if node == root:
            continue
        lineage = get_subtree(tree, node)
        cfl = lineage_cfl(cnv_df, lineage)
        results.append(pd.DataFrame({
            'cell': node,
            'depth': np.nan,  # to be filled if needed
            'subtreesize': len(lineage),
            'feature': cfl.index,
            'cfl': cfl.values
        }))
    return pd.concat(results, ignore_index=True)

def empirical_pvals(observed, permuted, alternative='greater'):
    """
    Compute empirical p-values for observed values vs. permutation background.
    observed: array-like, shape (n,)
    permuted: array-like, shape (n, n_perm)
    alternative: 'greater' or 'less'
    Returns: array of p-values
    """
    n_perm = permuted.shape[1]
    if alternative == 'greater':
        pvals = (np.sum(permuted >= observed[:, None], axis=1) + 1) / (n_perm + 1)
    else:
        pvals = (np.sum(permuted <= observed[:, None], axis=1) + 1) / (n_perm + 1)
    return pvals

def fdr_correction(pvals, alpha=0.05):
    """
    Benjamini-Hochberg FDR correction.
    If the input array is empty, returns two empty arrays.
    """
    pvals = np.asarray(pvals, dtype=float)
    if pvals.size == 0:
        return np.array([], dtype=bool), np.array([], dtype=float)
    rejected, pvals_corr, _, _ = multipletests(pvals, alpha=alpha, method='fdr_bh')
    return rejected, pvals_corr

def collect_significant_cnas(lsa_df, pval_col='pval', alpha=0.05):
    """
    Collect significant CNAs after FDR correction.
    lsa_df: DataFrame with columns including 'feature', 'cell', 'cfl', 'pval'
    Returns: DataFrame of significant CNAs
    """
    rejected, pvals_corr = fdr_correction(lsa_df[pval_col].values, alpha=alpha)
    lsa_df['fdr'] = pvals_corr
    lsa_df['significant'] = rejected
    return lsa_df[lsa_df['significant']].copy()

# Build cytoband order map for adjacency-aware merging
# (remove all code below) 