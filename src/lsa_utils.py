import numpy as np
import pandas as pd
from collections import defaultdict, deque
from statsmodels.stats.multitest import multipletests
import os

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
    Calculate cumulative fold level (CFL) for a set of cells (lineage) for each feature.
    Returns a pandas Series (index: feature, value: mean CNV across lineage).
    """
    return cnv_df.loc[list(lineage_cells)].mean(axis=0)

def lineage_cfl_per_region(cnv_df, lineage_cells, region):
    """
    Calculate CFL for a specific region in a lineage.
    Returns the mean CNV value minus 2 (deviation from diploid).
    """
    lineage_cnvs = cnv_df.loc[list(lineage_cells), region]
    return lineage_cnvs.mean() - 2

def get_depths(tree, root):
    """
    Compute depth of each node in the tree (distance from root).
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
    """
    children_count = {}
    for node in tree:
        children_count[node] = len(tree.get(node, {}))
    return children_count

def empirical_pvals_per_region(lsa_df, tree, real_cnv, permut_path, n_perms=500, root='root', debug=False):
    """
    Calculate empirical p-values for each (node, region) pair using permuted data.
    """
    pvals = []
    
    # Group by (cell, region) for efficient lookup
    obs_scores = {}
    for _, row in lsa_df.iterrows():
        key = (row['cell'], row['region'])
        obs_scores[key] = row['Score']
    
    # Collect permuted scores
    perm_scores = defaultdict(list)
    
    if debug:
        print(f"DEBUG: Loading {n_perms} permutation files...")
    
    found_perms = 0
    for j in range(1, n_perms + 1):
        perm_cnv_path = os.path.join(permut_path, f"permute.{j}.CNV.txt")
        if not os.path.exists(perm_cnv_path):
            continue
        
        found_perms += 1
        perm_cnv = pd.read_csv(perm_cnv_path, sep="\t", index_col=0)
        
        if debug and j == 1:  # Debug first permutation file
            print(f"DEBUG: Permutation columns: {list(perm_cnv.columns)}")
            print(f"DEBUG: Real columns: {list(real_cnv.columns)}")
            common_regions = set(perm_cnv.columns) & set(real_cnv.columns)
            print(f"DEBUG: Common regions: {list(common_regions)}")
        
        # For each node in the tree
        for node in tree:
            if node == root:
                continue
            
            lineage = get_subtree(tree, node)
            lineage = [cell for cell in lineage if cell != 'root']
            
            if len(lineage) < 5:
                continue
            
            # Calculate score for each region - use only regions that exist in both
            common_regions = set(perm_cnv.columns) & set(real_cnv.columns)
            for region in common_regions:
                lineage_cnvs = perm_cnv.loc[lineage, region]
                score = lineage_cnvs.mean() - 2
                perm_scores[(node, region)].append(score)
    
    if debug:
        print(f"DEBUG: Found {found_perms} permutation files")
        print(f"DEBUG: Collected permutation scores for {len(perm_scores)} (node, region) pairs")
    
    # Calculate p-values
    for _, row in lsa_df.iterrows():
        key = (row['cell'], row['region'])
        observed = row['Score']
        perm_vals = np.array(perm_scores.get(key, []))
        
        if len(perm_vals) == 0:
            pval = 1.0
        else:
            # Two-tailed test
            if observed > 0:  # Amplification
                pval = (np.sum(perm_vals >= observed) + 1) / (len(perm_vals) + 1)
            else:  # Deletion
                pval = (np.sum(perm_vals <= observed) + 1) / (len(perm_vals) + 1)
        
        pvals.append(pval)
    
    if debug:
        pval_array = np.array(pvals)
        print(f"DEBUG: P-value distribution: min={pval_array.min():.4f}, "
              f"median={np.median(pval_array):.4f}, max={pval_array.max():.4f}")
        print(f"DEBUG: P-values < 0.05: {np.sum(pval_array < 0.05)}")
    
    return pvals

def fdr_correction(pvals, alpha=0.05):
    """
    Benjamini-Hochberg FDR correction.
    """
    pvals = np.asarray(pvals, dtype=float)
    if pvals.size == 0:
        return np.array([], dtype=bool), np.array([], dtype=float)
    
    # Handle case where all p-values are 1
    if np.all(pvals == 1):
        return np.zeros(len(pvals), dtype=bool), pvals
    
    rejected, pvals_corr, _, _ = multipletests(pvals, alpha=alpha, method='fdr_bh')
    return rejected, pvals_corr

def refine_cna(res_df, tree, realcell_df):
    """
    Remove redundant lineages: if two lineages have the same CNA and one is a subset of another.
    """
    if res_df.empty:
        return res_df
    
    keep_idx = []
    
    # Group by region and CNA type
    for (region, cna_type), group in res_df.groupby(['region', 'CNA']):
        # Build lineage sets for this region/CNA
        lineage_sets = {}
        for _, row in group.iterrows():
            cell = row['cell']
            lineage = get_subtree(tree, cell)
            lineage_sets[cell] = lineage
        
        cells = list(lineage_sets.keys())
        redundant = set()
        
        # Find redundant (subset) lineages
        for i, c1 in enumerate(cells):
            if c1 in redundant:
                continue
            for j, c2 in enumerate(cells):
                if i == j or c2 in redundant:
                    continue
                # If c1 is a subset of c2, mark c1 as redundant
                if lineage_sets[c1].issubset(lineage_sets[c2]) and lineage_sets[c1] != lineage_sets[c2]:
                    redundant.add(c1)
                    break
        
        # Keep non-redundant cells
        keep_cells = [c for c in cells if c not in redundant]
        keep_idx.extend(group[group['cell'].isin(keep_cells)].index.tolist())
    
    return res_df.loc[keep_idx].reset_index(drop=True)

def lineage_score(cnv_df, lineage_cells, amp_cutoff=0.5, del_cutoff=-0.5):
    """
    Calculate AMP/DEL scores for a lineage using R's exact algorithm.
    This matches the lineageScore function in LSAfunction.R lines 709-733.
    """
    subcnv = cnv_df.loc[lineage_cells]
    total_cells = len(lineage_cells)
    
    # Initialize result arrays
    amp_scores = []
    del_scores = []
    
    # Process each region (column)
    for region in subcnv.columns:
        x = subcnv[region].values
        
        # Count cells with amplifications (>2) and deletions (<2)
        f1 = len(x[x > 2])  # amplification count
        f2 = len(x[x < 2])  # deletion count
        
        # Calculate amplification score (R logic)
        if f1 > 1:
            Gamp = f1 * np.mean(x[x > 2] - 2) / total_cells
        elif f1 == 1:
            Gamp = f1 * (x[x > 2] - 2)[0] / total_cells
        else:
            Gamp = 0
            
        # Calculate deletion score (R logic)
        if f2 > 1:
            Gdel = f2 * np.mean(x[x < 2] - 2) / total_cells
        elif f2 == 1:
            Gdel = f2 * (x[x < 2] - 2)[0] / total_cells
        else:
            Gdel = 0
            
        amp_scores.append(Gamp)
        del_scores.append(Gdel)
    
    return pd.DataFrame({
        'AMP': amp_scores, 
        'DEL': del_scores
    }, index=subcnv.columns)

def split_tree(tree, node):
    """
    Return all descendants (subtree) of a node, including itself.
    """
    return get_subtree(tree, node)

def collect_significant_cnas(lsa_df, pval_col='pvalue', alpha=0.05):
    """
    Collect significant CNAs after FDR correction.
    """
    rejected, pvals_corr = fdr_correction(lsa_df[pval_col].values, alpha=alpha)
    lsa_df['adjustp'] = pvals_corr
    lsa_df['significant'] = rejected
    return lsa_df[lsa_df['significant']].copy()

def find_nas(x):
    """
    Utility to check for NAs in a vector.
    """
    return 0 if pd.isnull(x).any() else 1