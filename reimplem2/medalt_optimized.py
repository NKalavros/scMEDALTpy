#!/usr/bin/env python
"""
MEDALT Optimized: Fast implementation for large-scale scRNA-seq CNV analysis
Optimized for 10k cells × 1k genes using inferCNV output
"""

import numpy as np
import pandas as pd
import numba
from numba import jit, prange
import sparse
from scipy import stats
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import networkx as nx
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import h5py
import zarr
from typing import Dict, List, Tuple, Optional, Set
import logging
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import multiprocessing as mp
from functools import lru_cache
import warnings
warnings.filterwarnings('ignore', category=numba.NumbaWarning)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class InferCNVReader:
    """Efficient reader for inferCNV output files"""
    
    @staticmethod
    def load_observations(filepath: str, top_k_genes: int = 1000, 
                         min_variance: float = 0.1) -> Tuple[np.ndarray, List[str], List[str], np.ndarray]:
        """
        Load inferCNV observations file with gene selection
        
        Args:
            filepath: Path to infercnv.observations.txt
            top_k_genes: Number of top variance genes to select
            min_variance: Minimum variance threshold for gene selection
            
        Returns:
            cnv_matrix: Cells × Genes matrix (integer copy numbers)
            cell_names: List of cell names
            gene_names: List of selected gene names
            gene_indices: Original indices of selected genes
        """
        logger.info(f"Loading inferCNV observations from {filepath}")
        
        # Read header to get dimensions
        with open(filepath, 'r') as f:
            header = f.readline().strip().split('\t')
            cell_names = header[1:]  # First column is gene names
            n_cells = len(cell_names)
        
        # Use chunks for memory efficiency
        chunk_size = 1000
        gene_data = []
        gene_names_all = []
        
        # Read in chunks and calculate variance on the fly
        for chunk in pd.read_csv(filepath, sep='\t', index_col=0, chunksize=chunk_size):
            # Calculate variance for each gene in chunk
            variances = chunk.var(axis=1)
            gene_data.append((chunk, variances))
            gene_names_all.extend(chunk.index.tolist())
        
        logger.info(f"Loaded {len(gene_names_all)} genes from {n_cells} cells")
        
        # Select top variance genes
        all_variances = pd.concat([gd[1] for gd in gene_data])
        
        # Filter by minimum variance and select top k
        high_var_mask = all_variances > min_variance
        high_var_genes = all_variances[high_var_mask]
        
        if len(high_var_genes) < top_k_genes:
            logger.warning(f"Only {len(high_var_genes)} genes pass variance threshold {min_variance}")
            top_k_genes = len(high_var_genes)
        
        top_genes = high_var_genes.nlargest(top_k_genes)
        selected_gene_names = top_genes.index.tolist()
        gene_indices = np.array([gene_names_all.index(g) for g in selected_gene_names])
        
        logger.info(f"Selected {len(selected_gene_names)} high-variance genes")
        logger.info(f"Variance range: {top_genes.min():.3f} - {top_genes.max():.3f}")
        
        # Build final matrix with only selected genes
        selected_data = []
        current_idx = 0
        
        for chunk_df, _ in gene_data:
            # Get genes from this chunk that are selected
            chunk_selected = chunk_df.loc[chunk_df.index.intersection(selected_gene_names)]
            if not chunk_selected.empty:
                selected_data.append(chunk_selected)
        
        # Combine and convert to matrix
        cnv_df = pd.concat(selected_data)
        cnv_df = cnv_df.loc[selected_gene_names]  # Ensure correct order
        
        # Convert to integer copy numbers
        # inferCNV uses centered log2 ratios, convert back
        cnv_matrix = np.round(2 ** cnv_df.values).astype(np.int8)  # int8 saves memory
        cnv_matrix = cnv_matrix.T  # Transpose to cells × genes
        
        # Clip extreme values
        cnv_matrix = np.clip(cnv_matrix, 0, 10)
        
        return cnv_matrix, cell_names, selected_gene_names, gene_indices
    
    @staticmethod
    def load_references(filepath: str, gene_indices: np.ndarray) -> np.ndarray:
        """Load reference cells with same gene selection"""
        ref_df = pd.read_csv(filepath, sep='\t', index_col=0)
        # Select same genes
        ref_matrix = ref_df.iloc[gene_indices].values
        ref_matrix = np.round(2 ** ref_matrix).astype(np.int8)
        return np.clip(ref_matrix.T, 0, 10)


@jit(nopython=True, parallel=False)
def compute_med_fast(a: np.ndarray, b: np.ndarray) -> int:
    """
    Numba-optimized MED calculation
    
    Args:
        a, b: Integer arrays representing copy number profiles
        
    Returns:
        Integer MED (uses large value instead of inf for numba compatibility)
    """
    n = len(a)
    if n != len(b):
        return 999999
    
    # Check homozygous deletion constraint
    for i in range(n):
        if a[i] == 0 and b[i] > 0:
            return 999999
    
    distance = 0
    i = 0
    
    while i < n:
        if a[i] == b[i]:
            i += 1
            continue
        
        # Find segment with consistent direction
        j = i
        direction = 1 if b[i] > a[i] else -1
        
        while j < n and a[j] != b[j]:
            curr_dir = 1 if b[j] > a[j] else -1
            if curr_dir != direction:
                break
            j += 1
        
        # Calculate max difference in segment
        max_diff = 0
        for k in range(i, j):
            diff = abs(b[k] - a[k])
            if diff > max_diff:
                max_diff = diff
        
        distance += max_diff
        i = j
    
    return distance


@jit(nopython=True, parallel=True)
def compute_distance_matrix_parallel(cnv_matrix: np.ndarray, 
                                   batch_start: int,
                                   batch_end: int) -> np.ndarray:
    """
    Compute a batch of pairwise distances in parallel
    
    Args:
        cnv_matrix: n_cells × n_genes matrix
        batch_start: Start index for this batch
        batch_end: End index for this batch
        
    Returns:
        Distance matrix for the batch
    """
    n_cells = cnv_matrix.shape[0]
    batch_size = batch_end - batch_start
    distances = np.full((batch_size, n_cells), 999999, dtype=np.int32)
    
    for i in prange(batch_size):
        cell_i = batch_start + i
        for cell_j in range(n_cells):
            if cell_i != cell_j:
                distances[i, cell_j] = compute_med_fast(
                    cnv_matrix[cell_i], cnv_matrix[cell_j]
                )
            else:
                distances[i, cell_j] = 0
    
    return distances


class OptimizedMEDALT:
    """Memory and speed optimized MEDALT implementation"""
    
    def __init__(self, cnv_matrix: np.ndarray, 
                 reference_profile: Optional[np.ndarray] = None,
                 n_components: int = 50,
                 batch_size: int = 100):
        """
        Initialize optimized MEDALT
        
        Args:
            cnv_matrix: n_cells × n_genes integer matrix
            reference_profile: Optional reference (normal) profile
            n_components: Number of PCA components for dimensionality reduction
            batch_size: Batch size for distance calculations
        """
        self.cnv_matrix = cnv_matrix.astype(np.int8)  # Memory efficient
        self.n_cells, self.n_genes = cnv_matrix.shape
        self.batch_size = batch_size
        
        # Set reference profile
        if reference_profile is None:
            self.reference_profile = np.full(self.n_genes, 2, dtype=np.int8)
        else:
            self.reference_profile = reference_profile.astype(np.int8)
        
        # For very large datasets, use PCA for initial clustering
        self.use_pca = self.n_cells > 5000
        if self.use_pca:
            logger.info(f"Using PCA with {n_components} components for initial clustering")
            self.pca_components = self._compute_pca(n_components)
        
        logger.info(f"Initialized OptimizedMEDALT with {self.n_cells} cells, {self.n_genes} genes")
    
    def _compute_pca(self, n_components: int) -> np.ndarray:
        """Compute PCA for dimensionality reduction"""
        # Use incremental PCA for large datasets
        from sklearn.decomposition import IncrementalPCA
        
        pca = IncrementalPCA(n_components=n_components, batch_size=1000)
        # Center data around diploid state
        centered_data = self.cnv_matrix - 2
        return pca.fit_transform(centered_data)
    
    def build_tree_fast(self, subsample_size: Optional[int] = None) -> nx.DiGraph:
        """
        Build MEDALT tree with optimizations
        
        Args:
            subsample_size: If provided, build tree on subsample first
            
        Returns:
            NetworkX directed graph
        """
        if subsample_size and self.n_cells > subsample_size:
            logger.info(f"Building initial tree on {subsample_size} subsampled cells")
            tree = self._build_subsampled_tree(subsample_size)
            logger.info("Assigning remaining cells to tree")
            tree = self._assign_remaining_cells(tree, subsample_size)
        else:
            tree = self._build_complete_tree()
        
        return tree
    
    def _build_subsampled_tree(self, subsample_size: int) -> nx.DiGraph:
        """Build tree on subsampled cells"""
        # Use PCA for intelligent subsampling if available
        if self.use_pca:
            # K-means clustering on PCA components
            from sklearn.cluster import MiniBatchKMeans
            kmeans = MiniBatchKMeans(n_clusters=subsample_size, batch_size=1000)
            labels = kmeans.fit_predict(self.pca_components)
            
            # Select one cell per cluster (closest to centroid)
            selected_indices = []
            for i in range(subsample_size):
                cluster_cells = np.where(labels == i)[0]
                if len(cluster_cells) > 0:
                    # Find cell closest to centroid
                    centroid = kmeans.cluster_centers_[i]
                    distances = np.sum((self.pca_components[cluster_cells] - centroid) ** 2, axis=1)
                    selected_indices.append(cluster_cells[np.argmin(distances)])
            
            selected_indices = np.array(selected_indices)
        else:
            # Random subsampling
            selected_indices = np.random.choice(self.n_cells, subsample_size, replace=False)
        
        # Build tree on subsample
        sub_matrix = self.cnv_matrix[selected_indices]
        sub_tree = self._build_tree_from_matrix(sub_matrix, selected_indices)
        
        return sub_tree
    
    def _build_tree_from_matrix(self, cnv_matrix: np.ndarray, 
                               cell_indices: Optional[np.ndarray] = None) -> nx.DiGraph:
        """Build tree from given matrix"""
        n_cells = cnv_matrix.shape[0]
        if cell_indices is None:
            cell_indices = np.arange(n_cells)
        
        # Compute distances to reference
        ref_distances = np.array([
            compute_med_fast(self.reference_profile, cnv_matrix[i])
            for i in range(n_cells)
        ])
        
        # Build sparse distance matrix for efficiency
        # Only compute k-nearest neighbors instead of all pairs
        k = min(50, n_cells - 1)  # Number of nearest neighbors
        
        if n_cells < 1000:
            # For small datasets, compute full matrix
            distances = self._compute_full_distances(cnv_matrix)
            return self._mst_from_full_distances(distances, ref_distances, cell_indices)
        else:
            # For large datasets, use approximate approach
            return self._build_approximate_tree(cnv_matrix, ref_distances, cell_indices, k)
    
    def _compute_full_distances(self, cnv_matrix: np.ndarray) -> np.ndarray:
        """Compute full distance matrix with batching"""
        n_cells = cnv_matrix.shape[0]
        distances = np.full((n_cells, n_cells), 999999, dtype=np.int32)
        
        # Process in batches for memory efficiency
        n_batches = (n_cells + self.batch_size - 1) // self.batch_size
        
        with ProcessPoolExecutor() as executor:
            futures = []
            
            for batch_idx in range(n_batches):
                batch_start = batch_idx * self.batch_size
                batch_end = min((batch_idx + 1) * self.batch_size, n_cells)
                
                future = executor.submit(
                    compute_distance_matrix_parallel,
                    cnv_matrix, batch_start, batch_end
                )
                futures.append((batch_start, batch_end, future))
            
            # Collect results
            for batch_start, batch_end, future in futures:
                batch_distances = future.result()
                distances[batch_start:batch_end] = batch_distances
        
        return distances
    
    def _mst_from_full_distances(self, distances: np.ndarray, 
                                ref_distances: np.ndarray,
                                cell_indices: np.ndarray) -> nx.DiGraph:
        """Build MST from full distance matrix"""
        n_cells = len(cell_indices)
        
        # Create graph with root
        G = nx.DiGraph()
        G.add_node('root')
        
        for i in range(n_cells):
            G.add_node(cell_indices[i])
            G.add_edge('root', cell_indices[i], weight=ref_distances[i])
        
        # Add edges between cells (convert to sparse for efficiency)
        # Only add edges below a threshold to keep graph sparse
        threshold = np.percentile(distances[distances < 999999], 90)
        
        for i in range(n_cells):
            for j in range(n_cells):
                if i != j and distances[i, j] < threshold:
                    G.add_edge(cell_indices[i], cell_indices[j], 
                             weight=distances[i, j])
        
        # Find MST
        try:
            mst = nx.minimum_spanning_arborescence(G)
        except:
            # Fallback to greedy
            mst = self._greedy_tree(G, cell_indices)
        
        return mst
    
    def _build_approximate_tree(self, cnv_matrix: np.ndarray,
                               ref_distances: np.ndarray,
                               cell_indices: np.ndarray,
                               k: int) -> nx.DiGraph:
        """Build approximate tree using k-NN approach"""
        from sklearn.neighbors import NearestNeighbors
        
        n_cells = cnv_matrix.shape[0]
        
        # Use PCA components if available for faster NN search
        if hasattr(self, 'pca_components'):
            # Map indices to PCA space
            pca_subset = self.pca_components[cell_indices]
            nn = NearestNeighbors(n_neighbors=k+1, algorithm='ball_tree')
            nn.fit(pca_subset)
            nn_indices = nn.kneighbors(pca_subset, return_distance=False)
        else:
            # Fall back to random sampling
            nn_indices = np.array([
                np.random.choice(n_cells, k+1, replace=False) 
                for _ in range(n_cells)
            ])
        
        # Build graph with k-NN edges
        G = nx.DiGraph()
        G.add_node('root')
        
        for i in range(n_cells):
            G.add_node(cell_indices[i])
            G.add_edge('root', cell_indices[i], weight=ref_distances[i])
        
        # Add k-NN edges
        for i in range(n_cells):
            for j_idx in nn_indices[i][1:]:  # Skip self
                if j_idx < n_cells:
                    dist = compute_med_fast(cnv_matrix[i], cnv_matrix[j_idx])
                    if dist < 999999:
                        G.add_edge(cell_indices[i], cell_indices[j_idx], weight=dist)
        
        # Find MST
        try:
            mst = nx.minimum_spanning_arborescence(G)
        except:
            mst = self._greedy_tree(G, cell_indices)
        
        return mst
    
    def _greedy_tree(self, G: nx.Graph, cell_indices: np.ndarray) -> nx.DiGraph:
        """Fallback greedy tree construction"""
        tree = nx.DiGraph()
        tree.add_node('root')
        
        # Add all nodes
        for idx in cell_indices:
            tree.add_node(idx)
        
        # Greedy assignment
        connected = {'root'}
        unconnected = set(cell_indices)
        
        while unconnected:
            best_edge = None
            best_weight = float('inf')
            
            # Find best edge from connected to unconnected
            for u in connected:
                for v in unconnected:
                    if G.has_edge(u, v):
                        weight = G[u][v]['weight']
                        if weight < best_weight:
                            best_weight = weight
                            best_edge = (u, v)
            
            if best_edge:
                tree.add_edge(best_edge[0], best_edge[1], weight=best_weight)
                connected.add(best_edge[1])
                unconnected.remove(best_edge[1])
            else:
                # Connect remaining directly to root
                for v in unconnected:
                    tree.add_edge('root', v, weight=999999)
                break
        
        return tree
    
    def _assign_remaining_cells(self, tree: nx.DiGraph, subsample_size: int) -> nx.DiGraph:
        """Assign remaining cells to existing tree"""
        assigned_cells = set(n for n in tree.nodes() if n != 'root')
        unassigned = set(range(self.n_cells)) - assigned_cells
        
        # Get tree cell profiles
        tree_cells = list(assigned_cells)
        tree_profiles = self.cnv_matrix[tree_cells]
        
        # Process unassigned cells in batches
        unassigned_list = list(unassigned)
        batch_size = 1000
        
        for i in range(0, len(unassigned_list), batch_size):
            batch = unassigned_list[i:i+batch_size]
            batch_profiles = self.cnv_matrix[batch]
            
            # Find nearest tree node for each cell
            for j, cell_idx in enumerate(batch):
                # Compute distances to all tree nodes
                distances = [
                    compute_med_fast(batch_profiles[j], tree_profiles[k])
                    for k in range(len(tree_cells))
                ]
                
                # Also check root
                ref_dist = compute_med_fast(self.reference_profile, batch_profiles[j])
                
                # Find minimum
                min_dist = ref_dist
                parent = 'root'
                
                for k, dist in enumerate(distances):
                    if dist < min_dist:
                        min_dist = dist
                        parent = tree_cells[k]
                
                # Add to tree
                tree.add_node(cell_idx)
                tree.add_edge(parent, cell_idx, weight=min_dist)
        
        return tree
    
    def _build_complete_tree(self) -> nx.DiGraph:
        """Build tree on complete dataset"""
        return self._build_tree_from_matrix(self.cnv_matrix)


class OptimizedLSA:
    """Optimized Lineage Speciation Analysis for large datasets"""
    
    def __init__(self, tree: nx.DiGraph, cnv_matrix: np.ndarray,
                 n_permutations: int = 500,
                 use_approximation: bool = True):
        """
        Initialize optimized LSA
        
        Args:
            tree: MEDALT tree
            cnv_matrix: n_cells × n_genes matrix
            n_permutations: Number of permutations
            use_approximation: Use faster approximation methods
        """
        self.tree = tree
        self.cnv_matrix = cnv_matrix.astype(np.int8)
        self.n_cells, self.n_genes = cnv_matrix.shape
        self.n_permutations = n_permutations
        self.use_approximation = use_approximation
        
        # Precompute tree structure
        self._precompute_tree_structure()
        
        logger.info(f"Initialized OptimizedLSA with {n_permutations} permutations")
    
    def _precompute_tree_structure(self):
        """Precompute tree properties for efficiency"""
        # Cache descendants for each node
        self.descendants_cache = {}
        self.subtree_sizes = {}
        
        for node in self.tree.nodes():
            if node != 'root':
                descendants = set(nx.descendants(self.tree, node))
                descendants.add(node)
                self.descendants_cache[node] = descendants
                self.subtree_sizes[node] = len(descendants)
    
    def run_analysis_fast(self, min_lineage_size: int = 20,
                         significance_threshold: float = 0.05) -> pd.DataFrame:
        """
        Run optimized LSA
        
        Args:
            min_lineage_size: Minimum lineage size
            significance_threshold: P-value threshold
            
        Returns:
            DataFrame with results
        """
        # Step 1: Find significant lineages
        logger.info("Finding significant lineages...")
        lineages = self._find_lineages_fast(min_lineage_size)
        
        if not lineages:
            logger.warning("No lineages found!")
            return pd.DataFrame()
        
        logger.info(f"Found {len(lineages)} lineages")
        
        # Step 2: Identify high-variance genes
        logger.info("Identifying high-variance genes...")
        target_genes = self._identify_target_genes()
        
        # Step 3: Calculate observed statistics
        logger.info("Calculating observed statistics...")
        observed_stats = self._calculate_observed_stats_fast(lineages, target_genes)
        
        # Step 4: Permutation testing
        logger.info(f"Running {self.n_permutations} permutations...")
        p_values = self._permutation_test_fast(observed_stats, lineages, target_genes)
        
        # Step 5: Format results
        results = self._format_results_fast(observed_stats, p_values, lineages)
        
        # Filter by significance
        results = results[results['pvalue'] < significance_threshold]
        
        return results
    
    def _find_lineages_fast(self, min_size: int) -> List[Dict]:
        """Fast lineage identification"""
        lineages = []
        
        # Sort nodes by subtree size
        candidate_nodes = [
            (node, size) for node, size in self.subtree_sizes.items()
            if size >= min_size
        ]
        candidate_nodes.sort(key=lambda x: x[1])
        
        # Remove redundant lineages
        selected_nodes = set()
        
        for node, size in candidate_nodes:
            # Check if this node is already covered
            is_covered = False
            for selected in selected_nodes:
                if node in self.descendants_cache.get(selected, set()):
                    is_covered = True
                    break
            
            if not is_covered:
                selected_nodes.add(node)
                
                # Get cells (excluding root)
                cells = [n for n in self.descendants_cache[node] if n != 'root']
                
                lineages.append({
                    'root': node,
                    'cells': np.array(cells, dtype=np.int32),
                    'size': len(cells)
                })
        
        return lineages
    
    def _identify_target_genes(self, top_k: int = 100) -> np.ndarray:
        """Identify genes with high variance across lineages"""
        # Calculate variance across all cells
        gene_vars = np.var(self.cnv_matrix, axis=0)
        
        # Also consider genes with bimodal distributions
        gene_bimodality = np.zeros(self.n_genes)
        
        for i in range(self.n_genes):
            values = self.cnv_matrix[:, i]
            unique_vals, counts = np.unique(values, return_counts=True)
            
            if len(unique_vals) >= 2:
                # Simple bimodality score
                sorted_counts = np.sort(counts)[::-1]
                if len(sorted_counts) >= 2:
                    gene_bimodality[i] = sorted_counts[1] / sorted_counts[0]
        
        # Combine variance and bimodality
        gene_scores = gene_vars + gene_bimodality
        
        # Select top genes
        target_genes = np.argsort(gene_scores)[-top_k:]
        
        return target_genes
    
    def _calculate_observed_stats_fast(self, lineages: List[Dict], 
                                     target_genes: np.ndarray) -> Dict:
        """Calculate observed statistics efficiently"""
        observed = {}
        
        for lineage in lineages:
            lineage_id = lineage['root']
            cells = lineage['cells']
            
            # Get CNV data for lineage cells
            lineage_cnvs = self.cnv_matrix[cells][:, target_genes]
            
            # Calculate statistics per gene
            for i, gene_idx in enumerate(target_genes):
                gene_cnvs = lineage_cnvs[:, i]
                
                # Mean copy number
                mean_cn = np.mean(gene_cnvs)
                
                # Skip if neutral
                if abs(mean_cn - 2.0) < 0.5:
                    continue
                
                # Determine direction
                direction = 'AMP' if mean_cn > 2.0 else 'DEL'
                
                # Calculate enrichment score
                # Compare to background (all other cells)
                other_cells = np.setdiff1d(np.arange(self.n_cells), cells)
                background_mean = np.mean(self.cnv_matrix[other_cells, gene_idx])
                
                enrichment = abs(mean_cn - background_mean)
                
                observed[(lineage_id, gene_idx, direction)] = {
                    'mean_cn': mean_cn,
                    'enrichment': enrichment,
                    'lineage_size': len(cells),
                    'cells': cells
                }
        
        return observed
    
    def _permutation_test_fast(self, observed_stats: Dict,
                             lineages: List[Dict],
                             target_genes: np.ndarray) -> Dict:
        """Fast permutation testing using parallel processing"""
        # Prepare data for parallel processing
        lineage_sizes = [lin['size'] for lin in lineages]
        
        # Use multiprocessing for permutations
        n_cores = mp.cpu_count()
        chunk_size = self.n_permutations // n_cores
        
        with ProcessPoolExecutor(max_workers=n_cores) as executor:
            futures = []
            
            for i in range(n_cores):
                start_idx = i * chunk_size
                end_idx = start_idx + chunk_size if i < n_cores - 1 else self.n_permutations
                
                future = executor.submit(
                    self._run_permutation_chunk,
                    start_idx, end_idx,
                    lineage_sizes, target_genes,
                    self.cnv_matrix, self.n_cells
                )
                futures.append(future)
            
            # Collect results
            all_null_stats = []
            for future in futures:
                all_null_stats.extend(future.result())
        
        # Calculate p-values
        p_values = {}
        
        for key, obs_data in observed_stats.items():
            lineage_id, gene_idx, direction = key
            obs_enrichment = obs_data['enrichment']
            
            # Find matching null distribution
            lineage_size = obs_data['lineage_size']
            null_enrichments = []
            
            for null_stats in all_null_stats:
                if (lineage_size, gene_idx, direction) in null_stats:
                    null_enrichments.append(null_stats[(lineage_size, gene_idx, direction)])
            
            if null_enrichments:
                null_enrichments = np.array(null_enrichments)
                p_value = np.mean(null_enrichments >= obs_enrichment)
                p_values[key] = max(p_value, 1.0 / len(null_enrichments))  # Avoid p=0
            else:
                p_values[key] = 1.0
        
        return p_values
    
    @staticmethod
    def _run_permutation_chunk(start_idx: int, end_idx: int,
                              lineage_sizes: List[int],
                              target_genes: np.ndarray,
                              cnv_matrix: np.ndarray,
                              n_cells: int) -> List[Dict]:
        """Run a chunk of permutations (for parallel processing)"""
        np.random.seed(start_idx)  # Reproducible randomness
        
        results = []
        
        for perm_idx in range(start_idx, end_idx):
            null_stats = {}
            
            # Permute cell labels
            permuted_indices = np.random.permutation(n_cells)
            
            # Calculate statistics for each lineage size
            current_idx = 0
            for lineage_size in lineage_sizes:
                # Get random cells of this size
                perm_cells = permuted_indices[current_idx:current_idx + lineage_size]
                current_idx += lineage_size
                
                # Calculate stats for target genes
                for gene_idx in target_genes:
                    gene_cnvs = cnv_matrix[perm_cells, gene_idx]
                    mean_cn = np.mean(gene_cnvs)
                    
                    if abs(mean_cn - 2.0) < 0.5:
                        continue
                    
                    direction = 'AMP' if mean_cn > 2.0 else 'DEL'
                    
                    # Enrichment vs background
                    other_cells = np.setdiff1d(permuted_indices, perm_cells)
                    background_mean = np.mean(cnv_matrix[other_cells, gene_idx])
                    enrichment = abs(mean_cn - background_mean)
                    
                    null_stats[(lineage_size, gene_idx, direction)] = enrichment
            
            results.append(null_stats)
        
        return results
    
    def _format_results_fast(self, observed_stats: Dict, 
                           p_values: Dict,
                           lineages: List[Dict]) -> pd.DataFrame:
        """Format results into DataFrame"""
        rows = []
        
        for key, obs_data in observed_stats.items():
            lineage_id, gene_idx, direction = key
            
            if key in p_values:
                rows.append({
                    'gene_idx': gene_idx,
                    'lineage_root': lineage_id,
                    'direction': direction,
                    'mean_cn': obs_data['mean_cn'],
                    'enrichment_score': obs_data['enrichment'],
                    'lineage_size': obs_data['lineage_size'],
                    'pvalue': p_values[key]
                })
        
        df = pd.DataFrame(rows)
        
        if not df.empty:
            # Add FDR correction
            df['qvalue'] = stats.false_discovery_control(df['pvalue'].values)
            
            # Sort by significance
            df = df.sort_values('pvalue')
        
        return df


class InferCNVPipeline:
    """Complete pipeline for inferCNV data analysis"""
    
    def __init__(self, output_dir: str = 'medalt_results'):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # Set up logging
        log_file = os.path.join(output_dir, 'medalt_analysis.log')
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        ))
        logger.addHandler(file_handler)
    
    def run_infercnv_analysis(self, 
                            observations_file: str,
                            reference_file: Optional[str] = None,
                            top_k_genes: int = 1000,
                            subsample_size: int = 2000,
                            n_permutations: int = 500,
                            min_lineage_size: int = 20):
        """
        Run complete analysis on inferCNV output
        
        Args:
            observations_file: Path to infercnv.observations.txt
            reference_file: Optional path to infercnv.references.txt
            top_k_genes: Number of high-variance genes to use
            subsample_size: Initial tree building subsample size
            n_permutations: Number of LSA permutations
            min_lineage_size: Minimum lineage size for LSA
        """
        # Load data
        logger.info("Loading inferCNV observations...")
        cnv_matrix, cell_names, gene_names, gene_indices = InferCNVReader.load_observations(
            observations_file, top_k_genes
        )
        
        # Save gene list
        gene_file = os.path.join(self.output_dir, 'selected_genes.txt')
        pd.DataFrame({'gene': gene_names, 'original_index': gene_indices}).to_csv(
            gene_file, sep='\t', index=False
        )
        
        # Build tree
        logger.info("Building MEDALT tree...")
        medalt = OptimizedMEDALT(cnv_matrix)
        tree = medalt.build_tree_fast(subsample_size=subsample_size)
        
        # Save tree
        self._save_tree(tree, cell_names)
        
        # Run LSA
        logger.info("Running Lineage Speciation Analysis...")
        lsa = OptimizedLSA(tree, cnv_matrix, n_permutations=n_permutations)
        results = lsa.run_analysis_fast(min_lineage_size=min_lineage_size)
        
        # Map gene indices back to names
        if not results.empty:
            results['gene_name'] = results['gene_idx'].apply(lambda x: gene_names[x])
            
            # Save results
            results_file = os.path.join(self.output_dir, 'lsa_results.txt')
            results.to_csv(results_file, sep='\t', index=False)
            
            logger.info(f"Found {len(results)} significant gene-lineage associations")
            logger.info(f"Results saved to {results_file}")
        
        # Generate summary
        self._save_summary(tree, results, cnv_matrix.shape)
        
        return tree, results
    
    def _save_tree(self, tree: nx.DiGraph, cell_names: List[str]):
        """Save tree in standard format"""
        tree_file = os.path.join(self.output_dir, 'medalt_tree.txt')
        
        with open(tree_file, 'w') as f:
            f.write("parent\tchild\tdistance\n")
            
            for parent, child, data in tree.edges(data=True):
                parent_name = 'root' if parent == 'root' else cell_names[parent]
                child_name = cell_names[child] if isinstance(child, int) else child
                
                f.write(f"{parent_name}\t{child_name}\t{data.get('weight', 0)}\n")
    
    def _save_summary(self, tree: nx.DiGraph, results: pd.DataFrame, data_shape: Tuple):
        """Save analysis summary"""
        summary = {
            'data_dimensions': {
                'n_cells': data_shape[0],
                'n_genes': data_shape[1]
            },
            'tree_statistics': {
                'n_nodes': tree.number_of_nodes() - 1,  # Exclude root
                'n_edges': tree.number_of_edges(),
                'max_depth': max(nx.shortest_path_length(tree, 'root').values())
            },
            'lsa_results': {
                'n_significant': len(results) if not results.empty else 0,
                'n_amplifications': len(results[results['direction'] == 'AMP']) if not results.empty else 0,
                'n_deletions': len(results[results['direction'] == 'DEL']) if not results.empty else 0
            }
        }
        
        import json
        summary_file = os.path.join(self.output_dir, 'analysis_summary.json')
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)


def main():
    """Example usage"""
    import argparse
    
    parser = argparse.ArgumentParser(description='MEDALT Optimized for inferCNV output')
    parser.add_argument('observations', help='Path to infercnv.observations.txt')
    parser.add_argument('--references', help='Path to infercnv.references.txt (optional)')
    parser.add_argument('--genes', type=int, default=1000, help='Number of top variance genes')
    parser.add_argument('--subsample', type=int, default=2000, help='Subsample size for tree building')
    parser.add_argument('--permutations', type=int, default=500, help='Number of permutations')
    parser.add_argument('--min-lineage', type=int, default=20, help='Minimum lineage size')
    parser.add_argument('--output', default='medalt_results', help='Output directory')
    
    args = parser.parse_args()
    
    # Run pipeline
    pipeline = InferCNVPipeline(args.output)
    tree, results = pipeline.run_infercnv_analysis(
        observations_file=args.observations,
        reference_file=args.references,
        top_k_genes=args.genes,
        subsample_size=args.subsample,
        n_permutations=args.permutations,
        min_lineage_size=args.min_lineage
    )
    
    print(f"\nAnalysis complete!")
    print(f"Results saved to: {args.output}")
    if not results.empty:
        print(f"\nTop 10 significant associations:")
        print(results[['gene_name', 'direction', 'enrichment_score', 'pvalue', 'qvalue']].head(10))


if __name__ == '__main__':
    main()