#!/usr/bin/env python
"""
MEDALT Memory-Optimized Implementation
Critical memory fixes for large-scale single-cell datasets (100k+ cells)
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Iterator
import networkx as nx
import logging
from collections import defaultdict
import psutil
import gc
import tempfile
import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import pickle
import warnings
warnings.filterwarnings('ignore')

# Optional imports with fallbacks
try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False

try:
    from scipy.sparse import csr_matrix, save_npz, load_npz
    from scipy import stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

try:
    from numba import jit, prange
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    # Fallback decorator that does nothing
    def jit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator
    
    def prange(*args, **kwargs):
        return range(*args, **kwargs)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MemoryMonitor:
    """Monitor and manage memory usage during analysis"""
    
    def __init__(self, max_memory_gb: float = 32.0):
        self.max_memory_bytes = max_memory_gb * 1024**3
        self.initial_memory = psutil.virtual_memory().used
        
    def get_memory_usage(self) -> float:
        """Get current memory usage in GB"""
        return psutil.virtual_memory().used / 1024**3
    
    def get_memory_increase(self) -> float:
        """Get memory increase since initialization in GB"""
        return (psutil.virtual_memory().used - self.initial_memory) / 1024**3
    
    def check_memory_limit(self) -> bool:
        """Check if memory usage exceeds limit"""
        return psutil.virtual_memory().used > self.max_memory_bytes
    
    def force_cleanup(self):
        """Force garbage collection and cleanup"""
        gc.collect()
        
    def log_memory_status(self, operation: str):
        """Log current memory status"""
        usage = self.get_memory_usage()
        increase = self.get_memory_increase()
        logger.info(f"{operation}: {usage:.2f}GB used (+{increase:.2f}GB)")


class SparseDistanceMatrix:
    """Memory-efficient sparse distance matrix using k-NN approach"""
    
    def __init__(self, n_cells: int, k_neighbors: int = 50,
                 temp_dir: Optional[str] = None):
        self.n_cells = n_cells
        self.k_neighbors = min(k_neighbors, n_cells - 1)
        self.temp_dir = temp_dir or tempfile.mkdtemp()
        
        if HAS_H5PY:
            # Use HDF5 for memory-mapped storage
            self.storage_file = os.path.join(self.temp_dir, 'distances.h5')
            self.h5_file = h5py.File(self.storage_file, 'w')
            
            # Create datasets for sparse storage
            self.indices = self.h5_file.create_dataset(
                'indices', (n_cells, self.k_neighbors), dtype=np.int32
            )
            self.distances = self.h5_file.create_dataset(
                'distances', (n_cells, self.k_neighbors), dtype=np.float32
            )
            self.use_h5py = True
        else:
            # Fallback to numpy memmap
            self.storage_file_idx = os.path.join(self.temp_dir, 'indices.npy')
            self.storage_file_dist = os.path.join(self.temp_dir, 'distances.npy')
            
            self.indices = np.memmap(self.storage_file_idx, dtype=np.int32, mode='w+',
                                   shape=(n_cells, self.k_neighbors))
            self.distances = np.memmap(self.storage_file_dist, dtype=np.float32, mode='w+',
                                     shape=(n_cells, self.k_neighbors))
            self.use_h5py = False
        
        logger.info(f"Initialized sparse distance matrix: {n_cells} cells, "
                   f"{self.k_neighbors} neighbors each")
    
    def set_neighbors(self, cell_id: int, neighbor_indices: np.ndarray, 
                     neighbor_distances: np.ndarray):
        """Set k-nearest neighbors for a cell"""
        k_actual = min(len(neighbor_indices), self.k_neighbors)
        
        # Sort by distance and take top k
        sorted_idx = np.argsort(neighbor_distances)[:k_actual]
        
        self.indices[cell_id, :k_actual] = neighbor_indices[sorted_idx]
        self.distances[cell_id, :k_actual] = neighbor_distances[sorted_idx]
        
        # Fill remaining with -1 (invalid)
        if k_actual < self.k_neighbors:
            self.indices[cell_id, k_actual:] = -1
            self.distances[cell_id, k_actual:] = np.inf
    
    def get_neighbors(self, cell_id: int) -> Tuple[np.ndarray, np.ndarray]:
        """Get k-nearest neighbors for a cell"""
        indices = self.indices[cell_id]
        distances = self.distances[cell_id]
        
        # Filter out invalid entries
        valid_mask = indices >= 0
        return indices[valid_mask], distances[valid_mask]
    
    def get_distance(self, cell_i: int, cell_j: int) -> float:
        """Get distance between two cells (if available)"""
        neighbors, distances = self.get_neighbors(cell_i)
        
        # Check if cell_j is in neighbors of cell_i
        idx = np.where(neighbors == cell_j)[0]
        if len(idx) > 0:
            return float(distances[idx[0]])
        
        # Check reverse direction
        neighbors, distances = self.get_neighbors(cell_j)
        idx = np.where(neighbors == cell_i)[0]
        if len(idx) > 0:
            return float(distances[idx[0]])
        
        return float('inf')
    
    def cleanup(self):
        """Clean up temporary files"""
        if hasattr(self, 'h5_file'):
            self.h5_file.close()
        if hasattr(self, 'storage_file') and os.path.exists(self.storage_file):
            os.remove(self.storage_file)
        if hasattr(self, 'storage_file_idx') and os.path.exists(self.storage_file_idx):
            os.remove(self.storage_file_idx)
        if hasattr(self, 'storage_file_dist') and os.path.exists(self.storage_file_dist):
            os.remove(self.storage_file_dist)
        if hasattr(self, 'temp_dir') and os.path.exists(self.temp_dir):
            try:
                os.rmdir(self.temp_dir)
            except OSError:
                pass  # Directory not empty


@jit(nopython=True, parallel=False)
def compute_med_vectorized(profile_a: np.ndarray, profile_b: np.ndarray) -> int:
    """
    Vectorized MED computation with proper bounds checking
    """
    n = len(profile_a)
    if n != len(profile_b):
        return 999999  # Invalid
    
    # Check homozygous deletion constraint
    for i in range(n):
        if profile_a[i] == 0 and profile_b[i] > 0:
            return 999999
    
    distance = 0
    i = 0
    
    while i < n:
        if profile_a[i] == profile_b[i]:
            i += 1
            continue
        
        # Find segment with consistent direction
        j = i
        direction = 1 if profile_b[i] > profile_a[i] else -1
        
        while j < n and profile_a[j] != profile_b[j]:
            curr_dir = 1 if profile_b[j] > profile_a[j] else -1
            if curr_dir != direction:
                break
            j += 1
        
        # Calculate max difference in segment
        max_diff = 0
        for k in range(i, j):
            diff = abs(profile_b[k] - profile_a[k])
            if diff > max_diff:
                max_diff = diff
        
        distance += max_diff
        i = j
    
    return min(distance, 999999)  # Prevent overflow


class ChunkedDistanceComputer:
    """Compute distances in chunks to manage memory"""
    
    def __init__(self, cnv_matrix: np.ndarray, chunk_size: int = 1000,
                 k_neighbors: int = 50, n_jobs: int = 4):
        self.cnv_matrix = cnv_matrix.astype(np.int8)  # Memory efficient
        self.n_cells, self.n_genes = cnv_matrix.shape
        self.chunk_size = chunk_size
        self.k_neighbors = k_neighbors
        self.n_jobs = n_jobs
        
        # Reference profile (diploid)
        self.ref_profile = np.full(self.n_genes, 2, dtype=np.int8)
        
        logger.info(f"Initialized chunked computer: {self.n_cells} cells, "
                   f"chunk size {chunk_size}, {k_neighbors} neighbors")
    
    def compute_chunk_distances(self, chunk_start: int, chunk_end: int) -> Dict:
        """Compute distances for a chunk of cells"""
        chunk_results = {}
        
        for i in range(chunk_start, chunk_end):
            if i >= self.n_cells:
                break
            
            # Compute distances to all other cells
            distances = []
            indices = []
            
            # Distance to reference
            ref_dist = compute_med_vectorized(self.ref_profile, self.cnv_matrix[i])
            distances.append(ref_dist)
            indices.append(-1)  # Special index for reference
            
            # Distances to other cells (sample if too many)
            if self.n_cells > 10000:
                # For very large datasets, sample cells
                sample_size = min(5000, self.n_cells - 1)
                other_cells = np.random.choice(
                    [j for j in range(self.n_cells) if j != i],
                    size=sample_size, replace=False
                )
            else:
                other_cells = [j for j in range(self.n_cells) if j != i]
            
            for j in other_cells:
                dist = compute_med_vectorized(self.cnv_matrix[i], self.cnv_matrix[j])
                distances.append(dist)
                indices.append(j)
            
            # Keep only k best neighbors
            distances = np.array(distances)
            indices = np.array(indices)
            
            # Sort and keep top k
            valid_mask = distances < 999999
            if np.any(valid_mask):
                valid_distances = distances[valid_mask]
                valid_indices = indices[valid_mask]
                
                sorted_idx = np.argsort(valid_distances)[:self.k_neighbors]
                
                chunk_results[i] = {
                    'indices': valid_indices[sorted_idx],
                    'distances': valid_distances[sorted_idx]
                }
            else:
                # No valid neighbors, connect to reference only
                chunk_results[i] = {
                    'indices': np.array([-1]),
                    'distances': np.array([ref_dist])
                }
        
        return chunk_results
    
    def compute_all_distances(self, memory_monitor: MemoryMonitor) -> SparseDistanceMatrix:
        """Compute all distances using chunked approach"""
        sparse_matrix = SparseDistanceMatrix(self.n_cells, self.k_neighbors)
        
        n_chunks = (self.n_cells + self.chunk_size - 1) // self.chunk_size
        logger.info(f"Computing distances in {n_chunks} chunks")
        
        # Process chunks in parallel
        with ProcessPoolExecutor(max_workers=self.n_jobs) as executor:
            futures = []
            
            for chunk_idx in range(n_chunks):
                chunk_start = chunk_idx * self.chunk_size
                chunk_end = min((chunk_idx + 1) * self.chunk_size, self.n_cells)
                
                future = executor.submit(
                    self.compute_chunk_distances, chunk_start, chunk_end
                )
                futures.append((chunk_start, chunk_end, future))
            
            # Collect results
            for chunk_start, chunk_end, future in futures:
                chunk_results = future.result()
                
                for cell_id, data in chunk_results.items():
                    sparse_matrix.set_neighbors(
                        cell_id, data['indices'], data['distances']
                    )
                
                # Memory cleanup after each chunk
                memory_monitor.force_cleanup()
                memory_monitor.log_memory_status(f"Processed chunk {chunk_start}-{chunk_end}")
                
                if memory_monitor.check_memory_limit():
                    logger.warning("Memory limit reached, forcing cleanup")
                    memory_monitor.force_cleanup()
        
        return sparse_matrix


class MemoryOptimizedMEDALT:
    """Memory-optimized MEDALT implementation for large datasets"""
    
    def __init__(self, cnv_matrix: np.ndarray, 
                 max_memory_gb: float = 32.0,
                 k_neighbors: int = 50,
                 chunk_size: int = 1000,
                 n_jobs: int = 4):
        """
        Initialize memory-optimized MEDALT
        
        Args:
            cnv_matrix: n_cells × n_genes matrix
            max_memory_gb: Maximum memory usage in GB
            k_neighbors: Number of nearest neighbors to store
            chunk_size: Size of processing chunks
            n_jobs: Number of parallel jobs
        """
        self.cnv_matrix = cnv_matrix.astype(np.int8)
        self.n_cells, self.n_genes = cnv_matrix.shape
        self.k_neighbors = k_neighbors
        self.chunk_size = chunk_size
        self.n_jobs = n_jobs
        
        # Initialize memory monitor
        self.memory_monitor = MemoryMonitor(max_memory_gb)
        
        # Initialize distance computer
        self.distance_computer = ChunkedDistanceComputer(
            cnv_matrix, chunk_size, k_neighbors, n_jobs
        )
        
        logger.info(f"Initialized MemoryOptimizedMEDALT: {self.n_cells} cells, "
                   f"{self.n_genes} genes, {k_neighbors} neighbors, "
                   f"max memory {max_memory_gb}GB")
    
    def build_tree(self) -> nx.DiGraph:
        """Build MEDALT tree using memory-optimized approach"""
        logger.info("Building tree with memory optimization...")
        
        # Step 1: Compute sparse distance matrix
        self.memory_monitor.log_memory_status("Starting distance computation")
        sparse_distances = self.distance_computer.compute_all_distances(self.memory_monitor)
        
        # Step 2: Build tree incrementally
        self.memory_monitor.log_memory_status("Building tree structure")
        tree = self._build_tree_from_sparse_distances(sparse_distances)
        
        # Step 3: Cleanup
        sparse_distances.cleanup()
        self.memory_monitor.force_cleanup()
        self.memory_monitor.log_memory_status("Tree construction complete")
        
        return tree
    
    def _build_tree_from_sparse_distances(self, sparse_distances: SparseDistanceMatrix) -> nx.DiGraph:
        """Build tree from sparse distance matrix using progressive approach"""
        
        # Create initial tree with root
        tree = nx.DiGraph()
        tree.add_node('root')
        
        # Progressive tree construction
        connected_cells = set()
        unconnected_cells = set(range(self.n_cells))
        
        # Start with cell closest to reference
        best_cell = None
        best_distance = float('inf')
        
        for cell_id in range(self.n_cells):
            neighbors, distances = sparse_distances.get_neighbors(cell_id)
            
            # Check if reference (-1) is in neighbors
            ref_idx = np.where(neighbors == -1)[0]
            if len(ref_idx) > 0:
                ref_dist = distances[ref_idx[0]]
                if ref_dist < best_distance:
                    best_distance = ref_dist
                    best_cell = cell_id
        
        # Connect best cell to root
        if best_cell is not None:
            tree.add_node(best_cell)
            tree.add_edge('root', best_cell, weight=best_distance)
            connected_cells.add(best_cell)
            unconnected_cells.remove(best_cell)
        
        # Progressively add remaining cells
        while unconnected_cells:
            best_edge = None
            best_weight = float('inf')
            
            # Find best connection from connected to unconnected
            for connected_cell in connected_cells:
                neighbors, distances = sparse_distances.get_neighbors(connected_cell)
                
                for neighbor_id, distance in zip(neighbors, distances):
                    if neighbor_id in unconnected_cells and distance < best_weight:
                        best_weight = distance
                        best_edge = (connected_cell, neighbor_id)
            
            # Also check connections from unconnected to connected
            for unconnected_cell in unconnected_cells:
                neighbors, distances = sparse_distances.get_neighbors(unconnected_cell)
                
                for neighbor_id, distance in zip(neighbors, distances):
                    if neighbor_id in connected_cells and distance < best_weight:
                        best_weight = distance
                        best_edge = (neighbor_id, unconnected_cell)
            
            # Add best edge
            if best_edge:
                parent, child = best_edge
                tree.add_node(child)
                tree.add_edge(parent, child, weight=best_weight)
                connected_cells.add(child)
                unconnected_cells.discard(child)
            else:
                # No valid connections found, connect remaining to root
                logger.warning(f"No valid connections found for {len(unconnected_cells)} cells")
                for cell in unconnected_cells:
                    tree.add_node(cell)
                    tree.add_edge('root', cell, weight=999999)
                break
            
            # Memory cleanup every 1000 cells
            if len(connected_cells) % 1000 == 0:
                self.memory_monitor.force_cleanup()
        
        logger.info(f"Tree built with {tree.number_of_nodes()} nodes and {tree.number_of_edges()} edges")
        return tree


class MemoryOptimizedLSA:
    """Memory-optimized Lineage Speciation Analysis"""
    
    def __init__(self, tree: nx.DiGraph, cnv_matrix: np.ndarray,
                 max_memory_gb: float = 32.0,
                 n_permutations: int = 500,
                 chunk_size: int = 1000):
        self.tree = tree
        self.cnv_matrix = cnv_matrix.astype(np.int8)
        self.n_cells, self.n_genes = cnv_matrix.shape
        self.n_permutations = n_permutations
        self.chunk_size = chunk_size
        
        self.memory_monitor = MemoryMonitor(max_memory_gb)
        
        # Precompute tree structure
        self._precompute_lineages()
        
        logger.info(f"Initialized MemoryOptimizedLSA: {n_permutations} permutations")
    
    def _precompute_lineages(self):
        """Precompute lineage information for efficiency"""
        self.lineages = {}
        
        for node in self.tree.nodes():
            if node == 'root':
                continue
            
            descendants = set(nx.descendants(self.tree, node))
            descendants.add(node)
            
            # Only keep actual cells (not root)
            cell_descendants = [n for n in descendants if n != 'root']
            
            if len(cell_descendants) >= 5:  # Minimum lineage size
                self.lineages[node] = {
                    'cells': np.array(cell_descendants, dtype=np.int32),
                    'size': len(cell_descendants)
                }
    
    def run_analysis(self, min_lineage_size: int = 20) -> pd.DataFrame:
        """Run memory-optimized LSA analysis"""
        logger.info("Starting memory-optimized LSA analysis...")
        
        # Filter lineages by size
        valid_lineages = {
            node: data for node, data in self.lineages.items()
            if data['size'] >= min_lineage_size
        }
        
        if not valid_lineages:
            logger.warning("No lineages found meeting minimum size requirement")
            return pd.DataFrame()
        
        logger.info(f"Found {len(valid_lineages)} valid lineages")
        
        # Calculate observed statistics
        self.memory_monitor.log_memory_status("Computing observed statistics")
        observed_stats = self._compute_observed_stats(valid_lineages)
        
        # Run permutation test in chunks
        self.memory_monitor.log_memory_status("Starting permutation testing")
        p_values = self._run_permutation_test_chunked(observed_stats, valid_lineages)
        
        # Format results
        results = self._format_results(observed_stats, p_values)
        
        self.memory_monitor.log_memory_status("LSA analysis complete")
        return results
    
    def _compute_observed_stats(self, lineages: Dict) -> Dict:
        """Compute observed statistics for lineages"""
        observed = {}
        
        for lineage_id, lineage_data in lineages.items():
            cells = lineage_data['cells']
            
            # Get CNV data for lineage
            lineage_cnvs = self.cnv_matrix[cells]
            
            for gene_idx in range(self.n_genes):
                gene_cnvs = lineage_cnvs[:, gene_idx]
                mean_cn = np.mean(gene_cnvs)
                
                # Skip neutral regions
                if abs(mean_cn - 2.0) < 0.3:
                    continue
                
                # Determine direction
                direction = 'AMP' if mean_cn > 2.0 else 'DEL'
                
                # Calculate enrichment vs background
                other_cells = np.setdiff1d(np.arange(self.n_cells), cells)
                if len(other_cells) > 0:
                    background_mean = np.mean(self.cnv_matrix[other_cells, gene_idx])
                    enrichment = abs(mean_cn - background_mean)
                else:
                    enrichment = abs(mean_cn - 2.0)
                
                observed[(lineage_id, gene_idx, direction)] = {
                    'mean_cn': mean_cn,
                    'enrichment': enrichment,
                    'lineage_size': len(cells)
                }
        
        return observed
    
    def _run_permutation_test_chunked(self, observed_stats: Dict, lineages: Dict) -> Dict:
        """Run permutation test in memory-efficient chunks"""
        
        # Process permutations in chunks
        chunk_size = min(50, self.n_permutations)  # Permutations per chunk
        n_chunks = (self.n_permutations + chunk_size - 1) // chunk_size
        
        all_null_stats = []
        
        for chunk_idx in range(n_chunks):
            chunk_start = chunk_idx * chunk_size
            chunk_end = min((chunk_idx + 1) * chunk_size, self.n_permutations)
            
            logger.info(f"Processing permutation chunk {chunk_idx + 1}/{n_chunks}")
            
            # Run chunk of permutations
            chunk_results = self._run_permutation_chunk(
                chunk_start, chunk_end, lineages
            )
            all_null_stats.extend(chunk_results)
            
            # Memory cleanup after each chunk
            self.memory_monitor.force_cleanup()
        
        # Calculate p-values
        p_values = {}
        for key, obs_data in observed_stats.items():
            lineage_id, gene_idx, direction = key
            obs_enrichment = obs_data['enrichment']
            lineage_size = obs_data['lineage_size']
            
            # Collect matching null statistics
            null_enrichments = []
            for null_stats in all_null_stats:
                if (lineage_size, gene_idx, direction) in null_stats:
                    null_enrichments.append(null_stats[(lineage_size, gene_idx, direction)])
            
            if null_enrichments:
                null_enrichments = np.array(null_enrichments)
                p_value = np.mean(null_enrichments >= obs_enrichment)
                p_values[key] = max(p_value, 1.0 / len(null_enrichments))
            else:
                p_values[key] = 1.0
        
        return p_values
    
    def _run_permutation_chunk(self, start_idx: int, end_idx: int, lineages: Dict) -> List[Dict]:
        """Run a chunk of permutations"""
        results = []
        
        for perm_idx in range(start_idx, end_idx):
            np.random.seed(perm_idx)  # Reproducible
            
            # Permute cell labels
            permuted_indices = np.random.permutation(self.n_cells)
            
            null_stats = {}
            
            for lineage_id, lineage_data in lineages.items():
                lineage_size = lineage_data['size']
                
                # Get random cells of same size
                perm_cells = permuted_indices[:lineage_size]
                
                # Calculate statistics
                for gene_idx in range(self.n_genes):
                    gene_cnvs = self.cnv_matrix[perm_cells, gene_idx]
                    mean_cn = np.mean(gene_cnvs)
                    
                    if abs(mean_cn - 2.0) < 0.3:
                        continue
                    
                    direction = 'AMP' if mean_cn > 2.0 else 'DEL'
                    
                    # Calculate enrichment
                    other_cells = permuted_indices[lineage_size:]
                    if len(other_cells) > 0:
                        background_mean = np.mean(self.cnv_matrix[other_cells, gene_idx])
                        enrichment = abs(mean_cn - background_mean)
                    else:
                        enrichment = abs(mean_cn - 2.0)
                    
                    null_stats[(lineage_size, gene_idx, direction)] = enrichment
            
            results.append(null_stats)
        
        return results
    
    def _format_results(self, observed_stats: Dict, p_values: Dict) -> pd.DataFrame:
        """Format results into DataFrame"""
        rows = []
        
        for key, obs_data in observed_stats.items():
            lineage_id, gene_idx, direction = key
            
            if key in p_values:
                rows.append({
                    'lineage_root': lineage_id,
                    'gene_idx': gene_idx,
                    'direction': direction,
                    'mean_cn': obs_data['mean_cn'],
                    'enrichment_score': obs_data['enrichment'],
                    'lineage_size': obs_data['lineage_size'],
                    'pvalue': p_values[key]
                })
        
        df = pd.DataFrame(rows)
        
        if not df.empty:
            # FDR correction
            from scipy import stats
            df['qvalue'] = stats.false_discovery_control(df['pvalue'].values)
            df = df.sort_values('pvalue')
        
        return df


def main():
    """Example usage of memory-optimized MEDALT"""
    
    # Create synthetic large dataset
    np.random.seed(42)
    n_cells = 5000  # Large dataset
    n_genes = 1000
    
    logger.info(f"Creating synthetic dataset: {n_cells} cells × {n_genes} genes")
    
    # Base diploid state
    cnv_matrix = np.random.poisson(2, (n_cells, n_genes)).astype(np.int8)
    
    # Add some structure
    # Lineage 1: amplification
    cnv_matrix[1000:2000, 200:300] = 4
    # Lineage 2: deletion
    cnv_matrix[2000:3000, 500:600] = 0
    
    # Run memory-optimized analysis
    logger.info("Running memory-optimized MEDALT analysis...")
    
    medalt = MemoryOptimizedMEDALT(
        cnv_matrix,
        max_memory_gb=16.0,  # 16GB limit
        k_neighbors=30,      # Reduced for speed
        chunk_size=500,      # Process in chunks
        n_jobs=4
    )
    
    # Build tree
    tree = medalt.build_tree()
    
    # Run LSA
    lsa = MemoryOptimizedLSA(
        tree, cnv_matrix,
        max_memory_gb=16.0,
        n_permutations=100,  # Reduced for demo
        chunk_size=500
    )
    
    results = lsa.run_analysis(min_lineage_size=50)
    
    logger.info(f"Analysis complete! Found {len(results)} significant events")
    
    if not results.empty:
        print("\nTop results:")
        print(results[['lineage_root', 'gene_idx', 'direction', 
                      'enrichment_score', 'pvalue', 'qvalue']].head(10))


if __name__ == '__main__':
    main()