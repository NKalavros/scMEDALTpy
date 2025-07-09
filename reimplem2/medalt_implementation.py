#!/usr/bin/env python
"""
MEDALT: Minimal Event Distance Aneuploidy Lineage Tree
Implementation of single-cell copy number lineage tracing
"""

import numpy as np
import pandas as pd
from collections import defaultdict
from scipy import stats
import networkx as nx
from typing import Dict, List, Tuple, Set
import warnings
from multiprocessing import Pool, cpu_count
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MED:
    """Minimal Event Distance calculator"""
    
    @staticmethod
    def compute(profile_a: np.ndarray, profile_b: np.ndarray) -> float:
        """
        Compute MED between two copy number profiles
        
        Args:
            profile_a, profile_b: Integer arrays representing copy number profiles
            
        Returns:
            Float representing minimal event distance (inf if impossible)
        """
        if len(profile_a) != len(profile_b):
            raise ValueError("Profiles must have same length")
        
        # Check homozygous deletion constraint
        for i in range(len(profile_a)):
            if profile_a[i] == 0 and profile_b[i] > 0:
                return float('inf')
        
        distance = 0
        i = 0
        n = len(profile_a)
        
        while i < n:
            if profile_a[i] == profile_b[i]:
                i += 1
                continue
            
            # Find maximal segment with consistent direction
            j = i
            direction = 1 if profile_b[i] > profile_a[i] else -1
            
            while j < n and profile_a[j] != profile_b[j]:
                curr_direction = 1 if profile_b[j] > profile_a[j] else -1
                if curr_direction != direction:
                    break
                j += 1
            
            # Calculate events needed for segment [i, j)
            max_diff = 0
            for k in range(i, j):
                max_diff = max(max_diff, abs(profile_b[k] - profile_a[k]))
            
            distance += max_diff
            i = j
        
        return distance


class MEDALT:
    """MEDALT tree construction using Edmonds' algorithm"""
    
    def __init__(self, cnv_matrix: np.ndarray, chromosome_boundaries: List[Tuple[int, int]] = None):
        """
        Initialize MEDALT tree builder
        
        Args:
            cnv_matrix: n_cells × n_bins integer matrix
            chromosome_boundaries: List of (start, end) tuples for chromosome regions
        """
        self.cnv_matrix = cnv_matrix.astype(int)
        self.n_cells, self.n_bins = cnv_matrix.shape
        self.chromosome_boundaries = chromosome_boundaries
        
        # Add virtual root (normal diploid)
        self.root_profile = np.full(self.n_bins, 2, dtype=int)
        self.root_id = 'root'
        
        # Cache for distance calculations
        self._distance_cache = {}
        
        logger.info(f"Initialized MEDALT with {self.n_cells} cells and {self.n_bins} bins")
    
    def build_tree(self) -> nx.DiGraph:
        """
        Construct RDMST using adapted Edmonds' algorithm
        
        Returns:
            NetworkX directed graph representing the lineage tree
        """
        logger.info("Computing distance matrix...")
        distances = self._compute_distance_matrix()
        
        logger.info("Building directed graph...")
        # Create directed graph with all edges
        G = nx.DiGraph()
        
        # Add nodes
        G.add_node(self.root_id, profile=self.root_profile)
        for i in range(self.n_cells):
            G.add_node(i, profile=self.cnv_matrix[i])
        
        # Add edges from root to all cells
        for i in range(self.n_cells):
            dist = distances[self.root_id][i]
            if dist < float('inf'):
                G.add_edge(self.root_id, i, weight=dist)
        
        # Add edges between cells
        for i in range(self.n_cells):
            for j in range(self.n_cells):
                if i != j and distances[i][j] < float('inf'):
                    G.add_edge(i, j, weight=distances[i][j])
        
        logger.info("Finding minimum spanning arborescence...")
        # Use NetworkX's implementation of Edmonds' algorithm
        try:
            mst = nx.minimum_spanning_arborescence(G, attr='weight')
        except nx.NetworkXException as e:
            logger.warning(f"Failed to find spanning arborescence: {e}")
            # Fall back to simple greedy approach
            mst = self._greedy_tree_construction(distances)
        
        logger.info(f"Tree constructed with {mst.number_of_nodes()} nodes and {mst.number_of_edges()} edges")
        return mst
    
    def _compute_distance_matrix(self) -> Dict:
        """Compute all pairwise MEDs (asymmetric)"""
        distances = defaultdict(dict)
        
        # Distances from root
        for i in range(self.n_cells):
            key = (self.root_id, i)
            if key not in self._distance_cache:
                self._distance_cache[key] = MED.compute(self.root_profile, self.cnv_matrix[i])
            distances[self.root_id][i] = self._distance_cache[key]
        
        # Pairwise distances between cells (asymmetric)
        for i in range(self.n_cells):
            for j in range(self.n_cells):
                if i != j:
                    key = (i, j)  # Keep directional - MED(i,j) != MED(j,i)
                    if key not in self._distance_cache:
                        self._distance_cache[key] = MED.compute(self.cnv_matrix[i], self.cnv_matrix[j])
                    distances[i][j] = self._distance_cache[key]
        
        return distances
    
    def _greedy_tree_construction(self, distances: Dict) -> nx.DiGraph:
        """Fallback greedy tree construction"""
        tree = nx.DiGraph()
        tree.add_node(self.root_id)
        
        # Add all cells
        for i in range(self.n_cells):
            tree.add_node(i)
        
        # Greedily connect each cell to its nearest ancestor
        connected = {self.root_id}
        unconnected = set(range(self.n_cells))
        
        while unconnected:
            best_edge = None
            best_dist = float('inf')
            
            for parent in connected:
                for child in unconnected:
                    if parent == self.root_id:
                        dist = distances[self.root_id][child]
                    else:
                        dist = distances[parent][child]
                    
                    if dist < best_dist:
                        best_dist = dist
                        best_edge = (parent, child)
            
            if best_edge:
                tree.add_edge(best_edge[0], best_edge[1], weight=best_dist)
                connected.add(best_edge[1])
                unconnected.remove(best_edge[1])
            else:
                # Connect remaining nodes directly to root
                for node in unconnected:
                    tree.add_edge(self.root_id, node, weight=float('inf'))
                break
        
        return tree


class LineageSpeciationAnalysis:
    """Identify CNAs associated with lineage expansion"""
    
    def __init__(self, tree: nx.DiGraph, cnv_matrix: np.ndarray, 
                 chromosome_boundaries: List[Tuple[int, int]] = None,
                 n_permutations: int = 500, data_type: str = 'RNA'):
        """
        Initialize LSA
        
        Args:
            tree: NetworkX directed graph from MEDALT
            cnv_matrix: n_cells × n_bins integer matrix
            chromosome_boundaries: List of (start, end) tuples for chromosomes
            n_permutations: Number of permutations for background distribution
            data_type: 'RNA' for scRNA-seq, 'DNA' for scDNA-seq
        """
        self.tree = tree
        self.cnv_matrix = cnv_matrix.astype(int)
        self.n_cells, self.n_bins = cnv_matrix.shape
        self.n_permutations = n_permutations
        self.data_type = data_type
        
        # Set up chromosome boundaries
        if chromosome_boundaries is None:
            # Assume single chromosome if not provided
            self.chromosome_boundaries = [(0, self.n_bins)]
        else:
            self.chromosome_boundaries = chromosome_boundaries
        
        logger.info(f"Initialized LSA with {n_permutations} permutations, data_type={data_type}")
    
    def run_analysis(self, min_lineage_size: int = 5, n_jobs: int = -1, gene_names: List[str] = None) -> pd.DataFrame:
        """
        Run complete LSA pipeline
        
        Args:
            min_lineage_size: Minimum size for lineage consideration
            n_jobs: Number of parallel jobs (-1 for all CPUs)
            
        Returns:
            DataFrame with LSA results
        """
        # Step 1: Dissect tree into lineages
        logger.info("Dissecting tree into lineages...")
        initial_lineages = self._dissect_tree(min_lineage_size)
        logger.info(f"Found {len(initial_lineages)} initial lineages")
        
        # Step 1.5: Refine lineages to remove nested/redundant ones
        lineages = self._refine_lineages(initial_lineages)
        logger.info(f"Using {len(lineages)} refined lineages")
        
        if not lineages:
            logger.warning("No lineages found!")
            return pd.DataFrame()
        
        # Step 2: Calculate observed CFLs
        logger.info("Calculating observed CFLs...")
        observed_cfls = self._calculate_cfls(lineages)
        
        # Step 3: Generate background distribution
        logger.info(f"Generating background distribution with {self.n_permutations} permutations...")
        background_cfls = self._generate_background_parallel(lineages, n_jobs)
        
        # Step 4: Calculate p-values
        logger.info("Calculating p-values...")
        results = self._calculate_pvalues(observed_cfls, background_cfls, lineages)
        
        # Step 5: Format output
        return self._format_results(results, gene_names)
    
    def _dissect_tree(self, min_size: int = 5) -> List[Dict]:
        """Partition tree into disjoint lineages"""
        lineages = []
        
        # Calculate subtree sizes
        subtree_sizes = {}
        for node in self.tree.nodes():
            descendants = nx.descendants(self.tree, node)
            descendants.add(node)  # Include the node itself
            # Only count non-root nodes (actual cells)
            cell_descendants = [n for n in descendants if n != 'root']
            subtree_sizes[node] = len(cell_descendants)
        
        # Select candidate nodes (excluding root)
        candidates = [(node, size) for node, size in subtree_sizes.items() 
                      if node != 'root' and size >= min_size]
        candidates.sort(key=lambda x: x[1])  # Sort by size
        
        # Remove redundant lineages
        selected = []
        for node, size in candidates:
            is_redundant = False
            
            # Check if this node is ancestor or descendant of already selected
            for selected_node in selected:
                if nx.has_path(self.tree, node, selected_node) or \
                   nx.has_path(self.tree, selected_node, node):
                    is_redundant = True
                    break
            
            if not is_redundant:
                selected.append(node)
                descendants = nx.descendants(self.tree, node)
                descendants.add(node)
                cell_descendants = [n for n in descendants if n != 'root']
                
                lineages.append({
                    'root': node,
                    'cells': cell_descendants,
                    'size': len(cell_descendants),
                    'depth': nx.shortest_path_length(self.tree, 'root', node)
                })
        
        return lineages
    
    def _calculate_cfls(self, lineages: List[Dict]) -> Dict:
        """Calculate Cumulative Fold Level for each CNA in each lineage"""
        cfls = {}
        
        for lineage in lineages:
            lineage_id = lineage['root']
            cells = lineage['cells']
            
            # Calculate CFLs for amplifications and deletions separately
            for bin_idx in range(self.n_bins):
                # Sum copy numbers across cells in lineage
                cn_values = [self.cnv_matrix[cell][bin_idx] for cell in cells]
                cfl = sum(cn_values)
                avg_cn = np.mean(cn_values)
                
                # Determine if amplification or deletion
                if avg_cn > 2:  # Amplification
                    direction = 'AMP'
                elif avg_cn < 2:  # Deletion
                    direction = 'DEL'
                else:
                    continue  # Skip neutral
                
                cfls[(lineage_id, bin_idx, direction)] = {
                    'cfl': cfl,
                    'avg_cn': avg_cn,
                    'cells': cells,
                    'size': len(cells)
                }
        
        return cfls
    
    def _generate_background_parallel(self, lineages: List[Dict], n_jobs: int = -1) -> Dict:
        """Generate background distribution using parallel processing"""
        if n_jobs == -1:
            n_jobs = cpu_count()
        
        # Prepare arguments for parallel processing
        perm_args = [(i, lineages) for i in range(self.n_permutations)]
        
        # Run permutations in parallel
        with Pool(n_jobs) as pool:
            all_results = pool.map(self._single_permutation, perm_args)
        
        # Combine results
        background = defaultdict(list)
        for perm_result in all_results:
            for key, values in perm_result.items():
                background[key].extend(values)
        
        return dict(background)
    
    def _single_permutation(self, args: Tuple) -> Dict:
        """Run single permutation (for parallel processing)"""
        perm_idx, original_lineages = args
        np.random.seed(perm_idx)  # For reproducibility
        
        # Permute matrix
        permuted_matrix = self._permute_by_chromosome()
        
        # Build tree with permuted data
        medalt = MEDALT(permuted_matrix, self.chromosome_boundaries)
        perm_tree = medalt.build_tree()
        
        # Create new LSA instance for permuted data
        perm_lsa = LineageSpeciationAnalysis(perm_tree, permuted_matrix, 
                                             self.chromosome_boundaries, 0, self.data_type)
        
        # Dissect permuted tree
        min_size = min(lin['size'] for lin in original_lineages)
        perm_initial_lineages = perm_lsa._dissect_tree(min_size)
        perm_lineages = perm_lsa._refine_lineages(perm_initial_lineages)
        
        # Calculate CFLs for permuted data
        background_cfls = defaultdict(list)
        
        for lineage in perm_lineages:
            size = lineage['size']
            cells = lineage['cells']
            
            for bin_idx in range(self.n_bins):
                cn_values = [permuted_matrix[cell][bin_idx] for cell in cells]
                cfl = sum(cn_values)
                avg_cn = np.mean(cn_values)
                
                if avg_cn > 2:
                    direction = 'AMP'
                elif avg_cn < 2:
                    direction = 'DEL'
                else:
                    continue
                
                background_cfls[(size, bin_idx, direction)].append(cfl)
        
        return dict(background_cfls)
    
    def _permute_by_chromosome(self) -> np.ndarray:
        """Permute CNV matrix by chromosome with data-type specific strategy"""
        permuted = self.cnv_matrix.copy()
        
        if self.data_type == 'DNA':
            # For DNA-seq: permute chromosomal segments across cells
            for chr_start, chr_end in self.chromosome_boundaries:
                # Shuffle cell assignments within chromosome
                shuffled_indices = np.random.permutation(self.n_cells)
                permuted[:, chr_start:chr_end] = self.cnv_matrix[shuffled_indices][:, chr_start:chr_end]
        
        elif self.data_type == 'RNA':
            # For RNA-seq: permute genes within same chromosome across cells
            for chr_start, chr_end in self.chromosome_boundaries:
                # Permute each gene position independently within chromosome
                for gene_pos in range(chr_start, chr_end):
                    shuffled_indices = np.random.permutation(self.n_cells)
                    permuted[:, gene_pos] = self.cnv_matrix[shuffled_indices, gene_pos]
        
        return permuted
    
    def _calculate_pvalues(self, observed_cfls: Dict, background_cfls: Dict, 
                          lineages: List[Dict]) -> Dict:
        """Calculate empirical p-values"""
        results = {}
        
        # Create lineage lookup
        lineage_lookup = {lin['root']: lin for lin in lineages}
        
        for (lineage_id, bin_idx, direction), obs_data in observed_cfls.items():
            observed_cfl = obs_data['cfl']
            lineage_size = obs_data['size']
            
            # Find matching background distribution
            background_key = (lineage_size, bin_idx, direction)
            
            if background_key in background_cfls:
                background = background_cfls[background_key]
            else:
                # Try nearest size
                background = self._find_nearest_background(lineage_size, bin_idx, 
                                                         direction, background_cfls)
            
            if background:
                # Calculate empirical p-value
                n_greater_equal = sum(1 for b in background if b >= observed_cfl)
                p_value = (n_greater_equal + 1) / (len(background) + 1)  # Add pseudocount
            else:
                p_value = 1.0
            
            results[(lineage_id, bin_idx, direction)] = {
                'cfl': observed_cfl,
                'avg_cn': obs_data['avg_cn'],
                'size': lineage_size,
                'pvalue': p_value,
                'n_background': len(background) if background else 0,
                'depth': lineage_lookup[lineage_id]['depth']
            }
        
        return results
    
    def _find_nearest_background(self, target_size: int, bin_idx: int, 
                                direction: str, background_cfls: Dict) -> List:
        """Find background distribution for nearest lineage size"""
        # Get all available sizes for this bin and direction
        available_sizes = [size for (size, b, d) in background_cfls.keys() 
                          if b == bin_idx and d == direction]
        
        if not available_sizes:
            return []
        
        # Find nearest size
        nearest_size = min(available_sizes, key=lambda x: abs(x - target_size))
        return background_cfls[(nearest_size, bin_idx, direction)]
    
    def _benjamini_hochberg_correction(self, pvalues: np.ndarray) -> np.ndarray:
        """Manual implementation of Benjamini-Hochberg FDR correction"""
        pvalues = np.array(pvalues)
        n = len(pvalues)
        
        # Sort p-values and keep track of original indices
        sorted_indices = np.argsort(pvalues)
        sorted_pvals = pvalues[sorted_indices]
        
        # Apply BH correction
        adjusted_pvals = np.zeros_like(sorted_pvals)
        for i in range(n-1, -1, -1):
            if i == n-1:
                adjusted_pvals[i] = sorted_pvals[i]
            else:
                adjusted_pvals[i] = min(adjusted_pvals[i+1], 
                                      sorted_pvals[i] * n / (i + 1))
        
        # Restore original order
        result = np.zeros_like(pvalues)
        result[sorted_indices] = adjusted_pvals
        
        return result
    
    def _refine_lineages(self, lineages: List[Dict]) -> List[Dict]:
        """Remove nested/redundant lineages similar to RefineCNA() in original R code"""
        refined_lineages = []
        
        # Sort lineages by size (largest first)
        sorted_lineages = sorted(lineages, key=lambda x: x['size'], reverse=True)
        
        for lineage in sorted_lineages:
            is_nested = False
            lineage_cells = set(lineage['cells'])
            
            # Check if this lineage is nested within any already selected lineage
            for selected in refined_lineages:
                selected_cells = set(selected['cells'])
                
                # If this lineage is a subset of a selected lineage, it's nested
                if lineage_cells.issubset(selected_cells):
                    is_nested = True
                    break
                
                # If overlap is too large, consider it redundant
                overlap = len(lineage_cells.intersection(selected_cells))
                if overlap > 0.7 * min(len(lineage_cells), len(selected_cells)):
                    is_nested = True
                    break
            
            if not is_nested:
                refined_lineages.append(lineage)
        
        logger.info(f"Refined {len(lineages)} lineages to {len(refined_lineages)} non-redundant lineages")
        return refined_lineages
    
    def _format_results(self, results: Dict, gene_names: List[str] = None) -> pd.DataFrame:
        """Format results as DataFrame"""
        rows = []
        
        for (lineage_id, bin_idx, direction), stats in results.items():
            # Get gene name
            if gene_names and bin_idx < len(gene_names):
                region_name = gene_names[bin_idx]
            else:
                region_name = f'bin_{bin_idx}'
            
            # Calculate score as normalized CFL
            score = (stats['cfl'] / stats['size']) - 2  # Subtract diploid baseline
            
            rows.append({
                'region': region_name,
                'Score': score,
                'pvalue': stats['pvalue'],
                'cell': f'{lineage_id}' if isinstance(lineage_id, str) else f'cell_{lineage_id}',
                'depth': stats['depth'],
                'subtreesize': stats['size'],
                'CNA': direction
            })
        
        df = pd.DataFrame(rows)
        
        if not df.empty:
            # Add Benjamini-Hochberg FDR correction
            try:
                from scipy.stats import false_discovery_control
                df['adjustp'] = false_discovery_control(df['pvalue'], method='bh')
            except ImportError:
                # Fallback to statsmodels if scipy version is too old
                try:
                    from statsmodels.stats.multitest import fdrcorrection
                    _, df['adjustp'] = fdrcorrection(df['pvalue'], method='indep')
                except ImportError:
                    # Manual Benjamini-Hochberg implementation
                    df['adjustp'] = self._benjamini_hochberg_correction(df['pvalue'])
            
            # Sort by adjusted p-value
            df = df.sort_values('adjustp')
        
        return df


def load_scRNA_data(filepath: str, window_size: int = 30) -> Tuple[np.ndarray, List[str], List[str]]:
    """
    Load scRNA-seq CNV data from inferCNV output
    
    Args:
        filepath: Path to CNV file
        window_size: Number of genes to average per window
        
    Returns:
        cnv_matrix: Integer copy number matrix (cells x genes)
        cell_names: List of cell names
        gene_names: List of gene names
    """
    logger.info(f"Loading data from {filepath}")
    
    # Read the data
    df = pd.read_csv(filepath, sep='\t', index_col=0)
    
    # Data is already in correct format: genes as rows, cells as columns
    # Convert from relative to absolute copy numbers
    # inferCNV outputs: 0.5 = CN1, 1.0 = CN2, 1.5 = CN3, etc.
    cn_matrix = np.round(df.values * 2).astype(int)
    
    # Transpose so rows are cells, columns are genes
    cn_matrix = cn_matrix.T
    
    n_cells, n_genes = cn_matrix.shape
    logger.info(f"Loaded {n_cells} cells and {n_genes} genes")
    
    cell_names = df.columns.tolist()
    gene_names = df.index.tolist()
    
    # If windowing is requested, apply it
    if window_size > 1 and window_size < n_genes:
        n_windows = n_genes // window_size
        windowed_matrix = np.zeros((n_cells, n_windows), dtype=int)
        windowed_genes = []
        
        for i in range(n_windows):
            start = i * window_size
            end = min((i + 1) * window_size, n_genes)
            # Average and round
            windowed_matrix[:, i] = np.round(np.mean(cn_matrix[:, start:end], axis=1))
            # Use first gene name in window as representative
            windowed_genes.append(gene_names[start])
        
        logger.info(f"Reduced to {n_windows} windows using window size {window_size}")
        return windowed_matrix, cell_names, windowed_genes
    
    return cn_matrix, cell_names, gene_names


def main():
    """Example usage with scRNA data"""
    import sys
    
    # Load example data
    filepath = 'scRNA.CNV.txt'
    cnv_matrix, cell_names, gene_names = load_scRNA_data(filepath, window_size=1)
    
    # Build MEDALT tree
    logger.info("Building MEDALT tree...")
    medalt = MEDALT(cnv_matrix)
    tree = medalt.build_tree()
    
    # Save tree
    with open('CNV.tree.txt', 'w') as f:
        f.write("from\tto\tdist\n")
        for parent, child, data in tree.edges(data=True):
            parent_name = parent if parent == 'root' else cell_names[parent]
            child_name = cell_names[child]
            f.write(f"{parent_name}\t{child_name}\t{int(data['weight'])}\n")
    
    logger.info("Tree saved to CNV.tree.txt")
    
    # Run LSA
    logger.info("Running Lineage Speciation Analysis...")
    lsa = LineageSpeciationAnalysis(tree, cnv_matrix, n_permutations=500)
    results = lsa.run_analysis(min_lineage_size=5, n_jobs=-1, gene_names=gene_names)
    
    # Save results
    if not results.empty:
        results.to_csv('gene.LSA.txt', sep='\t', index=False)
        logger.info(f"LSA results saved to gene.LSA.txt")
        logger.info(f"Found {len(results)} significant associations")
        print("\nTop 10 results:")
        print(results.head(10))
    else:
        logger.warning("No significant associations found")


if __name__ == '__main__':
    main()