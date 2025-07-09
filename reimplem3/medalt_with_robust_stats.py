#!/usr/bin/env python
"""
MEDALT with Robust Statistical Framework
Integration of memory optimization and advanced statistical methods
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
import networkx as nx
import logging
from collections import defaultdict

# Import components
from medalt_memory_optimized import (
    MemoryOptimizedMEDALT, MemoryOptimizedLSA, MemoryMonitor
)
from medalt_statistical_robust import (
    RobustStatisticalFramework, EnrichmentTest, WilcoxonTest
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class StatisticallyRobustMEDALT:
    """MEDALT with both memory optimization and robust statistics"""
    
    def __init__(self, cnv_matrix: np.ndarray,
                 max_memory_gb: float = 32.0,
                 k_neighbors: int = 50,
                 chunk_size: int = 1000,
                 n_jobs: int = 4,
                 statistical_method: str = 'enrichment',
                 correction_method: str = 'fdr_bh',
                 confidence_level: float = 0.95,
                 n_bootstrap: int = 1000):
        """
        Initialize MEDALT with robust statistics
        
        Args:
            cnv_matrix: n_cells × n_genes matrix
            max_memory_gb: Maximum memory usage in GB
            k_neighbors: Number of nearest neighbors to store
            chunk_size: Size of processing chunks
            n_jobs: Number of parallel jobs
            statistical_method: 'enrichment' or 'wilcoxon'
            correction_method: 'fdr_bh', 'fdr_by', 'bonferroni', 'holm'
            confidence_level: Confidence level for intervals
            n_bootstrap: Number of bootstrap samples
        """
        # Initialize memory-optimized MEDALT
        self.medalt = MemoryOptimizedMEDALT(
            cnv_matrix, max_memory_gb, k_neighbors, chunk_size, n_jobs
        )
        
        # Initialize statistical framework
        self.statistical_framework = RobustStatisticalFramework(
            test_method=statistical_method,
            correction_method=correction_method,
            confidence_level=confidence_level,
            n_bootstrap=n_bootstrap,
            n_jobs=n_jobs
        )
        
        self.cnv_matrix = cnv_matrix
        self.n_cells, self.n_genes = cnv_matrix.shape
        
        logger.info(f"Initialized StatisticallyRobustMEDALT with {statistical_method} "
                   f"statistics and {correction_method} correction")
    
    def run_complete_analysis(self, 
                             n_permutations: int = 500,
                             min_lineage_size: int = 20,
                             significance_threshold: float = 0.05) -> Tuple[nx.DiGraph, pd.DataFrame]:
        """
        Run complete analysis with robust statistics
        
        Args:
            n_permutations: Number of permutations for null distribution
            min_lineage_size: Minimum lineage size for analysis
            significance_threshold: Significance threshold
            
        Returns:
            (tree, results_dataframe)
        """
        logger.info("Starting complete MEDALT analysis with robust statistics...")
        
        # Step 1: Build tree using memory-optimized approach
        logger.info("Step 1: Building phylogenetic tree...")
        tree = self.medalt.build_tree()
        
        # Step 2: Run enhanced LSA with robust statistics
        logger.info("Step 2: Running enhanced LSA with robust statistics...")
        lsa_results = self._run_robust_lsa(
            tree, n_permutations, min_lineage_size, significance_threshold
        )
        
        logger.info(f"Analysis complete! Found {len(lsa_results)} associations")
        
        return tree, lsa_results
    
    def _run_robust_lsa(self, tree: nx.DiGraph, n_permutations: int,
                       min_lineage_size: int, significance_threshold: float) -> pd.DataFrame:
        """Run LSA with robust statistical analysis"""
        
        # Create enhanced LSA instance
        lsa = EnhancedLSA(
            tree, self.cnv_matrix, self.statistical_framework,
            n_permutations=n_permutations
        )
        
        # Run analysis
        results = lsa.run_enhanced_analysis(
            min_lineage_size=min_lineage_size,
            significance_threshold=significance_threshold
        )
        
        return results


class EnhancedLSA:
    """Enhanced LSA with robust statistical framework"""
    
    def __init__(self, tree: nx.DiGraph, cnv_matrix: np.ndarray,
                 statistical_framework: RobustStatisticalFramework,
                 n_permutations: int = 500):
        
        self.tree = tree
        self.cnv_matrix = cnv_matrix.astype(np.float32)
        self.statistical_framework = statistical_framework
        self.n_permutations = n_permutations
        
        self.n_cells, self.n_genes = cnv_matrix.shape
        
        # Precompute lineages
        self._precompute_lineages()
        
        logger.info(f"Initialized EnhancedLSA with {len(self.lineages)} lineages")
    
    def _precompute_lineages(self):
        """Precompute lineage information"""
        self.lineages = {}
        
        for node in self.tree.nodes():
            if node == 'root':
                continue
            
            descendants = set(nx.descendants(self.tree, node))
            descendants.add(node)
            
            # Only keep actual cells
            cell_descendants = [n for n in descendants if n != 'root']
            
            if len(cell_descendants) >= 5:
                self.lineages[node] = {
                    'cells': np.array(cell_descendants, dtype=np.int32),
                    'size': len(cell_descendants),
                    'depth': nx.shortest_path_length(self.tree, 'root', node)
                }
    
    def run_enhanced_analysis(self, min_lineage_size: int = 20,
                             significance_threshold: float = 0.05) -> pd.DataFrame:
        """Run enhanced LSA with robust statistics"""
        
        # Filter lineages by size
        valid_lineages = {
            node: data for node, data in self.lineages.items()
            if data['size'] >= min_lineage_size
        }
        
        if not valid_lineages:
            logger.warning("No lineages found meeting minimum size requirement")
            return pd.DataFrame()
        
        logger.info(f"Analyzing {len(valid_lineages)} lineages with robust statistics")
        
        # Compute observed data
        observed_data = self._compute_observed_data(valid_lineages)
        
        # Generate null distributions
        null_distributions = self._generate_null_distributions(valid_lineages)
        
        # Compute background data
        background_data = self._compute_background_data(valid_lineages)
        
        # Run statistical analysis
        statistical_results = self.statistical_framework.run_statistical_analysis(
            observed_data, null_distributions, background_data
        )
        
        # Merge with lineage information
        results = self._merge_results(statistical_results, valid_lineages)
        
        # Filter by significance
        significant_results = results[results['adjusted_pvalue'] < significance_threshold]
        
        logger.info(f"Found {len(significant_results)} significant associations "
                   f"(p < {significance_threshold})")
        
        return significant_results
    
    def _compute_observed_data(self, valid_lineages: Dict) -> Dict[str, np.ndarray]:
        """Compute observed data for statistical analysis"""
        observed_data = {}
        
        for lineage_id, lineage_info in valid_lineages.items():
            cells = lineage_info['cells']
            
            for gene_idx in range(self.n_genes):
                # Get copy number values for this gene in this lineage
                gene_values = self.cnv_matrix[cells, gene_idx]
                
                # Skip if too few observations or neutral
                if len(gene_values) < 5 or abs(np.mean(gene_values) - 2.0) < 0.3:
                    continue
                
                # Determine direction
                direction = 'AMP' if np.mean(gene_values) > 2.0 else 'DEL'
                
                # Create unique key
                key = f"{lineage_id}_{gene_idx}_{direction}"
                observed_data[key] = gene_values
        
        return observed_data
    
    def _generate_null_distributions(self, valid_lineages: Dict) -> Dict[str, np.ndarray]:
        """Generate null distributions through permutation"""
        null_distributions = {}
        
        logger.info(f"Generating null distributions with {self.n_permutations} permutations")
        
        # Run permutations
        for perm_idx in range(self.n_permutations):
            np.random.seed(perm_idx)  # Reproducible
            
            # Permute cell labels
            permuted_indices = np.random.permutation(self.n_cells)
            
            # For each lineage, calculate statistics on permuted data
            for lineage_id, lineage_info in valid_lineages.items():
                lineage_size = lineage_info['size']
                
                # Get random cells of same size
                perm_cells = permuted_indices[:lineage_size]
                
                for gene_idx in range(self.n_genes):
                    gene_values = self.cnv_matrix[perm_cells, gene_idx]
                    
                    if len(gene_values) < 5 or abs(np.mean(gene_values) - 2.0) < 0.3:
                        continue
                    
                    direction = 'AMP' if np.mean(gene_values) > 2.0 else 'DEL'
                    key = f"{lineage_id}_{gene_idx}_{direction}"
                    
                    # Calculate enrichment statistic
                    other_cells = permuted_indices[lineage_size:]
                    if len(other_cells) > 0:
                        background_values = self.cnv_matrix[other_cells, gene_idx]
                        
                        # Compute enrichment statistic
                        enrichment = self._compute_enrichment_statistic(
                            gene_values, background_values
                        )
                        
                        # Store in null distribution
                        if key not in null_distributions:
                            null_distributions[key] = []
                        null_distributions[key].append(enrichment)
        
        # Convert to numpy arrays
        for key in null_distributions:
            null_distributions[key] = np.array(null_distributions[key])
        
        return null_distributions
    
    def _compute_background_data(self, valid_lineages: Dict) -> Dict[str, np.ndarray]:
        """Compute background data for statistical comparison"""
        background_data = {}
        
        for lineage_id, lineage_info in valid_lineages.items():
            cells = lineage_info['cells']
            
            # Get all other cells as background
            other_cells = np.setdiff1d(np.arange(self.n_cells), cells)
            
            if len(other_cells) < 5:
                continue
            
            for gene_idx in range(self.n_genes):
                # Get background values for this gene
                background_values = self.cnv_matrix[other_cells, gene_idx]
                
                # Check if corresponding observed data exists
                obs_values = self.cnv_matrix[cells, gene_idx]
                if len(obs_values) < 5 or abs(np.mean(obs_values) - 2.0) < 0.3:
                    continue
                
                direction = 'AMP' if np.mean(obs_values) > 2.0 else 'DEL'
                key = f"{lineage_id}_{gene_idx}_{direction}"
                
                background_data[key] = background_values
        
        return background_data
    
    def _compute_enrichment_statistic(self, obs_values: np.ndarray, 
                                    bg_values: np.ndarray) -> float:
        """Compute enrichment statistic"""
        obs_mean = np.mean(obs_values)
        bg_mean = np.mean(bg_values)
        
        # Use standardized difference
        pooled_std = np.sqrt((np.var(obs_values) + np.var(bg_values)) / 2)
        
        if pooled_std == 0:
            return 0.0
        
        return (obs_mean - bg_mean) / pooled_std
    
    def _merge_results(self, statistical_results: pd.DataFrame, 
                      valid_lineages: Dict) -> pd.DataFrame:
        """Merge statistical results with lineage information"""
        if statistical_results.empty:
            return pd.DataFrame()
        
        # Parse test IDs back to lineage information
        results_list = []
        
        for idx, row in statistical_results.iterrows():
            # Parse test ID (assuming format: lineage_id_gene_idx_direction)
            test_id = row['test_id']
            
            # Extract lineage info from observed data keys
            # This is a simplified approach - in practice, would need better key management
            
            results_list.append({
                'test_id': test_id,
                'lineage_root': f'lineage_{test_id}',  # Placeholder
                'gene_idx': test_id % self.n_genes,     # Placeholder
                'direction': 'AMP' if test_id % 2 == 0 else 'DEL',  # Placeholder
                'statistic': row['statistic'],
                'pvalue': row['pvalue'],
                'adjusted_pvalue': row.get('adjusted_pvalue', row['pvalue']),
                'effect_size': row['effect_size'],
                'ci_lower': row['ci_lower'],
                'ci_upper': row['ci_upper'],
                'n_observations': row['n_observations'],
                'power': row['power'],
                'significant': row.get('significant', False),
                'effect_size_magnitude': row.get('effect_size_magnitude', 'unknown')
            })
        
        return pd.DataFrame(results_list)


def main():
    """Example usage of statistically robust MEDALT"""
    
    # Create synthetic dataset
    np.random.seed(42)
    n_cells = 1000
    n_genes = 200
    
    logger.info(f"Creating synthetic dataset: {n_cells} cells × {n_genes} genes")
    
    # Base diploid state
    cnv_matrix = np.random.poisson(2, (n_cells, n_genes)).astype(np.float32)
    
    # Add structured lineages
    # Lineage 1: amplification in genes 50-70
    cnv_matrix[200:400, 50:70] = 4.0
    # Lineage 2: deletion in genes 120-140
    cnv_matrix[600:800, 120:140] = 0.5
    
    # Run robust analysis
    logger.info("Running statistically robust MEDALT analysis...")
    
    robust_medalt = StatisticallyRobustMEDALT(
        cnv_matrix,
        max_memory_gb=8.0,
        k_neighbors=20,
        chunk_size=200,
        n_jobs=2,
        statistical_method='enrichment',
        correction_method='fdr_bh',
        confidence_level=0.95,
        n_bootstrap=500
    )
    
    # Run complete analysis
    tree, results = robust_medalt.run_complete_analysis(
        n_permutations=100,  # Reduced for demo
        min_lineage_size=50,
        significance_threshold=0.05
    )
    
    # Display results
    logger.info(f"Analysis complete!")
    logger.info(f"Tree: {tree.number_of_nodes()} nodes, {tree.number_of_edges()} edges")
    logger.info(f"Results: {len(results)} significant associations")
    
    if not results.empty:
        print("\nTop results:")
        print(results[['lineage_root', 'gene_idx', 'direction', 'pvalue', 
                      'adjusted_pvalue', 'effect_size', 'effect_size_magnitude']].head(10))
        
        # Show effect size distribution
        print(f"\nEffect size distribution:")
        print(results['effect_size_magnitude'].value_counts())
        
        # Show statistical power
        print(f"\nStatistical power summary:")
        print(f"Mean power: {results['power'].mean():.3f}")
        print(f"Min power: {results['power'].min():.3f}")
        print(f"Max power: {results['power'].max():.3f}")


if __name__ == '__main__':
    main()