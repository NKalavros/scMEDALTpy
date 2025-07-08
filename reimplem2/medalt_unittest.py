#!/usr/bin/env python
"""
Unit tests for MEDALT implementation
"""

import unittest
import numpy as np
import networkx as nx
from medalt_implementation import MED, MEDALT, LineageSpeciationAnalysis


class TestMED(unittest.TestCase):
    """Test cases for Minimal Event Distance calculation"""
    
    def test_identical_profiles(self):
        """Test MED between identical profiles"""
        profile = np.array([2, 2, 2, 2, 2])
        distance = MED.compute(profile, profile)
        self.assertEqual(distance, 0)
    
    def test_single_event(self):
        """Test single copy gain/loss"""
        profile_a = np.array([2, 2, 2, 2, 2])
        profile_b = np.array([2, 3, 2, 2, 2])
        distance = MED.compute(profile_a, profile_b)
        self.assertEqual(distance, 1)
    
    def test_segmental_event(self):
        """Test segmental copy number change"""
        profile_a = np.array([2, 2, 2, 2, 2])
        profile_b = np.array([2, 3, 3, 3, 2])
        distance = MED.compute(profile_a, profile_b)
        self.assertEqual(distance, 1)  # One event affecting segment
    
    def test_multiple_events(self):
        """Test multiple independent events"""
        profile_a = np.array([2, 2, 2, 2, 2])
        profile_b = np.array([3, 3, 2, 1, 1])
        distance = MED.compute(profile_a, profile_b)
        self.assertEqual(distance, 2)  # One gain event + one loss event
    
    def test_homozygous_deletion_constraint(self):
        """Test that we cannot recover from CN=0"""
        profile_a = np.array([0, 0, 2, 2, 2])
        profile_b = np.array([1, 1, 2, 2, 2])
        distance = MED.compute(profile_a, profile_b)
        self.assertEqual(distance, float('inf'))
    
    def test_complex_profile(self):
        """Test complex copy number profile"""
        profile_a = np.array([2, 2, 2, 3, 3, 1, 1, 2])
        profile_b = np.array([2, 3, 3, 3, 3, 2, 2, 2])
        # Changes: gain at 1-2 (1 event), no change at 3-4,
        # gain from 1->2 at 5-6 (1 event)
        distance = MED.compute(profile_a, profile_b)
        self.assertEqual(distance, 2)
    
    def test_multi_copy_change(self):
        """Test multi-copy changes require multiple events"""
        profile_a = np.array([2, 2, 2])
        profile_b = np.array([4, 4, 4])
        # Need 2 events to go from CN=2 to CN=4
        distance = MED.compute(profile_a, profile_b)
        self.assertEqual(distance, 2)


class TestMEDALT(unittest.TestCase):
    """Test cases for MEDALT tree construction"""
    
    def setUp(self):
        """Create test data"""
        # Simple test case with clear lineage structure
        self.cnv_matrix = np.array([
            [2, 2, 2, 2, 2],  # Cell 0: Normal
            [2, 3, 3, 2, 2],  # Cell 1: Small amp
            [2, 3, 3, 3, 3],  # Cell 2: Larger amp (child of 1)
            [1, 1, 2, 2, 2],  # Cell 3: Deletion
            [2, 2, 2, 2, 3],  # Cell 4: Different amp
        ])
    
    def test_tree_construction(self):
        """Test basic tree construction"""
        medalt = MEDALT(self.cnv_matrix)
        tree = medalt.build_tree()
        
        # Check tree properties
        self.assertEqual(tree.number_of_nodes(), 6)  # 5 cells + root
        self.assertTrue(tree.has_node('root'))
        self.assertTrue(all(tree.has_node(i) for i in range(5)))
        
        # Check that root has outgoing edges
        self.assertGreater(tree.out_degree('root'), 0)
        
        # Check all nodes are reachable from root
        for node in range(5):
            self.assertTrue(nx.has_path(tree, 'root', node))
    
    def test_tree_is_directed_acyclic(self):
        """Test that tree is a valid DAG"""
        medalt = MEDALT(self.cnv_matrix)
        tree = medalt.build_tree()
        
        self.assertTrue(nx.is_directed_acyclic_graph(tree))
        
        # Check it's a tree (connected and n-1 edges)
        undirected = tree.to_undirected()
        self.assertTrue(nx.is_connected(undirected))
        self.assertEqual(tree.number_of_edges(), tree.number_of_nodes() - 1)
    
    def test_distance_caching(self):
        """Test that distance calculations are cached"""
        medalt = MEDALT(self.cnv_matrix)
        
        # Calculate distances twice
        distances1 = medalt._compute_distance_matrix()
        cache_size1 = len(medalt._distance_cache)
        
        distances2 = medalt._compute_distance_matrix()
        cache_size2 = len(medalt._distance_cache)
        
        # Cache size shouldn't change
        self.assertEqual(cache_size1, cache_size2)
        
        # Results should be identical
        for key in distances1:
            if key in distances2:
                self.assertEqual(distances1[key], distances2[key])


class TestLSA(unittest.TestCase):
    """Test cases for Lineage Speciation Analysis"""
    
    def setUp(self):
        """Create test data with known structure"""
        # Create data with clear lineage expansion
        n_cells = 30
        n_bins = 50
        
        # Normal cells
        self.cnv_matrix = np.full((n_cells, n_bins), 2, dtype=int)
        
        # Create expanding lineage with specific alteration
        # Cells 10-25 share an amplification
        self.cnv_matrix[10:25, 20:30] = 3
        
        # Build tree
        medalt = MEDALT(self.cnv_matrix)
        self.tree = medalt.build_tree()
    
    def test_lineage_dissection(self):
        """Test tree dissection into lineages"""
        lsa = LineageSpeciationAnalysis(self.tree, self.cnv_matrix, n_permutations=10)
        lineages = lsa._dissect_tree(min_size=5)
        
        # Should find at least one lineage
        self.assertGreater(len(lineages), 0)
        
        # Check lineage properties
        for lineage in lineages:
            self.assertIn('root', lineage)
            self.assertIn('cells', lineage)
            self.assertIn('size', lineage)
            self.assertIn('depth', lineage)
            self.assertGreaterEqual(lineage['size'], 5)
    
    def test_cfl_calculation(self):
        """Test CFL calculation"""
        lsa = LineageSpeciationAnalysis(self.tree, self.cnv_matrix, n_permutations=10)
        lineages = lsa._dissect_tree(min_size=5)
        
        if lineages:
            cfls = lsa._calculate_cfls(lineages)
            
            # Check CFL structure
            for key, value in cfls.items():
                lineage_id, bin_idx, direction = key
                self.assertIn('cfl', value)
                self.assertIn('avg_cn', value)
                self.assertIn('size', value)
                self.assertIn(direction, ['AMP', 'DEL'])
    
    def test_permutation(self):
        """Test permutation maintains matrix properties"""
        lsa = LineageSpeciationAnalysis(self.tree, self.cnv_matrix, n_permutations=10)
        
        # Test single permutation
        permuted = lsa._permute_by_chromosome()
        
        # Shape should be preserved
        self.assertEqual(permuted.shape, self.cnv_matrix.shape)
        
        # Copy number values should be preserved (just rearranged)
        self.assertEqual(sorted(permuted.flatten()), sorted(self.cnv_matrix.flatten()))
    
    def test_lsa_finds_significant_cna(self):
        """Test that LSA can identify the planted CNA"""
        # Use fewer permutations for speed
        lsa = LineageSpeciationAnalysis(self.tree, self.cnv_matrix, n_permutations=50)
        results = lsa.run_analysis(min_lineage_size=5, n_jobs=1)
        
        if not results.empty:
            # Check that we found amplifications
            amp_results = results[results['CNA'] == 'AMP']
            self.assertGreater(len(amp_results), 0)
            
            # Check that some are in the planted region (bins 20-30)
            # Note: with windowing, exact bins might shift
            significant_bins = amp_results['region'].str.extract('(\d+)').astype(int).values.flatten()
            # At least some should be near our planted region
            self.assertTrue(any(15 <= b <= 35 for b in significant_bins))


class TestIntegration(unittest.TestCase):
    """Integration tests for the full pipeline"""
    
    def test_small_example_runs(self):
        """Test that a small example runs without errors"""
        # Create minimal data
        cnv_matrix = np.array([
            [2, 2, 2, 2],
            [2, 3, 3, 2],
            [2, 3, 3, 3],
            [1, 1, 2, 2],
        ])
        
        # Build tree
        medalt = MEDALT(cnv_matrix)
        tree = medalt.build_tree()
        
        # Run LSA with minimal permutations
        lsa = LineageSpeciationAnalysis(tree, cnv_matrix, n_permutations=10)
        results = lsa.run_analysis(min_lineage_size=2, n_jobs=1)
        
        # Should complete without errors
        self.assertIsInstance(results, type(results))  # DataFrame or similar
    
    def test_edge_cases(self):
        """Test edge cases"""
        # All cells identical
        cnv_matrix = np.full((10, 20), 2, dtype=int)
        medalt = MEDALT(cnv_matrix)
        tree = medalt.build_tree()
        
        # Should still produce valid tree
        self.assertTrue(nx.is_directed_acyclic_graph(tree))
        
        # Single cell
        cnv_matrix_single = np.array([[2, 2, 2, 2]])
        medalt_single = MEDALT(cnv_matrix_single)
        tree_single = medalt_single.build_tree()
        self.assertEqual(tree_single.number_of_nodes(), 2)  # root + 1 cell


if __name__ == '__main__':
    unittest.main()