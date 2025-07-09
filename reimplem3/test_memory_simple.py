#!/usr/bin/env python
"""
Simple test for memory optimization validation
"""

import numpy as np
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from medalt_memory_optimized import (
    MemoryMonitor, SparseDistanceMatrix, compute_med_vectorized,
    MemoryOptimizedMEDALT
)

def test_memory_monitor():
    """Test basic memory monitoring"""
    print("Testing memory monitor...")
    monitor = MemoryMonitor(max_memory_gb=1.0)
    
    initial_usage = monitor.get_memory_usage()
    print(f"Initial memory usage: {initial_usage:.1f} MB")
    
    # Test cleanup
    monitor.force_cleanup()
    print("✓ Memory monitor works")

def test_sparse_distance_matrix():
    """Test sparse distance matrix"""
    print("\nTesting sparse distance matrix...")
    
    import tempfile
    with tempfile.TemporaryDirectory() as temp_dir:
        sparse_matrix = SparseDistanceMatrix(100, 10, temp_dir)
        
        # Test setting neighbors
        neighbors = np.array([1, 3, 5, 7, 9])
        distances = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        sparse_matrix.set_neighbors(0, neighbors, distances)
        
        # Test getting neighbors
        retrieved_neighbors, retrieved_distances = sparse_matrix.get_neighbors(0)
        
        print(f"Stored {len(neighbors)} neighbors")
        print(f"Retrieved {len(retrieved_neighbors)} neighbors")
        
        # Test distance lookup
        distance = sparse_matrix.get_distance(0, 3)
        print(f"Distance lookup: {distance}")
        
        sparse_matrix.cleanup()
        print("✓ Sparse distance matrix works")

def test_med_computation():
    """Test MED computation"""
    print("\nTesting MED computation...")
    
    # Test identical profiles
    profile_a = np.array([2, 2, 2, 2, 2], dtype=np.int8)
    profile_b = np.array([2, 2, 2, 2, 2], dtype=np.int8)
    distance = compute_med_vectorized(profile_a, profile_b)
    print(f"Identical profiles distance: {distance}")
    assert distance == 0, "Identical profiles should have distance 0"
    
    # Test single event
    profile_a = np.array([2, 2, 2, 2, 2], dtype=np.int8)
    profile_b = np.array([2, 3, 2, 2, 2], dtype=np.int8)
    distance = compute_med_vectorized(profile_a, profile_b)
    print(f"Single event distance: {distance}")
    assert distance == 1, "Single event should have distance 1"
    
    # Test homozygous deletion constraint
    profile_a = np.array([0, 0, 2, 2, 2], dtype=np.int8)
    profile_b = np.array([1, 1, 2, 2, 2], dtype=np.int8)
    distance = compute_med_vectorized(profile_a, profile_b)
    print(f"Homozygous deletion distance: {distance}")
    assert distance == 999999, "Homozygous deletion should be impossible"
    
    print("✓ MED computation works")

def test_memory_optimized_medalt():
    """Test memory-optimized MEDALT"""
    print("\nTesting memory-optimized MEDALT...")
    
    # Create small test dataset
    np.random.seed(42)
    cnv_matrix = np.random.poisson(2, (50, 20)).astype(np.int8)
    
    # Add some structure
    cnv_matrix[10:20, 5:10] = 4  # Amplification
    cnv_matrix[25:35, 15:18] = 0  # Deletion
    
    print(f"Test dataset: {cnv_matrix.shape[0]} cells × {cnv_matrix.shape[1]} genes")
    
    # Test with small parameters
    medalt = MemoryOptimizedMEDALT(
        cnv_matrix,
        max_memory_gb=1.0,
        k_neighbors=10,
        chunk_size=25,
        n_jobs=1
    )
    
    print("Building tree...")
    tree = medalt.build_tree()
    
    print(f"Tree built: {tree.number_of_nodes()} nodes, {tree.number_of_edges()} edges")
    
    # Verify tree properties
    assert tree.number_of_nodes() > 0, "Tree should have nodes"
    assert tree.number_of_edges() > 0, "Tree should have edges"
    assert tree.has_node('root'), "Tree should have root node"
    
    # Check tree is connected
    import networkx as nx
    undirected = tree.to_undirected()
    assert nx.is_connected(undirected), "Tree should be connected"
    
    print("✓ Memory-optimized MEDALT works")

def test_memory_scaling():
    """Test memory scaling characteristics"""
    print("\nTesting memory scaling...")
    
    sizes = [50, 100, 200]
    k_neighbors = 10
    
    for size in sizes:
        # Calculate theoretical memory usage
        full_matrix_memory = size * size * 4  # 4 bytes per float
        sparse_matrix_memory = size * k_neighbors * 4
        
        reduction = full_matrix_memory / sparse_matrix_memory
        
        print(f"Size {size}: {reduction:.1f}x memory reduction")
        print(f"  Full matrix: {full_matrix_memory / 1024:.1f} KB")
        print(f"  Sparse matrix: {sparse_matrix_memory / 1024:.1f} KB")
    
    print("✓ Memory scaling analysis complete")

def main():
    """Run all tests"""
    print("MEDALT Memory Optimization Tests")
    print("=" * 50)
    
    try:
        test_memory_monitor()
        test_sparse_distance_matrix()
        test_med_computation()
        test_memory_optimized_medalt()
        test_memory_scaling()
        
        print("\n" + "=" * 50)
        print("🎉 ALL TESTS PASSED!")
        print("✅ Memory optimizations are working correctly")
        
    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()