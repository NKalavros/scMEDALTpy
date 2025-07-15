#!/usr/bin/env python2
"""
Test v2 optimizations against v1
"""

import time
import random

def create_test_data(num_cells=100, num_segments=5, segment_size=20):
    """Create test data for comparison"""
    random.seed(42)
    
    nodes = {}
    for i in range(num_cells):
        cell_name = "Cell_%03d" % i
        segments = []
        for j in range(num_segments):
            segment = []
            for k in range(segment_size):
                if random.random() < 0.05:
                    cn = 0
                elif random.random() < 0.7:
                    cn = 2
                else:
                    cn = random.choice([1, 3, 4])
                segment.append(cn)
            segments.append(segment)
        nodes[cell_name] = segments
    
    return nodes

def test_optimization_comparison():
    """Compare v1 vs v2 optimizations"""
    print "Testing v1 vs v2 Distance Calculation Optimizations"
    print "=" * 60
    
    # Test with increasing dataset sizes
    test_sizes = [50, 100, 200, 500]
    
    for num_cells in test_sizes:
        print "\\nTesting with %d cells..." % num_cells
        print "-" * 40
        
        # Create test data
        nodes = create_test_data(num_cells, 5, 20)
        
        # Test v1 (current optimized version)
        print "Testing v1 (current optimized)..."
        try:
            from ComputeDistance import matrixbuilder as matrixbuilder_v1
            
            start_time = time.time()
            keys_v1, matrix_v1 = matrixbuilder_v1(nodes)
            v1_time = time.time() - start_time
            
            print "  v1 time: %.3f seconds" % v1_time
            v1_success = True
            
        except Exception as e:
            print "  v1 failed: %s" % str(e)
            v1_time = float('inf')
            v1_success = False
        
        # Test v2 (new ultra-optimized version)
        print "Testing v2 (ultra-optimized)..."
        try:
            from ComputeDistance_v2 import matrixbuilder as matrixbuilder_v2
            
            start_time = time.time()
            keys_v2, matrix_v2 = matrixbuilder_v2(nodes)
            v2_time = time.time() - start_time
            
            print "  v2 time: %.3f seconds" % v2_time
            v2_success = True
            
        except Exception as e:
            print "  v2 failed: %s" % str(e)
            v2_time = float('inf')
            v2_success = False
        
        # Compare results
        if v1_success and v2_success:
            speedup = v1_time / v2_time if v2_time > 0 else 0
            print "  Speedup: %.2fx" % speedup
            
            # Quick consistency check (first few values)
            consistent = True
            for i in range(min(5, len(matrix_v1))):
                for j in range(min(5, len(matrix_v1[i]))):
                    diff = abs(matrix_v1[i][j] - matrix_v2[i][j])
                    if diff > 0.1:  # Allow for approximation differences
                        consistent = False
                        break
                if not consistent:
                    break
            
            print "  Results consistent: %s" % consistent
            
        elif v1_success:
            print "  v2 failed, v1 succeeded"
        elif v2_success:
            print "  v1 failed, v2 succeeded"
        else:
            print "  Both failed"
        
        # Memory usage estimate
        memory_mb = num_cells * num_cells * 4 / (1024.0 * 1024.0)  # float32
        print "  Estimated memory: %.1f MB" % memory_mb
    
    print "\\n" + "=" * 60
    print "RECOMMENDATION:"
    
    if num_cells >= 500:
        print "For datasets with 500+ cells, v2 optimizations are essential."
        print "v2 uses memory-efficient processing and approximations when needed."
    else:
        print "For smaller datasets, v1 optimizations are sufficient."
        print "v2 provides additional scalability for large datasets."

if __name__ == "__main__":
    test_optimization_comparison()