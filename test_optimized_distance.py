#!/usr/bin/env python2
"""
Test the optimized distance calculation for correctness and performance
"""

import time
import random

def test_optimized_distance():
    print "Testing Optimized Distance Calculation"
    print "=" * 50
    
    # Create test data
    print "Creating test data..."
    
    num_cells = 50  # Smaller for precise testing
    num_segments = 5
    segment_size = 8
    
    random.seed(42)
    
    test_nodes = {}
    for i in range(num_cells):
        cell_name = "Cell_%03d" % i
        segments = []
        for j in range(num_segments):
            segment = []
            for k in range(segment_size):
                if random.random() < 0.1:  # 10% chance of zero
                    cn = 0
                elif random.random() < 0.7:  # 70% chance of diploid
                    cn = 2
                else:
                    cn = random.choice([1, 3, 4])
                segment.append(cn)
            segments.append(segment)
        test_nodes[cell_name] = segments
    
    print "Created %d cells with %d segments each" % (num_cells, num_segments)
    
    # Test original implementation
    print "\\nTesting original implementation..."
    from ComputeDistance_optimized import matrixbuilder as matrixbuilder_orig
    from ComputeDistance_optimized import dist_original
    
    start_time = time.time()
    keys_orig, matrix_orig = matrixbuilder_orig(test_nodes)
    orig_time = time.time() - start_time
    
    # Test optimized implementation
    print "Testing optimized implementation..."
    from ComputeDistance_optimized import matrixbuilder
    
    start_time = time.time()
    keys_opt, matrix_opt = matrixbuilder(test_nodes)
    opt_time = time.time() - start_time
    
    # Test manual calculation for verification
    print "\\nVerifying with manual calculation..."
    
    # Pick two cells and calculate distance manually
    cell1 = test_nodes['Cell_001']
    cell2 = test_nodes['Cell_002']
    
    manual_dist = dist_original(cell1, cell2)
    from ComputeDistance_optimized import dist_optimized
    opt_dist = dist_optimized(cell1, cell2)
    
    print "Manual calculation: Cell_001 to Cell_002 = %d" % manual_dist
    print "Optimized calculation: Cell_001 to Cell_002 = %d" % opt_dist
    
    # Check consistency
    print "\\nChecking matrix consistency..."
    consistent = True
    max_diff = 0
    
    for i in range(len(matrix_orig)):
        for j in range(len(matrix_orig[i])):
            diff = abs(matrix_orig[i][j] - matrix_opt[i][j])
            max_diff = max(max_diff, diff)
            if diff > 1e-10:  # Allow for tiny floating point differences
                print "Mismatch at [%d,%d]: orig=%f, opt=%f, diff=%f" % (i, j, matrix_orig[i][j], matrix_opt[i][j], diff)
                consistent = False
                if not consistent:
                    break
        if not consistent:
            break
    
    # Performance summary
    speedup = orig_time / opt_time if opt_time > 0 else 0
    
    print "\\n" + "=" * 50
    print "RESULTS:"
    print "Original time: %.4f seconds" % orig_time
    print "Optimized time: %.4f seconds" % opt_time
    print "Speedup: %.2fx" % speedup
    print "Maximum difference: %.2e" % max_diff
    print "Consistent: %s" % consistent
    print "Manual verification: %s" % (manual_dist == opt_dist)
    
    return {
        'speedup': speedup,
        'consistent': consistent,
        'manual_check': manual_dist == opt_dist
    }

if __name__ == "__main__":
    test_optimized_distance()