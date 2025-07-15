#!/usr/bin/env python2
"""
Focused timing test for distance calculation optimizations
"""

import time
import sys
import os
sys.path.insert(0, '.')

# Test different distance calculation implementations
def test_distance_implementations():
    print "Distance Calculation Timing Test"
    print "=" * 50
    
    # Create test data
    print "Creating test data..."
    
    # Simulate 100 cells with 5 chromosomal segments each
    # Each segment has 10 bins (similar to real data after segmentation)
    num_cells = 100
    num_segments = 5
    segment_size = 10
    
    import random
    random.seed(42)
    
    test_nodes = {}
    for i in range(num_cells):
        cell_name = "Cell_%03d" % i
        segments = []
        for j in range(num_segments):
            segment = []
            for k in range(segment_size):
                # Generate copy numbers (mostly 2 with some variation)
                if random.random() < 0.8:
                    cn = 2
                else:
                    cn = random.choice([0, 1, 3, 4])
                segment.append(cn)
            segments.append(segment)
        test_nodes[cell_name] = segments
    
    print "Created %d cells with %d segments each" % (num_cells, num_segments)
    
    # Test original implementation
    print "\\nTesting original implementation..."
    
    # Import original functions
    from ComputeDistance import matrixbuilder as matrixbuilder_orig
    
    start_time = time.time()
    keys_orig, matrix_orig = matrixbuilder_orig(test_nodes)
    orig_time = time.time() - start_time
    
    print "Original implementation: %.3f seconds" % orig_time
    
    # Test numpy optimized implementation
    print "\\nTesting numpy optimized implementation..."
    
    from ComputeDistance_numpy import matrixbuilder as matrixbuilder_numpy
    
    start_time = time.time()
    keys_numpy, matrix_numpy = matrixbuilder_numpy(test_nodes)
    numpy_time = time.time() - start_time
    
    print "Numpy implementation: %.3f seconds" % numpy_time
    
    # Calculate speedup
    speedup = orig_time / numpy_time if numpy_time > 0 else 0
    print "\\nSpeedup: %.2fx" % speedup
    
    # Verify consistency (check first few values)
    print "\\nVerifying consistency..."
    consistent = True
    for i in range(min(5, len(matrix_orig))):
        for j in range(min(5, len(matrix_orig[i]))):
            if abs(matrix_orig[i][j] - matrix_numpy[i][j]) > 1e-6:
                print "Mismatch at [%d,%d]: orig=%.6f, numpy=%.6f" % (i, j, matrix_orig[i][j], matrix_numpy[i][j])
                consistent = False
                break
        if not consistent:
            break
    
    if consistent:
        print "Results are consistent!"
    else:
        print "WARNING: Results are inconsistent!"
    
    return {
        'original_time': orig_time,
        'numpy_time': numpy_time,
        'speedup': speedup,
        'consistent': consistent
    }

if __name__ == "__main__":
    results = test_distance_implementations()
    print "\\n" + "=" * 50
    print "FINAL RESULTS:"
    print "Original: %.3f seconds" % results['original_time']
    print "Numpy:    %.3f seconds" % results['numpy_time']
    print "Speedup:  %.2fx" % results['speedup']
    print "Consistent: %s" % results['consistent']