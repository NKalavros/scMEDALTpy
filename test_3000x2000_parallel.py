#!/usr/bin/env python2
"""
Test 3000x2000 scale with optimized parallel processing
"""

import time
import random

def create_large_scale_data(num_cells=2000, num_segments=3000, segment_size=1):
    """Create large scale test data"""
    random.seed(42)
    
    print "Creating %d cells with %d segments each..." % (num_cells, num_segments)
    print "This represents %d genes x %d cells" % (num_segments, num_cells)
    
    nodes = {}
    for i in range(num_cells):
        if i % 100 == 0:
            print "  Creating cell %d/%d (%.1f%%)" % (i, num_cells, 100.0 * i / num_cells)
            
        cell_name = "Cell_%04d" % i
        segments = []
        for j in range(num_segments):
            segment = []
            for k in range(segment_size):
                if random.random() < 0.02:  # 2% zeros (deletions)
                    cn = 0
                elif random.random() < 0.7:  # 70% diploid (normal)
                    cn = 2
                else:  # 28% aneuploid (gains/losses)
                    cn = random.choice([1, 3, 4])
                segment.append(cn)
            segments.append(segment)
        nodes[cell_name] = segments
    
    print "Dataset created successfully: %d cells x %d genes" % (len(nodes), len(nodes[nodes.keys()[0]]))
    return nodes

def test_3000x2000_scale():
    """Test 3000x2000 scale with optimized parallel processing"""
    print "3000x2000 Scale Test with Optimized Parallel Processing"
    print "=" * 70
    
    # Create the large dataset
    print "Step 1: Creating 3000 genes x 2000 cells dataset..."
    start_creation = time.time()
    nodes = create_large_scale_data(2000, 3000, 1)
    creation_time = time.time() - start_creation
    
    print "Dataset creation completed in %.1f seconds" % creation_time
    
    # Calculate expected metrics
    total_calcs = 2000 * 1999 / 2
    memory_mb = 2000 * 2000 * 4 / (1024.0 * 1024.0)
    
    print "\nDataset metrics:"
    print "Total distance calculations: %d" % total_calcs
    print "Estimated matrix memory: %.1f MB" % memory_mb
    print "Average genes per cell: %d" % len(nodes[nodes.keys()[0]])
    
    # Test optimized parallel processing
    print "\nStep 2: Running optimized parallel processing..."
    try:
        from ComputeDistance_parallel_optimized import matrixbuilder
        
        start_time = time.time()
        keys, matrix = matrixbuilder(nodes)
        elapsed = time.time() - start_time
        
        print "\nSUCCESS!"
        print "Processing time: %.1f seconds (%.2f minutes)" % (elapsed, elapsed / 60.0)
        print "Total time (including creation): %.1f seconds (%.2f minutes)" % (elapsed + creation_time, (elapsed + creation_time) / 60.0)
        print "Speed: %.1f calculations/second" % (total_calcs / elapsed)
        
        # Verify results
        print "\nResult verification:"
        print "Matrix dimensions: %dx%d" % (len(matrix), len(matrix[0]))
        
        # Check diagonal (should be zeros)
        print "Diagonal verification:"
        diagonal_correct = True
        for i in [0, 500, 1000, 1500, 1999]:
            val = matrix[i][i]
            print "  Matrix[%d][%d] = %.6f" % (i, i, val)
            if abs(val) > 1e-6:
                diagonal_correct = False
        print "Diagonal correct: %s" % diagonal_correct
        
        # Check symmetry
        print "Symmetry verification:"
        symmetric = True
        max_diff = 0
        test_pairs = [(0, 1), (0, 100), (100, 200), (500, 1000), (1000, 1500), (1500, 1999)]
        
        for i, j in test_pairs:
            val_ij = matrix[i][j]
            val_ji = matrix[j][i]
            diff = abs(val_ij - val_ji)
            max_diff = max(max_diff, diff)
            print "  Matrix[%d][%d] = %.6f, Matrix[%d][%d] = %.6f, diff = %.2e" % (i, j, val_ij, j, i, val_ji, diff)
            if diff > 1e-6:
                symmetric = False
        
        print "Symmetric: %s (max diff: %.2e)" % (symmetric, max_diff)
        
        # Sample distance values
        print "Sample distance values:"
        for i, j in test_pairs:
            print "  Distance[%d][%d] = %.3f" % (i, j, matrix[i][j])
        
        # Performance analysis
        print "\nPerformance analysis:"
        print "Processing efficiency: %.1f calculations/second" % (total_calcs / elapsed)
        print "Memory efficiency: %.1f MB/second processed" % (memory_mb / elapsed)
        
        # Estimate for even larger scales
        print "\nScaling estimates:"
        speed = total_calcs / elapsed
        
        # 3000x3000 estimate
        calcs_3000x3000 = 3000 * 2999 / 2
        time_3000x3000 = calcs_3000x3000 / speed
        print "Estimated time for 3000x3000: %.1f minutes" % (time_3000x3000 / 60.0)
        
        # 5000x2000 estimate  
        calcs_5000x2000 = 2000 * 1999 / 2  # Same as current
        time_5000x2000 = calcs_5000x2000 / speed * (5000.0 / 3000.0)  # Scale by gene count
        print "Estimated time for 5000x2000: %.1f minutes" % (time_5000x2000 / 60.0)
        
        # Save sample for validation
        print "\nSaving results sample..."
        sample_size = 100
        with open("3000x2000_sample.txt", "w") as f:
            f.write("# Sample from 3000x2000 parallel processing\n")
            f.write("# Format: i j distance\n")
            for i in range(min(sample_size, len(matrix))):
                for j in range(min(sample_size, len(matrix[0]))):
                    f.write("%d %d %.6f\n" % (i, j, matrix[i][j]))
        
        print "Sample saved to 3000x2000_sample.txt"
        
        return True, elapsed, matrix
        
    except Exception as e:
        print "ERROR: %s" % str(e)
        import traceback
        traceback.print_exc()
        return False, 0, None

if __name__ == "__main__":
    success, time_taken, matrix = test_3000x2000_scale()
    if success:
        print "\n" + "=" * 70
        print "SUCCESS: 3000x2000 scale processing completed!"
        print "Processing time: %.2f minutes" % (time_taken / 60.0)
        print "The optimized parallel processing can handle large-scale datasets!"
    else:
        print "\n3000x2000 scale test failed - needs further optimization."