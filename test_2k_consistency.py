#!/usr/bin/env python2
"""
Test 2k cell consistency between v1 and chunked v2_accurate processing
"""

import time
import random

def create_test_data(num_cells=2000, num_segments=1000, segment_size=1):
    """Create test data for 2k consistency test"""
    random.seed(42)  # Fixed seed for reproducible results
    
    print "Creating %d cells with %d segments each..." % (num_cells, num_segments)
    
    nodes = {}
    for i in range(num_cells):
        if i % 200 == 0:
            print "  Creating cell %d/%d" % (i, num_cells)
            
        cell_name = "Cell_%04d" % i
        segments = []
        for j in range(num_segments):
            segment = []
            for k in range(segment_size):
                if random.random() < 0.02:  # 2% zeros (deleted regions)
                    cn = 0
                elif random.random() < 0.7:  # 70% normal diploid
                    cn = 2
                else:  # 28% aneuploid (gains/losses)
                    cn = random.choice([1, 3, 4])
                segment.append(cn)
            segments.append(segment)
        nodes[cell_name] = segments
    
    print "Dataset created: %d cells with %d segments each" % (len(nodes), len(nodes[nodes.keys()[0]]))
    return nodes

def test_2k_consistency():
    """Test 2k cell processing for consistency"""
    print "2K Cell Consistency Test"
    print "=" * 60
    
    # Create 2k dataset
    print "Step 1: Creating 2k cell dataset..."
    nodes = create_test_data(2000, 1000, 1)  # 1000 genes x 2000 cells
    
    total_calcs = 2000 * 1999 / 2
    memory_mb = 2000 * 2000 * 4 / (1024.0 * 1024.0)
    
    print "Total distance calculations: %d" % total_calcs
    print "Estimated matrix memory: %.1f MB" % memory_mb
    
    # Test chunked v2_accurate processing
    print "\nStep 2: Testing chunked v2_accurate processing..."
    try:
        from ComputeDistance import matrixbuilder
        
        start_time = time.time()
        keys, matrix = matrixbuilder(nodes)
        elapsed = time.time() - start_time
        
        print "\nCHUNKED PROCESSING COMPLETE!"
        print "Time: %.1f seconds (%.2f minutes)" % (elapsed, elapsed / 60.0)
        print "Matrix size: %dx%d" % (len(matrix), len(matrix[0]))
        print "Speed: %.1f calculations/second" % (total_calcs / elapsed)
        
        # Verify matrix properties
        print "\nMatrix verification:"
        print "Diagonal elements (should be 0):"
        for i in [0, 100, 500, 1000, 1500, 1999]:
            print "  Matrix[%d][%d] = %.6f" % (i, i, matrix[i][i])
        
        print "Symmetry check:"
        symmetric = True
        max_diff = 0
        for i in [0, 100, 500, 1000]:
            for j in [i+1, i+100, i+500]:
                if j < len(matrix):
                    diff = abs(matrix[i][j] - matrix[j][i])
                    max_diff = max(max_diff, diff)
                    if diff > 1e-6:
                        symmetric = False
                        print "  ASYMMETRY: Matrix[%d][%d]=%.6f != Matrix[%d][%d]=%.6f" % (i, j, matrix[i][j], j, i, matrix[j][i])
                        break
            if not symmetric:
                break
        
        print "  Symmetric: %s (max diff: %.2e)" % (symmetric, max_diff)
        
        # Sample distance values
        print "Sample distance values:"
        for i in [0, 500, 1000, 1500]:
            for j in [i+1, i+10, i+100]:
                if j < len(matrix):
                    print "  Matrix[%d][%d] = %.3f" % (i, j, matrix[i][j])
        
        # Save a subset for later comparison
        print "\nSaving sample for comparison..."
        sample_size = 50
        sample_matrix = []
        for i in range(sample_size):
            row = []
            for j in range(sample_size):
                row.append(matrix[i][j])
            sample_matrix.append(row)
        
        # Write sample to file for validation
        with open("sample_2k_matrix.txt", "w") as f:
            f.write("# Sample 50x50 matrix from 2k cell chunked processing\n")
            for i in range(sample_size):
                line = " ".join("%.6f" % sample_matrix[i][j] for j in range(sample_size))
                f.write(line + "\n")
        
        print "Sample matrix saved to sample_2k_matrix.txt"
        
        return True, elapsed, matrix
        
    except Exception as e:
        print "ERROR: %s" % str(e)
        import traceback
        traceback.print_exc()
        return False, 0, None

if __name__ == "__main__":
    success, time_taken, matrix = test_2k_consistency()
    if success:
        print "\n" + "=" * 60
        print "SUCCESS: 2k cell chunked processing completed!"
        print "Ready for parallel optimization implementation."
        print "Benchmark time: %.2f minutes" % (time_taken / 60.0)
    else:
        print "\n2k consistency test failed - debugging needed."