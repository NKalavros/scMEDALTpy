#!/usr/bin/env python2
"""
Test 4k cells with optimized parallel processing for metrics
"""

import time
import random

def create_4k_test_data(num_cells=4000, num_segments=1000, segment_size=1):
    """Create 4k cell test data"""
    random.seed(42)
    
    print "Creating %d cells with %d segments each..." % (num_cells, num_segments)
    
    nodes = {}
    for i in range(num_cells):
        if i % 200 == 0:
            print "  Creating cell %d/%d (%.1f%%)" % (i, num_cells, 100.0 * i / num_cells)
            
        cell_name = "Cell_%04d" % i
        segments = []
        for j in range(num_segments):
            segment = []
            for k in range(segment_size):
                if random.random() < 0.02:  # 2% zeros
                    cn = 0
                elif random.random() < 0.7:  # 70% diploid
                    cn = 2
                else:  # 28% aneuploid
                    cn = random.choice([1, 3, 4])
                segment.append(cn)
            segments.append(segment)
        nodes[cell_name] = segments
    
    print "Dataset created: %d cells with %d segments each" % (len(nodes), len(nodes[nodes.keys()[0]]))
    return nodes

def test_4k_cells():
    """Test 4k cells with optimized parallel processing"""
    print "4K Cells Optimization Test"
    print "=" * 60
    
    # Create 4k dataset
    print "Step 1: Creating 4k cell dataset..."
    start_creation = time.time()
    nodes = create_4k_test_data(4000, 1000, 1)
    creation_time = time.time() - start_creation
    
    print "Dataset creation: %.1f seconds" % creation_time
    
    # Calculate metrics
    total_calcs = 4000 * 3999 / 2
    memory_mb = 4000 * 4000 * 4 / (1024.0 * 1024.0)
    
    print "\nDataset metrics:"
    print "Total cells: %d" % len(nodes)
    print "Segments per cell: %d" % len(nodes[nodes.keys()[0]])
    print "Total distance calculations: %d" % total_calcs
    print "Matrix memory estimate: %.1f MB" % memory_mb
    
    # Test optimized parallel processing
    print "\nStep 2: Running optimized parallel processing..."
    try:
        from ComputeDistance import matrixbuilder
        
        start_time = time.time()
        keys, matrix = matrixbuilder(nodes)
        elapsed = time.time() - start_time
        
        print "\nSUCCESS!"
        print "Processing time: %.1f seconds (%.2f minutes)" % (elapsed, elapsed / 60.0)
        print "Total time (with creation): %.1f seconds (%.2f minutes)" % (elapsed + creation_time, (elapsed + creation_time) / 60.0)
        
        # Performance metrics
        speed = total_calcs / elapsed
        print "\nPerformance metrics:"
        print "Speed: %.1f calculations/second" % speed
        print "Throughput: %.1f MB/second" % (memory_mb / elapsed)
        print "Efficiency: %.1f calculations/core/second" % (speed / 8)
        
        # Verify correctness
        print "\nCorrectness verification:"
        print "Matrix dimensions: %dx%d" % (len(matrix), len(matrix[0]))
        
        # Check diagonal
        diagonal_ok = True
        for i in [0, 1000, 2000, 3000, 3999]:
            if matrix[i][i] != 0.0:
                diagonal_ok = False
                break
        print "Diagonal zeros: %s" % diagonal_ok
        
        # Check symmetry
        symmetric_ok = True
        test_pairs = [(0, 1), (100, 200), (1000, 2000), (2000, 3000)]
        for i, j in test_pairs:
            if abs(matrix[i][j] - matrix[j][i]) > 1e-6:
                symmetric_ok = False
                break
        print "Matrix symmetric: %s" % symmetric_ok
        
        # Sample values
        print "\nSample distance values:"
        for i, j in [(0, 1), (100, 200), (1000, 2000), (2000, 3000)]:
            print "  Distance[%d][%d] = %.3f" % (i, j, matrix[i][j])
        
        # Scaling analysis
        print "\nScaling analysis:"
        print "4k cells performance: %.1f calculations/second" % speed
        
        # Compare to previous results
        prev_2k_speed = 7545  # From previous test
        efficiency_4k = speed / prev_2k_speed
        print "4k vs 2k efficiency: %.1f%%" % (efficiency_4k * 100)
        
        # Extrapolate to larger scales
        print "\nExtrapolation estimates:"
        
        # 5k cells
        calcs_5k = 5000 * 4999 / 2
        time_5k = calcs_5k / speed
        print "5k cells estimated time: %.1f minutes" % (time_5k / 60.0)
        
        # 8k cells
        calcs_8k = 8000 * 7999 / 2
        time_8k = calcs_8k / speed
        print "8k cells estimated time: %.1f minutes" % (time_8k / 60.0)
        
        # 10k cells
        calcs_10k = 10000 * 9999 / 2
        time_10k = calcs_10k / speed
        print "10k cells estimated time: %.1f minutes" % (time_10k / 60.0)
        
        # Resource utilization
        print "\nResource utilization:"
        print "CPU cores used: 8"
        print "Memory peak: ~%.1f MB" % memory_mb
        print "Processing efficiency: %.1f%%" % (100.0 * speed / 30000)  # Estimate vs theoretical max
        
        return True, elapsed, speed, matrix
        
    except Exception as e:
        print "ERROR: %s" % str(e)
        import traceback
        traceback.print_exc()
        return False, 0, 0, None

if __name__ == "__main__":
    success, time_taken, speed, matrix = test_4k_cells()
    if success:
        print "\n" + "=" * 60
        print "4K CELLS OPTIMIZATION METRICS:"
        print "Time: %.2f minutes" % (time_taken / 60.0)
        print "Speed: %.1f calculations/second" % speed
        print "Ready for full MEDALT pipeline testing"
        print "=" * 60
    else:
        print "\n4k cells test failed"