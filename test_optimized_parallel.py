#!/usr/bin/env python2
"""
Test optimized parallel processing
"""

import time
import random

def create_test_data(num_cells=1000, num_segments=1000, segment_size=1):
    """Create test data"""
    random.seed(42)
    
    print "Creating %d cells with %d segments each..." % (num_cells, num_segments)
    
    nodes = {}
    for i in range(num_cells):
        if i % 100 == 0:
            print "  Creating cell %d/%d" % (i, num_cells)
            
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
    
    return nodes

def test_optimized_parallel():
    """Test optimized parallel processing"""
    print "Optimized Parallel Processing Test"
    print "=" * 60
    
    # Test with realistic scale
    test_configs = [
        (500, 500, "500 genes x 500 cells"),
        (1000, 1000, "1000 genes x 1000 cells"),
        (2000, 1000, "2000 genes x 1000 cells"),
        (1000, 2000, "1000 genes x 2000 cells")
    ]
    
    results = []
    
    for num_genes, num_cells, description in test_configs:
        print "\n" + "=" * 60
        print "Testing: %s" % description
        
        # Create dataset
        nodes = create_test_data(num_cells, num_genes, 1)
        
        total_calcs = num_cells * (num_cells - 1) / 2
        memory_mb = num_cells * num_cells * 4 / (1024.0 * 1024.0)
        
        print "Distance calculations: %d" % total_calcs
        print "Estimated memory: %.1f MB" % memory_mb
        
        # Test optimized parallel
        print "\nTesting optimized parallel processing..."
        try:
            from ComputeDistance_parallel_optimized import matrixbuilder
            
            start_time = time.time()
            keys, matrix = matrixbuilder(nodes)
            elapsed = time.time() - start_time
            
            print "SUCCESS!"
            print "Time: %.1f seconds (%.2f minutes)" % (elapsed, elapsed / 60.0)
            print "Speed: %.1f calculations/second" % (total_calcs / elapsed)
            
            # Verify results
            print "\nResult verification:"
            print "Matrix size: %dx%d" % (len(matrix), len(matrix[0]))
            print "Matrix[0][0] = %.6f (should be 0)" % matrix[0][0]
            print "Matrix[0][1] = %.6f" % matrix[0][1]
            print "Matrix[1][0] = %.6f (should equal [0][1])" % matrix[1][0]
            print "Symmetric: %s" % (abs(matrix[0][1] - matrix[1][0]) < 1e-6)
            
            results.append((description, elapsed, total_calcs / elapsed))
            
            # Check if we can scale to 2000 cells
            if num_cells == 2000:
                print "\nSUCCESS: 2000 cells processed!"
                if elapsed < 1800:  # 30 minutes
                    print "Time is reasonable for production use"
                else:
                    print "Time is long but feasible with patience"
            
        except Exception as e:
            print "FAILED: %s" % str(e)
            import traceback
            traceback.print_exc()
            break
    
    # Summary
    print "\n" + "=" * 60
    print "OPTIMIZED PARALLEL PROCESSING RESULTS:"
    for desc, time_taken, speed in results:
        print "%s: %.1fs (%.1f calc/s)" % (desc, time_taken, speed)
    
    if results:
        print "\nPerformance scaling analysis:"
        for i, (desc, time_taken, speed) in enumerate(results):
            if i > 0:
                prev_speed = results[i-1][2]
                efficiency = speed / prev_speed
                print "%s: %.1f%% efficiency vs previous" % (desc, efficiency * 100)

if __name__ == "__main__":
    test_optimized_parallel()