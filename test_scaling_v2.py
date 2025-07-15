#!/usr/bin/env python2
"""
Test v2_accurate optimization scaling
"""

import time
import random

def create_test_data(num_cells=100, num_segments=500, segment_size=1):
    """Create test data for scaling analysis"""
    random.seed(42)
    
    print "Creating %d cells with %d segments each..." % (num_cells, num_segments)
    
    nodes = {}
    for i in range(num_cells):
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

def test_scaling_performance():
    """Test scaling with increasing dataset sizes"""
    print "v2_accurate Scaling Analysis"
    print "=" * 50
    
    # Test with increasing sizes
    test_configs = [
        (500, 500, "500 genes x 500 cells"),
        (1000, 500, "1000 genes x 500 cells"), 
        (1500, 500, "1500 genes x 500 cells"),
        (500, 1000, "500 genes x 1000 cells"),
        (1000, 1000, "1000 genes x 1000 cells"),
        (1500, 1000, "1500 genes x 1000 cells"),
        (2000, 1000, "2000 genes x 1000 cells")
    ]
    
    results = []
    
    for num_genes, num_cells, description in test_configs:
        print "\n" + "-" * 50
        print "Testing: %s" % description
        
        # Create dataset
        nodes = create_test_data(num_cells, num_genes, 1)
        
        total_calcs = num_cells * (num_cells - 1) / 2
        memory_mb = num_cells * num_cells * 4 / (1024.0 * 1024.0)
        
        print "Distance calculations: %d" % total_calcs
        print "Estimated memory: %.1f MB" % memory_mb
        
        # Test v2_accurate
        try:
            from ComputeDistance import matrixbuilder
            
            start_time = time.time()
            keys, matrix = matrixbuilder(nodes)
            elapsed = time.time() - start_time
            
            print "Time: %.1f seconds" % elapsed
            print "Speed: %.0f calculations/second" % (total_calcs / elapsed)
            
            # Check if we're hitting timeout issues
            if elapsed > 300:  # 5 minutes
                print "WARNING: Approaching timeout threshold"
                estimated_3000x2000 = elapsed * (3000.0 * 2000.0) / (num_genes * num_cells)
                print "Estimated time for 3000x2000: %.0f minutes" % (estimated_3000x2000 / 60.0)
                break
            
            results.append((description, elapsed, total_calcs / elapsed))
            
        except Exception as e:
            print "FAILED: %s" % str(e)
            break
    
    # Summary
    print "\n" + "=" * 50
    print "SCALING RESULTS:"
    for desc, time_taken, speed in results:
        print "%s: %.1fs (%.0f calc/s)" % (desc, time_taken, speed)
    
    if results:
        # Extrapolate to target size
        last_speed = results[-1][2]
        target_calcs = 2000 * 1999 / 2  # 2000 cells
        estimated_time = target_calcs / last_speed
        
        print "\nExtrapolation to 2000 cells:"
        print "Estimated time: %.0f seconds (%.1f minutes)" % (estimated_time, estimated_time / 60.0)
        
        if estimated_time < 600:  # 10 minutes
            print "FEASIBLE: Should complete within timeout"
        else:
            print "CHALLENGING: May need chunked processing"

if __name__ == "__main__":
    test_scaling_performance()