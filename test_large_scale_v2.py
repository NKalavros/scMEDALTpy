#!/usr/bin/env python2
"""
Test v2_accurate optimization with large-scale dataset (3000x2000)
"""

import time
import random

def create_large_test_data(num_cells=2000, num_segments=3000, segment_size=1):
    """Create large-scale test data"""
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
                if random.random() < 0.02:  # 2% zeros (deleted regions)
                    cn = 0
                elif random.random() < 0.7:  # 70% normal diploid
                    cn = 2
                else:  # 28% aneuploid
                    cn = random.choice([1, 3, 4])
                segment.append(cn)
            segments.append(segment)
        nodes[cell_name] = segments
    
    print "Created %d cells successfully" % len(nodes)
    return nodes

def test_large_scale():
    """Test v2_accurate with large dataset"""
    print "Large Scale Test: v2_accurate Optimization"
    print "=" * 60
    
    # Create large dataset
    print "Step 1: Creating large-scale dataset..."
    nodes = create_large_test_data(2000, 3000, 1)  # 3000 genes x 2000 cells
    
    print "\nDataset created: %d cells, %d genes per cell" % (len(nodes), len(nodes[nodes.keys()[0]]))
    
    # Estimate memory requirements
    total_distance_calcs = len(nodes) * (len(nodes) - 1) / 2
    memory_mb = len(nodes) * len(nodes) * 4 / (1024.0 * 1024.0)  # float32
    
    print "Estimated matrix memory: %.1f MB" % memory_mb
    print "Total distance calculations: %d" % total_distance_calcs
    
    # Test with v2_accurate optimization
    print "\nStep 2: Testing v2_accurate optimization..."
    try:
        from ComputeDistance import matrixbuilder
        
        start_time = time.time()
        keys, matrix = matrixbuilder(nodes)
        total_time = time.time() - start_time
        
        print "\nSUCCESS!"
        print "Total time: %.1f seconds (%.1f minutes)" % (total_time, total_time / 60.0)
        print "Matrix size: %dx%d" % (len(matrix), len(matrix[0]))
        print "Average time per distance calculation: %.2f ms" % (total_time * 1000.0 / total_distance_calcs)
        
        # Quick verification of results
        print "\nQuick result verification:"
        print "Matrix[0][0] = %.1f (should be 0)" % matrix[0][0]
        print "Matrix[0][1] = %.1f" % matrix[0][1]
        print "Matrix[1][0] = %.1f (should equal Matrix[0][1])" % matrix[1][0]
        
        # Performance analysis
        distances_per_second = total_distance_calcs / total_time
        print "\nPerformance: %.0f distance calculations per second" % distances_per_second
        
        if total_time < 300:  # Less than 5 minutes
            print "EXCELLENT: Fast enough for routine analysis"
        elif total_time < 1800:  # Less than 30 minutes  
            print "GOOD: Reasonable for batch processing"
        else:
            print "SLOW: May need further optimization"
            
        return True
        
    except Exception as e:
        print "ERROR: %s" % str(e)
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_large_scale()
    if success:
        print "\nv2_accurate optimization successfully scales to 3000x2000!"
    else:
        print "\nv2_accurate optimization needs further work for large scale."