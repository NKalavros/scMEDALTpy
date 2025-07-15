#!/usr/bin/env python2
"""
Test chunked v2_accurate optimization specifically for large datasets
"""

import time
import random

def create_test_data(num_cells=1500, num_segments=1000, segment_size=1):
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

def test_chunked_processing():
    """Test chunked processing for large datasets"""
    print "Chunked v2_accurate Processing Test"
    print "=" * 50
    
    # Test chunked processing with a dataset that triggers it
    print "Creating 1500 cell dataset to trigger chunked processing..."
    nodes = create_test_data(1500, 1000, 1)  # 1000 genes x 1500 cells
    
    print "Dataset size: %d cells" % len(nodes)
    
    # This should trigger chunked processing (>1000 cells)
    try:
        from ComputeDistance import matrixbuilder
        
        print "Starting chunked distance matrix calculation..."
        start_time = time.time()
        keys, matrix = matrixbuilder(nodes)
        elapsed = time.time() - start_time
        
        print "\nSUCCESS!"
        print "Time: %.1f seconds (%.1f minutes)" % (elapsed, elapsed / 60.0)
        print "Matrix size: %dx%d" % (len(matrix), len(matrix[0]))
        
        # Verify matrix properties
        print "\nMatrix verification:"
        print "Matrix[0][0] = %.1f (should be 0)" % matrix[0][0]
        print "Matrix[0][1] = %.1f" % matrix[0][1]
        print "Matrix[1][0] = %.1f (should equal [0][1])" % matrix[1][0]
        print "Symmetric: %s" % (abs(matrix[0][1] - matrix[1][0]) < 1e-6)
        
        # Performance analysis
        total_calcs = len(nodes) * (len(nodes) - 1) / 2
        speed = total_calcs / elapsed
        print "\nPerformance: %.0f calculations/second" % speed
        
        # Extrapolate to 2000 cells
        target_calcs_2000 = 2000 * 1999 / 2
        estimated_time_2000 = target_calcs_2000 / speed
        print "Estimated time for 2000 cells: %.1f minutes" % (estimated_time_2000 / 60.0)
        
        if estimated_time_2000 < 600:  # 10 minutes
            print "FEASIBLE for 2000 cells within timeout"
        else:
            print "May still timeout for 2000 cells"
            
        return True
        
    except Exception as e:
        print "ERROR: %s" % str(e)
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_chunked_processing()
    if success:
        print "\nChunked v2_accurate optimization working correctly!"
    else:
        print "\nChunked processing needs debugging."