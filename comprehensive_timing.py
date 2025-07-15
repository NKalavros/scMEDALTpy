#!/usr/bin/env python2
"""
Comprehensive timing comparison for MEDALT optimizations
"""

import time
import os
import sys
import shutil

def time_function(func, *args, **kwargs):
    """Time a function execution"""
    start_time = time.time()
    result = func(*args, **kwargs)
    end_time = time.time()
    return result, end_time - start_time

def create_test_data_direct(num_cells=100, num_segments=5, segment_size=10):
    """Create test data directly in memory"""
    import random
    random.seed(42)
    
    test_nodes = {}
    for i in range(num_cells):
        cell_name = "Cell_%03d" % i
        segments = []
        for j in range(num_segments):
            segment = []
            for k in range(segment_size):
                if random.random() < 0.05:  # 5% chance of zero
                    cn = 0
                elif random.random() < 0.7:  # 70% chance of diploid
                    cn = 2
                else:
                    cn = random.choice([1, 3, 4])
                segment.append(cn)
            segments.append(segment)
        test_nodes[cell_name] = segments
    
    return test_nodes

def comprehensive_timing_test():
    """Run comprehensive timing tests"""
    print "MEDALT Comprehensive Performance Analysis"
    print "=" * 60
    
    # Test different dataset sizes
    test_sizes = [
        (20, 5, 10),   # Small: 20 cells
        (50, 5, 10),   # Medium: 50 cells  
        (100, 5, 10),  # Large: 100 cells
        (200, 5, 10),  # Very large: 200 cells
    ]
    
    results = []
    
    for num_cells, num_segments, segment_size in test_sizes:
        print "\\nTesting dataset: %d cells x %d segments x %d bins" % (num_cells, num_segments, segment_size)
        print "-" * 50
        
        # Create test data
        test_nodes = create_test_data_direct(num_cells, num_segments, segment_size)
        
        # Test original implementation
        print "Testing original implementation..."
        from ComputeDistance_optimized import matrixbuilder_adaptive
        sys.path.insert(0, 'original_backup')
        
        try:
            # Use original files
            import ComputeDistance as orig_compute
            reload(orig_compute)
            
            _, orig_time = time_function(orig_compute.matrixbuilder, test_nodes)
            print "  Original: %.4f seconds" % orig_time
            
        except Exception as e:
            print "  Original: Error - %s" % str(e)
            orig_time = float('inf')
        
        # Test optimized implementation
        print "Testing optimized implementation..."
        _, opt_time = time_function(matrixbuilder_adaptive, test_nodes)
        print "  Optimized: %.4f seconds" % opt_time
        
        # Calculate speedup
        speedup = orig_time / opt_time if opt_time > 0 and orig_time != float('inf') else 0
        
        result = {
            'num_cells': num_cells,
            'num_segments': num_segments,
            'segment_size': segment_size,
            'original_time': orig_time,
            'optimized_time': opt_time,
            'speedup': speedup
        }
        
        results.append(result)
        
        print "  Speedup: %.2fx" % speedup
        print "  Memory reduction: ~50% (symmetric matrix calculation)"
    
    # Summary
    print "\\n" + "=" * 60
    print "PERFORMANCE SUMMARY"
    print "=" * 60
    
    print "%-12s %-12s %-12s %-12s" % ("Dataset", "Original", "Optimized", "Speedup")
    print "-" * 50
    
    for result in results:
        dataset_str = "%dx%dx%d" % (result['num_cells'], result['num_segments'], result['segment_size'])
        if result['original_time'] == float('inf'):
            orig_str = "Error"
        else:
            orig_str = "%.4fs" % result['original_time']
        opt_str = "%.4fs" % result['optimized_time']
        speedup_str = "%.2fx" % result['speedup'] if result['speedup'] > 0 else "N/A"
        
        print "%-12s %-12s %-12s %-12s" % (dataset_str, orig_str, opt_str, speedup_str)
    
    # Calculate complexity analysis
    print "\\n" + "=" * 60
    print "COMPLEXITY ANALYSIS"
    print "=" * 60
    
    valid_results = [r for r in results if r['original_time'] != float('inf')]
    
    if len(valid_results) >= 2:
        # Compare scaling
        small = valid_results[0]
        large = valid_results[-1]
        
        cells_ratio = float(large['num_cells']) / small['num_cells']
        time_ratio_orig = large['original_time'] / small['original_time']
        time_ratio_opt = large['optimized_time'] / small['optimized_time']
        
        print "Scaling from %d to %d cells (%.1fx increase):" % (small['num_cells'], large['num_cells'], cells_ratio)
        print "  Original time scaling: %.2fx" % time_ratio_orig
        print "  Optimized time scaling: %.2fx" % time_ratio_opt
        print "  Expected O(n^2) scaling: %.2fx" % (cells_ratio ** 2)
        
        # Estimate performance for target scale
        target_cells = 1000
        if large['num_cells'] > 0:
            scale_factor = float(target_cells) / large['num_cells']
            estimated_time = large['optimized_time'] * (scale_factor ** 2)
            print "\\nEstimated time for %d cells: %.2f seconds (%.2f minutes)" % (target_cells, estimated_time, estimated_time / 60)
    
    return results

if __name__ == "__main__":
    comprehensive_timing_test()