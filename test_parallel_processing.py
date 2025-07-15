#!/usr/bin/env python2
"""
Test parallel processing implementation
"""

import time
import random

def create_test_data(num_cells=1000, num_segments=500, segment_size=1):
    """Create test data for parallel processing test"""
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

def test_parallel_vs_sequential():
    """Test parallel vs sequential processing"""
    print "Parallel Processing Test"
    print "=" * 50
    
    # Test sizes
    test_sizes = [200, 500, 1000]
    
    for num_cells in test_sizes:
        print "\nTesting with %d cells..." % num_cells
        print "-" * 30
        
        # Create test data
        nodes = create_test_data(num_cells, 500, 1)
        
        total_calcs = num_cells * (num_cells - 1) / 2
        print "Total calculations: %d" % total_calcs
        
        # Test sequential (v2_accurate)
        print "Testing sequential v2_accurate..."
        try:
            from ComputeDistance import matrixbuilder as matrixbuilder_seq
            
            start_time = time.time()
            keys_seq, matrix_seq = matrixbuilder_seq(nodes)
            seq_time = time.time() - start_time
            
            print "  Sequential time: %.1f seconds" % seq_time
            seq_success = True
            
        except Exception as e:
            print "  Sequential failed: %s" % str(e)
            seq_time = float('inf')
            seq_success = False
            matrix_seq = None
        
        # Test parallel
        print "Testing parallel processing..."
        try:
            from ComputeDistance_parallel import matrixbuilder as matrixbuilder_par
            
            start_time = time.time()
            keys_par, matrix_par = matrixbuilder_par(nodes)
            par_time = time.time() - start_time
            
            print "  Parallel time: %.1f seconds" % par_time
            par_success = True
            
        except Exception as e:
            print "  Parallel failed: %s" % str(e)
            par_time = float('inf')
            par_success = False
            matrix_par = None
        
        # Compare results
        if seq_success and par_success:
            speedup = seq_time / par_time if par_time > 0 else 0
            print "  Speedup: %.2fx" % speedup
            
            # Check consistency
            consistent = True
            max_diff = 0
            
            for i in range(min(10, len(matrix_seq))):
                for j in range(min(10, len(matrix_seq[i]))):
                    diff = abs(matrix_seq[i][j] - matrix_par[i][j])
                    max_diff = max(max_diff, diff)
                    if diff > 1e-6:
                        consistent = False
                        print "    Mismatch at [%d,%d]: seq=%.6f, par=%.6f, diff=%.6f" % (i, j, matrix_seq[i][j], matrix_par[i][j], diff)
                        break
                if not consistent:
                    break
            
            print "  Results consistent: %s (max diff: %.2e)" % (consistent, max_diff)
            
            if consistent and speedup > 1.5:
                print "  SUCCESS: Parallel is faster and accurate!"
            elif consistent:
                print "  Accurate but similar speed"
            else:
                print "  Speed improvement but accuracy issues"
        
        # Memory usage
        memory_mb = num_cells * num_cells * 4 / (1024.0 * 1024.0)
        print "  Memory usage: %.1f MB" % memory_mb
        
        # Stop if sequential takes too long
        if seq_success and seq_time > 300:
            print "  Sequential too slow - parallel processing needed"
            break
    
    print "\n" + "=" * 50
    print "PARALLEL PROCESSING ANALYSIS COMPLETE"

if __name__ == "__main__":
    test_parallel_vs_sequential()