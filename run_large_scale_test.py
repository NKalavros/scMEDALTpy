#!/usr/bin/env python2
"""
Non-interactive large-scale test runner
"""

import time
import os
import sys
import subprocess

def run_large_scale_test():
    """Run the large-scale test without interaction"""
    
    data_file = "synthetic_data_3000g_2000c.txt"
    
    print "MEDALT Large-Scale Performance Test"
    print "=" * 50
    
    # Check dataset
    if not os.path.exists(data_file):
        print "ERROR: %s not found!" % data_file
        return False
    
    # Analyze dataset
    with open(data_file, 'r') as f:
        first_line = f.readline().strip()
        num_cells = len(first_line.split('\t')) - 1 if first_line else 0
        
        num_genes = 0
        for line in f:
            if line.strip():
                num_genes += 1
    
    print "Dataset: %d genes x %d cells" % (num_genes, num_cells)
    
    # Estimate complexity
    distance_ops = num_cells * (num_cells - 1) / 2
    estimated_memory_gb = (num_cells * num_cells * 8) / (1024.0 ** 3)
    time_estimate_seconds = (num_cells / 200.0) ** 2 * 0.55
    
    print "Estimates:"
    print "  Distance calculations: %d" % distance_ops
    print "  Memory needed: %.2f GB" % estimated_memory_gb
    print "  Distance calc time: %.1f seconds" % time_estimate_seconds
    
    # Setup output
    output_dir = "large_scale_test_output"
    if os.path.exists(output_dir):
        import shutil
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    
    # Run pipeline
    cmd = [
        "python2", "scTree.py",
        "-P", "./",
        "-I", data_file,
        "-D", "R", 
        "-G", "hg19",
        "-O", output_dir,
        "-R", "F",  # No permutations for speed
        "-W", "30"
    ]
    
    print "\\nRunning: %s" % " ".join(cmd)
    print "\\nStarting large-scale test..."
    print "-" * 50
    
    start_time = time.time()
    
    try:
        result = subprocess.call(cmd)
        success = (result == 0)
    except Exception as e:
        print "ERROR: %s" % str(e)
        success = False
    
    end_time = time.time()
    total_time = end_time - start_time
    
    print "-" * 50
    print "\\nResults:"
    print "Success: %s" % success
    print "Total time: %.2f seconds (%.2f minutes)" % (total_time, total_time / 60.0)
    
    if success:
        # Check outputs
        output_files = ["CNV.tree.txt", "segmental.LSA.txt", "gene.LSA.txt"]
        for filename in output_files:
            filepath = os.path.join(output_dir, filename)
            if os.path.exists(filepath):
                size = os.path.getsize(filepath)
                print "%s: %d bytes" % (filename, size)
    
    return success, total_time

if __name__ == "__main__":
    success, runtime = run_large_scale_test()
    
    if success:
        print "\\nLarge-scale test PASSED!"
        print "The optimized pipeline can handle ~2000 cells!"
        
        if runtime < 600:  # Less than 10 minutes
            print "Runtime is reasonable for this scale."
        else:
            print "Runtime is long - consider further optimizations."
    else:
        print "\\nLarge-scale test FAILED!"
        print "Need to investigate bottlenecks and optimize further."