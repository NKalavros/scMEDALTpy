#!/usr/bin/env python2
"""
Simple performance profiler for large-scale MEDALT testing
"""

import time
import os
import sys
import subprocess

def estimate_dataset_complexity(data_file):
    """Estimate computational complexity for the dataset"""
    print "Analyzing dataset: %s" % data_file
    
    # Count lines and cells
    with open(data_file, 'r') as f:
        first_line = f.readline().strip()
        num_cells = len(first_line.split('\t')) - 1 if first_line else 0
        
        num_genes = 0
        for line in f:
            if line.strip():
                num_genes += 1
    
    print "Dataset characteristics:"
    print "  Genes: %d" % num_genes
    print "  Cells: %d" % num_cells
    
    # Estimate complexity
    distance_ops = num_cells * (num_cells - 1) / 2  # Symmetric matrix
    estimated_segments = num_genes / 30  # Default binning
    estimated_memory_gb = (num_cells * num_cells * 8) / (1024.0 ** 3)
    
    print "Estimated computational requirements:"
    print "  Distance calculations: %d" % distance_ops
    print "  Genomic segments after binning: ~%d" % estimated_segments
    print "  Distance matrix memory: ~%.2f GB" % estimated_memory_gb
    
    # Time estimates based on our scaling data
    time_estimate_seconds = (num_cells / 200.0) ** 2 * 0.55  # Based on 200 cells = 0.55s
    time_estimate_minutes = time_estimate_seconds / 60.0
    
    print "Estimated distance calculation time: %.1f seconds (%.1f minutes)" % (time_estimate_seconds, time_estimate_minutes)
    
    return num_genes, num_cells, estimated_memory_gb, time_estimate_minutes

def run_timed_pipeline(data_file, output_dir, use_permutations=False):
    """Run pipeline with timing"""
    
    print "\\n" + "=" * 70
    print "RUNNING MEDALT PIPELINE"
    print "=" * 70
    print "Dataset: %s" % data_file
    print "Output: %s" % output_dir
    print "Permutations: %s" % ("Yes" if use_permutations else "No")
    print "=" * 70
    
    # Clean output directory
    if os.path.exists(output_dir):
        import shutil
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    
    # Build command
    perm_flag = "T" if use_permutations else "F"
    cmd = [
        "python2", "scTree.py",
        "-P", "./",
        "-I", data_file,
        "-D", "R",
        "-G", "hg19",
        "-O", output_dir,
        "-R", perm_flag,
        "-W", "30"
    ]
    
    print "Command: %s" % " ".join(cmd)
    print "\\nStarting pipeline..."
    
    # Run with timing
    start_time = time.time()
    
    try:
        # Use subprocess with real-time output
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        
        # Print output in real-time
        for line in iter(process.stdout.readline, ''):
            if line:
                print line.rstrip()
        
        process.wait()
        result = process.returncode
        success = (result == 0)
        
    except Exception as e:
        print "ERROR: %s" % str(e)
        success = False
    
    end_time = time.time()
    total_time = end_time - start_time
    
    # Results
    print "\\n" + "=" * 70
    print "PERFORMANCE RESULTS"
    print "=" * 70
    print "Success: %s" % success
    print "Total time: %.2f seconds (%.2f minutes)" % (total_time, total_time / 60.0)
    
    if success:
        print "Output files generated successfully!"
        
        # Check output files
        expected_files = ["CNV.tree.txt", "segmental.LSA.txt", "gene.LSA.txt"]
        for filename in expected_files:
            filepath = os.path.join(output_dir, filename)
            if os.path.exists(filepath):
                size = os.path.getsize(filepath)
                print "  %s: %d bytes" % (filename, size)
            else:
                print "  %s: MISSING!" % filename
    else:
        print "Pipeline failed!"
    
    return {
        'success': success,
        'total_time': total_time,
        'output_dir': output_dir
    }

def main():
    """Main function"""
    # Target dataset
    data_file = "synthetic_data_3000g_2000c.txt"
    
    if not os.path.exists(data_file):
        print "ERROR: %s not found!" % data_file
        print "Run generate_synthetic_data.py first."
        sys.exit(1)
    
    # Analyze dataset
    num_genes, num_cells, estimated_memory_gb, time_estimate_minutes = estimate_dataset_complexity(data_file)
    
    # Warning for large datasets
    if num_cells > 1000:
        print "\\n" + "!" * 50
        print "WARNING: This is a very large dataset!"
        print "Estimated memory usage: %.2f GB" % estimated_memory_gb
        print "Estimated time: %.1f minutes (just for distance calculation)" % time_estimate_minutes
        print "Total pipeline time could be much longer with permutations."
        print "!" * 50
        
        response = raw_input("\\nContinue with the test? (y/N): ")
        if response.lower() != 'y':
            print "Test aborted. Consider testing with smaller datasets first."
            sys.exit(0)
    
    # Run pipeline without permutations first
    print "\\nRunning pipeline WITHOUT permutations (faster test)..."
    results = run_timed_pipeline(
        data_file, 
        "large_scale_output_no_perm", 
        use_permutations=False
    )
    
    if results['success']:
        print "\\nSUCCESS! Large-scale test completed."
        print "The optimized pipeline can handle 3000 genes x 2000 cells!"
        
        # Ask about permutation test
        if results['total_time'] < 1800:  # Less than 30 minutes
            response = raw_input("\\nRun with permutations? This will take much longer (y/N): ")
            if response.lower() == 'y':
                print "\\nRunning with permutations..."
                perm_results = run_timed_pipeline(
                    data_file,
                    "large_scale_output_with_perm",
                    use_permutations=True
                )
                
                if perm_results['success']:
                    print "\\nFull pipeline with permutations also successful!"
                else:
                    print "\\nPermutation test failed, but base pipeline works."
        else:
            print "\\nSkipping permutation test due to long runtime."
            print "Base pipeline took %.1f minutes." % (results['total_time'] / 60.0)
    
    else:
        print "\\nLarge-scale test failed!"
        print "This indicates we need further optimizations for this scale."

if __name__ == "__main__":
    main()