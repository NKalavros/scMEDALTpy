#!/usr/bin/env python2
"""
Large-scale performance profiler for MEDALT pipeline
"""

import time
import os
import sys
import psutil
import subprocess
import threading

class PerformanceProfiler(object):
    def __init__(self):
        self.start_time = None
        self.memory_usage = []
        self.monitoring = False
        self.monitor_thread = None
    
    def start_monitoring(self):
        """Start performance monitoring"""
        self.start_time = time.time()
        self.memory_usage = []
        self.monitoring = True
        
        # Start memory monitoring thread
        self.monitor_thread = threading.Thread(target=self._monitor_memory)
        self.monitor_thread.daemon = True
        self.monitor_thread.start()
        
        print "Performance monitoring started..."
    
    def stop_monitoring(self):
        """Stop performance monitoring"""
        self.monitoring = False
        if self.monitor_thread:
            self.monitor_thread.join()
        
        end_time = time.time()
        total_time = end_time - self.start_time
        
        print "Performance monitoring stopped."
        return total_time
    
    def _monitor_memory(self):
        """Monitor memory usage in background thread"""
        process = psutil.Process(os.getpid())
        
        while self.monitoring:
            try:
                memory_info = process.memory_info()
                memory_mb = memory_info.rss / 1024.0 / 1024.0
                self.memory_usage.append(memory_mb)
                time.sleep(1)  # Sample every second
            except:
                break
    
    def get_peak_memory(self):
        """Get peak memory usage in MB"""
        return max(self.memory_usage) if self.memory_usage else 0
    
    def get_average_memory(self):
        """Get average memory usage in MB"""
        return sum(self.memory_usage) / len(self.memory_usage) if self.memory_usage else 0

def run_profiled_pipeline(data_file, output_dir, use_permutations=False, num_permutations=10):
    """Run pipeline with performance profiling"""
    
    print "=" * 70
    print "LARGE-SCALE MEDALT PERFORMANCE TEST"
    print "=" * 70
    print "Dataset: %s" % data_file
    print "Output: %s" % output_dir
    print "Permutations: %s (%d)" % (use_permutations, num_permutations if use_permutations else 0)
    print "=" * 70
    
    # Clean output directory
    if os.path.exists(output_dir):
        import shutil
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    
    # Start profiling
    profiler = PerformanceProfiler()
    profiler.start_monitoring()
    
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
    
    if use_permutations and hasattr(sys, 'argv'):
        # Add permutation count if the option exists
        try:
            cmd.extend(["-N", str(num_permutations)])
        except:
            pass
    
    print "Running command: %s" % " ".join(cmd)
    print ""
    
    # Run the pipeline
    try:
        result = subprocess.call(cmd)
        success = (result == 0)
    except Exception as e:
        print "ERROR: %s" % str(e)
        success = False
    
    # Stop profiling
    total_time = profiler.stop_monitoring()
    peak_memory = profiler.get_peak_memory()
    avg_memory = profiler.get_average_memory()
    
    # Results
    print ""
    print "=" * 70
    print "PERFORMANCE RESULTS"
    print "=" * 70
    print "Success: %s" % success
    print "Total time: %.2f seconds (%.2f minutes)" % (total_time, total_time / 60.0)
    print "Peak memory: %.1f MB (%.2f GB)" % (peak_memory, peak_memory / 1024.0)
    print "Average memory: %.1f MB" % avg_memory
    
    if success:
        print "Output files generated successfully!"
    else:
        print "Pipeline failed - check logs above"
    
    return {
        'success': success,
        'total_time': total_time,
        'peak_memory': peak_memory,
        'avg_memory': avg_memory
    }

def estimate_complexity(data_file):
    """Estimate the computational complexity for the dataset"""
    
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
    
    print "Estimated computational load:"
    print "  Distance calculations: %d" % distance_ops
    print "  Genomic segments: ~%d" % estimated_segments
    print "  Memory for distance matrix: ~%.1f MB" % (num_cells * num_cells * 8 / 1024.0 / 1024.0)
    
    return num_genes, num_cells

if __name__ == "__main__":
    # Test the large dataset
    data_file = "synthetic_data_3000g_2000c.txt"
    
    if not os.path.exists(data_file):
        print "ERROR: %s not found!" % data_file
        sys.exit(1)
    
    print "Analyzing dataset complexity..."
    num_genes, num_cells = estimate_complexity(data_file)
    
    if num_cells > 500:
        print ""
        print "WARNING: This is a very large dataset!"
        print "Estimated runtime could be several hours."
        print "Consider starting with smaller datasets first."
        
        response = raw_input("Continue? (y/N): ")
        if response.lower() != 'y':
            print "Aborting."
            sys.exit(0)
    
    # Run without permutations first
    print ""
    print "Testing without permutations (faster)..."
    results_no_perm = run_profiled_pipeline(
        data_file, 
        "large_scale_output_no_perm", 
        use_permutations=False
    )
    
    if results_no_perm['success']:
        print ""
        print "No-permutation test successful!"
        print "Consider running with permutations if time allows."
    else:
        print ""
        print "No-permutation test failed. Check for issues before proceeding."