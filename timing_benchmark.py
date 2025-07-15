#!/usr/bin/env python2
"""
Comprehensive timing benchmark for MEDALT pipeline optimization
"""

import time
import os
import sys
import subprocess
import shutil
import json

def run_command_with_timing(cmd, description):
    """Run a command and return the execution time and output"""
    print "=" * 60
    print "RUNNING: %s" % description
    print "COMMAND: %s" % cmd
    print "=" * 60
    
    start_time = time.time()
    try:
        result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        end_time = time.time()
        elapsed = end_time - start_time
        
        print "COMPLETED in %.2f seconds" % elapsed
        return elapsed, result, True
    except subprocess.CalledProcessError as e:
        end_time = time.time()
        elapsed = end_time - start_time
        
        print "FAILED after %.2f seconds" % elapsed
        print "ERROR OUTPUT:"
        print e.output
        return elapsed, e.output, False

def restore_original_files():
    """Restore original pipeline files"""
    print "Restoring original files..."
    for filename in ['scTree.py', 'ComputeDistance.py', 'Edmonds.py', 'Readfile.py', 'dataTransfer.R']:
        if os.path.exists('original_backup/%s' % filename):
            shutil.copy('original_backup/%s' % filename, filename)
    print "Original files restored."

def restore_optimized_files():
    """Restore optimized pipeline files"""
    print "Restoring optimized files..."
    subprocess.call(['git', 'stash', 'pop'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print "Optimized files restored."

def clean_output_dir(output_dir):
    """Clean output directory"""
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

def benchmark_pipeline(data_file, output_suffix, num_genes, num_cells, use_permutations=True, num_permutations=20):
    """Benchmark pipeline on given data"""
    results = {}
    
    # Test parameters
    perm_flag = "T" if use_permutations else "F"
    perm_str = "perm" if use_permutations else "noperm"
    
    # Original pipeline
    print "\\n" + "="*80
    print "TESTING ORIGINAL PIPELINE"
    print "="*80
    
    restore_original_files()
    
    output_dir = "benchmark_original_%s_%s" % (output_suffix, perm_str)
    clean_output_dir(output_dir)
    
    cmd = "python2 scTree.py -P ./ -I %s -D R -G hg19 -O %s -R %s -W 30" % (data_file, output_dir, perm_flag)
    if use_permutations:
        cmd += " -N %d" % num_permutations
    
    original_time, original_output, original_success = run_command_with_timing(cmd, "Original Pipeline")
    
    # Optimized pipeline
    print "\\n" + "="*80
    print "TESTING OPTIMIZED PIPELINE"
    print "="*80
    
    restore_optimized_files()
    
    output_dir = "benchmark_optimized_%s_%s" % (output_suffix, perm_str)
    clean_output_dir(output_dir)
    
    cmd = "python2 scTree.py -P ./ -I %s -D R -G hg19 -O %s -R %s -W 30" % (data_file, output_dir, perm_flag)
    if use_permutations:
        cmd += " -N %d" % num_permutations
    
    optimized_time, optimized_output, optimized_success = run_command_with_timing(cmd, "Optimized Pipeline")
    
    # Calculate speedup
    speedup = original_time / optimized_time if optimized_time > 0 else 0
    
    results = {
        'dataset': output_suffix,
        'num_genes': num_genes,
        'num_cells': num_cells,
        'permutations': use_permutations,
        'num_permutations': num_permutations if use_permutations else 0,
        'original_time': original_time,
        'optimized_time': optimized_time,
        'speedup': speedup,
        'original_success': original_success,
        'optimized_success': optimized_success
    }
    
    return results

def main():
    """Main benchmarking function"""
    print "MEDALT Pipeline Benchmarking"
    print "=" * 50
    
    # Test configurations
    test_configs = [
        {
            'data_file': 'synthetic_data_500g_100c.txt',
            'output_suffix': '500g_100c',
            'num_genes': 500,
            'num_cells': 100,
            'permutations': True,
            'num_permutations': 20
        },
        {
            'data_file': 'synthetic_data_500g_100c.txt',
            'output_suffix': '500g_100c',
            'num_genes': 500,
            'num_cells': 100,
            'permutations': False,
            'num_permutations': 0
        }
    ]
    
    all_results = []
    
    for config in test_configs:
        print "\\n" + "="*100
        print "BENCHMARK CONFIG: %s (permutations=%s)" % (config['output_suffix'], config['permutations'])
        print "="*100
        
        results = benchmark_pipeline(
            config['data_file'],
            config['output_suffix'],
            config['num_genes'],
            config['num_cells'],
            config['permutations'],
            config['num_permutations']
        )
        
        all_results.append(results)
        
        print "\\n" + "-"*60
        print "RESULTS FOR %s:" % config['output_suffix']
        print "  Original time: %.2f seconds" % results['original_time']
        print "  Optimized time: %.2f seconds" % results['optimized_time']
        print "  Speedup: %.2fx" % results['speedup']
        print "  Success: Original=%s, Optimized=%s" % (results['original_success'], results['optimized_success'])
        print "-"*60
    
    # Summary report
    print "\\n" + "="*80
    print "FINAL BENCHMARK RESULTS"
    print "="*80
    
    for result in all_results:
        perm_str = "with permutations" if result['permutations'] else "without permutations"
        print "\\nDataset: %d genes x %d cells (%s)" % (result['num_genes'], result['num_cells'], perm_str)
        print "  Original:  %.2f seconds" % result['original_time']
        print "  Optimized: %.2f seconds" % result['optimized_time']
        print "  Speedup:   %.2fx" % result['speedup']
        print "  Success:   Original=%s, Optimized=%s" % (result['original_success'], result['optimized_success'])
    
    # Save results to JSON
    with open('benchmark_results.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    
    print "\\nBenchmark results saved to benchmark_results.json"

if __name__ == "__main__":
    main()