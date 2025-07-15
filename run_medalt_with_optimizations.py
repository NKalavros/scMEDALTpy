#!/usr/bin/env python2
"""
Run MEDALT pipeline with both original and optimized versions for comparison
"""

import time
import os
import subprocess
import shutil

def backup_original():
    """Backup original ComputeDistance.py"""
    print "Backing up original ComputeDistance.py..."
    if os.path.exists("ComputeDistance_original_backup.py"):
        print "Backup already exists"
    else:
        shutil.copy("ComputeDistance.py", "ComputeDistance_original_backup.py")
        print "Backup created"

def restore_original():
    """Restore original ComputeDistance.py"""
    print "Restoring original ComputeDistance.py..."
    if os.path.exists("original_backup/ComputeDistance.py"):
        shutil.copy("original_backup/ComputeDistance.py", "ComputeDistance.py")
        print "Original restored"
    else:
        print "Original backup not found"

def run_medalt_original():
    """Run MEDALT with original distance calculation"""
    print "Running MEDALT with ORIGINAL distance calculation..."
    
    # Restore original
    restore_original()
    
    # Clean up any previous output
    if os.path.exists("example/medalt_original_output"):
        shutil.rmtree("example/medalt_original_output")
    
    start_time = time.time()
    
    # Run MEDALT
    cmd = ["python2", "scTree.py", "-P", "./", "-I", "./example/scRNA.CNV.txt", 
           "-O", "./example/medalt_original_output", "-D", "R", "-G", "hg19"]
    
    result = subprocess.call(cmd)
    
    elapsed = time.time() - start_time
    
    print "Original MEDALT completed in %.1f seconds" % elapsed
    
    if result == 0:
        print "SUCCESS: Original MEDALT ran successfully"
        return elapsed, True
    else:
        print "ERROR: Original MEDALT failed"
        return elapsed, False

def run_medalt_optimized():
    """Run MEDALT with optimized distance calculation"""
    print "Running MEDALT with OPTIMIZED distance calculation..."
    
    # Our current ComputeDistance.py is the optimized version
    # Clean up any previous output
    if os.path.exists("example/medalt_optimized_output"):
        shutil.rmtree("example/medalt_optimized_output")
    
    start_time = time.time()
    
    # Run MEDALT
    cmd = ["python2", "scTree.py", "-P", "./", "-I", "./example/scRNA.CNV.txt", 
           "-O", "./example/medalt_optimized_output", "-D", "R", "-G", "hg19"]
    
    result = subprocess.call(cmd)
    
    elapsed = time.time() - start_time
    
    print "Optimized MEDALT completed in %.1f seconds" % elapsed
    
    if result == 0:
        print "SUCCESS: Optimized MEDALT ran successfully"
        return elapsed, True
    else:
        print "ERROR: Optimized MEDALT failed"
        return elapsed, False

def compare_outputs():
    """Compare outputs between original and optimized versions"""
    print "Comparing outputs between original and optimized versions..."
    
    orig_dir = "example/medalt_original_output"
    opt_dir = "example/medalt_optimized_output"
    
    if not os.path.exists(orig_dir) or not os.path.exists(opt_dir):
        print "One or both output directories missing - cannot compare"
        return False
    
    # Compare CNV.tree.txt files
    orig_tree = os.path.join(orig_dir, "CNV.tree.txt")
    opt_tree = os.path.join(opt_dir, "CNV.tree.txt")
    
    if os.path.exists(orig_tree) and os.path.exists(opt_tree):
        # Read both files and compare
        with open(orig_tree, "r") as f:
            orig_lines = f.readlines()
        
        with open(opt_tree, "r") as f:
            opt_lines = f.readlines()
        
        if len(orig_lines) != len(opt_lines):
            print "WARNING: Different number of lines in tree files"
            return False
        
        differences = 0
        for i, (orig_line, opt_line) in enumerate(zip(orig_lines, opt_lines)):
            if orig_line.strip() != opt_line.strip():
                differences += 1
                if differences <= 5:  # Show first 5 differences
                    print "Line %d differs:" % (i+1)
                    print "  Original: %s" % orig_line.strip()
                    print "  Optimized: %s" % opt_line.strip()
        
        if differences == 0:
            print "SUCCESS: CNV.tree.txt files are identical!"
            return True
        else:
            print "WARNING: Found %d differences in tree files" % differences
            return False
    else:
        print "Tree files missing - cannot compare"
        return False

def main():
    """Main comparison function"""
    print "MEDALT Pipeline: Original vs Optimized Comparison"
    print "=" * 60
    
    # Backup current version
    backup_original()
    
    # Test 1: Run with original distance calculation
    print "\nTest 1: Original Distance Calculation"
    print "-" * 40
    orig_time, orig_success = run_medalt_original()
    
    # Test 2: Run with optimized distance calculation  
    print "\nTest 2: Optimized Distance Calculation"
    print "-" * 40
    opt_time, opt_success = run_medalt_optimized()
    
    # Compare results
    print "\nComparison Results"
    print "-" * 40
    
    if orig_success and opt_success:
        speedup = orig_time / opt_time if opt_time > 0 else 0
        print "Original time: %.1f seconds" % orig_time
        print "Optimized time: %.1f seconds" % opt_time
        print "Speedup: %.2fx" % speedup
        
        # Compare outputs
        outputs_match = compare_outputs()
        
        print "\nSUMMARY:"
        print "Speedup achieved: %.2fx" % speedup
        print "Output consistency: %s" % ("PASS" if outputs_match else "FAIL")
        
        if outputs_match and speedup > 1.1:
            print "SUCCESS: Optimization provides speedup with consistent results!"
        elif outputs_match:
            print "SUCCESS: Results consistent (minimal speedup on small dataset)"
        else:
            print "WARNING: Results differ - need to investigate"
            
    else:
        print "ERROR: One or both runs failed"
        if not orig_success:
            print "Original run failed"
        if not opt_success:
            print "Optimized run failed"

if __name__ == "__main__":
    main()