#!/usr/bin/env python2
"""
Test the optimized RDMST implementation with 1k x 1k real genes dataset
"""

import os
import sys
import time
from datetime import datetime as dt_

def test_rdmst_optimized():
    """Test the optimized RDMST implementation with 1k x 1k real genes"""
    
    print "Testing optimized RDMST with 1k x 1k real genes dataset..."
    
    # Check if input file exists
    input_file = "scRNA_1k_real_genes.CNV.txt"
    if not os.path.exists(input_file):
        print "Error: Input file %s not found" % input_file
        return False
    
    # Create output directory
    output_dir = "medalt_1k_rdmst_optimized"
    if os.path.exists(output_dir):
        os.system("rm -rf %s" % output_dir)
    os.makedirs(output_dir)
    
    # Run MEDALT with optimized rdmst
    start_time = time.time()
    
    print "Running MEDALT with optimized RDMST..."
    
    command = "python2 scTree.py -P ./ -I %s -O %s -D R -G hg19" % (input_file, output_dir)
    print "Command: %s" % command
    
    result = os.system(command)
    
    end_time = time.time()
    elapsed = end_time - start_time
    
    print ""
    print "="*60
    print "RDMST OPTIMIZATION TEST RESULTS"
    print "="*60
    print "Dataset: 1k cells x 1k real genes"
    print "Total time: %.2f seconds (%.2f minutes)" % (elapsed, elapsed/60.0)
    print "Exit code: %d" % result
    print ""
    
    # Check outputs
    tree_file = os.path.join(output_dir, "CNV.tree.txt")
    if os.path.exists(tree_file):
        print "Tree file created successfully: %s" % tree_file
        
        # Count lines in tree file
        with open(tree_file, "r") as f:
            lines = f.readlines()
            print "Tree file contains %d lines" % len(lines)
            
        # Show first few lines
        print "First 5 lines of tree file:"
        for i, line in enumerate(lines[:5]):
            print "  %d: %s" % (i+1, line.strip())
            
    else:
        print "Error: Tree file not created"
        return False
    
    # Check for other outputs
    gene_lsa_file = os.path.join(output_dir, "gene.LSA.txt")
    if os.path.exists(gene_lsa_file):
        print "Gene LSA file created successfully"
    
    segmental_lsa_file = os.path.join(output_dir, "segmental.LSA.txt")  
    if os.path.exists(segmental_lsa_file):
        print "Segmental LSA file created successfully"
    
    pdf_file = os.path.join(output_dir, "singlecell.tree.pdf")
    if os.path.exists(pdf_file):
        print "PDF visualization created successfully"
    
    print ""
    print "Optimized RDMST test completed successfully!"
    
    return True

if __name__ == "__main__":
    success = test_rdmst_optimized()
    sys.exit(0 if success else 1)