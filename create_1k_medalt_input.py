#!/usr/bin/env python2
"""
Create 1k cells x 1k genes MEDALT input file
"""

import random
import time

def create_1k_medalt_input():
    """Create 1k x 1k MEDALT-format input file"""
    random.seed(42)
    
    num_cells = 1000
    num_genes = 1000
    
    print "Creating MEDALT input file: %d cells x %d genes" % (num_cells, num_genes)
    
    output_file = "scRNA_1k_1k.CNV.txt"
    
    with open(output_file, "w") as f:
        # Write header row with cell names (no gene name column)
        cell_names = []
        for i in range(num_cells):
            cell_names.append("Cell_%04d" % i)
        
        f.write("\t".join(cell_names) + "\n")
        
        # Write gene rows with copy number data
        for gene_idx in range(num_genes):
            if gene_idx % 100 == 0:
                print "  Writing gene %d/%d" % (gene_idx, num_genes)
            
            gene_name = "Gene_%04d" % gene_idx
            row = [gene_name]
            
            # Create realistic copy number patterns with some biological structure
            for cell_idx in range(num_cells):
                # Add some cellular heterogeneity - some cells have more alterations
                if cell_idx < 200:  # First 200 cells: mostly normal
                    if random.random() < 0.02:  # 2% deletions
                        cn = random.uniform(0.2, 0.6)
                    elif random.random() < 0.85:  # 85% normal
                        cn = random.uniform(0.9, 1.1)
                    else:  # 13% gains
                        cn = random.uniform(1.2, 1.6)
                        
                elif cell_idx < 600:  # Middle 400 cells: intermediate alterations
                    if random.random() < 0.05:  # 5% deletions
                        cn = random.uniform(0.1, 0.7)
                    elif random.random() < 0.75:  # 75% normal
                        cn = random.uniform(0.8, 1.2)
                    elif random.random() < 0.15:  # 15% gains
                        cn = random.uniform(1.3, 1.8)
                    else:  # 5% amplifications
                        cn = random.uniform(1.9, 2.5)
                        
                else:  # Last 400 cells: highly altered
                    if random.random() < 0.08:  # 8% deletions
                        cn = random.uniform(0.05, 0.6)
                    elif random.random() < 0.65:  # 65% normal
                        cn = random.uniform(0.7, 1.3)
                    elif random.random() < 0.20:  # 20% gains
                        cn = random.uniform(1.4, 2.0)
                    else:  # 7% high amplifications
                        cn = random.uniform(2.1, 3.0)
                
                # Add some chromosomal-level structure
                if gene_idx < 200:  # First 200 genes: chr1-like
                    if cell_idx > 500 and random.random() < 0.1:  # Some cells have chr1 gain
                        cn = cn * 1.3
                elif gene_idx < 400:  # Next 200 genes: chr2-like  
                    if cell_idx > 700 and random.random() < 0.08:  # Some cells have chr2 loss
                        cn = cn * 0.7
                elif gene_idx < 600:  # Next 200 genes: chr3-like
                    if cell_idx > 800 and random.random() < 0.12:  # Some cells have chr3 amplification
                        cn = cn * 1.5
                
                # Ensure reasonable bounds
                cn = max(0.01, min(4.0, cn))
                row.append("%.6f" % cn)
            
            f.write("\t".join(row) + "\n")
    
    print "MEDALT input file created: %s" % output_file
    print "File size: %.1f MB" % (os.path.getsize(output_file) / (1024.0 * 1024.0))
    return output_file

def main():
    """Create 1k x 1k MEDALT input file"""
    import os
    
    start_time = time.time()
    input_file = create_1k_medalt_input()
    creation_time = time.time() - start_time
    
    print "File creation completed in %.1f seconds" % creation_time
    print "Ready to run MEDALT pipeline with: %s" % input_file

if __name__ == "__main__":
    main()