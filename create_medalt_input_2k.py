#!/usr/bin/env python2
"""
Create proper MEDALT input file with 2k cells x 1k genes
"""

import random
import time

def create_medalt_input(num_cells=2000, num_genes=1000):
    """Create MEDALT-format input file"""
    random.seed(42)
    
    print "Creating MEDALT input file: %d cells x %d genes" % (num_cells, num_genes)
    
    output_file = "scRNA_2k_1k.CNV.txt"
    
    with open(output_file, "w") as f:
        # Write header row with cell names
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
            
            for cell_idx in range(num_cells):
                # Create realistic copy number with some structure
                if random.random() < 0.05:  # 5% deletions
                    if random.random() < 0.4:  # Some homozygous deletions
                        cn = 0.0
                    else:  # Heterozygous deletions
                        cn = 0.5
                elif random.random() < 0.7:  # 70% normal diploid
                    cn = random.uniform(0.8, 1.2)  # Some noise around 1.0
                elif random.random() < 0.2:  # 20% gains
                    cn = random.uniform(1.3, 1.8)
                else:  # 5% amplifications
                    cn = random.uniform(2.0, 3.0)
                
                row.append("%.6f" % cn)
            
            f.write("\t".join(row) + "\n")
    
    print "MEDALT input file created: %s" % output_file
    return output_file

def main():
    """Create MEDALT input file"""
    input_file = create_medalt_input(2000, 1000)
    print "Ready to run MEDALT pipeline with: %s" % input_file

if __name__ == "__main__":
    main()