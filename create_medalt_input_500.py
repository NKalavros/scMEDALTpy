#!/usr/bin/env python2
"""
Create smaller MEDALT input file first (500 cells x 500 genes) to test the pipeline
"""

import random

def create_medalt_input_500():
    """Create smaller MEDALT-format input file"""
    random.seed(42)
    
    num_cells = 500
    num_genes = 500
    
    print "Creating MEDALT input file: %d cells x %d genes" % (num_cells, num_genes)
    
    output_file = "scRNA_500_500.CNV.txt"
    
    with open(output_file, "w") as f:
        # Write header row with cell names (similar to example format)
        cell_names = []
        for i in range(num_cells):
            cell_names.append("Cell_%04d" % i)
        
        f.write("\t".join(cell_names) + "\n")
        
        # Write gene rows with copy number data
        for gene_idx in range(num_genes):
            if gene_idx % 50 == 0:
                print "  Writing gene %d/%d" % (gene_idx, num_genes)
            
            gene_name = "Gene_%04d" % gene_idx
            row = [gene_name]
            
            for cell_idx in range(num_cells):
                # Create realistic copy number following example pattern
                if random.random() < 0.03:  # 3% low values (deletions)
                    cn = random.uniform(0.2, 0.7)
                elif random.random() < 0.75:  # 75% normal around 1.0
                    cn = random.uniform(0.8, 1.2)
                elif random.random() < 0.15:  # 15% moderate gains
                    cn = random.uniform(1.3, 1.8)
                else:  # 7% higher gains
                    cn = random.uniform(1.9, 2.5)
                
                row.append("%.6f" % cn)
            
            f.write("\t".join(row) + "\n")
    
    print "MEDALT input file created: %s" % output_file
    return output_file

def main():
    """Create smaller MEDALT input file"""
    input_file = create_medalt_input_500()
    print "Ready to run MEDALT pipeline with: %s" % input_file

if __name__ == "__main__":
    main()