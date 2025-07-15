#!/usr/bin/env python2
"""
Create 1k cells x 1k genes MEDALT input using real gene names
"""

import random

def create_1k_real_genes_input():
    """Create 1k x 1k MEDALT input with real gene names"""
    random.seed(42)
    
    num_cells = 1000
    num_genes = 1000
    
    print "Creating MEDALT input: %d cells x %d real genes" % (num_cells, num_genes)
    
    # Read real gene names from gencode file
    print "Reading gene names from gencode_v19_gene_pos.txt..."
    gene_names = []
    with open("gencode_v19_gene_pos.txt", "r") as f:
        for i, line in enumerate(f):
            if i >= num_genes:  # Take first 1000 genes
                break
            parts = line.strip().split("\t")
            gene_name = parts[0]
            gene_names.append(gene_name)
    
    print "Selected %d gene names" % len(gene_names)
    
    output_file = "scRNA_1k_real_genes.CNV.txt"
    
    with open(output_file, "w") as f:
        # Write header row with cell names
        cell_names = []
        for i in range(num_cells):
            cell_names.append("Cell_%04d" % i)
        
        f.write("\t".join(cell_names) + "\n")
        
        # Write gene rows with copy number data
        for gene_idx, gene_name in enumerate(gene_names):
            if gene_idx % 100 == 0:
                print "  Writing gene %d/%d (%s)" % (gene_idx, len(gene_names), gene_name)
            
            row = [gene_name]
            
            # Create realistic copy number patterns with cellular heterogeneity
            for cell_idx in range(num_cells):
                # Create three populations of cells with different alteration levels
                if cell_idx < 300:  # Population 1: mostly normal (diploid-like)
                    if random.random() < 0.03:  # 3% deletions
                        cn = random.uniform(0.3, 0.7)
                    elif random.random() < 0.85:  # 85% normal
                        cn = random.uniform(0.85, 1.15)
                    else:  # 12% gains
                        cn = random.uniform(1.2, 1.6)
                        
                elif cell_idx < 700:  # Population 2: intermediate alterations
                    if random.random() < 0.06:  # 6% deletions
                        cn = random.uniform(0.2, 0.8)
                    elif random.random() < 0.70:  # 70% normal
                        cn = random.uniform(0.75, 1.25)
                    elif random.random() < 0.18:  # 18% gains
                        cn = random.uniform(1.3, 1.9)
                    else:  # 6% amplifications
                        cn = random.uniform(2.0, 2.8)
                        
                else:  # Population 3: highly altered (aneuploid)
                    if random.random() < 0.10:  # 10% deletions
                        cn = random.uniform(0.1, 0.6)
                    elif random.random() < 0.60:  # 60% normal
                        cn = random.uniform(0.7, 1.3)
                    elif random.random() < 0.22:  # 22% gains
                        cn = random.uniform(1.4, 2.2)
                    else:  # 8% high amplifications
                        cn = random.uniform(2.3, 3.5)
                
                # Add some gene-specific patterns (simulate chromosomal alterations)
                if gene_idx < 100:  # First 100 genes: chr1-like region
                    if cell_idx > 600 and random.random() < 0.15:  # Some cells have this region amplified
                        cn = cn * 1.4
                elif gene_idx < 200:  # Next 100 genes: chr2-like region
                    if cell_idx > 800 and random.random() < 0.12:  # Some cells have this region deleted
                        cn = cn * 0.6
                elif gene_idx < 300:  # Next 100 genes: chr3-like region  
                    if cell_idx > 500 and cell_idx < 800 and random.random() < 0.18:  # Mid-population amplification
                        cn = cn * 1.6
                
                # Ensure reasonable bounds
                cn = max(0.05, min(4.0, cn))
                row.append("%.6f" % cn)
            
            f.write("\t".join(row) + "\n")
    
    print "MEDALT input file created: %s" % output_file
    
    # Show file info
    import os
    file_size_mb = os.path.getsize(output_file) / (1024.0 * 1024.0)
    print "File size: %.1f MB" % file_size_mb
    
    return output_file

def main():
    """Create 1k x 1k MEDALT input with real genes"""
    import time
    
    start_time = time.time()
    input_file = create_1k_real_genes_input()
    creation_time = time.time() - start_time
    
    print "File creation completed in %.1f seconds" % creation_time
    print "Ready to run MEDALT pipeline with real gene names!"
    print "Command: python2 scTree.py -P ./ -I ./%s -O ./medalt_1k_real_output -D R -G hg19" % input_file

if __name__ == "__main__":
    main()