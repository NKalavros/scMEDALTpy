#!/usr/bin/env python2
"""
Generate synthetic scRNA-seq data for MEDALT scalability testing
"""

import random
import os
import sys

def generate_synthetic_scrna_data(num_genes=500, num_cells=100, output_file="synthetic_data.txt"):
    """
    Generate realistic synthetic scRNA-seq copy number data
    
    Args:
        num_genes: Number of genes to generate
        num_cells: Number of cells to generate
        output_file: Output filename
    """
    
    # Set random seed for reproducibility
    random.seed(42)
    
    # Generate cell names
    cell_names = ["Cell_%03d" % i for i in range(1, num_cells + 1)]
    
    # Generate gene names (using pattern similar to real data)
    gene_names = ["GENE_%03d" % i for i in range(1, num_genes + 1)]
    
    # Generate realistic copy number data
    # Most cells should be near diploid (2), with some variation
    data = []
    
    for gene_idx, gene_name in enumerate(gene_names):
        row = [gene_name]
        
        for cell_idx, cell_name in enumerate(cell_names):
            # Base copy number around 1.0 (diploid relative to normal)
            base_cn = 1.0
            
            # Add some chromosomal-level variation
            # Create some regions with gains/losses
            if gene_idx < num_genes * 0.2:  # First 20% of genes - some amplification
                if cell_idx < num_cells * 0.3:  # 30% of cells have amplification
                    base_cn = random.uniform(1.2, 1.8)
                else:
                    base_cn = random.uniform(0.8, 1.2)
            elif gene_idx < num_genes * 0.4:  # Next 20% - some deletion
                if cell_idx < num_cells * 0.2:  # 20% of cells have deletion
                    base_cn = random.uniform(0.3, 0.8)
                else:
                    base_cn = random.uniform(0.9, 1.1)
            else:  # Remaining genes - mostly diploid with noise
                base_cn = random.uniform(0.8, 1.2)
            
            # Add some noise
            noise = random.uniform(-0.1, 0.1)
            final_cn = max(0.1, base_cn + noise)  # Ensure positive values
            
            row.append("%.3f" % final_cn)
        
        data.append(row)
    
    # Write to file
    with open(output_file, 'w') as f:
        # Write header
        f.write("\t".join([""] + cell_names) + "\n")
        
        # Write data
        for row in data:
            f.write("\t".join(row) + "\n")
    
    print "Generated synthetic data: %d genes x %d cells" % (num_genes, num_cells)
    print "Output file: %s" % output_file
    
    return output_file

if __name__ == "__main__":
    # Generate different sizes for testing
    sizes = [
        (500, 100),   # Target test size
        (1000, 200),  # Medium size
        (2000, 500),  # Large size
    ]
    
    for num_genes, num_cells in sizes:
        output_file = "synthetic_data_%dg_%dc.txt" % (num_genes, num_cells)
        generate_synthetic_scrna_data(num_genes, num_cells, output_file)