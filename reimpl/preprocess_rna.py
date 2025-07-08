#!/usr/bin/env python3
"""
RNA Data Preprocessing for MEDALT

Python equivalent of the R dataTransfer.R script for processing scRNA.CNV.txt files.
Converts gene-by-cell expression matrices into chromosomal segment format for MEDALT analysis.

Usage:
    python3 preprocess_rna.py <input_file> [genes_per_bin] [gene_positions_file]

Example:
    python3 preprocess_rna.py scRNA.CNV.txt 30 ../MEDALT/gencode_v19_gene_pos.txt
"""

import pandas as pd
import numpy as np
import sys
import os
from typing import Tuple, Dict, List

def load_gene_positions(gene_pos_file: str) -> pd.DataFrame:
    """
    Load gene position information from gencode file
    Format: gene_name, chromosome, start, end
    """
    try:
        gene_pos = pd.read_csv(gene_pos_file, sep='\t', header=None, 
                              names=['gene', 'chromosome', 'start', 'end'])
        
        # Clean chromosome names (remove 'chr' prefix if present)
        gene_pos['chromosome'] = gene_pos['chromosome'].astype(str)
        gene_pos['chromosome'] = gene_pos['chromosome'].str.replace('chr', '', regex=False)
        
        # Convert X to 23, remove Y and M
        gene_pos['chromosome'] = gene_pos['chromosome'].replace('X', '23')
        gene_pos = gene_pos[~gene_pos['chromosome'].isin(['Y', 'M'])].copy()
        
        # Convert chromosome to numeric
        gene_pos['chromosome'] = pd.to_numeric(gene_pos['chromosome'], errors='coerce')
        gene_pos = gene_pos.dropna(subset=['chromosome']).copy()
        gene_pos['chromosome'] = gene_pos['chromosome'].astype(int)
        
        return gene_pos
        
    except Exception as e:
        print(f"Error loading gene positions file: {e}")
        return None

def process_rna_data(input_file: str, gene_pos_file: str, genes_per_bin: int = 30) -> str:
    """
    Process scRNA.CNV.txt file into chromosomal segments
    
    Equivalent to R RNAinput function:
    1. Read gene expression data
    2. Round and scale to integer copy numbers (* 2)
    3. Match genes with chromosomal positions
    4. Group genes into chromosomal bins
    5. Average expression within each bin
    6. Output in MEDALT format
    """
    
    print(f"Processing RNA data: {input_file}")
    print(f"Genes per bin: {genes_per_bin}")
    print(f"Gene positions: {gene_pos_file}")
    
    # Step 1: Load gene expression data
    try:
        # Read the data - first row is cell names, first column is gene names
        data = pd.read_csv(input_file, sep='\t', index_col=0)
        print(f"Loaded expression data: {data.shape[0]} genes × {data.shape[1]} cells")
    except Exception as e:
        print(f"Error reading input file: {e}")
        return None
    
    # Step 2: Load gene position information
    gene_pos = load_gene_positions(gene_pos_file)
    if gene_pos is None:
        return None
    
    print(f"Loaded gene positions: {len(gene_pos)} genes")
    
    # Step 3: Round and scale expression data (equivalent to R: round(data*2))
    data = np.round(data * 2).astype(int)
    
    # Step 4: Match genes with positions
    common_genes = set(data.index) & set(gene_pos['gene'])
    print(f"Matching genes found: {len(common_genes)}")
    
    if len(common_genes) == 0:
        print("ERROR: No matching genes found between expression data and gene positions!")
        return None
    
    # Filter data to common genes
    data_filtered = data.loc[list(common_genes)].copy()
    gene_pos_filtered = gene_pos[gene_pos['gene'].isin(common_genes)].copy()
    
    # Merge data with gene positions
    gene_pos_filtered.set_index('gene', inplace=True)
    data_with_pos = data_filtered.join(gene_pos_filtered, how='inner')
    
    print(f"Data with positions: {data_with_pos.shape[0]} genes")
    
    # Step 5: Group genes into chromosomal bins
    chromosomes = sorted(data_with_pos['chromosome'].unique())
    print(f"Chromosomes found: {chromosomes}")
    
    segmented_data = []
    segment_names = []
    
    for chrom in chromosomes:
        # Get genes for this chromosome
        chrom_data = data_with_pos[data_with_pos['chromosome'] == chrom].copy()
        
        # Sort by start position
        chrom_data = chrom_data.sort_values('start')
        
        # Get expression values (exclude position columns)
        expr_cols = [col for col in chrom_data.columns if col not in ['chromosome', 'start', 'end']]
        expr_data = chrom_data[expr_cols]
        
        # Calculate number of bins
        num_genes = len(chrom_data)
        num_bins = max(1, round(num_genes / genes_per_bin))
        
        print(f"  Chromosome {chrom}: {num_genes} genes → {num_bins} bins")
        
        if num_bins > 1:
            # Create multiple bins
            for bin_idx in range(num_bins - 1):
                start_idx = bin_idx * genes_per_bin
                end_idx = (bin_idx + 1) * genes_per_bin
                
                bin_data = expr_data.iloc[start_idx:end_idx]
                bin_mean = bin_data.mean(axis=0)
                
                segmented_data.append(bin_mean)
                segment_names.append(f"chr{chrom}_{bin_idx + 1}")
            
            # Last bin (may have different size)
            last_start = (num_bins - 1) * genes_per_bin
            last_bin_data = expr_data.iloc[last_start:]
            last_bin_mean = last_bin_data.mean(axis=0)
            
            segmented_data.append(last_bin_mean)
            segment_names.append(f"chr{chrom}_{num_bins}")
            
        else:
            # Single bin for this chromosome
            bin_mean = expr_data.mean(axis=0)
            segmented_data.append(bin_mean)
            segment_names.append(f"chr{chrom}_1")
    
    # Step 6: Create final dataframe
    segmented_df = pd.DataFrame(segmented_data, index=segment_names)
    segmented_df = segmented_df.round().astype(int)  # Round to integers
    
    # Step 7: Transpose to match MEDALT expected format (cells as rows, segments as columns)
    segmented_df = segmented_df.T
    
    print(f"Final segmented data: {segmented_df.shape[0]} cells × {segmented_df.shape[1]} segments")
    
    # Step 8: Save output
    output_file = f"{input_file}.CNV.txt"
    segmented_df.to_csv(output_file, sep='\t')
    
    print(f"✓ Preprocessed data saved to: {output_file}")
    return output_file

def main():
    """Main entry point"""
    
    if len(sys.argv) < 2:
        print("Usage: python3 preprocess_rna.py <input_file> [genes_per_bin] [gene_positions_file]")
        print()
        print("Example:")
        print("  python3 preprocess_rna.py scRNA.CNV.txt")
        print("  python3 preprocess_rna.py scRNA.CNV.txt 30")
        print("  python3 preprocess_rna.py scRNA.CNV.txt 30 ../MEDALT/gencode_v19_gene_pos.txt")
        sys.exit(1)
    
    input_file = sys.argv[1]
    genes_per_bin = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    
    # Default gene positions file
    default_gene_pos = "../MEDALT/gencode_v19_gene_pos.txt"
    gene_pos_file = sys.argv[3] if len(sys.argv) > 3 else default_gene_pos
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found!")
        sys.exit(1)
    
    # Check if gene positions file exists
    if not os.path.exists(gene_pos_file):
        print(f"Error: Gene positions file '{gene_pos_file}' not found!")
        print("Please provide the path to gencode_v19_gene_pos.txt or similar file")
        sys.exit(1)
    
    print("="*60)
    print("MEDALT RNA DATA PREPROCESSING")
    print("="*60)
    print(f"Input file: {input_file}")
    print(f"Genes per bin: {genes_per_bin}")
    print(f"Gene positions: {gene_pos_file}")
    print()
    
    try:
        output_file = process_rna_data(input_file, gene_pos_file, genes_per_bin)
        
        if output_file:
            print()
            print("="*60)
            print("PREPROCESSING COMPLETED SUCCESSFULLY!")
            print("="*60)
            print(f"✓ Output file: {output_file}")
            print(f"✓ Ready for MEDALT analysis")
            print()
            print("Next steps:")
            print(f"  python3 run_medalt.py {output_file}")
            
        else:
            print("❌ Preprocessing failed!")
            sys.exit(1)
            
    except Exception as e:
        print(f"Error during preprocessing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()