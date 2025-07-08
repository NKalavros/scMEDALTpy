#!/usr/bin/env python3
"""
Create R-style binning that exactly matches the R RNAinput function
"""

import pandas as pd
import numpy as np

def r_style_rna_binning(inputfile, reference_file, delt=30):
    """
    Exactly replicate the R RNAinput function from dataTransfer.R
    """
    # Step 1: data=read.csv(inputfile,sep="\t")
    data = pd.read_csv(inputfile, sep="\t", index_col=0)
    
    # Step 2: data=round(data*2)
    data = np.round(data * 2)
    
    # Step 3: geneInfo=read.csv(reference,sep="\t",header=F)
    gene_info = pd.read_csv(reference_file, sep="\t", header=None)
    
    # Step 4-6: Match genes and create position data
    gene_names = data.index.tolist()
    gene_list = gene_info.iloc[:, 0].astype(str).tolist()
    
    # Match genes to reference
    matched_indices = []
    matched_genes = []
    for gene in gene_names:
        try:
            idx = gene_list.index(gene)
            matched_indices.append(idx)
            matched_genes.append(gene)
        except ValueError:
            continue
    
    if not matched_indices:
        raise ValueError("No genes found in reference")
    
    # Create newdata: [chr, start, end, gene_data...]
    newdata = []
    for i, gene_idx in enumerate(matched_indices):
        gene_row = gene_info.iloc[gene_idx]
        gene_name = matched_genes[i]
        gene_data = data.loc[gene_name].values
        
        # chr, start, end, then all cell values
        row = [gene_row.iloc[1], gene_row.iloc[2], gene_row.iloc[3]] + list(gene_data)
        newdata.append(row)
    
    newdata = pd.DataFrame(newdata)
    
    # Step 7-9: Process chromosome names (lines 29-35 in R)
    newdata.iloc[:, 0] = newdata.iloc[:, 0].astype(str)
    
    # Extract chromosome numbers (remove "chr" prefix)
    chromo_nums = []
    for chr_str in newdata.iloc[:, 0]:
        if chr_str.startswith('chr'):
            chrom_num = chr_str[3:]  # Remove 'chr'
        else:
            chrom_num = chr_str
        
        if chrom_num == "X":
            chrom_num = "23"
        chromo_nums.append(chrom_num)
    
    newdata.iloc[:, 0] = chromo_nums
    
    # Step 10: Filter out M and Y chromosomes (line 36)
    newdata = newdata[~newdata.iloc[:, 0].isin(["M", "Y"])]
    
    # Step 11: Get chromosome list (lines 37-38)
    chrom_list = list(range(1, 24))  # 1 to 23
    present_chroms = [int(c) for c in newdata.iloc[:, 0].unique() if c.isdigit()]
    chrom = [c for c in chrom_list if c in present_chroms]
    
    print(f"DEBUG R-style: Present chromosomes: {present_chroms}")
    print(f"DEBUG R-style: Final chromosome list: {chrom}")
    
    # Step 12-15: Binning process (lines 39-65)
    segdata = []
    chrregion = []
    
    for i in chrom:
        print(f"DEBUG R-style: Processing chromosome {i}")
        
        # Get data for this chromosome
        subdata = newdata[newdata.iloc[:, 0] == str(i)]
        
        # Sort by start position (line 44)
        subdata = subdata.sort_values(by=subdata.columns[1])
        
        # Calculate number of bins
        kk = len(subdata) / delt
        intekk = round(kk)
        
        print(f"DEBUG R-style: Chromosome {i} has {len(subdata)} genes, kk={kk}, intekk={intekk}")
        
        if intekk > 1:
            # Multiple bins
            for j in range(1, intekk):  # 1 to intekk-1
                start_idx = (j-1) * delt
                end_idx = j * delt
                bin_data = subdata.iloc[start_idx:end_idx, 3:]  # Skip chr, start, end columns
                
                # Calculate mean for each cell across genes in this bin
                bin_means = bin_data.mean(axis=0)
                segdata.append(bin_means.values)
                chrregion.append(f"{i}_{j}")
                
                print(f"DEBUG R-style: Bin {i}_{j}: genes {start_idx}-{end_idx}")
            
            # Last bin (gets remaining genes)
            start_idx = (intekk-1) * delt
            bin_data = subdata.iloc[start_idx:, 3:]
            bin_means = bin_data.mean(axis=0)
            segdata.append(bin_means.values)
            chrregion.append(f"{i}_{intekk}")
            
            print(f"DEBUG R-style: Bin {i}_{intekk}: genes {start_idx}-{len(subdata)}")
        else:
            # Single bin for whole chromosome
            bin_data = subdata.iloc[:, 3:]
            bin_means = bin_data.mean(axis=0)
            segdata.append(bin_means.values)
            chrregion.append(f"{i}_1")
            
            print(f"DEBUG R-style: Single bin {i}_1: all {len(subdata)} genes")
    
    # Create final segdata matrix
    segdata = np.array(segdata)
    
    # Create row names (line 61)
    row_names = [f"chr{region}" for region in chrregion]
    
    # Create column names (cell names)
    col_names = data.columns.tolist()
    
    # Transpose and round (line 62)
    segdata = np.round(segdata.T).astype(int)
    
    # Create DataFrame
    result_df = pd.DataFrame(segdata, index=col_names, columns=row_names)
    
    print(f"DEBUG R-style: Final matrix shape: {result_df.shape}")
    print(f"DEBUG R-style: Column names: {list(result_df.columns)}")
    
    return result_df

if __name__ == "__main__":
    # Test the R-style binning
    input_file = "/Users/nikolas/Desktop/Projects/MEDALT_new/example/scRNA.CNV.txt"
    ref_file = "/Users/nikolas/Desktop/Projects/MEDALT_new/gencode_v19_gene_pos.txt"
    
    result = r_style_rna_binning(input_file, ref_file, delt=30)
    
    # Save result
    output_file = "/Users/nikolas/Desktop/Projects/MEDALT_new/r_style_binning_test.csv"
    result.to_csv(output_file, sep='\t')
    
    print(f"R-style binning saved to: {output_file}")
    print(f"Shape: {result.shape}")
    print(f"First few columns: {list(result.columns[:5])}")