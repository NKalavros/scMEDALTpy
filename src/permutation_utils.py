import os
import numpy as np
import pandas as pd
from typing import Optional

def load_cytobands(band_file):
    bands = {}
    with open(band_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                chrom, start, end, band = parts[:4]
                if chrom not in bands:
                    bands[chrom] = []
                bands[chrom].append((int(start), int(end), band))
    return bands

def get_band_name(bands, chrom, start, end):
    if chrom not in bands:
        return f"{chrom}:{start}-{end}"
    
    overlaps = []
    for b in bands[chrom]:
        if not (end < b[0] or start > b[1]):
            overlaps.append(b)
    
    if not overlaps:
        return f"{chrom}:{start}-{end}"
    
    if len(overlaps) == 1:
        return f"{chrom}:{overlaps[0][2]}"
    else:
        # Keep individual segments, don't merge
        return f"{chrom}:{overlaps[0][2]}"

def permute_gene_cnv(cnv_df: pd.DataFrame, reference_df: pd.DataFrame, 
                     n_permutations: int = 500, bin_size: int = 30, 
                     outdir: str = "./permutation", prefix: str = "permute", 
                     band_file: str = None, seed: int = 42, debug: bool = False):
    """
    Generate permuted CNV profiles by shuffling genes within chromosomes.
    Matches the R implementation behavior.
    """
    os.makedirs(outdir, exist_ok=True)
    
    # Set random seed for reproducibility
    np.random.seed(seed)
    
    # Prepare chromosome info for each gene
    reference_df.columns = ["gene", "chr", "start", "end"]
    gene_chr = reference_df.set_index("gene")["chr"].astype(str)
    gene_chr = gene_chr.reindex(cnv_df.index)
    gene_chr = gene_chr.fillna('NA')
    
    # Load cytoband info if available
    bands = load_cytobands(band_file) if band_file else None
    
    # Get valid chromosomes
    valid_chroms = sorted(set(gene_chr) - {'NA', 'M', 'Y'}, key=lambda x: (len(x), x))
    
    if debug:
        print(f"DEBUG: Generating {n_permutations} permutations")
        print(f"DEBUG: Valid chromosomes: {valid_chroms}")
        print(f"DEBUG: Total genes: {len(cnv_df)}")
    
    # Function to bin permuted CNV
    def bin_cnv(permuted_cnv, reference, bin_size, bands=None):
        # Match genes between permuted CNV and reference
        gene_info = reference.set_index("gene")
        valid_genes = permuted_cnv.index.intersection(gene_info.index)
        gene_info = gene_info.loc[valid_genes].reset_index()
        permuted_cnv = permuted_cnv.loc[valid_genes]
        
        # Remove sex chromosomes
        gene_info = gene_info[~gene_info["chr"].isin(["M", "Y"])]
        chroms = sorted(set(gene_info["chr"]) - {"M", "Y"}, key=lambda x: (len(x), x))
        
        segdata = []
        chrregion = []
        
        for chrom in chroms:
            chrom_genes = gene_info[gene_info["chr"] == chrom]
            chrom_genes = chrom_genes.sort_values("start")
            
            if len(chrom_genes) == 0:
                continue
            
            chrom_cnv = permuted_cnv.loc[chrom_genes["gene"]]
            n_genes = len(chrom_genes)
            n_bins = int(np.ceil(n_genes / bin_size))
            
            for j in range(n_bins):
                start_idx = j * bin_size
                end_idx = min((j + 1) * bin_size, n_genes)
                
                # Get bin data
                bin_genes = chrom_genes.iloc[start_idx:end_idx]
                bin_data = chrom_cnv.loc[bin_genes["gene"]]
                bin_mean = bin_data.mean(axis=0)
                
                segdata.append(bin_mean.values)
                
                # Get region name
                bin_start = bin_genes["start"].iloc[0]
                bin_end = bin_genes["end"].iloc[-1]
                
                if bands:
                    region = get_band_name(bands, chrom, bin_start, bin_end)
                else:
                    region = f"chr{chrom}:bin{j+1}"
                
                chrregion.append(region)
        
        # Convert to dataframe
        segdata = np.vstack(segdata)
        segdata = np.round(segdata)
        segdata = pd.DataFrame(segdata, columns=permuted_cnv.columns, index=chrregion)
        return segdata.T
    
    # Generate permutations
    for perm_idx in range(1, n_permutations + 1):
        if debug and perm_idx % 100 == 0:
            print(f"DEBUG: Generated {perm_idx} permutations...")
        
        # Permute gene order within each chromosome
        permuted_idx = []
        
        for chrom in valid_chroms:
            chrom_genes = gene_chr[gene_chr == chrom].index.tolist()
            if len(chrom_genes) > 0:
                # Shuffle genes within chromosome
                shuffled = np.random.permutation(chrom_genes)
                permuted_idx.extend(shuffled)
        
        # Add any NA chromosome genes at the end (unshuffled)
        na_genes = gene_chr[gene_chr == 'NA'].index.tolist()
        permuted_idx.extend(na_genes)
        
        # Create permuted CNV matrix
        permuted_cnv = cnv_df.loc[permuted_idx]
        
        # Write gene-level permuted CNV
        gene_cnv_path = os.path.join(outdir, f"{prefix}.{perm_idx}.gene.CNV.txt")
        permuted_cnv.to_csv(gene_cnv_path, sep="\t", header=True, index=True, quoting=3)
        
        # Bin the permuted CNV
        binned_cnv = bin_cnv(permuted_cnv, reference_df, bin_size, bands)
        
        # Write binned CNV
        binned_cnv_path = os.path.join(outdir, f"{prefix}.{perm_idx}.CNV.txt")
        binned_cnv.to_csv(binned_cnv_path, sep="\t", header=True, index=True, quoting=3)
    
    if debug:
        print(f"DEBUG: Successfully generated {n_permutations} permutations")
        print(f"DEBUG: Files saved to {outdir}")
    
    return True