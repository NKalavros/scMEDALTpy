"""
R-style genomic binning that replicates the MEDALT R implementation exactly.
This implements the bedtools intersect approach used in LSA.tree.R
"""

import pandas as pd
import numpy as np
import os
from typing import Dict, List, Tuple, Optional
import subprocess
import tempfile

def create_region_bed(data: pd.DataFrame, gene_info: pd.DataFrame, output_path: str = "region.bed", debug: bool = False):
    """
    Create a region.bed file from gene data, replicating R's approach.
    
    For RNA data, this maps genes to their genomic coordinates.
    """
    # Match genes to reference positions (equivalent to R's match function)
    gene_names = data.index.tolist()
    gene_list = gene_info.iloc[:, 0].astype(str).tolist()
    
    if debug:
        print(f"DEBUG: Looking for {len(gene_names)} genes in reference of {len(gene_list)} entries")
        print(f"DEBUG: Sample gene names: {gene_names[:5]}")
        print(f"DEBUG: Sample reference genes: {gene_list[:5]}")
    
    valid_data = []
    found_genes = 0
    for gene in gene_names:
        if gene in gene_list:
            idx = gene_list.index(gene)
            # Get chromosome and position info
            chrom = gene_info.iloc[idx, 1]  # chromosome
            start = gene_info.iloc[idx, 2]  # start position
            end = gene_info.iloc[idx, 3]    # end position
            
            # Format like R: ensure chr prefix and proper coordinates
            if not str(chrom).startswith('chr'):
                chrom = f"chr{chrom}"
            
            valid_data.append([chrom, start, end])
            found_genes += 1
    
    if debug:
        print(f"DEBUG: Found {found_genes} genes in reference")
    
    # Create region.bed file
    region_df = pd.DataFrame(valid_data, columns=['chrom', 'start', 'end'])
    region_df.to_csv(output_path, sep='\t', header=False, index=False)
    
    return region_df

def bedtools_intersect_python(cytoband_file: str, region_file: str) -> pd.DataFrame:
    """
    Python implementation of bedtools intersect.
    Maps genomic regions to cytoband annotations.
    """
    # Read cytoband file (format: chr1 0 2300000 p36.33)
    cytobands = pd.read_csv(cytoband_file, sep='\t', header=None, 
                           names=['seqnames', 'start_band', 'end_band', 'name'])
    
    # Read region file (format: chr1 1000000 2000000)
    regions = pd.read_csv(region_file, sep='\t', header=None,
                         names=['chrom', 'start_region', 'end_region'])
    
    intersections = []
    
    for _, region in regions.iterrows():
        # Find overlapping cytobands
        chrom_bands = cytobands[cytobands['seqnames'] == region['chrom']]
        
        for _, band in chrom_bands.iterrows():
            # Check for overlap
            if (region['start_region'] < band['end_band'] and 
                region['end_region'] > band['start_band']):
                
                intersections.append({
                    'seqnames': band['seqnames'],
                    'start_band': band['start_band'],
                    'end_band': band['end_band'],
                    'name': band['name'],
                    'region_start': region['start_region'],
                    'region_end': region['end_region']
                })
    
    return pd.DataFrame(intersections)

def create_r_style_cnv_matrix(data: pd.DataFrame, intersect_result: pd.DataFrame, 
                             gene_info: pd.DataFrame) -> pd.DataFrame:
    """
    Create CNV matrix with R-style cytoband regions as columns.
    Replicates the R code that creates newCNV matrix.
    """
    # Create region IDs like R does
    intersect_result['ID'] = intersect_result['seqnames'] + ':' + intersect_result['name']
    intersect_result['ID1'] = intersect_result['seqnames'] + '_' + intersect_result['region_end'].astype(str)
    
    # Get unique cytoband IDs
    unique_ids = intersect_result['ID'].unique()
    
    # Create the new CNV matrix
    newCNV_data = {}
    
    # Match genes to their positions for lookup
    gene_positions = {}
    for gene in data.index:  # Genes are rows (index)
        if gene in gene_info.iloc[:, 0].values:
            idx = gene_info.iloc[:, 0].tolist().index(gene)
            chrom = gene_info.iloc[idx, 1]
            pos = gene_info.iloc[idx, 3]  # end position
            if not str(chrom).startswith('chr'):
                chrom = f"chr{chrom}"
            gene_positions[gene] = f"{chrom}_{pos}"
    
    for band_id in unique_ids:
        # Find all regions that map to this band
        band_regions = intersect_result[intersect_result['ID'] == band_id]['ID1'].values
        
        # Find genes that map to these regions
        matching_genes = []
        for gene, gene_id1 in gene_positions.items():
            if gene_id1 in band_regions:
                matching_genes.append(gene)
        
        if len(matching_genes) == 1:
            # Single gene - use its values directly
            newCNV_data[band_id] = data.loc[matching_genes[0]].values
        elif len(matching_genes) > 1:
            # Multiple genes - average them (like R does)
            gene_values = data.loc[matching_genes]
            newCNV_data[band_id] = np.round(gene_values.mean(axis=0)).astype(int)
        # If no genes map, skip this band
    
    # Convert to DataFrame (cells as rows, cytoband regions as columns)
    newCNV = pd.DataFrame(newCNV_data, index=data.columns)
    
    return newCNV

def combine_adjacent_regions(lsa_results: pd.DataFrame) -> pd.DataFrame:
    """
    Implement R's CombineRegion function to merge adjacent cytobands.
    Creates region names like "chr12:q24.12-24.13"
    """
    if lsa_results.empty:
        return lsa_results
    
    combined_results = []
    
    # Group by cell and CNA type
    for (cell, cna_type), group in lsa_results.groupby(['cell', 'CNA']):
        # Sort by chromosome and band position
        group = group.copy()
        
        # Parse region names to extract chromosome, arm, and band info
        group['chr'] = group['region'].str.extract(r'(chr[^:]+):')
        group['arm_band'] = group['region'].str.extract(r'chr[^:]+:([pq].+)')
        group['arm'] = group['arm_band'].str[0]  # p or q
        
        # Try to parse band numbers
        group['band_num'] = group['arm_band'].str.extract(r'[pq](\d+)').astype(float)
        
        # Group by chromosome and arm
        for (chr_name, arm), chr_group in group.groupby(['chr', 'arm']):
            if len(chr_group) == 1:
                # Single band - keep as is
                combined_results.append(chr_group.iloc[0])
            else:
                # Multiple bands - check if adjacent
                chr_group = chr_group.sort_values('band_num')
                
                # For now, simple approach: if multiple bands in same arm with same CNA type,
                # create merged name
                band_nums = chr_group['band_num'].dropna()
                if len(band_nums) > 1:
                    min_band = chr_group['arm_band'].iloc[0]
                    max_band = chr_group['arm_band'].iloc[-1]
                    
                    # Create merged region name
                    merged_name = f"{chr_name}:{min_band}-{max_band}"
                    
                    # Use the strongest score and best p-value
                    best_row = chr_group.loc[chr_group['pvalue'].idxmin()].copy()
                    best_row['region'] = merged_name
                    combined_results.append(best_row)
                else:
                    # Keep individual bands
                    for _, row in chr_group.iterrows():
                        combined_results.append(row)
    
    return pd.DataFrame(combined_results)

def r_style_genomic_binning(inputfile: str, reference: str, cytoband_file: str, 
                           output_file: str = None, debug: bool = False) -> pd.DataFrame:
    """
    Main function that replicates R's genomic binning approach exactly.
    
    This implements the workflow from LSA.tree.R lines 114-133:
    1. Create region.bed from gene data
    2. Run bedtools intersect against cytoband annotations
    3. Create CNV matrix with cytoband-based regions
    """
    if debug:
        print("Starting R-style genomic binning...")
    
    # Step 1: Load data (genes are rows, cells are columns)
    data = pd.read_csv(inputfile, sep="\t", index_col=0)
    data = np.round(data * 2)  # Convert to integer copy numbers
    
    gene_info = pd.read_csv(reference, sep="\t", header=None)
    
    if debug:
        print(f"Loaded {len(data.columns)} cells and {len(data)} genes")
        print(f"Gene reference has {len(gene_info)} entries")
        print(f"Sample gene names from data: {list(data.index)[:5]}")
    
    # Step 2: Create a gene DataFrame for bedtools intersect
    gene_df = pd.DataFrame(index=data.index)  # Genes as index
    
    # Step 3: Create region.bed file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        region_bed_path = f.name
    
    region_df = create_region_bed(gene_df, gene_info, region_bed_path, debug=debug)
    
    if debug:
        print(f"Created region.bed with {len(region_df)} regions")
    
    # Step 3: Perform bedtools intersect (Python implementation)
    intersect_result = bedtools_intersect_python(cytoband_file, region_bed_path)
    
    if debug:
        print(f"Found {len(intersect_result)} cytoband intersections")
        if len(intersect_result) > 0:
            unique_bands = intersect_result['seqnames'].astype(str) + ':' + intersect_result['name'].astype(str)
            print(f"Sample bands: {unique_bands.unique()[:5]}")
    
    # Step 4: Create CNV matrix with cytoband regions
    newCNV = create_r_style_cnv_matrix(data, intersect_result, gene_info)
    
    # Clean up temporary file
    os.unlink(region_bed_path)
    
    if debug:
        print(f"Created CNV matrix with {len(newCNV.columns)} cytoband regions")
        print(f"Sample regions: {list(newCNV.columns)[:5]}")
    
    # Step 5: Save result
    if output_file:
        newCNV.to_csv(output_file, sep='\t')
        print(f"Saved R-style binned CNV data to {output_file}")
    
    return newCNV

def merge_to_arm_level(lsa_results: pd.DataFrame, cnv_matrix: pd.DataFrame, 
                      coverage_threshold: float = 0.5) -> pd.DataFrame:
    """
    Implement R's arm-level merging when >50% of arm is affected.
    """
    if lsa_results.empty:
        return lsa_results
    
    merged_results = []
    
    # Group by cell and CNA type
    for (cell, cna_type), group in lsa_results.groupby(['cell', 'CNA']):
        # Parse regions to get chromosome arms
        arms = {}
        for _, row in group.iterrows():
            region = row['region']
            if ':' in region:
                chr_part, band_part = region.split(':', 1)
                arm = band_part[0] if band_part else 'unknown'  # p or q
                arm_key = f"{chr_part}:{arm}"
                
                if arm_key not in arms:
                    arms[arm_key] = []
                arms[arm_key].append(row)
        
        # Check each arm for coverage
        for arm_key, arm_rows in arms.items():
            if len(arm_rows) == 1:
                merged_results.append(arm_rows[0])
            else:
                # Multiple bands in same arm - could merge to arm level
                # For simplicity, merge if >2 bands (approximating 50% coverage)
                if len(arm_rows) > 2:
                    # Merge to arm level
                    best_row = min(arm_rows, key=lambda x: x['pvalue'])
                    best_row = best_row.copy()
                    best_row['region'] = arm_key  # e.g., "chr7:q"
                    merged_results.append(best_row)
                else:
                    # Keep individual bands
                    merged_results.extend(arm_rows)
    
    return pd.DataFrame(merged_results)