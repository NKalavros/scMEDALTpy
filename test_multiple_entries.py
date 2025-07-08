#!/usr/bin/env python3
"""
Test script for multiple R-style entries functionality
"""

import sys
import os
sys.path.append('/Users/nikolas/Desktop/Projects/MEDALT_new')

import pandas as pd
from src.region_mapping import create_multiple_r_entries

def test_multiple_entries():
    # Load existing results
    lsa_file = 'example/outputRNA_fixed/segmental.LSA.txt'
    if not os.path.exists(lsa_file):
        print(f"Error: {lsa_file} not found")
        return
    
    print("Testing multiple R-style entries functionality...")
    
    # Read current results
    df = pd.read_csv(lsa_file, sep='\t')
    print(f"Original results: {len(df)} entries")
    print("Original cells:", df['cell'].tolist())
    
    # Apply multiple entries function
    enhanced_df = create_multiple_r_entries(df)
    print(f"Enhanced results: {len(enhanced_df)} entries")
    
    # Show entries per cell
    cell_counts = enhanced_df.groupby('cell').size().to_dict()
    print("Entries per cell:")
    for cell, count in cell_counts.items():
        print(f"  {cell}: {count}")
    
    # Show all entries
    print("\nAll enhanced entries:")
    for _, row in enhanced_df.iterrows():
        print(f"  {row['cell']}: {row['region']} ({row['CNA']}, Score={row['Score']:.2f})")
    
    # Save enhanced results for comparison
    output_file = 'example/test_multiple_entries/segmental.LSA.txt'
    enhanced_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nEnhanced results saved to: {output_file}")

if __name__ == "__main__":
    test_multiple_entries()