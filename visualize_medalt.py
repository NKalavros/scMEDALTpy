#!/usr/bin/env python3
"""
MEDALT Visualization Tool
Standalone script for creating visualizations from MEDALT results
"""

import argparse
import os
import sys
from src.visualization import create_medalt_visualizations

def main():
    parser = argparse.ArgumentParser(
        description="Create visualizations for MEDALT analysis results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "tree_file", 
        help="Path to the CNV tree file (CNV.tree.txt)"
    )
    
    parser.add_argument(
        "--lsa-file", 
        help="Path to the LSA results file (segmental.LSA.txt)"
    )
    
    parser.add_argument(
        "--output-path", 
        default=".", 
        help="Output directory for visualization files"
    )
    
    parser.add_argument(
        "--format", 
        choices=['pdf', 'png', 'svg'], 
        default='pdf',
        help="Output format for visualizations"
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.tree_file):
        print(f"Error: Tree file not found: {args.tree_file}")
        sys.exit(1)
    
    if args.lsa_file and not os.path.exists(args.lsa_file):
        print(f"Error: LSA file not found: {args.lsa_file}")
        sys.exit(1)
    
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)
        print(f"Created output directory: {args.output_path}")
    
    # Create visualizations
    print("Creating MEDALT visualizations...")
    print(f"Tree file: {args.tree_file}")
    print(f"LSA file: {args.lsa_file or 'None'}")
    print(f"Output path: {args.output_path}")
    
    try:
        create_medalt_visualizations(
            tree_file=args.tree_file,
            lsa_file=args.lsa_file,
            output_path=args.output_path
        )
        print("\nVisualization complete!")
        print(f"Check {args.output_path} for generated files:")
        print("- singlecell_tree.pdf: Single cell phylogenetic tree")
        if args.lsa_file:
            print("- LSA_tree.pdf: LSA tree with genomic annotations")
        print("- MEDALT_visualization_report.pdf: Comprehensive report")
        
    except Exception as e:
        print(f"Error creating visualizations: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()