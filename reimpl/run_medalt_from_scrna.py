#!/usr/bin/env python3
"""
Complete MEDALT Pipeline from scRNA.CNV.txt

An integrated pipeline that handles the complete workflow from raw scRNA.CNV.txt files
through preprocessing, analysis, and visualization.

Usage:
    python3 run_medalt_from_scrna.py <scRNA_file> [output_file] [genes_per_bin]

Example:
    python3 run_medalt_from_scrna.py scRNA.CNV.txt
    python3 run_medalt_from_scrna.py scRNA.CNV.txt results.txt 30
"""

import sys
import os
from datetime import datetime

def main():
    """Main entry point for complete scRNA MEDALT pipeline"""
    
    if len(sys.argv) < 2:
        print("Usage: python3 run_medalt_from_scrna.py <scRNA_file> [output_file] [genes_per_bin]")
        print()
        print("Example:")
        print("  python3 run_medalt_from_scrna.py scRNA.CNV.txt")
        print("  python3 run_medalt_from_scrna.py scRNA.CNV.txt results.txt")
        print("  python3 run_medalt_from_scrna.py scRNA.CNV.txt results.txt 30")
        sys.exit(1)
    
    scrna_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "medalt_scrna_results.txt"
    genes_per_bin = int(sys.argv[3]) if len(sys.argv) > 3 else 30
    
    if not os.path.exists(scrna_file):
        print(f"Error: Input file '{scrna_file}' not found!")
        sys.exit(1)
    
    print("="*70)
    print("MEDALT COMPLETE PIPELINE - FROM scRNA.CNV.txt TO RESULTS")
    print("="*70)
    print(f"Input scRNA file: {scrna_file}")
    print(f"Output file: {output_file}")
    print(f"Genes per bin: {genes_per_bin}")
    print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    try:
        # Step 1: Preprocess scRNA data
        print("STEP 1: Preprocessing scRNA data...")
        print("-" * 40)
        
        from preprocess_rna import process_rna_data
        
        # Default gene positions file
        gene_pos_file = "../MEDALT/gencode_v19_gene_pos.txt"
        if not os.path.exists(gene_pos_file):
            print(f"Error: Gene positions file '{gene_pos_file}' not found!")
            print("Please ensure the MEDALT directory with gencode_v19_gene_pos.txt is available")
            sys.exit(1)
        
        processed_file = process_rna_data(scrna_file, gene_pos_file, genes_per_bin)
        
        if not processed_file:
            print("❌ Preprocessing failed!")
            sys.exit(1)
        
        print(f"✓ Preprocessing completed: {processed_file}")
        print()
        
        # Step 2: Run MEDALT analysis
        print("STEP 2: Running MEDALT analysis...")
        print("-" * 40)
        
        from run_medalt import run_medalt_pipeline
        
        success = run_medalt_pipeline(processed_file, output_file)
        
        if not success:
            print("❌ MEDALT analysis failed!")
            sys.exit(1)
        
        print()
        
        # Step 3: Run LSA analysis
        print("STEP 3: Running Lineage Speciation Analysis (LSA)...")
        print("-" * 40)
        
        try:
            from lsa_analysis import LSAAnalyzer
            
            lsa_dir = "./lsa_results"
            analyzer = LSAAnalyzer(lsa_dir)
            
            # Run LSA with fewer permutations for speed (can be increased for production)
            lsa_results = analyzer.run_lsa_analysis(output_file, processed_file, n_permutations=500)
            
            if lsa_results is not None and len(lsa_results) > 0:
                print(f"✓ LSA analysis completed with {len(lsa_results)} significant events")
                lsa_file = os.path.join(lsa_dir, "segmental.LSA.txt")
            else:
                print("✓ LSA analysis completed (no significant events found)")
                lsa_file = None
            
        except Exception as e:
            print(f"Note: LSA analysis failed: {e}")
            lsa_file = None
        
        print()
        
        # Step 4: Generate comprehensive visualizations
        print("STEP 4: Generating comprehensive visualizations...")
        print("-" * 40)
        
        try:
            from visualize_medalt import MEDALTVisualizer
            from lsa_visualization import LSAVisualizer
            
            viz_dir = "./scrna_plots"
            visualizer = MEDALTVisualizer(viz_dir)
            
            # Generate basic visualizations
            visualizer.plot_single_cell_tree(output_file)
            visualizer.plot_distance_heatmap(output_file)
            visualizer.plot_tree_statistics(output_file)
            visualizer.create_comprehensive_report(output_file, "scRNA_MEDALT_report.pdf")
            
            print(f"✓ Basic visualizations saved to {viz_dir}/")
            
            # Generate LSA visualizations if available
            if lsa_file and os.path.exists(lsa_file):
                lsa_visualizer = LSAVisualizer(viz_dir)
                lsa_plots = lsa_visualizer.create_lsa_visualizations(lsa_file, output_file)
                if lsa_plots:
                    print(f"✓ LSA visualizations saved to {viz_dir}/")
            
        except Exception as e:
            print(f"Note: Visualization generation failed: {e}")
            print("Install matplotlib and networkx for full visualization support")
        
        print()
        print("="*70)
        print("COMPLETE PIPELINE SUCCESSFUL!")
        print("="*70)
        
        print("Generated files:")
        print(f"  ✓ Preprocessed data: {processed_file}")
        print(f"  ✓ MEDALT results: {output_file}")
        if os.path.exists("./scrna_plots"):
            print(f"  ✓ Visualizations: ./scrna_plots/")
        
        print()
        print("Pipeline summary:")
        print(f"  • Input: {scrna_file}")
        print(f"  • Genes per bin: {genes_per_bin}")
        print(f"  • Final results: {output_file}")
        
        # Show key results
        try:
            import pandas as pd
            results_df = pd.read_csv(output_file, sep='\t', comment='#')
            print(f"  • Tree edges: {len(results_df)}")
            
            # Find root from results
            from_nodes = set(results_df.iloc[:, 0])
            to_nodes = set(results_df.iloc[:, 1])
            root_nodes = from_nodes - to_nodes
            
            if root_nodes:
                root = list(root_nodes)[0]
                print(f"  • Root cell: {root}")
                
                # Check if G05 is in results
                all_nodes = from_nodes | to_nodes
                g05_cells = [node for node in all_nodes if 'G05' in str(node)]
                if g05_cells:
                    print(f"  • G05 cell found: {g05_cells[0]}")
        
        except Exception:
            pass
        
        print()
        print("🎉 Complete scRNA MEDALT analysis finished successfully!")
        
    except Exception as e:
        print(f"ERROR: Pipeline failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
