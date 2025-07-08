#!/usr/bin/env python3
"""
MEDALT Pipeline Wrapper Script

A complete pipeline wrapper for the MEDALT reimplementation that processes
single-cell CNV data and generates minimum spanning trees for lineage analysis.

Usage:
    python3 run_medalt.py <input_file> [output_file]

Example:
    python3 run_medalt.py ../MEDALT/example/outputRNA_pytest/2_scRNA.CNV_bin_30.csv results.txt
"""

import sys
import os
from datetime import datetime

def run_medalt_pipeline(input_file, output_file=None):
    """Run the complete MEDALT pipeline"""
    
    # Import MEDALT modules
    from Readfile import read
    from Edmonds import create_tree
    from mdmst import compute_rdmst
    
    print("="*60)
    print("MEDALT Python Reimplementation Pipeline")
    print("="*60)
    print(f"Input file: {input_file}")
    print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Validate input file
    if not os.path.exists(input_file):
        print(f"ERROR: Input file '{input_file}' not found!")
        return False
    
    try:
        # Step 1: Read and parse data
        print("Step 1: Reading CNV data...")
        nodes, root = read(input_file)
        print(f"✓ Successfully read {len(nodes)} cells")
        print(f"✓ Selected root cell: {root}")
        
        # Check for key cells
        if 'HNSCC5_p7_HNSCC5_P7_G05' in nodes:
            print("✓ Target cell G05 found in dataset")
        
        print()
        
        # Step 2: Create distance graph
        print("Step 2: Computing pairwise distances...")
        node_list = list(nodes.keys())
        tree_dict = create_tree(nodes, node_list, root)
        print(f"✓ Created distance matrix for {len(node_list)} cells")
        
        # Show key distances
        if root in tree_dict and 'HNSCC5_p7_HNSCC5_P7_G05' in tree_dict[root]:
            g05_dist = tree_dict[root]['HNSCC5_p7_HNSCC5_P7_G05']
            print(f"✓ Root → G05 distance: {g05_dist}")
        
        print()
        
        # Step 3: Compute minimum spanning tree
        print("Step 3: Computing minimum spanning tree...")
        result = compute_rdmst(tree_dict, root)
        
        if not result:
            print("ERROR: Failed to compute minimum spanning tree!")
            return False
        
        mst, weight = result
        print(f"✓ MST computed successfully")
        print(f"✓ Total tree weight: {weight}")
        print(f"✓ Tree contains {len(mst)} nodes")
        
        # Show root connections
        if root in mst and mst[root]:
            root_connections = list(mst[root].items())
            print(f"✓ Root connects to {len(root_connections)} direct children:")
            for cell, dist in root_connections[:3]:  # Show first 3
                print(f"   - {cell} (distance: {dist})")
            if len(root_connections) > 3:
                print(f"   ... and {len(root_connections) - 3} more")
        
        print()
        
        # Step 4: Save results
        if output_file is None:
            output_file = "medalt_results.txt"
        
        print(f"Step 4: Saving results to {output_file}...")
        
        with open(output_file, 'w') as f:
            f.write("# MEDALT Python Reimplementation Results\n")
            f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# Input: {input_file}\n")
            f.write(f"# Root: {root}\n")
            f.write(f"# Nodes: {len(nodes)}\n")
            f.write(f"# MST Weight: {weight}\n")
            f.write("# Format: from_cell\tto_cell\tdistance\n")
            f.write("#\n")
            f.write("from\tto\tdist\n")
            
            # Write tree edges
            edge_count = 0
            for from_cell in mst:
                for to_cell, distance in mst[from_cell].items():
                    f.write(f"{from_cell}\t{to_cell}\t{distance}\n")
                    edge_count += 1
            
            f.write(f"#\n")
            f.write(f"# Total edges: {edge_count}\n")
        
        print(f"✓ Results saved to {output_file}")
        print(f"✓ Wrote {edge_count} tree edges")
        
        print()
        print("="*60)
        print("PIPELINE COMPLETED SUCCESSFULLY!")
        print("="*60)
        
        # Summary statistics
        print("Summary:")
        print(f"  • Input cells: {len(nodes)}")
        print(f"  • Root cell: {root}")
        print(f"  • MST weight: {weight}")
        print(f"  • Output file: {output_file}")
        
        # Optional: Generate visualizations
        try:
            print()
            print("Generating visualizations...")
            from visualize_medalt import MEDALTVisualizer
            
            viz_dir = os.path.join(os.path.dirname(output_file), "plots")
            visualizer = MEDALTVisualizer(viz_dir)
            
            visualizer.plot_single_cell_tree(output_file)
            visualizer.plot_tree_statistics(output_file)
            
            print(f"✓ Visualizations saved to {viz_dir}/")
            
        except ImportError:
            print("Note: Install matplotlib and networkx for visualizations")
        except Exception as e:
            print(f"Note: Visualization generation failed: {e}")
        
        return True
        
    except Exception as e:
        print(f"ERROR: Pipeline failed with exception: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Main entry point"""
    
    if len(sys.argv) < 2:
        print("Usage: python3 run_medalt.py <input_file> [output_file]")
        print()
        print("Example:")
        print("  python3 run_medalt.py ../MEDALT/example/outputRNA_pytest/2_scRNA.CNV_bin_30.csv")
        print("  python3 run_medalt.py data.csv results.txt")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    success = run_medalt_pipeline(input_file, output_file)
    
    if success:
        print("\n🎉 MEDALT pipeline completed successfully!")
        sys.exit(0)
    else:
        print("\n❌ MEDALT pipeline failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()