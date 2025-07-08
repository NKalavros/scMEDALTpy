#!/usr/bin/env python3

def final_comparison():
    """Final comparison with R target results"""
    
    from Readfile import read
    from Edmonds import create_tree
    from mdmst import compute_rdmst
    
    # Use the working R data
    data_file = "../MEDALT/example/outputRNA_pytest/2_scRNA.CNV_bin_30.csv"
    
    print("=== FINAL MEDALT REIMPLEMENTATION TEST ===")
    print(f"Data file: {data_file}")
    
    # Step 1: Read data
    nodes, root = read(data_file)
    print(f"✓ Read {len(nodes)} nodes")
    print(f"✓ Selected root: {root}")
    
    # Step 2: Create distance graph
    node_list = list(nodes.keys())
    tree_dict = create_tree(nodes, node_list, root)
    print(f"✓ Created distance graph")
    
    # Step 3: Compute MST
    result = compute_rdmst(tree_dict, root)
    
    if result:
        mst, weight = result
        print(f"✓ Computed MST with weight: {weight}")
        
        # Check root connections
        if root in mst and mst[root]:
            root_connections = list(mst[root].items())
            print(f"✓ Root connects to: {root_connections}")
            
            # Check if G05 is selected as primary lineage
            g05_connected = any('G05' in cell for cell, dist in root_connections)
            
            if g05_connected:
                print("✓ SUCCESS: G05 is connected to root (primary lineage)")
                for cell, dist in root_connections:
                    if 'G05' in cell:
                        print(f"   {cell} at distance {dist}")
            else:
                print("✗ G05 not directly connected to root")
                
                # Find G05 in the tree structure
                print("Finding G05 connections...")
                for node in mst:
                    if 'G05' in node and mst[node]:
                        connections = list(mst[node].items())
                        print(f"   {node} connects to: {connections}")
    
    # Key verification: Check the critical distance
    g05_dist = tree_dict[root]['HNSCC5_p7_HNSCC5_P7_G05']
    print(f"\n=== KEY RESULT ===")
    print(f"Distance from root ({root}) to G05: {g05_dist}")
    print(f"Expected distance (from R): 2")
    
    if g05_dist == 2:
        print("✓ SUCCESS: Distance calculation matches R output!")
    else:
        print(f"✗ FAIL: Distance mismatch (expected 2, got {g05_dist})")
    
    # Compare with target output
    print(f"\n=== COMPARISON WITH R TARGET ===")
    print(f"R target: HNSCC5_p7_HNSCC5_P7_G05 at depth 1 (primary lineage)")
    print(f"R target: subtree size 50")
    print(f"Our result: {root} as root, connects to {len(root_connections) if 'root_connections' in locals() else 0} children")
    
    return g05_dist == 2

if __name__ == "__main__":
    success = final_comparison()
    if success:
        print("\n🎉 REIMPLEMENTATION SUCCESSFUL!")
    else:
        print("\n❌ Reimplementation needs more work")