#!/usr/bin/env python3

def verify_tree_distances():
    """Verify distances from the R tree output"""
    
    from Readfile import read
    from ComputeDistance import dist
    
    # Use the working R data
    data_file = "../MEDALT/example/outputRNA_pytest/2_scRNA.CNV_bin_30.csv"
    
    print("Verifying R tree distances...")
    
    # Step 1: Read data
    nodes, root = read(data_file)
    
    # Force E11 as root to match R behavior
    root = 'HNSCC5_p5_P5_E11'
    
    # Test specific connections from R tree
    connections_to_test = [
        (root, 'HNSCC5_p7_HNSCC5_P7_G05', 2),
        ('HNSCC5_p7_HNSCC5_P7_G05', 'HNSCC5_p13_P13_E03', 1),
        ('HNSCC5_p7_HNSCC5_P7_G05', 'HNSCC5_p14_HNSCC5_P14_LN_C09', 1),
        ('HNSCC5_p7_HNSCC5_P7_G05', 'HNSCC5_p5_P5_F11', 1),
        (root, root, 0),  # Self distance should be 0
    ]
    
    for from_cell, to_cell, expected_dist in connections_to_test:
        if from_cell in nodes and to_cell in nodes:
            calculated_dist = dist(nodes[from_cell], nodes[to_cell])
            status = "✓" if calculated_dist == expected_dist else "✗"
            print(f"{status} {from_cell} → {to_cell}: expected {expected_dist}, calculated {calculated_dist}")
            
            if calculated_dist != expected_dist and expected_dist is not None:
                print(f"    {from_cell}: {nodes[from_cell]}")
                print(f"    {to_cell}: {nodes[to_cell]}")
        else:
            print(f"✗ Missing cell: {from_cell} or {to_cell}")
    
    # Let's check if all E11 distances are also wrong
    print(f"\nChecking E11 distances:")
    e11_name = 'HNSCC5_p5_P5_E11'
    if e11_name in nodes:
        e11_dist_from_root = dist(nodes[root], nodes[e11_name])
        print(f"E11 distance from root: {e11_dist_from_root}")
        
        # Check what R says E11 connects to
        # From tree: HNSCC5_p3_HNSCC5_P3_H01 → HNSCC5_p5_P5_E11 with distance 0
        h01_name = 'HNSCC5_p3_HNSCC5_P3_H01'
        if h01_name in nodes:
            h01_e11_dist = dist(nodes[h01_name], nodes[e11_name])
            print(f"H01 → E11 distance: {h01_e11_dist} (R says 0)")

if __name__ == "__main__":
    verify_tree_distances()