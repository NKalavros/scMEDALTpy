#!/usr/bin/env python3

import sys
import os

# Add MEDALT path to use original modules
sys.path.insert(0, "../MEDALT")

# Import original MEDALT modules
from Readfile import read
from Edmonds import create_tree
from mdmst import compute_rdmst

def test_original():
    """Test original MEDALT code"""
    
    # Use the working R data
    data_file = "../MEDALT/example/outputRNA_pytest/2_scRNA.CNV_bin_30.csv"
    
    print("Testing original MEDALT...")
    print(f"Reading data from: {data_file}")
    
    # Step 1: Read data
    nodes, root = read(data_file)
    print(f"Read {len(nodes)} nodes, root = {root}")
    
    # Check specific cells
    if 'HNSCC5_p7_HNSCC5_P7_G05' in nodes:
        print("G05 node:", nodes['HNSCC5_p7_HNSCC5_P7_G05'])
    if root in nodes:
        print("Root node:", nodes[root])
    
    # Step 2: Create distance graph
    node_list = list(nodes.keys())
    print(f"Creating tree with {len(node_list)} nodes...")
    tree_dict = create_tree(nodes, node_list, root)
    
    # Check root distances
    if root in tree_dict:
        print("Root distances (first 10):")
        sorted_distances = sorted(tree_dict[root].items(), key=lambda x: x[1])
        for cell, dist in sorted_distances[:10]:
            print(f"  {cell}: {dist}")
        
        # Check G05 specifically
        if 'HNSCC5_p7_HNSCC5_P7_G05' in tree_dict[root]:
            g05_dist = tree_dict[root]['HNSCC5_p7_HNSCC5_P7_G05']
            print(f"G05 distance from root: {g05_dist}")
    
    # Step 3: Compute MST
    print("Computing minimum spanning tree...")
    result = compute_rdmst(tree_dict, root)
    
    if result:
        mst, weight = result
        print(f"MST weight: {weight}")
        print(f"MST nodes: {len(mst)}")
        
        # Check root connection
        if root in mst and mst[root]:
            root_connections = list(mst[root].items())
            print(f"Root connects to: {root_connections}")
    else:
        print("Failed to compute MST")

if __name__ == "__main__":
    test_original()