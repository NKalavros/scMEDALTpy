#!/usr/bin/env python3

import sys
import os

def compare_data_structures():
    """Compare data structures between reimpl and expected R output"""
    
    from Readfile import read
    from ComputeDistance import dist
    
    # Use the working R data
    data_file = "../MEDALT/example/outputRNA_pytest/2_scRNA.CNV_bin_30.csv"
    
    print("Comparing data structures...")
    print(f"Reading data from: {data_file}")
    
    # Step 1: Read data
    nodes, root = read(data_file)
    print(f"Read {len(nodes)} nodes, root = {root}")
    
    # Check the data structure of key cells
    key_cells = ['HNSCC5_p7_HNSCC5_P7_G05', 'HNSCC5_p5_P5_E11']
    
    for cell in key_cells:
        if cell in nodes:
            print(f"\n{cell}:")
            print(f"  Structure: {nodes[cell]}")
            print(f"  Lengths: {[len(seg) for seg in nodes[cell]]}")
            total_segments = sum(len(seg) for seg in nodes[cell])
            print(f"  Total segments: {total_segments}")
    
    print(f"\nRoot structure: {nodes[root]}")
    print(f"Root lengths: {[len(seg) for seg in nodes[root]]}")
    
    # Calculate distances manually
    for cell in key_cells:
        if cell in nodes:
            distance = dist(nodes[root], nodes[cell])
            print(f"\nDistance from root to {cell}: {distance}")
            
            # Break down distance by chromosome
            for i, (root_seg, cell_seg) in enumerate(zip(nodes[root], nodes[cell])):
                from ComputeDistance import disthelper
                seg_dist = disthelper(root_seg, cell_seg)
                print(f"  Chromosome {i+1}: {seg_dist} (root: {root_seg}, cell: {cell_seg})")
    
    # Let's also check what the distance calculation produces for every cell
    print(f"\nAll distances from root:")
    all_distances = []
    for cell_name in nodes.keys():
        if cell_name != root:
            distance = dist(nodes[root], nodes[cell_name])
            all_distances.append((cell_name, distance))
    
    # Sort by distance
    all_distances.sort(key=lambda x: x[1])
    
    print("Top 10 closest cells to root:")
    for i, (cell, distance) in enumerate(all_distances[:10]):
        print(f"  {i+1}. {cell}: {distance}")
    
    # Find G05 position in the list
    for i, (cell, distance) in enumerate(all_distances):
        if 'G05' in cell:
            print(f"\nG05 position in distance ranking: {i+1} with distance {distance}")
            break

if __name__ == "__main__":
    compare_data_structures()