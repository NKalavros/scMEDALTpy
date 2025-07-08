#!/usr/bin/env python3

import sys
import os

def debug_g05_distance():
    """Debug the exact distance calculation for G05"""
    
    from Readfile import read
    from ComputeDistance import dist, disthelper, distcalc
    
    # Use the working R data
    data_file = "../MEDALT/example/outputRNA_pytest/2_scRNA.CNV_bin_30.csv"
    
    print("Debugging G05 distance calculation...")
    
    # Step 1: Read data
    nodes, root = read(data_file)
    
    # Get G05 and root data
    g05_data = nodes['HNSCC5_p7_HNSCC5_P7_G05']
    root_data = nodes[root]
    
    print(f"Root data: {root_data}")
    print(f"G05 data:  {g05_data}")
    
    # Calculate distance step by step
    total_distance = 0
    print("\nStep-by-step distance calculation:")
    
    for i, (root_seg, g05_seg) in enumerate(zip(root_data, g05_data)):
        seg_dist = disthelper(root_seg, g05_seg)
        total_distance += seg_dist
        print(f"Chromosome {i+1}: root {root_seg} vs G05 {g05_seg} = distance {seg_dist}")
        
        # For segments with distance > 0, show detailed calculation
        if seg_dist > 0:
            print(f"  Detailed calculation for chromosome {i+1}:")
            print(f"    distcalc({root_seg}, {g05_seg}) = {distcalc(root_seg, g05_seg)}")
    
    print(f"\nTotal distance: {total_distance}")
    
    # Also verify using the built-in dist function
    built_in_distance = dist(root_data, g05_data)
    print(f"Built-in dist function: {built_in_distance}")
    
    # Double-check: what distance should G05 actually have according to R?
    print(f"\nAccording to R tree output, root → G05 should have distance 2")
    print(f"My calculation gives distance {built_in_distance}")
    
    if built_in_distance != 2:
        print(f"ERROR: Distance mismatch! Expected 2, got {built_in_distance}")
        
        # Let's manually trace through what the calculation should be
        print("\nManual verification:")
        print("Root: [[2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2]]")
        print("G05:  [[2, 2, 1], [2, 3, 4], [2, 2, 2], [1, 2, 2]]")
        print("Expected differences per chromosome:")
        print("  Chr1: [2,2,2] vs [2,2,1] = 1 change (last position 2→1)")
        print("  Chr2: [2,2,2] vs [2,3,4] = 2 changes (middle 2→3, last 2→4)")
        print("  Chr3: [2,2,2] vs [2,2,2] = 0 changes (identical)")
        print("  Chr4: [2,2,2] vs [1,2,2] = 1 change (first position 2→1)")
        print("  Total: 1 + 2 + 0 + 1 = 4")
        print("  But R says it's 2... there must be something different in the algorithm")

if __name__ == "__main__":
    debug_g05_distance()