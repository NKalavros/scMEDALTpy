#!/usr/bin/env python3

import sys
import os

def debug_readfile():
    """Debug the readfile parsing"""
    
    data_file = "../MEDALT/example/outputRNA_pytest/2_scRNA.CNV_bin_30.csv"
    
    print("Debugging Readfile parsing...")
    print(f"File: {data_file}")
    
    nodes = {}
    charlist=[]
    chromosome=[]
    CNV={}
    
    # Build chromosome list
    for ele in range(1,23):
        chromosome.append("chr"+str(ele))
    chromosome.append("chrX")
    chromosome.append("chrY")
    
    # Read header
    data=open(data_file)
    line=next(data)
    line=line[0:-1].split("\t")
    print("Header:", line)
    
    segDist={}
    k=0
    for i, ele in enumerate(line):
        print(f"Processing header element {i}: '{ele}'")
        if i == 0:  # Skip first column (row names)
            print("  Skipping first column")
            continue
        if ":" in ele:  # Handle cytoband format chr10:p15.3-q22.1
            ele = ele.split(":")[0]
        else:
            ele=ele.split("_")[0]
        print(f"  Extracted chromosome: '{ele}' -> position {k}")
        segDist.setdefault(ele,[]).append(k)
        k=k+1
    
    print("segDist:", segDist)
    
    # Build charlist
    for ele in chromosome:
        if ele in segDist:
            range_tuple = (min(segDist[ele]),max(segDist[ele])+1)
            charlist.append(range_tuple)
            print(f"Chromosome {ele}: {range_tuple} = {range_tuple[1] - range_tuple[0]} segments")
    
    print("Final charlist:", charlist)
    
    # Read first few cells
    cell_count = 0
    for line in data:
        if cell_count >= 3:  # Only process first 3 cells
            break
            
        array = line.split()
        name = array.pop(0)
        print(f"\nProcessing cell: {name}")
        print(f"Array length: {len(array)}")
        print(f"Array: {array}")
        
        snip = []
        CNVvalue = []
        for i, (a,b) in enumerate(charlist):
            print(f"  Chromosome {i}: positions {a}:{b}")
            if b <= len(array):
                values = list(map(int, [float(x) for x in array[a:b]]))
                print(f"    Values: {values} (length {len(values)})")
                snip.append(values)
                CNVvalue.extend(values)
            else:
                print(f"    ERROR: Range {a}:{b} exceeds array length {len(array)}")
        
        nodes[name] = snip
        CNV[name]=list(set(CNVvalue))
        print(f"  Final snip: {snip}")
        print(f"  CNV values: {sorted(CNV[name])}")
        cell_count += 1
    
    data.close()
    
    # Find root
    root = 'NA'
    for ele in CNV.keys():
        if CNV[ele] == [2]:
            root=ele
            print(f"Found diploid root: {root}")
            break
    
    if root == "NA":
        print("No diploid cell found, creating artificial root")
        snip=[]
        for (a,b) in charlist:
            root_segment = [2]*(b-a)
            print(f"Root segment {a}:{b} = {root_segment} (length {len(root_segment)})")
            snip.append(root_segment)
        nodes['root']=snip
        root='root'
        print(f"Artificial root: {snip}")
    
    return nodes, root

if __name__ == "__main__":
    debug_readfile()