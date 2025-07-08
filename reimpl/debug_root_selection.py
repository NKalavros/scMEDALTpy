#!/usr/bin/env python3

def debug_root_selection():
    """Debug which cell should be selected as root"""
    
    from Readfile import read
    
    # Use the working R data
    data_file = "../MEDALT/example/outputRNA_pytest/2_scRNA.CNV_bin_30.csv"
    
    print("Debugging root selection...")
    
    # Read the data to get CNV values
    data=open(data_file)
    line=next(data)
    line=line[0:-1].split("\t")
    
    # Build chromosome structure
    chromosome = []
    for ele in range(1,23):
        chromosome.append("chr"+str(ele))
    chromosome.append("chrX")
    chromosome.append("chrY")
    
    segDist = {}
    k = 0
    for i, ele in enumerate(line):
        if i == 0:  # Skip first column (row names)
            continue
        if ":" in ele:  # Handle cytoband format chr10:p15.3-q22.1
            ele = ele.split(":")[0]
        else:
            ele = ele.split("_")[0]
        segDist.setdefault(ele,[]).append(k)
        k = k + 1
    
    charlist = []
    for ele in chromosome:
        if ele in segDist:
            charlist.append((min(segDist[ele]),max(segDist[ele])+1))
    
    # Process all cells and calculate deviations
    CNV = {}
    cell_deviations = {}
    
    print("Analyzing all cells:")
    
    for line in data:
        array = line.split()
        name = array.pop(0)
        
        CNVvalue = []
        for (a,b) in charlist:
            if b <= len(array):
                values = list(map(int, [float(x) for x in array[a:b]]))
                CNVvalue.extend(values)
        
        CNV[name] = list(set(CNVvalue))
        
        # Calculate deviations from diploid (2)
        total_deviations = sum(abs(val - 2) for val in CNV[name])
        cell_deviations[name] = total_deviations
        
        # Check if perfectly diploid
        if CNV[name] == [2]:
            print(f"DIPLOID CELL FOUND: {name}")
    
    data.close()
    
    # Find the most diploid cells
    sorted_cells = sorted(cell_deviations.items(), key=lambda x: x[1])
    
    print(f"\nTop 10 most diploid cells:")
    for i, (cell, deviations) in enumerate(sorted_cells[:10]):
        cnv_values = sorted(CNV[cell])
        print(f"  {i+1}. {cell}: {deviations} deviations, CNV values: {cnv_values}")
    
    # Check what R actually used
    print(f"\nAccording to R tree, E11 should be root.")
    e11_name = 'HNSCC5_p5_P5_E11'
    if e11_name in cell_deviations:
        e11_deviations = cell_deviations[e11_name]
        e11_cnv = sorted(CNV[e11_name])
        e11_rank = sorted_cells.index((e11_name, e11_deviations)) + 1
        print(f"E11 deviations: {e11_deviations}, CNV values: {e11_cnv}, rank: {e11_rank}")

if __name__ == "__main__":
    debug_root_selection()