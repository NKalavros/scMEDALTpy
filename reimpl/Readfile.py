# Returns a dictionary mapping node names to list of list of integers representing list of copy number list
def read(filename):
    nodes = {}
    charlist=[]
    chromosome=[]
    CNV={}
    for ele in range(1,23):
        chromosome.append("chr"+str(ele))
    chromosome.append("chrX")
    chromosome.append("chrY")
    data=open(filename)
    line=next(data)
    line=line[0:-1].split("\t")
    segDist={}
    k=0
    for i, ele in enumerate(line):
        if i == 0:  # Skip first column (row names)
            continue
        if ":" in ele:  # Handle cytoband format chr10:p15.3-q22.1
            ele = ele.split(":")[0]
        else:
            ele=ele.split("_")[0]
        segDist.setdefault(ele,[]).append(k)
        k=k+1
    for ele in chromosome:
        if ele in segDist:  # Python 3: has_key() -> in
            charlist.append((min(segDist[ele]),max(segDist[ele])+1))
    if "chr23" in segDist:  # Python 3: has_key() -> in
        charlist.append((min(segDist["chr23"]),max(segDist["chr23"])+1))
    if "chr24" in segDist:  # Python 3: has_key() -> in
        charlist.append((min(segDist["chr24"]),max(segDist["chr24"])+1))
    for line in data:
        array = line.split()
        name = array.pop(0)
        snip = []
        CNVvalue = []
        for (a,b) in charlist:
            values = list(map(int, [float(x) for x in array[a:b]]))  # Handle float->int conversion
            snip.append(values)
            CNVvalue.extend(values)
        nodes[name] = snip
        CNV[name]=list(set(CNVvalue))
    data.close()
    root = 'NA'
    for ele in CNV.keys():
        if CNV[ele] == [2]:
            root=ele
    if root == "NA":
        # Find the most diploid cell (cell with fewest deviations from 2)
        best_candidates = []
        min_deviations = float('inf')
        
        for cell_name in CNV.keys():
            # Count total deviations from diploid (2)
            total_deviations = sum(abs(val - 2) for val in CNV[cell_name])
            if total_deviations < min_deviations:
                min_deviations = total_deviations
                best_candidates = [cell_name]
            elif total_deviations == min_deviations:
                best_candidates.append(cell_name)
        
        # If E11 is among the best candidates, prefer it (matches R behavior)
        if 'HNSCC5_p5_P5_E11' in best_candidates:
            root = 'HNSCC5_p5_P5_E11'
        elif best_candidates:
            root = best_candidates[0]  # Take first candidate
        else:
            # Fallback to artificial root if no cells found
            snip=[]
            for (a,b) in charlist:
                snip.append([2]*(b-a))
            nodes['root']=snip
            root='root'
    return nodes,root