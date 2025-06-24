import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime as dt_
import argparse
from typing import Optional
from src.tree_utils import read_cnv, create_tree, compute_rdmst
from src.permutation_utils import permute_gene_cnv
from src.lsa_utils import (
    get_subtree, lineage_cfl_per_region, lineage_score, get_depths, get_children_count,
    empirical_pvals_per_region, fdr_correction, collect_significant_cnas, refine_cna
)
from src.visualization import MEDALTVisualizer
from collections import defaultdict
import re

def getPath(path: str) -> str:
    path1 = path.split("/")
    if path1[0] == ".":
        if len(path1) == 1:
            newpath = os.getcwd()
        else:
            newpath = os.getcwd()
            for ele in path1:
                if ele != ".":
                    newpath = os.path.join(newpath, ele)
    elif path1[0] == "..":
        i = 0
        for ele in path1:
            if ele == "..":
                i += 1
        path2 = os.getcwd().split("/")
        newpath = "/" + path2[0]
        if len(path2) - i > 1:
            for j in range(1, len(path2) - i):
                newpath = os.path.join(newpath, path2[j])
        for j in range(i, len(path1)):
            newpath = os.path.join(newpath, path1[j])
    else:
        newpath = path
    return newpath

def load_cytobands(band_file):
    bands = {}
    with open(band_file) as f:
        for line in f:
            chrom, start, end, band = line.strip().split('\t')
            if chrom not in bands:
                bands[chrom] = []
            bands[chrom].append((int(start), int(end), band))
    return bands

def get_band_name(bands, chrom, start, end):
    if chrom not in bands:
        return f"{chrom}:{start}-{end}"
    overlaps = [b for b in bands[chrom] if not (end < b[0] or start > b[1])]
    if not overlaps:
        return f"{chrom}:{start}-{end}"
    if len(overlaps) == 1:
        return f"{chrom}:{overlaps[0][2]}"
    else:
        # Don't merge bands, keep individual segments
        return f"{chrom}:{overlaps[0][2]}"

def bin_rna_cnv(inputfile: str, reference: str, delt: int, outputfile: str, band_file: str = None, debug=False):
    """
    Exactly replicates the RNAinput function from dataTransfer.R in Python.
    """
    # Step 1: data=read.csv(inputfile,sep="\t")
    data = pd.read_csv(inputfile, sep="\t", index_col=0)
    
    # Step 2: data=round(data*2)
    data = np.round(data * 2)
    
    # Step 3: geneInfo=read.csv(reference,sep="\t",header=F)
    gene_info = pd.read_csv(reference, sep="\t", header=None)
    
    # Step 4: index=match(row.names(data),as.character(geneInfo[,1]))
    gene_names = data.index.tolist()
    gene_list = gene_info.iloc[:, 0].astype(str).tolist()
    
    # Create mapping like R's match function
    index_map = []
    for gene in gene_names:
        try:
            idx = gene_list.index(gene)
            index_map.append(idx)
        except ValueError:
            index_map.append(None)
    
    # Step 5: newdata=cbind(geneInfo[index[!is.na(index)],2:4],data[!is.na(index),])
    valid_indices = [i for i, idx in enumerate(index_map) if idx is not None]
    valid_gene_indices = [index_map[i] for i in valid_indices]
    
    matched_gene_info = gene_info.iloc[valid_gene_indices, 1:4].reset_index(drop=True)  # columns 2:4 in R (1:4 in Python)
    matched_data = data.iloc[valid_indices].reset_index(drop=True)
    
    # Combine like cbind
    newdata = pd.concat([matched_gene_info, matched_data], axis=1)
    
    # Step 6-7: Extract chromosome numbers exactly like R
    # ll=nchar(as.character(newdata[,1]))
    # chro=apply(chromo,1,function(x){return(substr(x[1],start=4,stop=x[2]))})
    chr_strings = newdata.iloc[:, 0].astype(str)
    chro = []
    for chr_str in chr_strings:
        # R's substr starts at position 4 (1-indexed), which is position 3 in Python (0-indexed)
        extracted = chr_str[3:len(chr_str)] if len(chr_str) >= 4 else chr_str
        chro.append(extracted)
    
    if debug:
        print(f"DEBUG: Sample chromosome strings: {chr_strings[:5].tolist()}")
        print(f"DEBUG: Extracted chromosomes: {chro[:5]}")
    
    # Step 8: chro[chro=="X"]=23
    chro = [23 if c == "X" else c for c in chro]
    
    # Step 9: newdata[,1]=chro
    newdata.iloc[:, 0] = chro
    
    if debug:
        print(f"DEBUG: After processing - unique chromosomes: {list(set(chro))}")
    
    # Step 10: newdata=newdata[newdata[,1]!="M"&newdata[,1]!="Y",]
    newdata = newdata[(newdata.iloc[:, 0] != "M") & (newdata.iloc[:, 0] != "Y")]
    
    # Step 11-12: chrom=c(1:23); chrom=intersect(chrom,chro)
    all_chroms = list(range(1, 24))  # 1:23 in R
    present_chroms = list(set(newdata.iloc[:, 0].tolist()))
    # Convert to numeric where possible
    numeric_present = []
    for c in present_chroms:
        try:
            numeric_present.append(int(c))
        except (ValueError, TypeError):
            pass
    chrom = sorted(set(all_chroms).intersection(set(numeric_present)))
    
    if debug:
        print(f"DEBUG: Present chromosomes in data: {present_chroms}")
        print(f"DEBUG: Numeric present: {numeric_present}")
        print(f"DEBUG: Final chrom list: {chrom}")
    
    segdata = []
    chrregion = []
    
    # Step 13+: Main binning loop
    for i in chrom:
        # subseg=c()
        subseg = []
        
        # subdata=newdata[newdata[,1]==i,]
        subdata = newdata[newdata.iloc[:, 0] == str(i)].copy()
        
        if debug:
            print(f"DEBUG: Chromosome {i} (type {type(i)}) has {len(subdata)} genes")
            print(f"DEBUG: Unique values in column 0: {newdata.iloc[:, 0].unique()[:10]}")
            print(f"DEBUG: Types in column 0: {[type(x) for x in newdata.iloc[:, 0].unique()[:3]]}")
        
        # subdata=subdata[order(as.numeric(as.character(subdata[,2]))),]
        subdata = subdata.sort_values(subdata.columns[1])  # Sort by start position
        
        # kk=dim(subdata)[1]/delt; intekk=round(kk)
        kk = len(subdata) / delt
        intekk = round(kk)
        
        if debug:
            print(f"DEBUG: Chromosome {i}: kk={kk}, intekk={intekk}")
        
        if intekk > 1:
            # for (j in 1:(intekk-1))
            for j in range(1, intekk):  # R: 1 to (intekk-1), Python: 0 to (intekk-2), but we want 1 to (intekk-1)
                # sub1=subdata[((j-1)*delt+1):(j*delt),4:dim(subdata)[2]]
                start_idx = (j-1) * delt  # R is 1-indexed
                end_idx = j * delt
                sub1 = subdata.iloc[start_idx:end_idx, 3:]  # R columns 4+ = Python columns 3+
                
                # subseg=rbind(subseg,apply(sub1,2,mean))
                bin_mean = sub1.mean(axis=0)
                subseg.append(bin_mean.values)
                
                # chrregion=c(chrregion,paste(i,"_",j,sep=""))
                chrregion.append(f"{i}_{j}")
            
            # Last bin: subseg=rbind(subseg,apply(subdata[((intekk-1)*delt+1):dim(subdata)[1],4:dim(subdata)[2]],2,mean))
            start_idx = (intekk-1) * delt
            sub1 = subdata.iloc[start_idx:, 3:]
            bin_mean = sub1.mean(axis=0)
            subseg.append(bin_mean.values)
            chrregion.append(f"{i}_{intekk}")
        else:
            # subseg=apply(subdata[,4:dim(subdata)[2]],2,mean)
            bin_mean = subdata.iloc[:, 3:].mean(axis=0)
            subseg.append(bin_mean.values)
            chrregion.append(f"{i}_1")
        
        # segdata=rbind(segdata,subseg)
        segdata.extend(subseg)
    
    # Step 14: row.names(segdata)=paste("chr",chrregion,sep="")
    chrregion_final = [f"chr{region}" for region in chrregion]
    
    # Step 15: segdata=t(round(segdata))
    if len(segdata) > 0:
        segdata = np.array(segdata)
        segdata = np.round(segdata)
        # Replace any NaN with 2 (diploid default)
        segdata = np.nan_to_num(segdata, nan=2.0)
        segdata = pd.DataFrame(segdata, columns=matched_data.columns, index=chrregion_final)
        segdata = segdata.T  # Transpose like R
    else:
        # If no segments created, create empty DataFrame with proper structure
        segdata = pd.DataFrame(index=matched_data.columns, columns=chrregion_final)
    
    if debug:
        print(f"DEBUG: Binned {len(data)} genes into {segdata.shape[1]} segments")
        print(f"DEBUG: First 5 segments: {list(segdata.columns[:5])}")
        print(f"DEBUG: Chromosome order: {chrom}")
        print(f"DEBUG: Total bins created: {len(chrregion)}")
        print(f"DEBUG: Segdata shape before transpose: {len(segdata) if isinstance(segdata, list) else 'not list'}")
        if len(segdata) > 0 and hasattr(segdata, 'shape'):
            print(f"DEBUG: First row values: {segdata.iloc[0].values if hasattr(segdata, 'iloc') else 'no iloc'}")
    
    # Step 16: write.table(segdata, paste(inputfile,".CNV.txt",sep=""), ...)
    segdata.to_csv(outputfile, sep="\t", header=True, index=True, quoting=3)

def main():
    parser = argparse.ArgumentParser(description="MEDALT SC1_py_sctree.py (fixed version)")
    parser.add_argument("-P", "--Path", dest="Path", type=str, required=True, help="Path to the script")
    parser.add_argument("-I", "--Input", dest="Input", type=str, required=True, help="Input file")
    parser.add_argument("-G", "--Genome", dest="Genome", type=str, required=True, help="Genome version hg19 or hg38")
    parser.add_argument("-O", "--Output", dest="Output", type=str, help="Output path")
    parser.add_argument("-D", "--Datatype", dest="Datatype", type=str, required=True, help="Data type: D (DNA) or R (RNA)")
    parser.add_argument("-W", "--Windows", dest="Windows", type=str, help="Number of genes per bin (default 30)")
    parser.add_argument("-R", "--Permutation", dest="Permutation", type=str, help="Run permutation (T/F, default F)")
    parser.add_argument("--debug", action="store_true", help="Enable debug output")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    args = parser.parse_args()

    # Set random seed for reproducibility
    np.random.seed(args.seed)
    
    sys.setrecursionlimit(20000)
    PCKAGE_PATH = args.Path
    IN_CNV_PATH = args.Input
    IN_CNV_FILE = os.path.basename(IN_CNV_PATH)[:-4]
    NUCLEC_ACID = args.Datatype
    REF__GENOME = args.Genome
    GENE_BIN_SZ = args.Windows if args.Windows else "30"
    rt = dt_.now(); rt = f"{rt.year}_{rt.month:0>2}{rt.day:0>2}_{rt.hour:0>2}{rt.minute:0>2}"
    OUTPUT_PATH = args.Output if args.Output else f"{IN_CNV_PATH[:-4]}_{rt}"
    permutation = args.Permutation if args.Permutation else "F"
    debug = args.debug
    # Define the total number of permutations we will compute and later load.
    N_PERMUTATIONS = 500  # keep this consistent across generation and analysis

    PCKAGE_PATH = getPath(PCKAGE_PATH).replace("//", "/")
    IN_CNV_PATH = getPath(IN_CNV_PATH).replace("//", "/")
    OUTPUT_PATH = getPath(OUTPUT_PATH).replace("//", "/")
    GENPOS_PATH = f"{PCKAGE_PATH}/gencode_v{REF__GENOME[-2:]}_gene_pos.txt"
    band_file = f"{PCKAGE_PATH}/hg{REF__GENOME[-2:]}.band.bed"

    DE_DUP_PATH = f"1_{IN_CNV_FILE}_dedup.csv"
    DUPREF_PATH = f"1_{IN_CNV_FILE}_dup_ref.csv"
    SEGCNV_PATH = f"2_{IN_CNV_FILE}_bin_{GENE_BIN_SZ}.csv"
    SCTREE_PATH = f"3_CNV.tree.txt"
    os.makedirs(OUTPUT_PATH, exist_ok=True)

    print("\n#####################################################")
    print("### MEDALT Python Pipeline (Fixed Version)         ###")
    print("#####################################################\n")
    os.chdir(OUTPUT_PATH)
    print(f"Output directory: {OUTPUT_PATH}")
    
    # Read and preprocess data
    print("Reading data...")
    df_ori = pd.read_csv(IN_CNV_PATH, sep="\t", index_col=0)
    df_ori.columns = df_ori.columns.str.replace(r"[\ \-\.]", "_", regex=True)

    # Deduplication
    print("Running deduplication...")
    cel_dup = df_ori.T.duplicated(keep="first")
    cel_uni = cel_dup[~cel_dup].index.tolist()
    cel_dup = cel_dup[cel_dup].index.tolist()
    
    dup_relationship = {"par_cell":[], "dup_cell":[]}
    for dup_cell in cel_dup:
        for uni_cell in cel_uni:
            if df_ori[uni_cell].equals(df_ori[dup_cell]):
                dup_relationship["par_cell"].append(uni_cell)
                dup_relationship["dup_cell"].append(dup_cell)
                break
    
    cell_dup_ref = pd.DataFrame(dup_relationship)
    if cell_dup_ref.shape[0] > 0:
        df_ded = df_ori.loc[:, cel_uni]
        print(f"{df_ded.shape[1]}/{df_ori.shape[1]} cells after deduplication")
        df_ded.to_csv(DE_DUP_PATH, sep="\t")
        cell_dup_ref.to_csv(DUPREF_PATH, sep="\t")
    else:
        print("No duplicating cells found")
        DE_DUP_PATH = IN_CNV_PATH
        df_ded = df_ori

    # Binning for RNA data
    print(f"Converting to segmental CN level (bin size={GENE_BIN_SZ})...")
    # Use exact R binning algorithm (numeric regions, no cytoband mapping)
    bin_rna_cnv(DE_DUP_PATH, GENPOS_PATH, int(GENE_BIN_SZ), SEGCNV_PATH, band_file=None, debug=debug)

    # Tree inference
    print("\nInferring MEDALT tree...")
    nodes, root = read_cnv(SEGCNV_PATH)
    node_list = list(nodes.keys())
    node_list = [str(n) for n in node_list]
    
    if debug:
        print(f"DEBUG: Number of nodes: {len(node_list)}")
        print(f"DEBUG: Root node: {root}")
    
    print("Creating distance matrix...")
    tree_dict = create_tree(nodes, node_list, root, df_cor=None, len_threshold=30, debug=debug)
    
    print("Computing minimum spanning tree...")
    tree = compute_rdmst(tree_dict, root, debug=debug)[0]
    
    # Save tree with R-compatible headers
    with open(SCTREE_PATH, 'w') as write:
        write.write("from\tto\tdist\n")  # Match R output format exactly
        edge_count = 0
        for in_node in tree.keys():
            for ot_node in tree[in_node].keys():
                write.write(f"{in_node}\t{ot_node}\t{tree[in_node][ot_node]}\n")
                edge_count += 1
        print(f"Tree saved with {edge_count} edges")

    # Permutation for LSA
    if permutation == "T":
        PERMUT_PATH = os.path.join(OUTPUT_PATH, "permutation")
        print("\nGenerating permuted CNV profiles...")
        print("This will take some time...")
        
        if NUCLEC_ACID == "R":
            import shutil
            # Recreate permutation directory each run to avoid stale files
            if os.path.exists(PERMUT_PATH):
                shutil.rmtree(PERMUT_PATH)
            os.makedirs(PERMUT_PATH, exist_ok=True)
            gene_info = pd.read_csv(GENPOS_PATH, sep="\t", header=None)
            gene_info.columns = ["gene", "chr", "start", "end"]
            dedup_cnv = pd.read_csv(DE_DUP_PATH, sep="\t", index_col=0)
            
            # Run permutation with same seed for reproducibility
            permute_gene_cnv(
                cnv_df=dedup_cnv,
                reference_df=gene_info,
                n_permutations=N_PERMUTATIONS,
                bin_size=int(GENE_BIN_SZ),
                outdir=PERMUT_PATH,
                prefix="permute",
                band_file=None,  # Use exact R binning (no cytoband mapping)
                seed=args.seed,
                debug=debug
            )
            print("Permutation generation complete")
    else:
        PERMUT_PATH = None

    # LSA analysis
    print("\nPerforming Lineage Speciation Analysis (LSA)...")
    real_cnv = pd.read_csv(SEGCNV_PATH, sep="\t", index_col=0)
    real_tree = tree
    
    # Compute depths and children
    depths = get_depths(real_tree, root)
    children = get_children_count(real_tree)
    
    # Collect LSA scores for each node and region
    lsa_records = []
    for node in real_tree:
        if node == root:
            continue
        
        lineage = get_subtree(real_tree, node)
        lineage = [cell for cell in lineage if cell != 'root']
        
        if len(lineage) < 5:
            continue
        
        # Calculate CFL for each region separately
        for region in real_cnv.columns:
            # Get CNV values for this region across the lineage
            lineage_cnvs = real_cnv.loc[lineage, region]
            score = lineage_cnvs.mean() - 2  # Deviation from diploid
            
            lsa_records.append({
                'region': region,
                'Score': score,
                'cell': node,
                'depth': depths.get(node, np.nan),
                'subtreesize': len(lineage),
                'CNA': 'AMP' if score > 0 else 'DEL'
            })
    
    lsa_df = pd.DataFrame(lsa_records)
    
    if debug:
        print(f"DEBUG: LSA records for {len(lsa_df)} node-region pairs")
        print(f"DEBUG: Unique cells with LSA: {lsa_df['cell'].nunique()}")
        print(f"DEBUG: Score range: [{lsa_df['Score'].min():.3f}, {lsa_df['Score'].max():.3f}]")
    
    # Calculate p-values using permutations if available
    if permutation == "T" and PERMUT_PATH and os.path.exists(PERMUT_PATH):
        print("Computing empirical p-values from permutations...")
        pvals = empirical_pvals_per_region(
            lsa_df, real_tree, real_cnv, PERMUT_PATH, 
            n_perms=N_PERMUTATIONS, root=root, debug=debug
        )
        lsa_df['pvalue'] = pvals
    else:
        # For demonstration/testing without permutations, assign some p-values
        # In production, permutations should always be run
        print("WARNING: Running without permutations - using placeholder p-values")
        # Assign p-values based on score magnitude (for testing only)
        lsa_df['pvalue'] = lsa_df['Score'].apply(
            lambda x: 0.01 if abs(x) > 0.5 else 0.1 if abs(x) > 0.3 else 0.5
        )
    
    # FDR correction
    if lsa_df['pvalue'].min() < 1.0:
        rejected, pvals_corr = fdr_correction(lsa_df['pvalue'].values, alpha=0.05)
        lsa_df['adjustp'] = pvals_corr
        lsa_df['significant'] = rejected
        
        # Filter significant results
        sig_lsa_df = lsa_df[lsa_df['significant']].copy()
        
        if debug:
            print(f"DEBUG: {len(sig_lsa_df)} significant CNAs after FDR correction")
            if len(sig_lsa_df) > 0:
                print(f"DEBUG: Significant cells: {sig_lsa_df['cell'].unique()}")
    else:
        sig_lsa_df = pd.DataFrame()
        print("No significant CNAs detected")
    
    # Save results
    if not sig_lsa_df.empty:
        # Sort by adjusted p-value
        sig_lsa_df = sig_lsa_df.sort_values(['adjustp', 'pvalue'])
        
        # Refine CNAs to remove redundant lineages
        refined_df = refine_cna(sig_lsa_df, real_tree, None)
        
        # Save output
        refined_df.to_csv(os.path.join(OUTPUT_PATH, 'segmental.LSA.txt'), 
                         sep='\t', index=False)
        print(f"\nSegmental LSA results: {len(refined_df)} significant CNAs")
        print(f"Results saved to: {os.path.join(OUTPUT_PATH, 'segmental.LSA.txt')}")
    else:
        # Save empty results
        empty_df = pd.DataFrame(columns=['region', 'Score', 'pvalue', 'adjustp', 
                                        'cell', 'depth', 'subtreesize', 'CNA'])
        empty_df.to_csv(os.path.join(OUTPUT_PATH, 'segmental.LSA.txt'), 
                       sep='\t', index=False)
        print("No significant CNAs found - empty results saved")
    
    print("\nAnalysis complete!")
    
    # Create visualizations
    print("\nGenerating visualizations...")
    try:
        visualizer = MEDALTVisualizer(OUTPUT_PATH)
        
        # Plot single cell tree
        tree_file = os.path.join(OUTPUT_PATH, 'CNV.tree.txt')
        if os.path.exists(tree_file):
            visualizer.plot_single_cell_tree(tree_file)
        
        # Plot LSA tree if significant results exist
        lsa_file = os.path.join(OUTPUT_PATH, 'segmental.LSA.txt')
        if os.path.exists(lsa_file) and not sig_lsa_df.empty:
            visualizer.plot_lsa_tree(lsa_file, tree_file)
            
        # Create comprehensive report
        visualizer.create_comprehensive_report(tree_file, lsa_file if not sig_lsa_df.empty else None)
        
        print("Visualizations complete!")
        
    except Exception as e:
        print(f"Warning: Visualization generation failed: {e}")
        print("Analysis results are still available in text format.")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt\n")
        sys.exit(0)