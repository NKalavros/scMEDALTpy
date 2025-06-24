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
    Replicates the RNAinput function from dataTransfer.R in Python.
    """
    data = pd.read_csv(inputfile, sep="\t", index_col=0)
    data = np.round(data * 2)
    gene_info = pd.read_csv(reference, sep="\t", header=None)
    gene_info.columns = ["gene", "chr", "start", "end"]
    
    # Match genes
    idx = data.index.to_series().map(lambda g: gene_info["gene"].tolist().index(g) if g in gene_info["gene"].tolist() else np.nan)
    valid = ~idx.isna()
    matched_gene_info = gene_info.iloc[idx[valid].astype(int)]
    matched_data = data.loc[valid]
    
    newdata = pd.concat([matched_gene_info.reset_index(drop=True), matched_data.reset_index(drop=True)], axis=1)
    newdata = newdata[~newdata["chr"].isin(["M", "Y"])]
    
    chroms = sorted(set(newdata["chr"]) - set(["M", "Y"]), key=lambda x: (len(x), x))
    segdata = []
    chrregion = []
    bands = load_cytobands(band_file) if band_file else None
    
    for chrom in chroms:
        subdata = newdata[newdata["chr"] == chrom].sort_values("start")
        n_genes = subdata.shape[0]
        n_bins = max(1, int(np.round(n_genes / delt)))
        starts = subdata["start"].values
        ends = subdata["end"].values
        
        for j in range(n_bins):
            start_idx = j * delt
            end_idx = min((j + 1) * delt, n_genes)
            bin_data = subdata.iloc[start_idx:end_idx, 4:]
            bin_mean = bin_data.mean(axis=0)
            segdata.append(bin_mean.values)
            bin_start = starts[start_idx]
            bin_end = ends[end_idx-1]
            
            if bands:
                region = get_band_name(bands, chrom, bin_start, bin_end)
            else:
                region = f"chr{chrom}_{j+1}"
            chrregion.append(region)
    
    segdata = np.vstack(segdata)
    segdata = np.round(segdata)
    segdata = pd.DataFrame(segdata, columns=matched_data.columns, index=chrregion)
    segdata = segdata.T
    
    if debug:
        print(f"DEBUG: Binned {data.shape[0]} genes into {segdata.shape[1]} segments")
        print(f"DEBUG: First 5 segments: {list(segdata.columns[:5])}")
    
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
<<<<<<< HEAD
    debug = args.debug
=======
    # Define the total number of permutations we will compute and later load.
    N_PERMUTATIONS = 500  # keep this consistent across generation and analysis
>>>>>>> main

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
    # Use numeric bins like the original R implementation (no cytoband mapping)
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
            dedup_cnv = pd.read_csv(DE_DUP_PATH, sep="\t", index_col=0)
            
            # Run permutation with same seed for reproducibility
            permute_gene_cnv(
                cnv_df=dedup_cnv,
                reference_df=gene_info,
<<<<<<< HEAD
                n_permutations=500,  # Match R version
=======
                n_permutations=N_PERMUTATIONS,
>>>>>>> main
                bin_size=int(GENE_BIN_SZ),
                outdir=PERMUT_PATH,
                prefix="permute",
                band_file=band_file,
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
<<<<<<< HEAD
    
    if debug:
        print(f"DEBUG: LSA records for {len(lsa_df)} node-region pairs")
        print(f"DEBUG: Unique cells with LSA: {lsa_df['cell'].nunique()}")
        print(f"DEBUG: Score range: [{lsa_df['Score'].min():.3f}, {lsa_df['Score'].max():.3f}]")
    
    # Calculate p-values using permutations if available
    if permutation == "T" and PERMUT_PATH and os.path.exists(PERMUT_PATH):
        print("Computing empirical p-values from permutations...")
        pvals = empirical_pvals_per_region(
            lsa_df, real_tree, real_cnv, PERMUT_PATH, 
            n_perms=500, root=root, debug=debug
=======
    # Load permutation results (binned CNVs)
    perm_cfls = defaultdict(list)  # (node, feature) -> list of permuted CFLs
    for j in range(1, N_PERMUTATIONS + 1):
        perm_cnv_path = os.path.join(PERMUT_PATH, f"permute.{j}.CNV.txt")
        if not os.path.exists(perm_cnv_path):
            continue
        perm_cnv = pd.read_csv(perm_cnv_path, sep="\t", index_col=0)
        for node in real_tree:
            if node == root:
                continue
            lineage = get_subtree(real_tree, node)
            lineage = [cell for cell in lineage if cell != 'root']
            cfl = lineage_cfl(perm_cnv, lineage)
            for feature, score in cfl.items():
                perm_cfls[(node, feature)].append(score)
    # Compute empirical p-values for each (node, feature)
    pvals = []
    for idx, row in lsa_df.iterrows():
        key = (row['cell'], row['region'])
        observed = row['Score']
        perm_scores = np.array(perm_cfls.get(key, [0]*N_PERMUTATIONS))
        if len(perm_scores) == 0:
            pval = 1.0
        else:
            if row['CNA'] == 'AMP':
                pval = (np.sum(perm_scores >= observed) + 1) / (len(perm_scores) + 1)
            else:
                pval = (np.sum(perm_scores <= observed) + 1) / (len(perm_scores) + 1)
        pvals.append(pval)
    lsa_df['pvalue'] = pvals
    # FDR correction
    rejected, pvals_corr = fdr_correction(lsa_df['pvalue'].values, alpha=0.05)
    lsa_df['adjustp'] = pvals_corr
    lsa_df['significant'] = rejected
    # Output segmental LSA results (significant only)
    sig_lsa_df = lsa_df[lsa_df['significant']].copy()

    # If no significant CNAs, gracefully skip downstream merging/pruning
    if sig_lsa_df.empty:
        print("No significant CNAs detected after FDR correction — skipping merging and pruning steps.")
        merged_out = pd.DataFrame(columns=['region', 'Score', 'pvalue', 'adjustp',
                                           'cell', 'depth', 'subtreesize', 'CNA'])
    else:
        # Build cytoband order map for adjacency-aware merging
        cytoband_order = {}
        with open(band_file) as f:
            for line in f:
                chrom, start, end, band = line.strip().split('\t')
                if chrom not in cytoband_order:
                    cytoband_order[chrom] = []
                cytoband_order[chrom].append(band)

        # Helper to parse band region names (e.g., chr7:q22.3-q31.2)
        def parse_band(region):
            m = re.match(r'(chr\w+):(\w+\.?\d*)(?:-(\w+\.?\d*))?', region)
            if not m:
                return None, None, None
            chrom = m.group(1)
            start_band = m.group(2)
            end_band = m.group(3) if m.group(3) else m.group(2)
            return chrom, start_band, end_band

        # Sort and merge adjacent regions with same cell, CNA, and contiguous bands (using cytoband order)
        def merge_regions(df):
            merged = []
            for (cell, cna), group in df.groupby(['cell', 'CNA']):
                parsed = group['region'].apply(parse_band)
                group = group.assign(
                    chrom=parsed.apply(lambda x: x[0]),
                    start_band=parsed.apply(lambda x: x[1]),
                    end_band=parsed.apply(lambda x: x[2])
                )
                group = group.sort_values(['chrom', 'start_band'], key=lambda x: x.map(lambda b: cytoband_order.get(group['chrom'].iloc[0], []).index(b) if b in cytoband_order.get(group['chrom'].iloc[0], []) else -1))
                prev = None
                for idx, row in group.iterrows():
                    if prev is None:
                        prev = row.copy()
                    else:
                        # Check adjacency in cytoband order
                        chrom = row['chrom']
                        band_list = cytoband_order.get(chrom, [])
                        try:
                            prev_end_idx = band_list.index(prev['end_band'])
                            curr_start_idx = band_list.index(row['start_band'])
                        except ValueError:
                            merged.append(prev)
                            prev = row.copy()
                            continue
                        if (chrom == prev['chrom'] and curr_start_idx == prev_end_idx + 1):
                            prev['end_band'] = row['end_band']
                            prev['Score'] = (prev['Score'] + row['Score']) / 2
                            prev['pvalue'] = min(prev['pvalue'], row['pvalue'])
                            prev['adjustp'] = min(prev['adjustp'], row['adjustp'])
                            prev['subtreesize'] = max(prev['subtreesize'], row['subtreesize'])
                        else:
                            merged.append(prev)
                            prev = row.copy()
                if prev is not None:
                    merged.append(prev)
            out = []
            for row in merged:
                region = f"{row['chrom']}:{row['start_band']}" if row['start_band'] == row['end_band'] else f"{row['chrom']}:{row['start_band']}-{row['end_band']}"
                out.append({
                    'region': region,
                    'Score': row['Score'],
                    'pvalue': row['pvalue'],
                    'adjustp': row['adjustp'],
                    'cell': row['cell'],
                    'depth': row['depth'],
                    'subtreesize': row['subtreesize'],
                    'CNA': row['CNA']
                })
            return pd.DataFrame(out)

        merged_out = merge_regions(sig_lsa_df)
        # Build realcell_df (cell, depth, subtreesize)
        realcell_df = pd.DataFrame({
            'cell': list(depths.keys()),
            'depth': [depths[c] for c in depths]
        })
        realcell_df['subtreesize'] = realcell_df['cell'].map(
            lambda n: len(get_subtree(real_tree, n))
>>>>>>> main
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

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt\n")
        sys.exit(0)