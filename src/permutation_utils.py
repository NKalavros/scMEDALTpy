import os
import numpy as np
import pandas as pd
from typing import Optional
import networkx as nx

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
        return f"{chrom}:{overlaps[0][2]}-{overlaps[-1][2]}"

def permute_gene_cnv(cnv_df: pd.DataFrame, reference_df: pd.DataFrame, n_permutations: int = 500, bin_size: int = 30, outdir: str = "./permutation", prefix: str = "permute", band_file: str = None): 
    """
    For each permutation:
      - For each chromosome, shuffle the gene rows within that chromosome independently across all cells.
      - Bin the permuted gene-level CNV.
      - Write both the permuted gene-level CNV and the binned CNV to the output directory.
    """
    os.makedirs(outdir, exist_ok=True)
    # Prepare chromosome info for each gene
    gene_chr = reference_df.set_index(0)[1].astype(str)
    gene_chr = gene_chr.reindex(cnv_df.index)
    gene_chr = gene_chr.fillna('NA')
    bands = load_cytobands(band_file) if band_file else None
    # For binning
    def bin_cnv(permuted_cnv, reference, bin_size, bands=None):
        gene_info = reference.set_index(0).loc[permuted_cnv.index]
        gene_info = gene_info.reset_index()
        gene_info.columns = ["gene", "chr", "start", "end"]
        gene_info = gene_info[~gene_info["chr"].isin(["M", "Y"])]
        chroms = sorted(set(gene_info["chr"]) - set(["M", "Y"]))
        segdata = []
        chrregion = []
        for chrom in chroms:
            chrom_gene_names = gene_info[gene_info["chr"] == chrom]["gene"]
            subdata = permuted_cnv.loc[chrom_gene_names]
            n_genes = subdata.shape[0]
            n_bins = int(np.ceil(n_genes / bin_size))
            starts = gene_info[gene_info["chr"] == chrom]["start"].values
            ends = gene_info[gene_info["chr"] == chrom]["end"].values
            for j in range(n_bins):
                start_idx = j * bin_size
                end_idx = min((j + 1) * bin_size, n_genes)
                bin_data = subdata.iloc[start_idx:end_idx, :]
                bin_mean = bin_data.mean(axis=0)
                segdata.append(bin_mean.values)
                bin_start = starts[start_idx]
                bin_end = ends[end_idx-1]
                if bands:
                    region = get_band_name(bands, chrom, bin_start, bin_end)
                else:
                    region = f"{chrom}_{j+1}"
                chrregion.append(region)
        segdata = np.vstack(segdata)
        segdata = np.round(segdata)
        segdata = pd.DataFrame(segdata, columns=permuted_cnv.columns, index=chrregion)
        segdata = segdata.T
        return segdata
    for j in range(1, n_permutations + 1):
        # Permute gene order within each chromosome
        permuted_idx = []
        for chrom in gene_chr.unique():
            if chrom == 'NA':
                continue
            chrom_genes = gene_chr[gene_chr == chrom].index.tolist()
            permuted_genes = np.random.permutation(chrom_genes)
            permuted_idx.extend(permuted_genes)
        permuted_cnv = cnv_df.loc[permuted_idx]
        # Write permuted gene-level CNV
        gene_cnv_path = os.path.join(outdir, f"{prefix}.{j}.gene.CNV.txt")
        permuted_cnv.to_csv(gene_cnv_path, sep="\t", header=True, index=True, quoting=3)
        # Bin and write binned CNV
        binned_cnv = bin_cnv(permuted_cnv, reference_df, bin_size, bands=bands)
        binned_cnv_path = os.path.join(outdir, f"{prefix}.{j}.CNV.txt")
        binned_cnv.to_csv(binned_cnv_path, sep="\t", header=True, index=True, quoting=3) 