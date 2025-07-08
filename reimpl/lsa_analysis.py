#!/usr/bin/env python3
"""
Lineage Speciation Analysis (LSA) for MEDALT

Python implementation of the R LSA.tree.R and LSAfunction.R scripts.
Performs permutation-based statistical testing to identify copy number alterations
associated with cellular lineage expansion.

Usage:
    python3 lsa_analysis.py <tree_file> <cnv_file> [output_dir] [n_permutations]

Example:
    python3 lsa_analysis.py tree.txt cnv.txt ./lsa_results/ 100
"""

import pandas as pd
import numpy as np
import sys
import os
from typing import Dict, List, Tuple, Optional
import networkx as nx
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

class LSAAnalyzer:
    """
    Lineage Speciation Analysis implementation
    """
    
    def __init__(self, output_dir: str = "./lsa_results"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # Load chromosome band information
        self.chromosome_bands = self._load_chromosome_bands()
    
    def _load_chromosome_bands(self) -> Dict:
        """
        Load chromosome band information (simplified version)
        In practice, this would load from hg19.band.bed or hg38.band.bed
        """
        # Simplified chromosome bands for demonstration
        bands = {}
        
        # Define major chromosome arms and bands
        for chrom in range(1, 23):
            chr_name = f"chr{chrom}"
            bands[chr_name] = {
                'p': [f"{chr_name}:p11", f"{chr_name}:p12", f"{chr_name}:p13"],
                'q': [f"{chr_name}:q11", f"{chr_name}:q12", f"{chr_name}:q13", f"{chr_name}:q21", f"{chr_name}:q22"]
            }
        
        # Add X chromosome
        bands['chrX'] = {
            'p': ['chrX:p11', 'chrX:p21', 'chrX:p22'],
            'q': ['chrX:q11', 'chrX:q12', 'chrX:q21', 'chrX:q22']
        }
        
        return bands
    
    def _map_segments_to_bands(self, cnv_data: pd.DataFrame) -> Dict:
        """
        Map CNV segments to chromosome bands
        """
        segment_to_band = {}
        
        for segment in cnv_data.columns:
            if '_' in segment:
                # Parse segment like chr7_1, chr8_2
                parts = segment.split('_')
                if len(parts) >= 2:
                    chrom = parts[0]
                    segment_num = int(parts[1])
                    
                    if chrom in self.chromosome_bands:
                        # Map to bands based on segment number
                        all_bands = self.chromosome_bands[chrom]['p'] + self.chromosome_bands[chrom]['q']
                        if segment_num <= len(all_bands):
                            band = all_bands[min(segment_num - 1, len(all_bands) - 1)]
                            segment_to_band[segment] = band
                        else:
                            # Default mapping
                            segment_to_band[segment] = f"{chrom}:q"
                    else:
                        segment_to_band[segment] = segment
            else:
                segment_to_band[segment] = segment
        
        return segment_to_band
    
    def _calculate_tree_depth_and_subtree_size(self, tree_df: pd.DataFrame, root: str) -> Dict:
        """
        Calculate tree depth and subtree size for each cell
        """
        # Create networkx graph
        G = nx.from_pandas_edgelist(tree_df, 
                                   source=tree_df.columns[0], 
                                   target=tree_df.columns[1], 
                                   create_using=nx.DiGraph())
        
        cell_info = {}
        
        # Calculate depth (distance from root)
        try:
            depths = nx.single_source_shortest_path_length(G, root)
        except:
            # If root not found, use first node
            root = list(G.nodes())[0]
            depths = nx.single_source_shortest_path_length(G, root)
        
        # Calculate subtree size (number of descendants)
        def get_subtree_size(node):
            descendants = nx.descendants(G, node)
            return len(descendants) + 1  # Include the node itself
        
        for node in G.nodes():
            depth = depths.get(node, 0)
            subtree_size = get_subtree_size(node)
            
            cell_info[node] = {
                'depth': depth,
                'subtreesize': subtree_size
            }
        
        return cell_info
    
    def _calculate_lineage_scores(self, cnv_data: pd.DataFrame, cell_info: Dict) -> pd.DataFrame:
        """
        Calculate lineage speciation scores for each cell-segment combination
        """
        scores = []
        
        for cell in cnv_data.index:
            if cell not in cell_info:
                continue
                
            depth = cell_info[cell]['depth']
            subtree_size = cell_info[cell]['subtreesize']
            
            for segment in cnv_data.columns:
                cnv_value = cnv_data.loc[cell, segment]
                
                # Calculate deviation from diploid (2)
                deviation = cnv_value - 2
                
                if deviation != 0:  # Only consider non-diploid segments
                    # Score based on deviation magnitude and subtree size
                    score = deviation * np.log(subtree_size + 1)
                    
                    scores.append({
                        'cell': cell,
                        'segment': segment,
                        'cnv_value': cnv_value,
                        'deviation': deviation,
                        'score': score,
                        'depth': depth,
                        'subtreesize': subtree_size
                    })
        
        return pd.DataFrame(scores)
    
    def _permute_cnv_data(self, cnv_data: pd.DataFrame) -> pd.DataFrame:
        """
        Permute CNV data by chromosome (shuffle cells within each chromosome)
        """
        permuted_data = cnv_data.copy()
        
        # Group segments by chromosome
        chr_groups = {}
        for segment in cnv_data.columns:
            if '_' in segment:
                chrom = segment.split('_')[0]
                if chrom not in chr_groups:
                    chr_groups[chrom] = []
                chr_groups[chrom].append(segment)
        
        # Permute within each chromosome
        for chrom, segments in chr_groups.items():
            if len(segments) > 1:
                # Shuffle cell order for this chromosome's segments
                shuffled_indices = np.random.permutation(len(cnv_data))
                for segment in segments:
                    permuted_data[segment] = cnv_data[segment].iloc[shuffled_indices].values
        
        return permuted_data
    
    def _run_permutation_test(self, cnv_data: pd.DataFrame, tree_df: pd.DataFrame, 
                             root: str, n_permutations: int = 100) -> Dict:
        """
        Run permutation test to calculate empirical p-values
        """
        print(f"Running {n_permutations} permutations for p-value calculation...")
        
        # Calculate cell info once
        cell_info = self._calculate_tree_depth_and_subtree_size(tree_df, root)
        
        # Calculate real scores
        real_scores = self._calculate_lineage_scores(cnv_data, cell_info)
        
        # Run permutations
        permutation_scores = []
        
        for i in range(n_permutations):
            if i % 20 == 0:
                print(f"  Permutation {i+1}/{n_permutations}")
            
            # Permute CNV data
            permuted_cnv = self._permute_cnv_data(cnv_data)
            
            # Calculate scores for permuted data
            perm_scores = self._calculate_lineage_scores(permuted_cnv, cell_info)
            permutation_scores.append(perm_scores)
        
        # Calculate p-values
        print("Calculating empirical p-values...")
        pvalue_results = []
        
        for _, real_row in real_scores.iterrows():
            real_score = abs(real_row['score'])
            cell = real_row['cell']
            segment = real_row['segment']
            
            # Count how many permutation scores are >= real score
            extreme_count = 0
            total_count = 0
            
            for perm_df in permutation_scores:
                # Find matching cell-segment combinations in permuted data
                matching_rows = perm_df[
                    (perm_df['cell'] == cell) & 
                    (perm_df['segment'] == segment)
                ]
                
                for _, perm_row in matching_rows.iterrows():
                    perm_score = abs(perm_row['score'])
                    if perm_score >= real_score:
                        extreme_count += 1
                    total_count += 1
            
            # Calculate empirical p-value
            if total_count > 0:
                pvalue = (extreme_count + 1) / (total_count + 1)  # Add 1 for continuity
            else:
                pvalue = 1.0
            
            pvalue_results.append({
                'cell': cell,
                'segment': segment,
                'score': real_row['score'],
                'deviation': real_row['deviation'],
                'depth': real_row['depth'],
                'subtreesize': real_row['subtreesize'],
                'pvalue': pvalue
            })
        
        return pd.DataFrame(pvalue_results)
    
    def _filter_significant_results(self, pvalue_results: pd.DataFrame, 
                                   cutoff: float = 0.1) -> pd.DataFrame:
        """
        Filter significant results and map to chromosome bands
        """
        # Filter by p-value
        significant = pvalue_results[pvalue_results['pvalue'] <= cutoff].copy()
        
        if len(significant) == 0:
            return pd.DataFrame()
        
        # Map segments to chromosome bands
        segment_to_band = self._map_segments_to_bands(
            pd.DataFrame(columns=significant['segment'].unique())
        )
        
        significant['region'] = significant['segment'].map(segment_to_band)
        
        # Determine CNA type
        significant['CNA'] = 'AMP'
        significant.loc[significant['score'] < 0, 'CNA'] = 'DEL'
        
        # Calculate adjusted p-values (Benjamini-Hochberg)
        from scipy.stats import false_discovery_control
        try:
            significant['adjustp'] = false_discovery_control(significant['pvalue'])
        except:
            # Fallback if scipy not available
            significant['adjustp'] = significant['pvalue'] * len(significant)
        
        # Rename columns to match R output format
        result = significant[['region', 'score', 'pvalue', 'adjustp', 'cell', 'depth', 'subtreesize', 'CNA']].copy()
        result.columns = ['region', 'Score', 'pvalue', 'adjustp', 'cell', 'depth', 'subtreesize', 'CNA']
        
        return result
    
    def run_lsa_analysis(self, tree_file: str, cnv_file: str, 
                        n_permutations: int = 100) -> Optional[pd.DataFrame]:
        """
        Run complete LSA analysis
        """
        print("="*60)
        print("LINEAGE SPECIATION ANALYSIS (LSA)")
        print("="*60)
        print(f"Tree file: {tree_file}")
        print(f"CNV file: {cnv_file}")
        print(f"Permutations: {n_permutations}")
        print(f"Output directory: {self.output_dir}")
        print()
        
        # Load data
        try:
            tree_df = pd.read_csv(tree_file, sep='\t', comment='#')
            cnv_data = pd.read_csv(cnv_file, sep='\t', index_col=0)
            
            print(f"✓ Loaded tree: {len(tree_df)} edges")
            print(f"✓ Loaded CNV data: {cnv_data.shape[0]} cells × {cnv_data.shape[1]} segments")
            
        except Exception as e:
            print(f"Error loading data: {e}")
            return None
        
        # Find root
        from_nodes = set(tree_df.iloc[:, 0])
        to_nodes = set(tree_df.iloc[:, 1])
        root_nodes = from_nodes - to_nodes
        root = list(root_nodes)[0] if root_nodes else list(from_nodes)[0]
        
        print(f"✓ Root cell: {root}")
        print()
        
        # Run permutation test
        try:
            pvalue_results = self._run_permutation_test(cnv_data, tree_df, root, n_permutations)
            print(f"✓ Permutation test completed")
            
            # Show distribution of p-values first
            print(f"P-value distribution:")
            print(f"  • Min p-value: {pvalue_results['pvalue'].min():.4f}")
            print(f"  • Max p-value: {pvalue_results['pvalue'].max():.4f}")
            print(f"  • Median p-value: {pvalue_results['pvalue'].median():.4f}")
            print(f"  • P-values < 0.1: {len(pvalue_results[pvalue_results['pvalue'] < 0.1])}")
            print(f"  • P-values < 0.05: {len(pvalue_results[pvalue_results['pvalue'] < 0.05])}")
            
            # Filter significant results
            significant_results = self._filter_significant_results(pvalue_results)
            
            if len(significant_results) > 0:
                print(f"✓ Found {len(significant_results)} significant LSA events")
                
                # Save results
                output_file = os.path.join(self.output_dir, "segmental.LSA.txt")
                significant_results.to_csv(output_file, sep='\t', index=False)
                print(f"✓ Results saved to {output_file}")
                
                # Print summary
                print()
                print("LSA Results Summary:")
                print(f"  • Total significant events: {len(significant_results)}")
                print(f"  • Cells with LSA: {len(significant_results['cell'].unique())}")
                print(f"  • Amplifications: {len(significant_results[significant_results['CNA'] == 'AMP'])}")
                print(f"  • Deletions: {len(significant_results[significant_results['CNA'] == 'DEL'])}")
                
                # Show top events
                print("\nTop 5 most significant events:")
                top_events = significant_results.sort_values('pvalue').head(5)
                for _, row in top_events.iterrows():
                    print(f"  • {row['cell']}: {row['region']} ({row['CNA']}) p={row['pvalue']:.4f}")
                
                return significant_results
                
            else:
                print("No significant LSA events found")
                
                # Save all results for debugging
                debug_file = os.path.join(self.output_dir, "all_lsa_scores.txt")
                pvalue_results.to_csv(debug_file, sep='\t', index=False)
                print(f"All scores saved to {debug_file} for debugging")
                
                return pd.DataFrame()
                
        except Exception as e:
            print(f"Error during LSA analysis: {e}")
            import traceback
            traceback.print_exc()
            return None

def main():
    """Main entry point"""
    
    if len(sys.argv) < 3:
        print("Usage: python3 lsa_analysis.py <tree_file> <cnv_file> [output_dir] [n_permutations]")
        print()
        print("Example:")
        print("  python3 lsa_analysis.py tree.txt cnv.txt")
        print("  python3 lsa_analysis.py tree.txt cnv.txt ./lsa_results/ 100")
        sys.exit(1)
    
    tree_file = sys.argv[1]
    cnv_file = sys.argv[2]
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "./lsa_results"
    n_permutations = int(sys.argv[4]) if len(sys.argv) > 4 else 100
    
    # Check if files exist
    for file_path in [tree_file, cnv_file]:
        if not os.path.exists(file_path):
            print(f"Error: File '{file_path}' not found!")
            sys.exit(1)
    
    # Run LSA analysis
    analyzer = LSAAnalyzer(output_dir)
    results = analyzer.run_lsa_analysis(tree_file, cnv_file, n_permutations)
    
    if results is not None and len(results) > 0:
        print(f"\n🎉 LSA analysis completed successfully!")
        print(f"Results saved to: {output_dir}/segmental.LSA.txt")
    else:
        print(f"\n⚠️  LSA analysis completed but no significant events found")

if __name__ == "__main__":
    main()