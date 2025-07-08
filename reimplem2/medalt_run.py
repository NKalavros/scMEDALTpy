from medalt_optimized import *
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # import necessary modules


pipeline = InferCNVPipeline('output_dir')
tree, results = pipeline.run_infercnv_analysis(
    observations_file='infercnv.observations.txt',
    top_k_genes=1000,
    subsample_size=2000,  # Automatic for large datasets
    n_permutations=500,
    min_lineage_size=20
)