# MEDALT: Minimal Event Distance Aneuploidy Lineage Tree

A Python implementation of the MEDALT algorithm for single-cell copy number variation (CNV) analysis and phylogenetic tree reconstruction.

## Overview

MEDALT reconstructs phylogenetic trees from single-cell genomic data using minimal event distance calculations and performs Lineage Speciation Analysis (LSA) to identify copy number alterations associated with lineage expansion.

## Installation

### Prerequisites

- Python 3.8 or higher
- Required packages: numpy, pandas, scipy, networkx, matplotlib, seaborn

### Setup

```bash
# Create conda environment
conda env create -f environment.yml
conda activate medalt-py312

# Or install dependencies manually
pip install numpy pandas scipy networkx matplotlib seaborn
```

## Input Format

The input file should be a tab-separated CNV matrix where:
- **Rows**: Genes/genomic regions
- **Columns**: Single cells
- **Values**: Copy number estimates (e.g., 0.5 = CN1, 1.0 = CN2, 1.5 = CN3)

Example format:
```
	CELL1	CELL2	CELL3	...
GENE1	1.0	1.5	0.5	...
GENE2	1.2	1.0	2.0	...
...
```

## Usage

### Quick Start

For the provided example data:

```bash
python run_medalt.py
```

This will:
1. Process `scRNA.CNV.txt` in the current directory
2. Generate `CNV.tree.txt` (phylogenetic tree)
3. Generate `gene.LSA.txt` (LSA results)

### Advanced Usage

#### 1. Core Implementation

```python
from medalt_implementation import MEDALT, LineageSpeciationAnalysis, load_scRNA_data

# Load data
cnv_matrix, cell_names, gene_names = load_scRNA_data('scRNA.CNV.txt', window_size=1)

# Build tree
medalt = MEDALT(cnv_matrix)
tree = medalt.build_tree()

# Run LSA
lsa = LineageSpeciationAnalysis(tree, cnv_matrix, n_permutations=500)
results = lsa.run_analysis(min_lineage_size=5, gene_names=gene_names)
```

#### 2. Complete Pipeline

```python
from medalt_utils import MEDALTPipeline

# Run full pipeline
pipeline = MEDALTPipeline(output_dir='results')
tree, lsa_results = pipeline.run_analysis(
    input_file='scRNA.CNV.txt',
    data_type='R',
    genome_version='hg19',
    window_size=30,
    n_permutations=500,
    min_lineage_size=5
)
```

#### 3. Command Line Interface

```bash
python medalt_utils.py -I scRNA.CNV.txt -D R -G hg19 -O output_dir
```

**Arguments:**
- `-I`: Input CNV file
- `-D`: Data type ('R' for scRNA-seq, 'D' for scDNA-seq)
- `-G`: Genome version ('hg19' or 'hg38')
- `-O`: Output directory
- `-W`: Window size for scRNA-seq (default: 30)
- `-R`: Number of permutations (default: 500)
- `-M`: Minimum lineage size (default: 5)
- `-J`: Number of parallel jobs (default: all CPUs)

## Output Files

### 1. CNV.tree.txt
Phylogenetic tree in edge list format:
```
from	to	dist
root	CELL1	2
CELL1	CELL2	1
...
```

### 2. gene.LSA.txt
Lineage Speciation Analysis results:
```
region	Score	pvalue	adjustp	cell	depth	subtreesize	CNA
GENE1	1.42	0.0099	0.0186	CELL1	1	50	AMP
GENE2	-0.96	0.0099	0.0495	CELL1	1	50	DEL
...
```

**Columns:**
- `region`: Gene or genomic region name
- `Score`: Normalized CFL score (deviation from diploid)
- `pvalue`: Empirical p-value from permutation test
- `adjustp`: FDR-adjusted p-value
- `cell`: Root cell of the lineage
- `depth`: Tree depth of the lineage
- `subtreesize`: Number of cells in the lineage
- `CNA`: Copy number alteration type (AMP/DEL)

### 3. segmental.LSA.txt
Chromosome-level aggregated results (if using windowed data).

### 4. Additional Files
- `analysis_summary.json`: Summary statistics
- `singlecell.tree.pdf`: Tree visualization
- `cnv_heatmap.png`: CNV heatmap ordered by tree structure
- `LSA.tree.pdf`: LSA results visualization

## Algorithm Details

### 1. Minimal Event Distance (MED)
Calculates the minimum number of copy number events needed to transform one profile into another, considering:
- Homozygous deletion constraint (CN=0 cannot increase)
- Segmental nature of copy number alterations
- Directionality of changes

### 2. Tree Construction
Uses Edmonds' minimum spanning arborescence algorithm to find the optimal phylogenetic tree that minimizes total evolutionary distance.

### 3. Lineage Speciation Analysis (LSA)
Identifies copy number alterations associated with lineage expansion through:
- Tree dissection into non-overlapping lineages
- Cumulative Fold Level (CFL) calculation
- Permutation-based background distribution
- Multiple testing correction

## Performance Options

### For Large Datasets (>5k cells)
Use the optimized version:
```python
from medalt_optimized import OptimizedMEDALT, OptimizedLSA

# Memory-efficient processing
medalt = OptimizedMEDALT(cnv_matrix, chunk_size=1000)
tree = medalt.build_tree_chunked()
```

### Parameter Tuning
- **window_size**: Larger values reduce noise but lose resolution
- **n_permutations**: More permutations improve statistical power
- **min_lineage_size**: Larger values focus on major lineages
- **n_jobs**: Parallel processing for permutation tests

## Testing

Run the test suite:
```bash
python -m unittest medalt_unittest.py
```

Test specific components:
```bash
python medalt_unittest.py TestMED.test_single_event
python medalt_unittest.py TestMEDALT.test_tree_construction
python medalt_unittest.py TestLSA.test_permutation_test
```

## Troubleshooting

### Common Issues

1. **Memory errors with large datasets**
   - Use `medalt_optimized.py` instead
   - Reduce `window_size` or use subsampling

2. **No significant results**
   - Increase `n_permutations`
   - Decrease `min_lineage_size`
   - Check data quality and preprocessing

3. **Slow performance**
   - Increase `n_jobs` for parallel processing
   - Use larger `window_size` for scRNA-seq data
   - Consider subsampling very large datasets

### Dependencies
If you encounter import errors:
```bash
pip install --upgrade numpy pandas scipy networkx matplotlib seaborn
```

## Citation

If you use this implementation, please cite the original MEDALT paper:
> Chen, K. et al. MEDALT: single-cell copy number lineage tracing enabling gene discovery. Genome Biology (2021).

## License

This implementation is provided for research purposes. Please refer to the original MEDALT repository for licensing information.