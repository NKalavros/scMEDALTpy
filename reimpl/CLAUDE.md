# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Python 3 reimplementation of the MEDALT (Minimum Edit Distance Algorithm for Lineage Tracing) algorithm. The project provides a complete pipeline for analyzing single-cell RNA sequencing data to construct phylogenetic trees and perform lineage speciation analysis.

## Key Commands

### Main Pipeline Commands
```bash
# Complete pipeline from scRNA data (recommended)
python3 run_medalt_from_scrna.py scRNA.CNV.txt

# Alternative step-by-step workflow
python3 preprocess_rna.py scRNA.CNV.txt 30
python3 run_medalt.py scRNA.CNV.txt.CNV.txt
python3 lsa_analysis.py results.txt scRNA.CNV.txt.CNV.txt
```

### Testing and Validation
```bash
# Basic functionality test
python3 test_reimpl.py

# Compare with R implementation
python3 compare_with_r.py

# Verify tree distances
python3 verify_tree_distances.py
```

### Visualization Commands
```bash
# Generate all visualizations
python3 visualize_medalt.py results.txt

# LSA-specific plots
python3 lsa_visualization.py segmental.LSA.txt results.txt
```

### Dependencies Installation
```bash
pip install pandas numpy matplotlib networkx seaborn scipy
```

## Architecture Overview

### Core Algorithm Components
- **Readfile.py**: Data reading and root selection (fixed to use most diploid actual cell)
- **ComputeDistance.py**: Distance calculations between cells (exact match to R implementation)
- **mdmst.py**: Minimum spanning tree algorithm using Edmonds' algorithm
- **Edmonds.py**: Tree construction from distance matrix

### Pipeline Components
- **run_medalt_from_scrna.py**: Main integrated pipeline script
- **preprocess_rna.py**: Converts gene expression data to chromosomal segments
- **run_medalt.py**: Core MEDALT analysis with visualization
- **lsa_analysis.py**: Lineage Speciation Analysis with permutation testing
- **visualize_medalt.py**: Comprehensive visualization module
- **lsa_visualization.py**: LSA-specific plotting

### Key Design Decisions
1. **Root Selection**: Uses most diploid actual cell (typically E11-type) rather than artificial root
2. **Distance Calculation**: Exact consistency with R implementation for Root→G05 distance = 2
3. **Data Format**: Supports both raw scRNA.CNV.txt and pre-processed CNV matrices
4. **Modular Design**: Each component can be run independently or as part of integrated pipeline

### Input/Output Flow
1. **Input**: scRNA.CNV.txt (gene × cell expression matrix)
2. **Preprocessing**: Gene expression → chromosomal segments (configurable bins)
3. **Distance Calculation**: Pairwise distances between cells
4. **Tree Construction**: Minimum spanning tree using Edmonds' algorithm
5. **LSA Analysis**: Permutation testing for lineage-associated copy number alterations
6. **Visualization**: Publication-quality PDFs and comprehensive reports

### Critical Implementation Notes
- All distance calculations must match R implementation exactly
- Root selection logic prefers E11-type cells when multiple equally diploid cells exist
- Permutation testing in LSA uses empirical p-values based on null distributions
- Visualization maintains R plot styles while adding Python enhancements

### Testing Strategy
- **test_reimpl.py**: Basic functionality validation
- **compare_with_r.py**: Ensures consistency with original R implementation
- **verify_tree_distances.py**: Validates distance calculations
- Expected key result: Root→G05 distance = 2 (matches R output)

### Performance Considerations
- Handles 50+ cells and 300+ genes efficiently
- Configurable LSA permutations (balance speed vs accuracy)
- Memory-efficient processing for large datasets
- Progress tracking for long-running analyses