# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MEDALT (Minimal Event Distance Aneuploidy Lineage Tree) is a bioinformatics tool for single-cell copy number variation (CNV) analysis. It reconstructs phylogenetic trees from single-cell genomic data using minimal event distance calculations and performs Lineage Speciation Analysis (LSA) to identify copy number alterations associated with lineage expansion.

## Development Commands

### Environment Setup
```bash
# Create and activate conda environment
conda env create -f ../environment.yml
conda activate medalt-py312
```

### Running the Pipeline
```bash
# Basic implementation
python medalt_implementation.py

# Optimized for large datasets (10k+ cells)
python medalt_optimized.py observations.txt --genes 1000 --subsample 2000

# Full pipeline with utilities
python medalt_utils.py -I scRNA.CNV.txt -D R -G hg19 -O output_dir

# InferCNV integration
python medalt_infercnv.py
```

### Testing
```bash
# Run comprehensive unit tests
python -m unittest medalt_unittest.py

# Run specific test methods
python medalt_unittest.py TestMED.test_single_event
python medalt_unittest.py TestMEDALT.test_tree_construction
python medalt_unittest.py TestLSA.test_permutation_test
```

## Code Architecture

### Core Components
- `MED` class - Minimal event distance calculations with Numba optimization
- `MEDALT` class - Tree construction using Edmonds' minimum spanning arborescence algorithm
- `LineageSpeciationAnalysis` class - Statistical analysis for significant CNAs
- `OptimizedMEDALT` and `OptimizedLSA` classes - High-performance versions for large datasets

### Key Files
- `medalt_implementation.py` - Main algorithm implementation
- `medalt_optimized.py` - Memory-optimized version for large datasets
- `medalt_utils.py` - Pipeline utilities and visualization tools
- `medalt_unittest.py` - Comprehensive test suite
- `medalt_infercnv.py` - InferCNV integration

### Architecture Patterns
- Object-oriented design with clear separation of concerns
- Modular components usable independently
- Performance optimization through Numba JIT compilation
- Sparse matrix operations for memory efficiency
- Parallel processing for permutation tests

## Data Formats

### Input
- scRNA-seq: Gene expression matrix (inferCNV format)
- scDNA-seq: Integer copy number profiles
- Processed CNV: Tab-separated chromosomal segment data

### Output
- `CNV.tree.txt` - Phylogenetic tree structure
- `gene.LSA.txt` - Gene-level LSA results
- `segmental.LSA.txt` - Chromosomal segment LSA results
- `parallel.LSA.txt` - Parallel evolution events

## Performance Considerations

### Memory Optimization
- Use `medalt_optimized.py` for datasets >5k cells
- Implement chunked processing for distance calculations
- Utilize sparse matrices for large CNV profiles
- Consider PCA-based dimensionality reduction for extremely large datasets

### Statistical Robustness
- Multiple testing corrections (FDR-BH, Bonferroni) implemented
- Bootstrap confidence intervals available
- Permutation tests with configurable iterations
- Effect size calculations for biological significance

## Testing Strategy

The test suite covers:
- Unit tests for core algorithm components
- Integration tests for full pipeline functionality
- Performance validation for optimization features
- Comparison tests against R implementation for consistency

Tests use both synthetic and real biological data to ensure robustness across different scenarios.