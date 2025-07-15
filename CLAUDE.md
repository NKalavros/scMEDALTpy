# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

MEDALT (Minimal Event Distance Aneuploidy Lineage Tree) is a bioinformatics tool for performing lineage tracing using copy number profiles from single cell sequencing technology. It infers rooted directed minimal spanning trees (RDMST) to represent aneuploidy evolution of tumor cells and identifies copy number alterations associated with lineage expansion.

## System Requirements

- **Python 2.7** (legacy codebase)
- **R 3.5+** with required packages:
  - `igraph` - for graph operations and visualization
  - `HelloRanges` - for genomic range operations

## Common Commands

### Basic Usage
```bash
# Run MEDALT on scDNA-seq data
python scTree.py -P ./ -I input_file.txt -D D -G hg19 -O output_directory

# Run MEDALT on scRNA-seq data  
python scTree.py -P ./ -I input_file.txt -D R -G hg19 -O output_directory

# Run with permutation tree reconstruction (slower but more accurate)
python scTree.py -P ./ -I input_file.txt -D D -G hg19 -O output_directory -R T

# Run with custom smoothing window for scRNA-seq data
python scTree.py -P ./ -I input_file.txt -D R -G hg19 -O output_directory -W 50
```

### Testing Examples
```bash
# Run provided examples
cd example/
bash run.example.sh

# Individual example commands
python scTree.py -P ./ -I ./example/scDNA.CNV.txt -D D -G hg19 -O ./example/outputDNA
python scTree.py -P ./ -I ./example/scRNA.CNV.txt -D R -G hg19 -O ./example/outputRNA
```

## Code Architecture

### Main Pipeline (scTree.py)
The entry point that orchestrates the entire MEDALT workflow:
1. **Data preprocessing**: Converts input formats to internal representation
2. **Root identification**: Finds diploid cells or creates artificial root
3. **Tree inference**: Computes minimal spanning tree using MED distance
4. **Permutation analysis**: Generates background distribution for significance testing
5. **LSA analysis**: Identifies CNAs associated with lineage expansion

### Core Components

**Distance Calculation**:
- `ComputeDistance.py`: Implements MED (Minimal Event Distance) algorithm
- `distcalc()`: Core distance function between copy number profiles
- `zerodisthelper()`: Handles missing data (0 values) with penalty

**Tree Construction**:
- `Edmonds.py`: Tree class and graph creation utilities
- `mdmst.py`: Minimal directed spanning tree implementation
- `Readfile.py`: Input file parsing and root node identification

**R Components**:
- `dataTransfer.R`: Converts input formats to segmental copy number profiles
- `LSA.tree.R`: Main LSA analysis and visualization
- `LSAfunction.R`: Statistical functions for lineage speciation analysis
- `permutationCNA.R`: Generates permuted datasets for background estimation

### Data Flow
1. Input → `dataTransfer.R` → Segmental CNV matrix
2. CNV matrix → `Readfile.py` → Node dictionary + root identification
3. Nodes → `Edmonds.py` + `ComputeDistance.py` → Distance graph
4. Distance graph → `mdmst.py` → Minimal spanning tree
5. Tree + CNV → `LSA.tree.R` → Significance analysis + visualization

## Input/Output Formats

### Input Files
- **scDNA-seq**: Tab-separated with columns: chr, pos, cell1, cell2, ...
- **scRNA-seq**: Tab-separated with genes as rows, cells as columns (relative copy numbers)

### Output Files
- `CNV.tree.txt`: Tree structure (parent, child, distance)
- `segmental.LSA.txt`: Broad CNAs with significance scores
- `gene.LSA.txt`: Focal CNAs with significance scores  
- `singlecell.tree.pdf`: Tree visualization
- `LSA.tree.pdf`: CNA association visualization

## Key Parameters

- `-D`: Data type ("D" for DNA-seq, "R" for RNA-seq)
- `-G`: Genome version ("hg19" or "hg38")
- `-W`: Smoothing window size for RNA-seq (default: 30 genes)
- `-R`: Permutation tree reconstruction ("T" for full analysis, "F" for faster background estimation)

## Reference Files

The package includes genome annotation files:
- `gencode_v19_gene_pos.txt` / `gencode_v38_gene_pos.txt`: Gene position annotations
- `hg19.band.bed` / `hg38.band.bed`: Chromosomal band annotations
- `pathwaygene.txt`: Pathway gene annotations