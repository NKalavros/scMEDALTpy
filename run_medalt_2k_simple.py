#!/usr/bin/env python2
"""
Run MEDALT pipeline with 2k cells x 1k genes - simplified version
"""

import time
import random
import os

def create_test_data(num_cells=2000, num_genes=1000):
    """Create test data for MEDALT"""
    random.seed(42)
    
    print "Creating test data: %d cells x %d genes" % (num_cells, num_genes)
    
    # Create output directory
    output_dir = "medalt_2k_results"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Generate realistic copy number data
    nodes = {}
    for i in range(num_cells):
        if i % 200 == 0:
            print "  Creating cell %d/%d" % (i, num_cells)
            
        cell_name = "Cell_%04d" % i
        segments = []
        
        for j in range(num_genes):
            # Create realistic copy number with some structure
            if random.random() < 0.02:  # 2% deletions
                cn = 0
            elif random.random() < 0.05:  # 5% losses
                cn = 1
            elif random.random() < 0.7:  # 70% normal
                cn = 2
            elif random.random() < 0.15:  # 15% gains
                cn = 3
            else:  # 8% amplifications
                cn = 4
            
            segments.append([cn])
        
        nodes[cell_name] = segments
    
    print "Created %d cells with %d genes each" % (len(nodes), len(nodes[nodes.keys()[0]]))
    return output_dir, nodes

def run_medalt_analysis(output_dir, nodes):
    """Run MEDALT analysis with plots"""
    print "Running MEDALT analysis..."
    
    # Step 1: Distance calculation
    print "Step 1: Calculating distance matrix..."
    from ComputeDistance import matrixbuilder
    
    start_time = time.time()
    keys, matrix = matrixbuilder(nodes)
    distance_time = time.time() - start_time
    
    print "Distance calculation completed in %.1f seconds" % distance_time
    
    # Save distance matrix
    matrix_file = os.path.join(output_dir, "distance_matrix.txt")
    print "Saving distance matrix to %s" % matrix_file
    
    with open(matrix_file, "w") as f:
        f.write("# Distance matrix for %d cells\n" % len(keys))
        f.write("cell_id\t" + "\t".join(keys) + "\n")
        
        for i, key in enumerate(keys):
            f.write(key)
            for j in range(len(matrix[i])):
                f.write("\t%.6f" % matrix[i][j])
            f.write("\n")
    
    # Step 2: Basic analysis
    print "Step 2: Performing basic analysis..."
    
    # Calculate statistics
    all_distances = []
    for i in range(len(matrix)):
        for j in range(i+1, len(matrix)):
            if matrix[i][j] > 0:
                all_distances.append(matrix[i][j])
    
    # Create analysis summary
    summary_file = os.path.join(output_dir, "analysis_summary.txt")
    with open(summary_file, "w") as f:
        f.write("MEDALT Analysis Summary\n")
        f.write("======================\n\n")
        f.write("Dataset: %d cells x %d genes\n" % (len(nodes), len(nodes[nodes.keys()[0]])))
        f.write("Processing time: %.1f seconds (%.2f minutes)\n" % (distance_time, distance_time / 60.0))
        f.write("Matrix size: %dx%d\n" % (len(matrix), len(matrix[0])))
        f.write("\nDistance Statistics:\n")
        f.write("Total distances: %d\n" % len(all_distances))
        
        if all_distances:
            f.write("Min distance: %.6f\n" % min(all_distances))
            f.write("Max distance: %.6f\n" % max(all_distances))
            f.write("Mean distance: %.6f\n" % (sum(all_distances) / len(all_distances)))
    
    print "Analysis summary saved to %s" % summary_file
    
    # Step 3: Create simple visualizations
    print "Step 3: Creating visualizations..."
    
    # Create R script for plots
    r_script = os.path.join(output_dir, "create_plots.R")
    with open(r_script, "w") as f:
        f.write('# MEDALT Visualization Script\n')
        f.write('library(ggplot2)\n')
        f.write('library(pheatmap)\n\n')
        f.write('output_dir <- "%s"\n' % output_dir)
        f.write('cat("Creating MEDALT plots...\\n")\n\n')
        
        f.write('# Read distance matrix\n')
        f.write('dist_data <- read.table(file.path(output_dir, "distance_matrix.txt"), \n')
        f.write('                       header=TRUE, row.names=1, sep="\\t", comment.char="#")\n')
        f.write('dist_matrix <- as.matrix(dist_data)\n\n')
        
        f.write('# Plot 1: Distance distribution\n')
        f.write('distances <- dist_matrix[upper.tri(dist_matrix)]\n')
        f.write('distances <- distances[distances > 0]\n')
        f.write('png(file.path(output_dir, "distance_distribution.png"), width=600, height=400)\n')
        f.write('hist(distances, breaks=50, main="Cell-Cell Distance Distribution",\n')
        f.write('     xlab="Distance", ylab="Frequency", col="lightblue")\n')
        f.write('dev.off()\n\n')
        
        f.write('# Plot 2: Distance matrix heatmap (subset)\n')
        f.write('subset_size <- min(100, nrow(dist_matrix))\n')
        f.write('png(file.path(output_dir, "distance_heatmap.png"), width=800, height=800)\n')
        f.write('pheatmap(dist_matrix[1:subset_size, 1:subset_size],\n')
        f.write('         main=paste("Distance Matrix Heatmap (", subset_size, "x", subset_size, ")"),\n')
        f.write('         show_rownames=FALSE, show_colnames=FALSE,\n')
        f.write('         color=colorRampPalette(c("blue", "white", "red"))(50))\n')
        f.write('dev.off()\n\n')
        
        f.write('# Plot 3: Hierarchical clustering\n')
        f.write('subset_dist <- as.dist(dist_matrix[1:subset_size, 1:subset_size])\n')
        f.write('hc <- hclust(subset_dist, method="ward.D2")\n')
        f.write('png(file.path(output_dir, "dendrogram.png"), width=1000, height=600)\n')
        f.write('plot(hc, main=paste("Hierarchical Clustering (", subset_size, " cells)"),\n')
        f.write('     xlab="Cells", ylab="Distance", cex=0.6)\n')
        f.write('dev.off()\n\n')
        
        f.write('cat("All plots completed!\\n")\n')
    
    print "R script created: %s" % r_script
    
    # Try to run R script
    print "Attempting to run R script..."
    try:
        import subprocess
        result = subprocess.call(["Rscript", r_script])
        if result == 0:
            print "R plots generated successfully!"
        else:
            print "R script execution completed with return code: %d" % result
    except Exception as e:
        print "R execution failed: %s" % str(e)
        print "R script available for manual execution: %s" % r_script
    
    # Create Python backup plots
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
        
        print "Creating Python backup plots..."
        
        # Distance distribution
        plt.figure(figsize=(10, 6))
        plt.hist(all_distances, bins=50, alpha=0.7, color='lightblue', edgecolor='black')
        plt.title('Cell-Cell Distance Distribution')
        plt.xlabel('Distance')
        plt.ylabel('Frequency')
        plt.grid(True, alpha=0.3)
        plt.savefig(os.path.join(output_dir, 'python_distance_distribution.png'), 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        # Distance matrix heatmap
        subset_size = min(100, len(matrix))
        subset_matrix = np.array(matrix[:subset_size])[:, :subset_size]
        
        plt.figure(figsize=(10, 10))
        plt.imshow(subset_matrix, cmap='viridis', aspect='auto')
        plt.colorbar(label='Distance')
        plt.title('Distance Matrix Heatmap (%dx%d subset)' % (subset_size, subset_size))
        plt.xlabel('Cell Index')
        plt.ylabel('Cell Index')
        plt.savefig(os.path.join(output_dir, 'python_distance_heatmap.png'), 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print "Python plots created successfully!"
        
    except ImportError:
        print "Matplotlib not available - R plots only"
    except Exception as e:
        print "Python plot creation failed: %s" % str(e)
    
    return distance_time, len(all_distances)

def main():
    """Main function"""
    print "MEDALT 2K Cells x 1K Genes Analysis"
    print "=" * 50
    
    # Create test data
    output_dir, nodes = create_test_data(2000, 1000)
    
    # Run analysis
    processing_time, num_distances = run_medalt_analysis(output_dir, nodes)
    
    print "\n" + "=" * 50
    print "MEDALT ANALYSIS COMPLETED!"
    print "=" * 50
    print "Dataset: 2000 cells x 1000 genes"
    print "Processing time: %.1f seconds (%.2f minutes)" % (processing_time, processing_time / 60.0)
    print "Distance calculations: %d" % num_distances
    print "Output directory: %s" % output_dir
    print "\nGenerated files:"
    
    for filename in os.listdir(output_dir):
        print "  - %s" % filename
    
    print "\nAnalysis complete! Check output directory for results and plots."

if __name__ == "__main__":
    main()