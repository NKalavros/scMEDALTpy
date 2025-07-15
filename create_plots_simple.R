# Simple MEDALT Visualization Script
# Install required packages if needed
if (!require(ggplot2, quietly = TRUE)) {
  install.packages("ggplot2", repos="https://cran.rstudio.com/")
  library(ggplot2)
}

output_dir <- "medalt_2k_results"
cat("Creating MEDALT plots...\n")

# Read distance matrix (first few lines to understand structure)
cat("Reading distance matrix...\n")
dist_data <- read.table(file.path(output_dir, "distance_matrix.txt"), 
                       header=TRUE, row.names=1, sep="\t", comment.char="#")
dist_matrix <- as.matrix(dist_data)

cat("Matrix dimensions:", dim(dist_matrix), "\n")

# Plot 1: Basic distance distribution
distances <- dist_matrix[upper.tri(dist_matrix)]
distances <- distances[distances > 0]

cat("Creating distance distribution plot...\n")
png(file.path(output_dir, "distance_distribution.png"), width=600, height=400)
hist(distances, breaks=50, main="Cell-Cell Distance Distribution",
     xlab="Distance", ylab="Frequency", col="lightblue")
dev.off()

# Plot 2: Distance matrix heatmap (small subset for visualization)
subset_size <- min(50, nrow(dist_matrix))
cat("Creating distance heatmap for", subset_size, "x", subset_size, "subset...\n")

png(file.path(output_dir, "distance_heatmap.png"), width=600, height=600)
image(1:subset_size, 1:subset_size, 
      dist_matrix[1:subset_size, 1:subset_size],
      main=paste("Distance Matrix Heatmap (", subset_size, "x", subset_size, ")"),
      xlab="Cell Index", ylab="Cell Index",
      col=heat.colors(100))
dev.off()

# Plot 3: Basic clustering dendrogram
cat("Creating hierarchical clustering dendrogram...\n")
subset_dist <- as.dist(dist_matrix[1:subset_size, 1:subset_size])
hc <- hclust(subset_dist, method="complete")

png(file.path(output_dir, "dendrogram.png"), width=800, height=600)
plot(hc, main=paste("Hierarchical Clustering (", subset_size, " cells)"),
     xlab="Cells", ylab="Distance", cex=0.5)
dev.off()

# Plot 4: Summary statistics
cat("Creating summary statistics plot...\n")
png(file.path(output_dir, "summary_stats.png"), width=600, height=400)
par(mfrow=c(2,1))

# Distance summary
dist_summary <- summary(distances)
barplot(dist_summary, main="Distance Summary Statistics", 
        col="lightgreen", cex.names=0.8)

# Matrix properties
matrix_props <- c(nrow(dist_matrix), ncol(dist_matrix), length(distances))
names(matrix_props) <- c("Rows", "Cols", "Distances")
barplot(matrix_props, main="Matrix Properties", 
        col="lightcoral", log="y")

dev.off()

cat("All plots completed successfully!\n")
cat("Generated files:\n")
cat("- distance_distribution.png\n")
cat("- distance_heatmap.png\n") 
cat("- dendrogram.png\n")
cat("- summary_stats.png\n")