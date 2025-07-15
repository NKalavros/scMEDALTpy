# MEDALT Visualization Script
library(ggplot2)
library(pheatmap)

output_dir <- "medalt_2k_results"
cat("Creating MEDALT plots...\n")

# Read distance matrix
dist_data <- read.table(file.path(output_dir, "distance_matrix.txt"), 
                       header=TRUE, row.names=1, sep="\t", comment.char="#")
dist_matrix <- as.matrix(dist_data)

# Plot 1: Distance distribution
distances <- dist_matrix[upper.tri(dist_matrix)]
distances <- distances[distances > 0]
png(file.path(output_dir, "distance_distribution.png"), width=600, height=400)
hist(distances, breaks=50, main="Cell-Cell Distance Distribution",
     xlab="Distance", ylab="Frequency", col="lightblue")
dev.off()

# Plot 2: Distance matrix heatmap (subset)
subset_size <- min(100, nrow(dist_matrix))
png(file.path(output_dir, "distance_heatmap.png"), width=800, height=800)
pheatmap(dist_matrix[1:subset_size, 1:subset_size],
         main=paste("Distance Matrix Heatmap (", subset_size, "x", subset_size, ")"),
         show_rownames=FALSE, show_colnames=FALSE,
         color=colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

# Plot 3: Hierarchical clustering
subset_dist <- as.dist(dist_matrix[1:subset_size, 1:subset_size])
hc <- hclust(subset_dist, method="ward.D2")
png(file.path(output_dir, "dendrogram.png"), width=1000, height=600)
plot(hc, main=paste("Hierarchical Clustering (", subset_size, " cells)"),
     xlab="Cells", ylab="Distance", cex=0.6)
dev.off()

cat("All plots completed!\n")
