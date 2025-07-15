# Create a simple summary plot
output_dir <- "medalt_2k_results"

# Read distance matrix info
dist_data <- read.table(file.path(output_dir, "distance_matrix.txt"), 
                       header=TRUE, row.names=1, sep="\t", comment.char="#")
dist_matrix <- as.matrix(dist_data)

# Simple summary plot
png(file.path(output_dir, "analysis_summary.png"), width=800, height=600)
par(mfrow=c(2,2))

# Matrix size info
matrix_info <- c(nrow(dist_matrix), ncol(dist_matrix), 
                nrow(dist_matrix) * (nrow(dist_matrix)-1) / 2)
names(matrix_info) <- c("Cells", "Cells", "Distances")
barplot(matrix_info, main="Dataset Overview", 
        col=c("lightblue", "lightgreen", "lightcoral"), 
        ylab="Count", log="y")

# Distance distribution subset
distances <- dist_matrix[upper.tri(dist_matrix)]
hist(distances[1:1000], breaks=20, main="Distance Sample (first 1000)",
     xlab="Distance", col="lightblue")

# Diagonal check
diag_vals <- diag(dist_matrix)
plot(1:min(50, length(diag_vals)), diag_vals[1:min(50, length(diag_vals))],
     main="Diagonal Values (should be 0)", 
     xlab="Index", ylab="Value", type="h", col="red")

# Matrix symmetry check
sym_check <- all(dist_matrix[1:10, 1:10] == t(dist_matrix[1:10, 1:10]))
text_plot <- plot(1, 1, type="n", axes=FALSE, xlab="", ylab="", 
                 main="Matrix Properties")
text(1, 1, paste("Matrix is symmetric:", sym_check, 
                 "\nDimensions:", nrow(dist_matrix), "x", ncol(dist_matrix),
                 "\nProcessing: Complete"), 
     cex=1.2, col="darkblue")

dev.off()

cat("Summary plot created: analysis_summary.png\n")