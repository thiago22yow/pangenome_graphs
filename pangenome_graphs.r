############################################
#heatmap of gene intersections among groups#
############################################

# Open a plot screen divided into 2 lines and 2 columns
par(mfrow = c(3, 1))

# load intersection matrix (wide matrix with organisms in first row and first column, values correspond to intersection sizes)
data1 <- read.table('intersection_matrix.txt', sep = '\t', header = TRUE, row.names = 1)

# Function to calculate quantiles and obtain breaks
quantile_breaks <- function(xs, n = 11) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

# Apply the function to the matrix quantiles
data_breaks <- quantile_breaks(as.matrix(data1), n = 11)

# Cluster dendograms
data_cluster_cols <- hclust(dist(t(data1)))
data_cluster_rows <- hclust(dist(data1))

# Order dendograms
library(dendsort)
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
data_cluster_cols <- sort_hclust(data_cluster_cols)
data_cluster_rows <- sort_hclust(data_cluster_rows)

#color
library(RColorBrewer)
cols <- brewer.pal(3, "YlOrRd")
pal <- colorRampPalette(cols)

# Generate heatmap based on quantile breaks
library(pheatmap)
library(viridis)
pheatmap(data1, 
         border_color = 'gray74', 
         cluster_cols = data_cluster_cols, 
         cluster_rows = data_cluster_rows,
         cellwidth = 20,
         color = pal(8),
         breaks = data_breaks,
         main = "Plot 1")

###################################################
# calculate heaps alpha and plot rarefaction curve#
###################################################

library(micropan)
library(vegan)

# load matrix absence presence matrix (genes in columns and organisms in rows, 0 for absence and 1 for presence)
data2 <- read.table('binary_presence_absence_matrix.txt', sep = '\t', header = TRUE, row.names = 1)


# calculate coefficients 'k' and 'α' setting the number of random permutation in 1000
heap <- heaps(data2, n.perm = 1000)

# compute rarefaction curve with 1000 permutations
rf <- specaccum(data2, "random", permutations = 1000)

# plot the accumulation curve with the shaded area representing confidence intervals
plot(rf, ci.type = "poly", col = "darkblue", lwd = 2,
     ci.lty = 0, ci.col = "lightblue3",
     xlab = "Number of genomes",
     ylab = "Number of gene families",
     ylim = c(0, 2000),
     main = "Plot 2")

# set the legend with the value of the 'α' coefficient
legend(x = "bottomright", legend = paste("\u03B1 =", round(heap[2], 2)), fill = "blue")

################
#gene frequency#
################

#transpose the table absence presence (rows must be genes)
df_t <- t(data2)

# get proper breaks
breaks <- seq(min(rowSums(df_t)), max(rowSums(df_t)) + 1, by = 1)

# Plot histogram
hist(rowSums(df_t), xlab = "Number of genomes containing a gene",
     ylab = "Number of genes", main = "Gene frequency",
     ylim = c(0, 800), xlim = c(0, 40),
     breaks = breaks,
     main = "Plot 3")

# Adjust numeric labels on x-axis (every 5 units)
axis(side = 1, at = seq(0, ncol(data) + 1, by = 5))