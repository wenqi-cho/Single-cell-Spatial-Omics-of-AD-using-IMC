mat <- t(assay(spe_filtered, "exprs")[rowData(spe_filtered)$use_channel,])

library(dplyr)
library(pheatmap)
# k=20
clusters2 <- spe_filtered$nn_clusters_corrected
trem2  <- spe_filtered$trem2
# Compute average expression per marker and cluster
# use `mat`
avg_expr <- mat %>%
  data.frame() %>%
  mutate(cluster = clusters2) %>%
  group_by(cluster) %>%
  summarise_all(mean)


# Extract marker names from row names of 'avg_expr'
markers <- rownames(avg_expr)

# Remove the 'cluster' column from 'avg_expr' for plotting
avg_expr_plot <- avg_expr[, !names(avg_expr) %in% "cluster"]

# Scale the marker intensity alongside clusters
heatmap_scale <- function(mat) {
  mat_scaled <- scale(mat, scale = FALSE)
  return(mat_scaled)
}

# Create the heatmap using pheatmap
pheatmap(avg_expr_plot, 
         scale = "column",                          # Disable scaling by rows or columns
         color=colorRampPalette(c("navy", "white", "red"))(50),      
         cluster_rows = FALSE,                
         cluster_cols = FALSE,                     # Enable column clustering
         show_rownames = TRUE,                    # Show row names (cluster IDs)
         show_colnames = TRUE,                    # Show column names (marker names)
         main = "Marker Expression Heatmap - SNN",      
         fontsize_row = 12,                        # Adjust row label font size
         fontsize_col = 12,                        # Adjust column label font size
         #scale_row = heatmap_scale                 # Apply custom scaling function to rows
)

## k=16
clusters16 <- spe_filtered$nn_clusters_corrected16
trem2  <- spe_filtered$trem2

# Compute average expression per marker and cluster
# use `mat`
avg_expr <- mat %>%
  data.frame() %>%
  mutate(cluster = clusters16) %>%
  group_by(cluster) %>%
  summarise_all(mean)


# Extract marker names from row names of 'avg_expr'
markers <- rownames(avg_expr)

# Remove the 'cluster' column from 'avg_expr' for plotting
avg_expr_plot <- avg_expr[, !names(avg_expr) %in% "cluster"]

# Scale the marker intensity alongside clusters
heatmap_scale <- function(mat) {
  mat_scaled <- scale(mat, scale = FALSE)
  return(mat_scaled)
}

# Create the heatmap using pheatmap
pheatmap(avg_expr_plot, 
         scale = "column",                          # Disable scaling by rows or columns
         color=colorRampPalette(c("navy", "white", "red"))(50),      
         cluster_rows = FALSE,                
         cluster_cols = FALSE,                     # Enable column clustering
         show_rownames = TRUE,                    # Show row names (cluster IDs)
         show_colnames = TRUE,                    # Show column names (marker names)
         main = "Marker Expression Heatmap - SNN",      
         fontsize_row = 8,                        # Adjust row label font size
         fontsize_col = 8,                        # Adjust column label font size
         #scale_row = heatmap_scale                 # Apply custom scaling function to rows
)

## k=14
clusters14 <- spe_filtered$nn_clusters_corrected14
trem2  <- spe_filtered$trem2

# Compute average expression per marker and cluster
# use `mat`
avg_expr <- mat %>%
  data.frame() %>%
  mutate(cluster = clusters14) %>%
  group_by(cluster) %>%
  summarise_all(mean)


# Extract marker names from row names of 'avg_expr'
markers <- rownames(avg_expr)

# Remove the 'cluster' column from 'avg_expr' for plotting
avg_expr_plot <- avg_expr[, !names(avg_expr) %in% "cluster"]

# Scale the marker intensity alongside clusters
heatmap_scale <- function(mat) {
  mat_scaled <- scale(mat, scale = FALSE)
  return(mat_scaled)
}

# Create the heatmap using pheatmap
pheatmap(avg_expr_plot, 
         scale = "column",                          # Disable scaling by rows or columns
         color=colorRampPalette(c("navy", "white", "red"))(50),      
         cluster_rows = FALSE,                
         cluster_cols = FALSE,                     # Enable column clustering
         show_rownames = TRUE,                    # Show row names (cluster IDs)
         show_colnames = TRUE,                    # Show column names (marker names)
         main = "Marker Expression Heatmap - SNN",      
         fontsize_row = 12,                        # Adjust row label font size
         fontsize_col = 12,                        # Adjust column label font size
         #scale_row = heatmap_scale                 # Apply custom scaling function to rows
)

###########
# Compute average expression per marker and cluster - SOM
clusters1 <- spe_filtered$som_clusters_corrected
avg_expr1 <- mat %>%
  data.frame() %>%
  mutate(cluster = clusters1) %>%
  group_by(cluster) %>%
  summarise_all(mean)


# Extract marker names from row names of 'avg_expr'
markers <- rownames(avg_expr1)

# Remove the 'cluster' column from 'avg_expr' for plotting
avg_expr_plot1 <- avg_expr1[, !names(avg_expr1) %in% "cluster"]

# Scale the marker intensity alongside clusters
heatmap_scale <- function(mat1) {
  mat_scaled <- scale(mat1, scale = FALSE)
  return(mat_scaled)
}

# Create the heatmap using pheatmap
pheatmap(avg_expr_plot1, 
         scale = "column",                          # Disable scaling by rows or columns
         color=colorRampPalette(c("navy", "white", "red"))(50),      
         cluster_rows = FALSE,                
         cluster_cols = FALSE,                     # Enable column clustering
         show_rownames = TRUE,                    # Show row names (cluster IDs)
         show_colnames = TRUE,                    # Show column names (marker names)
         main = "Marker Expression Heatmap - FlowSOM",      
         fontsize_row = 12,                        # Adjust row label font size
         fontsize_col = 12,                        # Adjust column label font size
         #scale_row = heatmap_scale                 # Apply custom scaling function to rows
)
########### Phenograph
# Compute average expression per marker and cluster - SOM
clusters1 <- spe_filtered$pg_clusters_corrected
avg_expr1 <- mat %>%
  data.frame() %>%
  mutate(cluster = clusters1) %>%
  group_by(cluster) %>%
  summarise_all(mean)


# Extract marker names from row names of 'avg_expr'
markers <- rownames(avg_expr1)

# Remove the 'cluster' column from 'avg_expr' for plotting
avg_expr_plot1 <- avg_expr1[, !names(avg_expr1) %in% "cluster"]

# Scale the marker intensity alongside clusters
heatmap_scale <- function(mat1) {
  mat_scaled <- scale(mat1, scale = FALSE)
  return(mat_scaled)
}

# Create the heatmap using pheatmap
pheatmap(avg_expr_plot1, 
         scale = "column",                          # Disable scaling by rows or columns
         color=colorRampPalette(c("navy", "white", "red"))(50),      
         cluster_rows = FALSE,                
         cluster_cols = FALSE,                     # Enable column clustering
         show_rownames = TRUE,                    # Show row names (cluster IDs)
         show_colnames = TRUE,                    # Show column names (marker names)
         main = "Marker Expression Heatmap - Phenograph",      
         fontsize_row = 12,                        # Adjust row label font size
         fontsize_col = 12,                        # Adjust column label font size
         #scale_row = heatmap_scale                 # Apply custom scaling function to rows
)
###########
# Compute average expression per marker and cluster - SNN
avg_expr2 <- mat %>%
  data.frame() %>%
  mutate(cluster = clusters2) %>%
  group_by(cluster) %>%
  summarise_all(mean)


# Extract marker names from row names of 'avg_expr'
clusters <- rownames(avg_expr2)

# Remove the 'cluster' column from 'avg_expr' for plotting
avg_expr_plot2 <- avg_expr2[, !names(avg_expr2) %in% "cluster"]



range01 <- function(x){(x-min(x))/(max(x)-min(x))}

mat_01= sapply(avg_expr_plot2,range01)
rownames(mat_01)=rownames(avg_expr_plot2)

# Scale the marker intensity alongside clusters
heatmap_scale <- function(mat) {
  mat_scaled <- scale(mat, scale = FALSE)
  return(mat_scaled)
}

# Create the heatmap using pheatmap
pheatmap(avg_expr_plot2, 
         scale = "column",                          # Disable scaling by rows or columns
         color=colorRampPalette(c("navy", "white", "red"))(50),      
         cluster_rows = FALSE,                
         cluster_cols = FALSE,                     # Enable column clustering
         show_rownames = TRUE,                    # Show row names (cluster IDs)
         show_colnames = TRUE,                    # Show column names (marker names)
         main = "Marker Expression Heatmap - SNN",      
         fontsize_row = 8,                        # Adjust row label font size
         fontsize_col = 8,                        # Adjust column label font size
         #scale_row = heatmap_scale                 # Apply custom scaling function to rows
)

pheatmap(mat_01, 
         #scale = "column",                          # Disable scaling by rows or columns
         color=colorRampPalette(c("navy", "white", "red"))(50),      
         cluster_rows = FALSE,                
         cluster_cols = FALSE,                     # Enable column clustering
         show_rownames = TRUE,                    # Show row names (cluster IDs)
         show_colnames = TRUE,                    # Show column names (marker names)
         main = "Marker Expression Heatmap - SNN",      
         fontsize_row = 8,                        # Adjust row label font size
         fontsize_col = 8,                        # Adjust column label font size
         #scale_row = heatmap_scale                 # Apply custom scaling function to rows
)
