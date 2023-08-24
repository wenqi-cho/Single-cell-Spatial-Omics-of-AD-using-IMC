library(CATALYST)
library(SpatialExperiment)
library(SummarizedExperiment)
library(scater)
library(imcRtools)
library(dittoSeq)

# spe_filtered <- readRDS("spe_filtered.rds")


## Run FlowSOM and ConsensusClusterPlus clustering

####
## ConsensusClusterPlus: perform hierarchical consensus clustering 
library(kohonen)
library(ConsensusClusterPlus)

# Select integrated cells
mat <- reducedDim(spe_filtered, "fastMNN")
library(scran)
library(bluster)
set.seed(1234)
# Perform SOM clustering
som.out <- clusterRows(mat, SomParam(100), full = TRUE)

# Cluster the 100 SOM codes into larger clusters
ccp <- ConsensusClusterPlus(t(som.out$objects$som$codes[[1]]),
                            maxK = 20,
                            reps = 50, 
                            distance = "euclidean", 
                            seed = 220410, 
                            plot = NULL)

# Visualize delta area plot
CATALYST:::.plot_delta_area(ccp)

# Link ConsensusClusterPlus clusters with SOM codes and save in object
som.cluster <- ccp[[12]][["consensusClass"]][som.out$clusters]
spe_filtered$som_clusters_corrected <- as.factor(som.cluster)

dittoDimPlot(spe_filtered, var = "som_clusters_corrected", 
             reduction.use = "UMAP_mnnCorrected2", size = 0.4,
             do.label = TRUE) +
  ggtitle("FlowSOM clusters expression on UMAP, integrated cells")

# Sample cells
set.seed(220619)
cur_cells <- sample(seq_len(ncol(spe_filtered)), 10000)
library(viridis)
dittoHeatmap(spe_filtered[,cur_cells], 
             genes = rownames(spe_filtered)[rowData(spe_filtered)$use_channel],
             assay = "exprs", scale = "none",
             heatmap.colors = viridis(100), 
             annot.by = c("som_clusters_corrected", "diagnosis"))

dittoHeatmap(spe_filtered[,cur_cells], 
             genes = rownames(spe_filtered)[rowData(spe_filtered)$use_channel],
             assay = "exprs", scale = "none",
             heatmap.colors = viridis(100), 
             annot.by = c("som_clusters_corrected", "diagnosis","trem2"))


table(colData(spe_filtered)$som_clusters_corrected,colData(spe_filtered)$diagnosis)
df=as.data.frame.matrix(table(colData(spe_filtered)$som_clusters_corrected,colData(spe_filtered)$OnPlaqueObject))

df$perc=df$OnPlaque/(df$OnPlaque+df$OutsidePlaque)
meta=colData(spe_filtered)

data.frame(meta) %>%
  group_by(som_clusters_corrected) %>%
  summarise(n=n()/71165)


table(colData(spe_filtered)$som_clusters_corrected,colData(spe_filtered)$distance)

###
## PhenoGraph clustering approach: constructs a graph by detecting the k nearest neighbours based on
## euclidean distance in expression space. Jaccard index is used to quantify the overlap in shared neighbour sets. 
## Louvain modularity optimisation is used to detect connected communities and partition the graph into clusters of cells.
library(Rphenoannoy)
library(igraph)
library(dittoSeq)
library(viridis)

# mat1 <- t(assay(spe_filtered, "exprs")[rowData(spe_filtered)$use_channel,])
# using integrated cells: mat <- reducedDim(spe_filtered, "fastMNN")
out <- Rphenoannoy(mat, k = 60)

clusters <- factor(membership(out[[2]]))

spe_filtered$pg_clusters_corrected <- clusters

dittoDimPlot(spe_filtered, var = "pg_clusters_corrected", 
             reduction.use = "UMAP_mnnCorrected2", size = 0.4,
             do.label = TRUE) +
  ggtitle("Phenograph clusters expression on UMAP, integrated cells")

# dittoHeatmap(spe_filtered[,cur_cells], 
#              genes = rownames(spe_filtered)[rowData(spe_filtered)$use_channel],
#              assay = "exprs", scale = "none",
#              heatmap.colors = viridis(100), 
#              annot.by = c("pg_clusters_corrected","diagnosis", "trem2"))


###
## Compare between clustering approaches
library(patchwork)
library(pheatmap)
library(gridExtra)

tab1 <- table(paste("Rphenograph", spe_filtered$pg_clusters_corrected),
              paste("SOM", spe_filtered$som_clusters_corrected))

pheatmap(log10(tab1+10), color=viridis(100))
