library(imcRtools)

spe=read_steinbock(path="~/summer-project/combined_data/mtg_dna/deepcell")

spe # dim 20 73054
counts(spe)[1:5,1:5]
head(colData(spe))
head(spatialCoords(spe))
colPair(spe, "neighborhood")
head(rowData(spe)) ## information extracted from the panel.csv file 
# rownames(spe@assays@data@listData[["counts"]])

# # Reading custom files
# ### When not using steinbock, the single-cell information has to be read in from custom files. Generate a SpatialExperiment object from single-cell data contained in individual files.
# library(readr)
# cur_intensities <- read_csv("~/mtg_dna/deepcell/intensities/19891053_001.csv")
# cur_regionprops <- read_csv("~/mtg_dna/deepcell/regionprops/19891053_001.csv")
# dim(cur_intensities) # 442  21
# dim(cur_regionprops) # 442   7
# colnames(cur_intensities)
# # [1] "Object" "80ArAr" "127I"   "131Xe"  "134Xe"  "138Ba"  "Ab"     "ApoE"   "CD163"  "TSPO"   "NfL"    "CD68"   "MAP2"  
# # [14] "GFAP"   "pTau"   "Iba1"   "PLP"    "HLA-DR" "DNA 1"  "DNA 2"  "208Pb"
# colnames(cur_regionprops)
# # [1] "Object"            "area"              "centroid-0"        "centroid-1"        "axis_major_length" "axis_minor_length"
# # [7] "eccentricity"
# counts <- cur_intensities[,-1]
# meta <- cur_regionprops[,c("Object", "area", "axis_major_length", "axis_minor_length")]
# coords <- cur_regionprops[,c("centroid-1", "centroid-0")]
# # construct the SpatialExperiment object.
# library(SpatialExperiment)
# spe2 <- SpatialExperiment(assays = list(counts = t(counts)),
#                           colData = meta, 
#                           sample_id = "19891053_001",
#                           spatialCoords = as.matrix(coords))
# 
# cur_pairs <- read_csv("~/mtg_dna/deepcell/neighbors/19891053_001.csv")
# 
# edgelist <- SelfHits(from = match(cur_pairs$Object, spe2$Object),
#                      to = match(cur_pairs$Neighbor, spe2$Object),
#                      nnode = ncol(spe2))
# 
# colPair(spe2, "neighborhood") <- edgelist

# cur_cells <- sample(seq_len(ncol(spe)), 2000)
# dittoHeatmap(spe[,cur_cells], genes = rownames(spe)[rowData(spe)$use_channel],
#              assay = "exprs", cluster_cols = TRUE, scale = "none",
#              heatmap.colors = viridis(100), annot.by = "indication",
#              annotation_colors = list(indication = metadata(spe)$color_vectors$indication))

# Single-cell processing
## Add additional metadata
### set the colnames of the object to generate unique identifiers per cell:
colnames(spe) <- paste0(spe$sample_id, "_", spe$ObjectNumber)
library(stringr)
library(readr)
# metadata <- read_csv("~/mtg_dna/deepcell/metadata.csv")
# metadata_new <- readRDS("metadata_new.rds")
# metadata_filtered <- readRDS("metadata_filtered.rds")
# metadata_filtered2 <- readRDS("metadata_filtered2.rds")

# I want to extract only the sample ids, without the trailing '_001', '_002', etc. 
spe$patient_id <- str_extract(spe$sample_id, "^[^_]+")

metadata_filtered=metadata[which(metadata$CaseID %in% unique(spe$patient_id)),]
table(metadata_filtered$diagnosis)
# table(metadata$diagnosis)
dim(metadata_filtered)
# dim(metadata)

## Group Braak column into bins
# Convert "Braak" column to numeric, keeping NAs as is
data_clean <- as.numeric(as.character(metadata_filtered$Braak))

# Define the bin boundaries
bins <- c(0, 2, 4, 6)

# Create a new variable "BraakGroup" based on the bins (ignoring "N/A" values)
metadata_filtered$BraakGroup <- cut(data_clean, bins, labels = c("Braak 0-2", "Braak 3-4", "Braak 5-6"), include.lowest = TRUE)

### ADDING METADATA
spe$patient_id <- as.vector(spe$patient_id)
spe$ROI <- as.vector(str_extract(spe$sample_id, "00[1-4]"))
spe$trem2 <- as.factor(metadata_filtered$trem2[match(spe$patient_id, metadata_filtered$CaseID)])
spe$diagnosis <- as.factor(metadata_filtered$diagnosis[match(spe$patient_id, metadata_filtered$CaseID)])
spe$sex <- as.factor(metadata_filtered$sex[match(spe$patient_id, metadata_filtered$CaseID)])
spe$age <- as.numeric(metadata_filtered$age[match(spe$patient_id, metadata_filtered$CaseID)])
spe$PMD <- as.numeric(metadata_filtered$PMD[match(spe$patient_id, metadata_filtered$CaseID)])
spe$trem2_all <- as.factor(metadata_filtered$trem2_all[match(spe$patient_id, metadata_filtered$CaseID)])
spe$apoe <- as.factor(metadata_filtered$apoe_group[match(spe$patient_id, metadata_filtered$CaseID)])
spe$amyloid <- as.numeric(metadata_filtered$amyloid_beta[match(spe$patient_id, metadata_filtered$CaseID)])
spe$pTau <- as.numeric(metadata_filtered$pTau[match(spe$patient_id, metadata_filtered$CaseID)])
spe$pct_pTau_pos <- as.numeric(metadata_filtered$pct_pTau_pos[match(spe$patient_id, metadata_filtered$CaseID)])
spe$PHF1 <- as.numeric(metadata_filtered$PHF1[match(spe$patient_id, metadata_filtered$CaseID)])
spe$diagnosis=as.factor(spe$diagnosis)
spe$Braak <- as.numeric(metadata_filtered$Braak[match(spe$patient_id, metadata_filtered$CaseID)])
spe$BraakGroup <- as.factor(metadata_filtered$BraakGroup[match(spe$patient_id, metadata_filtered$CaseID)])

# saveRDS(spe,'spe.rds')


### Transform counts
## Ridge plots by AD/Control; variants (AD only ); by variants (Control only);
## and also grouping by Braak (create a new variable to bin Braak 0-II , Braak 3-4 and braak 5-6)
library(dittoSeq)
# colnames(spe) <- paste0(spe$sample_id, "_", spe$ObjectNumber)
# dittoRidgePlot(spe, var = "", group.by = "sample_id", assay = "counts")
assay(spe, "exprs") <- asinh(counts(spe)/1)
dittoRidgePlot(spe, var = "MAP2", group.by = "patient_id", assay = "exprs") +
  ggtitle("MAP2 - after transformation")
dittoRidgePlot(spe, var = "GFAP", group.by = "patient_id", assay = "exprs") +
  ggtitle("GFAP - after transformation")
dittoRidgePlot(spe, var = "CD68", group.by = "patient_id", assay = "exprs") +
  ggtitle("CD68 - after transformation")

## Ab
dittoRidgePlot(spe, var = "Ab", group.by = "diagnosis", assay = "exprs") +
  ggtitle("Ab - after transformation")
dittoRidgePlot(spe, var = "Ab", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis') +
  ggtitle("Ab - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe, var = "Ab", group.by = "trem2", assay = "exprs", split.by='diagnosis') +
  ggtitle("Ab - after transformation")
dittoRidgePlot(spe, var = "Ab", group.by = "trem2_all", assay = "exprs", split.by='diagnosis') +
  ggtitle("Ab - after transformation")


## CD68
dittoRidgePlot(spe, var = "CD68", group.by = "diagnosis", assay = "exprs") +
  ggtitle("CD68 - after transformation")
dittoRidgePlot(spe, var = "CD68", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis') +
  ggtitle("CD68 - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe, var = "CD68", group.by = "trem2", assay = "exprs", split.by='diagnosis') +
  ggtitle("CD68 - after transformation")
dittoRidgePlot(spe, var = "CD68", group.by = "trem2_all", assay = "exprs", split.by='diagnosis') +
  ggtitle("CD68 - after transformation")

## "ApoE"   
dittoRidgePlot(spe, var = "ApoE", group.by = "diagnosis", assay = "exprs") +
  ggtitle("ApoE - after transformation")
dittoRidgePlot(spe, var = "ApoE", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis') +
  ggtitle("ApoE - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe, var = "ApoE", group.by = "trem2", assay = "exprs", split.by='diagnosis') +
  ggtitle("ApoE - after transformation")
dittoRidgePlot(spe, var = "ApoE", group.by = "trem2_all", assay = "exprs", split.by='diagnosis') +
  ggtitle("ApoE - after transformation")

## "NfL"
dittoRidgePlot(spe, var = "NfL", group.by = "diagnosis", assay = "exprs") +
  ggtitle("NfL - after transformation")
dittoRidgePlot(spe, var = "NfL", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis') +
  ggtitle("NfL - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe, var = "NfL", group.by = "trem2", assay = "exprs", split.by='diagnosis') +
  ggtitle("NfL - after transformation")
dittoRidgePlot(spe, var = "NfL", group.by = "trem2_all", assay = "exprs", split.by='diagnosis') +
  ggtitle("NfL - after transformation")

## "GFAP"
dittoRidgePlot(spe, var = "GFAP", group.by = "diagnosis", assay = "exprs") +
  ggtitle("GFAP - after transformation")
dittoRidgePlot(spe, var = "GFAP", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis') +
  ggtitle("GFAP - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe, var = "GFAP", group.by = "trem2", assay = "exprs", split.by='diagnosis') +
  ggtitle("GFAP - after transformation")
dittoRidgePlot(spe, var = "GFAP", group.by = "trem2_all", assay = "exprs", split.by='diagnosis') +
  ggtitle("GFAP - after transformation")

## "MAP2"
dittoRidgePlot(spe, var = "MAP2", group.by = "diagnosis", assay = "exprs") +
  ggtitle("MAP2 - after transformation")
dittoRidgePlot(spe, var = "MAP2", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis') +
  ggtitle("MAP2 - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe, var = "MAP2", group.by = "trem2", assay = "exprs", split.by='diagnosis') +
  ggtitle("MAP2 - after transformation")
dittoRidgePlot(spe, var = "MAP2", group.by = "trem2_all", assay = "exprs", split.by='diagnosis') +
  ggtitle("MAP2 - after transformation")

## "HLA-DR"
dittoRidgePlot(spe, var = "HLA-DR", group.by = "diagnosis", assay = "exprs") +
  ggtitle("HLA-DR - after transformation")
dittoRidgePlot(spe, var = "HLA-DR", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis') +
  ggtitle("HLA-DR - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe, var = "HLA-DR", group.by = "trem2", assay = "exprs", split.by='diagnosis') +
  ggtitle("HLA-DR - after transformation")
dittoRidgePlot(spe, var = "HLA-DR", group.by = "trem2_all", assay = "exprs", split.by='diagnosis') +
  ggtitle("HLA-DR - after transformation")

## "Iba1"
dittoRidgePlot(spe, var = "Iba1", group.by = "diagnosis", assay = "exprs") +
  ggtitle("Iba1 - after transformation")
dittoRidgePlot(spe, var = "Iba1", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis') +
  ggtitle("Iba1 - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe, var = "Iba1", group.by = "trem2", assay = "exprs", split.by='diagnosis') +
  ggtitle("Iba1 - after transformation")
dittoRidgePlot(spe, var = "Iba1", group.by = "trem2_all", assay = "exprs", split.by='diagnosis') +
  ggtitle("Iba1 - after transformation")

## "CD163"
dittoRidgePlot(spe, var = "CD163", group.by = "diagnosis", assay = "exprs") +
  ggtitle("CD163 - after transformation")
dittoRidgePlot(spe, var = "CD163", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis') +
  ggtitle("CD163 - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe, var = "CD163", group.by = "trem2", assay = "exprs", split.by='diagnosis') +
  ggtitle("CD163 - after transformation")
dittoRidgePlot(spe, var = "CD163", group.by = "trem2_all", assay = "exprs", split.by='diagnosis') +
  ggtitle("CD163 - after transformation")

## "TSPO"
dittoRidgePlot(spe, var = "TSPO", group.by = "diagnosis", assay = "exprs") +
  ggtitle("TSPO - after transformation")
dittoRidgePlot(spe, var = "TSPO", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis') +
  ggtitle("TSPO - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe, var = "TSPO", group.by = "trem2", assay = "exprs", split.by='diagnosis') +
  ggtitle("TSPO - after transformation")
dittoRidgePlot(spe, var = "TSPO", group.by = "trem2_all", assay = "exprs", split.by='diagnosis') +
  ggtitle("TSPO - after transformation")

### Define interesting channels
#### Store an additional entry in the rowData slot that specifies the markers of interest. Deselect nuclear markers which were used for cell segmentation, and keep all other biological targets.
rowData(spe)$use_channel <- !grepl("DNA 1|DNA 2", rownames(spe))

### Define color schemes
# library(RColorBrewer)
# color_vectors <- list()
# ROI <- setNames(brewer.pal(length(unique(spe$ROI)), name = "BrBG"), 
#                 unique(spe$ROI))
# patient_id <- setNames(brewer.pal(length(unique(spe$patient_id)), name = "Set1"), 
#                        unique(spe$patient_id))
# sample_id <- setNames(c(brewer.pal(6, "YlOrRd")[3:5],
#                         brewer.pal(6, "PuBu")[3:6],
#                         brewer.pal(6, "YlGn")[3:5],
#                         brewer.pal(6, "BuPu")[3:6]),
#                       unique(spe$sample_id))
# indication <- setNames(brewer.pal(length(unique(spe$indication)), name = "Set2"), 
#                        unique(spe$indication))
# 
# color_vectors$ROI <- ROI
# color_vectors$patient_id <- patient_id
# color_vectors$sample_id <- sample_id
# color_vectors$indication <- indication
# 
# metadata(spe)$color_vectors <- color_vectors

## Read in images
library(cytomapper)
# images <- loadImages("~/mtg_dna/deepcell/img/")
images <- readRDS("saved_images.rds")
# masks <- loadImages("~/mtg_dna/deepcell/masks/", as.is = TRUE)
masks <- readRDS("masks.rds")
# masks_plaques <- loadImages("masks_plaques/", pattern = ".tiff", as.is = TRUE)
# saveRDS(masks_plaques, "masks_plaques.rds")
masks_plaques <- readRDS("masks_plaques.rds")
channelNames(images) <- rownames(spe)
images
all.equal(names(images), names(masks))

patient_id <- str_extract(names(images), "^[^_]+")
diagnosis <- metadata_filtered$diagnosis[match(patient_id, metadata_filtered$CaseID)]
trem2 <- metadata_filtered$trem2[match(patient_id, metadata_filtered$CaseID)]
# Save the matched sample_id, patient_id, diagnosis, and trem2 information within the elementMetadata slot of
# of the multi-channel images and segmentation masks objects.
mcols(images) <- mcols(masks) <- mcols(masks_plaques) <- DataFrame(sample_id = names(images),
                                                                  patient_id = patient_id,
                                                                  diagnosis = diagnosis,
                                                                  trem2 = trem2)

## Generate single-cell data from images
sce <- measureObjects(masks, image = images, img_id = "sample_id")
# saveRDS(sce, "sce.rds")

### Single-cell visualisation
#### Dimensionality reduction
library(scater)

set.seed(220225)
spe <- runUMAP(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs") 

reducedDims(spe)
spe@int_colData$reducedDims$UMAP
library(ggplot2)
library(patchwork)

library(RColorBrewer)
diagnosis_color <- setNames(brewer.pal(length(unique(spe$diagnosis)), name='Set2'),
                            unique(spe$diagnosis))

spe$diag_color=ifelse(spe$diagnosis=='AD',"#66C2A5","#FC8D62")

metadata(spe)
color_vectors=list()
color_vectors$diagnosis <- list(AD="#66C2A5",
                                Control= "#FC8D62")

metadata(spe)$color_vectors <- color_vectors
# visualise patient id
dittoDimPlot(spe,var='patient_id',reduction.use='UMAP')
#dittoDimPlot(spe,var='patient_id',reduction.use='UMAP')+scale_color_manual(values=metadata(spe)$color_vectors$diagnosis)

# visualise diagnosis
dittoDimPlot(spe, var = "diagnosis", reduction.use = "UMAP", size = 0.2) + 
  scale_color_manual(values = metadata(spe)$color_vectors$diagnosis) +
  ggtitle("Diagnosis on UMAP")
library(viridis)
# visualise marker expression
p1 <- dittoDimPlot(spe, var = "CD68", reduction.use = "UMAP", 
                   assay = "exprs", size = 0.2) +
  scale_color_viridis(name = "CD68") +
  ggtitle("CD68 expression on UMAP")
p2 <- dittoDimPlot(spe, var = "GFAP", reduction.use = "UMAP", 
                   assay = "exprs", size = 0.2) +
  scale_color_viridis(name = "GFAP") +
  ggtitle("GFAP expression on UMAP")
p3 <- dittoDimPlot(spe, var = "MAP2", reduction.use = "UMAP", 
                   assay = "exprs", size = 0.2) +
  scale_color_viridis(name = "MAP2") +
  ggtitle("MAP2 expression on UMAP")
p4 <- dittoDimPlot(spe, var = "Iba1", reduction.use = "UMAP", 
                   assay = "exprs", size = 0.2) +
  scale_color_viridis(name = "Iba1") +
  ggtitle("Iba1 expression on UMAP")
p5 <- dittoDimPlot(spe, var = "HLA-DR", reduction.use = "UMAP", 
                   assay = "exprs", size = 0.2) +
  scale_color_viridis(name = "HLA-DR") +
  ggtitle("HLA-DR expression on UMAP")
p6 <- dittoDimPlot(spe, var = "TSPO", reduction.use = "UMAP", 
                   assay = "exprs", size = 0.2) +
  scale_color_viridis(name = "TSPO") +
  ggtitle("TSPO expression on UMAP")


(p1 + p2)/ (p3 + p4)/ (p5 + p6)
# Visualise single-cell marker expression in form of heatmap.
# per cell
cur_cells <- sample(seq_len(ncol(spe)),2000)
dittoHeatmap(spe[,cur_cells], assay='exprs', cluster_cols=TRUE,annot.by = 'trem2_all',
             heatmap.colors = viridis(100))
dittoHeatmap(spe[,cur_cells], assay='exprs', cluster_cols=TRUE,annot.by = 'trem2',
             heatmap.colors = viridis(100))
dittoHeatmap(spe[,cur_cells], assay='exprs', cluster_cols=TRUE,annot.by = 'diagnosis',
             heatmap.colors = viridis(100))
# table(spe$diag_color)
# #66C2A5 #FC8D62 
# 50421   22633

# Image area covered by cells
library(dplyr)

colData(spe) %>%
  as.data.frame() %>%
  group_by(sample_id) %>%
  summarize(cell_area = sum(area),
            no_pixels = mean(width_px) * mean(height_px)) %>%
  mutate(covered_area = cell_area / no_pixels) %>%
  ggplot() +
  geom_point(aes(reorder(sample_id,covered_area), covered_area)) + 
  theme_minimal(base_size = 15) +
  ylim(c(0, 1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  ylab("% covered area") + xlab("")

# % covered area by cells, only showing the top 50 here
# colData(spe) %>%
#   as.data.frame() %>%
#   group_by(sample_id) %>%
#   summarize(cell_area = sum(area),
#             no_pixels = mean(width_px) * mean(height_px)) %>%
#   mutate(covered_area = cell_area / no_pixels) %>%
#   top_n(50, covered_area) %>%
#   ggplot() +
#   geom_point(aes(reorder(sample_id, covered_area), covered_area)) + 
#   theme_minimal(base_size = 15) +
#   ylim(c(0, 1)) + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
#   ylab("% covered area") + xlab("")

# Cell size
colData(spe) %>%
  as.data.frame() %>%
  group_by(sample_id) %>%
  ggplot() +
  geom_boxplot(aes(sample_id, area)) +
  theme_minimal(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  ylab("Cell area") + xlab("")

summary(spe$area)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00   22.00   53.00   68.36   98.00  875.00 

sum(spe$area < 5)
# [1] 1445

spe <- spe[, spe$area >= 5]

## Spatial visualisation
library(ggplot2)
library(viridis)

patient_ids <- unique(metadata_filtered$CaseID)
# Steinbock interaction graph 
# Custom width and height for the plots
plot_width <- 10 
plot_height <- 10
## MAP2
# Loop through each patient ID
for (patient_id in patient_ids) {
  # Generate the plot
  plot <- plotSpatial(spe[, grepl(patient_id, spe$sample_id)], 
                      node_color_by = "MAP2", 
                      assay_type = "exprs",
                      img_id = "sample_id", 
                      draw_edges = TRUE, 
                      colPairName = "neighborhood", 
                      nodes_first = FALSE, 
                      node_size_by = "area", 
                      directed = FALSE,
                      edge_color_fix = "grey") + 
    scale_size_continuous(range = c(0.1, 2)) +
    ggtitle("MAP2 expression")
  
  # Export the plot as a PDF file
  filename <- paste0("patient_", patient_id, ".pdf")
  ggsave(filename, plot, width = plot_width, height = plot_height, device = "pdf")
}

## GFAP
# Loop through each patient ID
for (patient_id in patient_ids) {
  # Generate the plot
  plot <- plotSpatial(spe[, grepl(patient_id, spe$sample_id)], 
                      node_color_by = "GFAP", 
                      assay_type = "exprs",
                      img_id = "sample_id", 
                      draw_edges = TRUE, 
                      colPairName = "neighborhood", 
                      nodes_first = FALSE, 
                      node_size_by = "area", 
                      directed = FALSE,
                      edge_color_fix = "grey") + 
    scale_size_continuous(range = c(0.1, 2)) +
    ggtitle("GFAP expression")
  
  # Export the plot as a PDF file
  filename <- paste0("patient_", patient_id, ".pdf")
  ggsave(filename, plot, width = plot_width, height = plot_height, device = "pdf")
}

## Iba1
# Loop through each patient ID
for (patient_id in patient_ids) {
  # Generate the plot
  plot <- plotSpatial(spe[, grepl(patient_id, spe$sample_id)], 
                      node_color_by = "Iba1", 
                      assay_type = "exprs",
                      img_id = "sample_id", 
                      draw_edges = TRUE, 
                      colPairName = "neighborhood", 
                      nodes_first = FALSE, 
                      node_size_by = "area", 
                      directed = FALSE,
                      edge_color_fix = "grey") + 
    scale_size_continuous(range = c(0.1, 2)) +
    ggtitle("Iba1 expression")
  
  # Export the plot as a PDF file
  filename <- paste0("patient_", patient_id, ".pdf")
  ggsave(filename, plot, width = plot_width, height = plot_height, device = "pdf")
}

## CD68
# Loop through each patient ID
for (patient_id in patient_ids) {
  # Generate the plot
  plot <- plotSpatial(spe[, grepl(patient_id, spe$sample_id)], 
                      node_color_by = "CD68", 
                      assay_type = "exprs",
                      img_id = "sample_id", 
                      draw_edges = TRUE, 
                      colPairName = "neighborhood", 
                      nodes_first = FALSE, 
                      node_size_by = "area", 
                      directed = FALSE,
                      edge_color_fix = "grey") + 
    scale_size_continuous(range = c(0.1, 2)) +
    ggtitle("CD68 expression")
  
  # Export the plot as a PDF file
  filename <- paste0("patient_", patient_id, ".pdf")
  ggsave(filename, plot, width = plot_width, height = plot_height, device = "pdf")
}

## CD163
# Loop through each patient ID
for (patient_id in patient_ids) {
  # Generate the plot
  plot <- plotSpatial(spe[, grepl(patient_id, spe$sample_id)], 
                      node_color_by = "CD163", 
                      assay_type = "exprs",
                      img_id = "sample_id", 
                      draw_edges = TRUE, 
                      colPairName = "neighborhood", 
                      nodes_first = FALSE, 
                      node_size_by = "area", 
                      directed = FALSE,
                      edge_color_fix = "grey") + 
    scale_size_continuous(range = c(0.1, 2)) +
    ggtitle("CD163 expression")
  
  # Export the plot as a PDF file
  filename <- paste0("patient_", patient_id, ".pdf")
  ggsave(filename, plot, width = plot_width, height = plot_height, device = "pdf")
}

## HLA-DR
# Loop through each patient ID
for (patient_id in patient_ids) {
  # Generate the plot
  plot <- plotSpatial(spe[, grepl(patient_id, spe$sample_id)], 
                      node_color_by = "HLA-DR", 
                      assay_type = "exprs",
                      img_id = "sample_id", 
                      draw_edges = TRUE, 
                      colPairName = "neighborhood", 
                      nodes_first = FALSE, 
                      node_size_by = "area", 
                      directed = FALSE,
                      edge_color_fix = "grey") + 
    scale_size_continuous(range = c(0.1, 2)) +
    ggtitle("HLA-DR expression")
  
  # Export the plot as a PDF file
  filename <- paste0("patient_", patient_id, ".pdf")
  ggsave(filename, plot, width = plot_width, height = plot_height, device = "pdf")
}

## Ab
# Loop through each patient ID
for (patient_id in patient_ids) {
  # Generate the plot
  plot <- plotSpatial(spe[, grepl(patient_id, spe$sample_id)], 
                      node_color_by = "Ab", 
                      assay_type = "exprs",
                      img_id = "sample_id", 
                      draw_edges = TRUE, 
                      colPairName = "neighborhood", 
                      nodes_first = FALSE, 
                      node_size_by = "area", 
                      directed = FALSE,
                      edge_color_fix = "grey") + 
    scale_size_continuous(range = c(0.1, 2)) +
    ggtitle("Ab expression")
  
  # Export the plot as a PDF file
  filename <- paste0("patient_", patient_id, ".pdf")
  ggsave(filename, plot, width = plot_width, height = plot_height, device = "pdf")
}

## ApoE
# Loop through each patient ID
for (patient_id in patient_ids) {
  # Generate the plot
  plot <- plotSpatial(spe[, grepl(patient_id, spe$sample_id)], 
                      node_color_by = "ApoE", 
                      assay_type = "exprs",
                      img_id = "sample_id", 
                      draw_edges = TRUE, 
                      colPairName = "neighborhood", 
                      nodes_first = FALSE, 
                      node_size_by = "area", 
                      directed = FALSE,
                      edge_color_fix = "grey") + 
    scale_size_continuous(range = c(0.1, 2)) +
    ggtitle("ApoE expression")
  
  # Export the plot as a PDF file
  filename <- paste0("patient_", patient_id, ".pdf")
  ggsave(filename, plot, width = plot_width, height = plot_height, device = "pdf")
}

## TSPO
# Loop through each patient ID
for (patient_id in patient_ids) {
  # Generate the plot
  plot <- plotSpatial(spe[, grepl(patient_id, spe$sample_id)], 
                      node_color_by = "TSPO", 
                      assay_type = "exprs",
                      img_id = "sample_id", 
                      draw_edges = TRUE, 
                      colPairName = "neighborhood", 
                      nodes_first = FALSE, 
                      node_size_by = "area", 
                      directed = FALSE,
                      edge_color_fix = "grey") + 
    scale_size_continuous(range = c(0.1, 2)) +
    ggtitle("TSPO expression")
  
  # Export the plot as a PDF file
  filename <- paste0("patient_", patient_id, ".pdf")
  ggsave(filename, plot, width = plot_width, height = plot_height, device = "pdf")
}

## NfL
# Loop through each patient ID
for (patient_id in patient_ids) {
  # Generate the plot
  plot <- plotSpatial(spe[, grepl(patient_id, spe$sample_id)], 
                      node_color_by = "NfL", 
                      assay_type = "exprs",
                      img_id = "sample_id", 
                      draw_edges = TRUE, 
                      colPairName = "neighborhood", 
                      nodes_first = FALSE, 
                      node_size_by = "area", 
                      directed = FALSE,
                      edge_color_fix = "grey") + 
    scale_size_continuous(range = c(0.1, 2)) +
    ggtitle("NfL expression")
  
  # Export the plot as a PDF file
  filename <- paste0("patient_", patient_id, ".pdf")
  ggsave(filename, plot, width = plot_width, height = plot_height, device = "pdf")
}

library(cytomapper)
# Image visualisation with six markers: Ab, Iba1, Cd68, MAP2, GFAP, APOE
for(i in seq(1,176,by=4)){
  # normalisation means to scale the pixel intensities per channel between 0 and 1
  norm_images <- cytomapper::normalize(images[i:(i+3)])
  
  
  pdf(file = paste0(names(norm_images)[1],"_pixelPlot.pdf"),
      paper='a4')# The height of the plot in inches
  
  plotPixels(norm_images,
             colour_by=c('Ab','Iba1','CD68','MAP2','GFAP',
                         'ApoE'),
             bcg = list(Ab = c(0, 15, 1),
                        Iba1 = c(0, 10, 1),
                        CD68 = c(0, 20, 1),
                        MAP2 = c(0, 10, 1),
                        GFAP = c(0, 20, 1),
                        APOE = c(0, 10, 1)))
  
  dev.off()
  
}

channelNames(norm_images)

# Image visualisation with only DNA: DNA1; DNA2
for(i in seq(1,176,by=4)){
  
  norm_images <- cytomapper::normalize(images[i:(i+3)])
  
  
  pdf(file = paste0(names(norm_images)[1],"_DNA_pixelPlot.pdf"),
      paper='a4')# The height of the plot in inches
  
  plotPixels(norm_images,
             colour_by=c('DNA 1','DNA 2'),
             bcg = list(`DNA 1` = c(0, 20, 1),
                        `DNA 2` = c(0, 15, 1)))
  
  dev.off()
  
}

# Image visualisation with only Ab
for(i in seq(1,176,by=4)){
  
  norm_images <- cytomapper::normalize(images[i:(i+3)])
  
  
  pdf(file = paste0(names(norm_images)[1],"_Ab_pixelPlot.pdf"),
      paper='a4')# The height of the plot in inches
  
  plotPixels(norm_images,
             colour_by=c('Ab'),
             bcg = list(`Ab` = c(0, 20, 1)))
  
  dev.off()
  
}


# Segmentation mask visualisation
# plotCells(masks, object = spe,
#           img_id = "sample_id", cell_id = "ObjectNumber",
#           colour_by = c('Ab','Iba1','CD68','MAP2','GFAP',
#                         'ApoE'),
#           display='single')

# library(gridExtra)
# num_images_per_page <- 4
# # Loop through the images and create separate PDF files
# for (i in seq(1, 176, by = num_images_per_page)) {
#   # Select a subset of images for the current page
#   subset_masks <- masks[i:(i + 3)]
#   
#   # Create a PDF file for the current page
#   pdf(file = paste0("Page_", i, "-", i + 3, "_mask_visualisation.pdf"),
#       paper = 'a4') # Adjust the paper size if needed
#   
#   # Loop through the subset of images and plot them
#   for (j in 1:length(subset_masks)) {
#     plotPixels(subset_masks[j],
#                colour_by = c('Ab', 'Iba1', 'CD68', 'MAP2', 'GFAP', 'ApoE'))
#   }
#   
#   # Close the current PDF file
#   dev.off()
# }
# Cell visualisation
# https://bodenmillergroup.github.io/IMCDataAnalysis/image-visualization.html
# Ab
plotCells(masks,
          object = spe, 
          cell_id = "ObjectNumber", img_id = "sample_id",
          colour_by = "Ab",
          exprs_values = "exprs")
plotCells(masks_plaques,
          object = spe, 
          cell_id = "ObjectNumber", img_id = "sample_id",
          colour_by = "Ab",
          exprs_values = "exprs")

# Iba1, CD68, MAP2, GFAP
plotCells(masks,
          object = spe, 
          cell_id = "ObjectNumber", img_id = "sample_id",
          colour_by = c("Iba1", "CD68", "MAP2", "GFAP"))
plotCells(masks,
          object = spe, 
          cell_id = "ObjectNumber", img_id = "sample_id",
          colour_by = c("Iba1", "CD68", "MAP2", "GFAP"),
          exprs_values = "exprs")
plotCells(masks_plaques,
          object = spe, 
          cell_id = "ObjectNumber", img_id = "sample_id",
          colour_by = c("Iba1", "CD68", "MAP2", "GFAP"),
          exprs_values = "exprs")
plotCells(masks_plaques,
          object = spe, 
          cell_id = "ObjectNumber", img_id = "sample_id",
          colour_by = c("DNA 1", "DNA 2"),
          exprs_values = "exprs")
# Outlining cells on image
# plotPixels(image = images,
#            mask = masks_plaques,
#            object = spe, 
#            cell_id = "ObjectNumber", img_id = "sample_id",
#            colour_by = c("Ab", "Iba1", "CD68", "MAP2", "GFAP"),
#            outline_by = "mask",
#            bcg = list(Ab = c(0, 200, 1),
#                       Iba1 = c(0, 10, 1),
#                       CD68 = c(0, 15, 1),
#                       MAP2 = c(0, 10, 1),
#                       GFAP = c(0, 15, 1)),
#            thick = TRUE)

sce$mask=ifelse(sce$s.area<100,'small','big')
sce$mask=as.factor(sce$mask)
#sce2$mask=ifelse(sce2$s.area<100,‘small’,‘big’)
#sce2$mask=as.factor(sce2$mask)
set.seed(123)
for(i in seq(1,176,by=4)){
  norm_images <- cytomapper::normalize(images[i:(i+3)])
  masks_plaques1 <- masks_plaques[i:(i+3)]
  pdf(file = paste0(names(norm_images)[1],"_outline_plaques.pdf"),
      paper='a4')# The height of the plot in inches
  plotPixels(norm_images, mask = masks_plaques1,
           object = sce, img_id = "sample_id",
           cell_id = "object_id",
           colour_by = c("Ab","GFAP","Iba1","CD68"),
           bcg = list(Ab = c(0, 30, 5),
                      GFAP = c(0, 10, 1),
                      Iba1 = c(0, 20, 1),
                      CD68 = c(0, 15, 1)),
           outline_by="mask",
           scale=T,thick = TRUE) ## original output in: pixel_plots_outline_plaques/19891053_001_pixelPlot.pdf
  dev.off()
  
}
plotPixels(images, mask = masks,
           object = sce, img_id = "sample_id",
           cell_id = "object_id",
           colour_by = c("Ab","GFAP","Iba1","CD68","MAP2"),bcg = list(Ab = c(0, 30, 1),
                                                                      GFAP = c(0, 10, 1),
                                                                      Iba1 = c(0, 20, 1),
                                                                      CD68 = c(0, 15, 1),
                                                                      MAP2 = c(0, 10, 1)),
           outline_by="mask",
           scale=T,thick = TRUE)

########## On Plaque objects
setwd("~/mtg_dna/deepcell/matched_objects")
file_names=list.files( recursive = TRUE, pattern = "*.csv")

read_csv_column <- function(x){ data=data.frame(read_csv(x))
if(nrow(data)>0){
  data$name =sub(x = x, pattern = ".csv", replacement = "")
}
return(data)}
all_df <- do.call(rbind,lapply(file_names,read_csv_column))

all_df$DataFrame=paste0(all_df$name,'_',all_df$masks)

spe$OnPlaqueObject=ifelse(rownames(colData(spe)) %in% all_df$DataFrame,'OnPlaque','OutsidePlaque')
spe$PlaqueObj=all_df$masks_plaques[match(rownames(colData(spe)),all_df$DataFrame)]

#sanity check

all(is.na(spe$PlaqueObj) == (spe$OnPlaqueObject=='OutsidePlaque'))
spe$PlaqueObj[rownames(colData(spe))=="19891053_001_2"]

###### Perform sample correction
# # https://bodenmillergroup.github.io/IMCDataAnalysis/batch-effects.html
# library(batchelor)
# set.seed(220228)
# out <- fastMNN(spe, batch = spe$patient_id,
#                auto.merge = TRUE,
#                subset.row = rowData(spe)$use_channel,
#                assay.type = "exprs")
# 
# # Transfer the correction results to the main spe object
# reducedDim(spe, "fastMNN") <- reducedDim(out, "corrected")
# 
# 
# 
# ###### Session 4 single-cell analysis
# assay(sce, "exprs") <- asinh(counts(sce)/5)
# plotSpotHeatmap(sce)

# ##### check masks