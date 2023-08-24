### Creating spe_filtered
# Cell size
colData(spe) %>%
  as.data.frame() %>%
  group_by(sample_id) %>%
  ggplot() +
  geom_boxplot(aes(sample_id, area)) +
  theme_minimal(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  ylab("Cell area") + xlab("")

test=data.frame(colData(spe))
test$unique_id=rownames(test)

# Removing cells that are >Q3 and <Q1. 
Q1 <- quantile(test$area, .05)
Q3 <- quantile(test$area, .95)
IQR <- IQR(test$area)
no_outliers <- subset(test, test$area > (Q1 - 1.5*IQR) &  test$area < (Q3 + 1.5*IQR))

dim(test)
dim(no_outliers)

no_outliers %>%
  as.data.frame() %>%
  group_by(sample_id) %>%
  ggplot() +
  geom_boxplot(aes(sample_id, area)) +
  theme_minimal(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  ylab("Cell area") + xlab("")


`%nin%` = Negate(`%in%`)

test$unique_id[test$unique_id %nin% no_outliers$unique_id]
unique(test$sample_id[test$unique_id %nin% no_outliers$unique_id])
unique(test$patient_id[test$unique_id %nin% no_outliers$unique_id])

dim(test)
spe
spe_filtered=spe[,test$unique_id %in% no_outliers$unique_id]
dim(spe_filtered)

# saveRDS(spe_filtered,"spe_filtered.rds")

# spe_filtered <- readRDS("spe_filtered.rds")

####### Using spe_filtered
# Image area covered by cells
library(dplyr)

colData(spe_filtered) %>%
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

# Cell size
colData(spe_filtered) %>%
  as.data.frame() %>%
  group_by(sample_id) %>%
  ggplot() +
  geom_boxplot(aes(sample_id, area)) +
  theme_minimal(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  ylab("Cell area") + xlab("")

# Calculate how many cells there are per ROI
colData(spe_filtered) %>%
  as.data.frame() %>%
  group_by(sample_id) %>%
  summarize(num_cells = n())

#### Transform counts
## Ridge plots by AD/Control; variants (AD only ); by variants (Control only);
## and also grouping by Braak (create a new variable to bin Braak 0-II , Braak 3-4 and braak 5-6)

library(dittoSeq)
# # colnames(spe) <- paste0(spe$sample_id, "_", spe$ObjectNumber)
# # dittoRidgePlot(spe, var = "", group.by = "sample_id", assay = "counts")
# assay(spe, "exprs") <- asinh(counts(spe)/1)
# dittoRidgePlot(spe, var = "MAP2", group.by = "patient_id", assay = "exprs") +
#   ggtitle("MAP2 - after transformation")
# dittoRidgePlot(spe, var = "GFAP", group.by = "patient_id", assay = "exprs") +
#   ggtitle("GFAP - after transformation")
# dittoRidgePlot(spe, var = "CD68", group.by = "patient_id", assay = "exprs") +
#   ggtitle("CD68 - after transformation")

# color.panel = dittoColors()
## Ab
dittoRidgePlot(spe_filtered, var = "Ab", group.by = "diagnosis", assay = "exprs",
               colors = seq_along(color.panel)) +
  ggtitle("Ab - after transformation")
dittoRidgePlot(spe_filtered, var = "Ab", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis', color.panel = c("#009E73", "#F0E442", "#0072B2")) +
  ggtitle("Ab - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe_filtered, var = "Ab", group.by = "trem2", assay = "exprs", split.by='diagnosis', color.panel = c("salmon", "#2CFFC6", "purple")) +
  ggtitle("Ab - after transformation")
dittoRidgePlot(spe_filtered, var = "Ab", group.by = "trem2_all", assay = "exprs", split.by='diagnosis', color.panel = c("#1C91D4", "red")) +
  ggtitle("Ab - after transformation")


## CD68
dittoRidgePlot(spe_filtered, var = "CD68", group.by = "diagnosis", assay = "exprs") +
  ggtitle("CD68 - after transformation")
dittoRidgePlot(spe_filtered, var = "CD68", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis', color.panel = c("#009E73", "#F0E442", "#0072B2")) +
  ggtitle("CD68 - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe_filtered, var = "CD68", group.by = "trem2", assay = "exprs", split.by='diagnosis', color.panel = c("salmon", "#2CFFC6", "purple")) +
  ggtitle("CD68 - after transformation")
dittoRidgePlot(spe_filtered, var = "CD68", group.by = "trem2_all", assay = "exprs", split.by='diagnosis', color.panel = c("#1C91D4", "red")) +
  ggtitle("CD68 - after transformation")

## "ApoE"   
dittoRidgePlot(spe_filtered, var = "ApoE", group.by = "diagnosis", assay = "exprs") +
  ggtitle("ApoE - after transformation")
dittoRidgePlot(spe_filtered, var = "ApoE", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis', color.panel = c("#009E73", "#F0E442", "#0072B2")) +
  ggtitle("ApoE - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe_filtered, var = "ApoE", group.by = "trem2", assay = "exprs", split.by='diagnosis', color.panel = c("salmon", "#2CFFC6", "purple")) +
  ggtitle("ApoE - after transformation")
dittoRidgePlot(spe_filtered, var = "ApoE", group.by = "trem2_all", assay = "exprs", split.by='diagnosis', color.panel = c("#1C91D4", "red")) +
  ggtitle("ApoE - after transformation")

## "NfL"
dittoRidgePlot(spe_filtered, var = "NfL", group.by = "diagnosis", assay = "exprs") +
  ggtitle("NfL - after transformation")
dittoRidgePlot(spe_filtered, var = "NfL", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis', color.panel = c("#009E73", "#F0E442", "#0072B2")) +
  ggtitle("NfL - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe_filtered, var = "NfL", group.by = "trem2", assay = "exprs", split.by='diagnosis', color.panel = c("salmon", "#2CFFC6", "purple")) +
  ggtitle("NfL - after transformation")
dittoRidgePlot(spe_filtered, var = "NfL", group.by = "trem2_all", assay = "exprs", split.by='diagnosis', color.panel = c("#1C91D4", "red")) +
  ggtitle("NfL - after transformation")

## "GFAP"
dittoRidgePlot(spe_filtered, var = "GFAP", group.by = "diagnosis", assay = "exprs") +
  ggtitle("GFAP - after transformation")
dittoRidgePlot(spe_filtered, var = "GFAP", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis', color.panel = c("#009E73", "#F0E442", "#0072B2")) +
  ggtitle("GFAP - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe_filtered, var = "GFAP", group.by = "trem2", assay = "exprs", split.by='diagnosis', color.panel = c("salmon", "#2CFFC6", "purple")) +
  ggtitle("GFAP - after transformation")
dittoRidgePlot(spe_filtered, var = "GFAP", group.by = "trem2_all", assay = "exprs", split.by='diagnosis', color.panel = c("#1C91D4", "red")) +
  ggtitle("GFAP - after transformation")

## "MAP2"
dittoRidgePlot(spe_filtered, var = "MAP2", group.by = "diagnosis", assay = "exprs") +
  ggtitle("MAP2 - after transformation")
dittoRidgePlot(spe_filtered, var = "MAP2", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis', color.panel = c("#009E73", "#F0E442", "#0072B2")) +
  ggtitle("MAP2 - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe_filtered, var = "MAP2", group.by = "trem2", assay = "exprs", split.by='diagnosis', color.panel = c("salmon", "#2CFFC6", "purple")) +
  ggtitle("MAP2 - after transformation")
dittoRidgePlot(spe_filtered, var = "MAP2", group.by = "trem2_all", assay = "exprs", split.by='diagnosis', color.panel = c("#1C91D4", "red")) +
  ggtitle("MAP2 - after transformation")

## "HLA-DR"
dittoRidgePlot(spe_filtered, var = "HLA-DR", group.by = "diagnosis", assay = "exprs") +
  ggtitle("HLA-DR - after transformation")
dittoRidgePlot(spe_filtered, var = "HLA-DR", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis', color.panel = c("#009E73", "#F0E442", "#0072B2")) +
  ggtitle("HLA-DR - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe_filtered, var = "HLA-DR", group.by = "trem2", assay = "exprs", split.by='diagnosis', color.panel = c("salmon", "#2CFFC6", "purple")) +
  ggtitle("HLA-DR - after transformation")
dittoRidgePlot(spe_filtered, var = "HLA-DR", group.by = "trem2_all", assay = "exprs", split.by='diagnosis', color.panel = c("#1C91D4", "red")) +
  ggtitle("HLA-DR - after transformation")

## "Iba1"
dittoRidgePlot(spe_filtered, var = "Iba1", group.by = "diagnosis", assay = "exprs") +
  ggtitle("Iba1 - after transformation")
dittoRidgePlot(spe_filtered, var = "Iba1", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis', color.panel = c("#009E73", "#F0E442", "#0072B2")) +
  ggtitle("Iba1 - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe_filtered, var = "Iba1", group.by = "trem2", assay = "exprs", split.by='diagnosis', color.panel = c("salmon", "#2CFFC6", "purple")) +
  ggtitle("Iba1 - after transformation")
dittoRidgePlot(spe_filtered, var = "Iba1", group.by = "trem2_all", assay = "exprs", split.by='diagnosis', color.panel = c("#1C91D4", "red")) +
  ggtitle("Iba1 - after transformation")

## "CD163"
dittoRidgePlot(spe_filtered, var = "CD163", group.by = "diagnosis", assay = "exprs") +
  ggtitle("CD163 - after transformation")
dittoRidgePlot(spe_filtered, var = "CD163", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis', color.panel = c("#009E73", "#F0E442", "#0072B2")) +
  ggtitle("CD163 - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe_filtered, var = "CD163", group.by = "trem2", assay = "exprs", split.by='diagnosis', color.panel = c("salmon", "#2CFFC6", "purple")) +
  ggtitle("CD163 - after transformation")
dittoRidgePlot(spe_filtered, var = "CD163", group.by = "trem2_all", assay = "exprs", split.by='diagnosis', color.panel = c("#1C91D4", "red")) +
  ggtitle("CD163 - after transformation")

## "TSPO"
dittoRidgePlot(spe_filtered, var = "TSPO", group.by = "diagnosis", assay = "exprs") +
  ggtitle("TSPO - after transformation")
dittoRidgePlot(spe_filtered, var = "TSPO", group.by = "BraakGroup", assay = "exprs", split.by='diagnosis', color.panel = c("#009E73", "#F0E442", "#0072B2")) +
  ggtitle("TSPO - after transformation (grouped by Braak, split by diagnosis)")
dittoRidgePlot(spe_filtered, var = "TSPO", group.by = "trem2", assay = "exprs", split.by='diagnosis', color.panel = c("salmon", "#2CFFC6", "purple")) +
  ggtitle("TSPO - after transformation")
dittoRidgePlot(spe_filtered, var = "TSPO", group.by = "trem2_all", assay = "exprs", split.by='diagnosis', color.panel = c("#1C91D4", "red")) +
  ggtitle("TSPO - after transformation")



## Read in images
library(cytomapper)
# images <- loadImages("~/mtg_dna/deepcell/img/")
# saveRDS(images, "saved_images.rds")
# images <- readRDS("saved_images.rds")
# masks <- loadImages("~/mtg_dna/deepcell/masks/", as.is = TRUE)
# saveRDS(masks, "masks.rds")
# masks <- readRDS('masks.rds')
# masks_plaques <- loadImages("masks_plaques/", pattern = ".tiff", as.is = TRUE)
# new_masks <- loadImages("new_masks/", pattern = ".tiff", as.is = TRUE)
# saveRDS(new_masks, "new_masks.rds")
# saveRDS(masks_plaques, "masks_plaques.rds")
# new_plaques_global <- loadImages("new_plaques_global/", pattern = ".tiff", as.is = TRUE)
# saveRDS(new_plaques_global, "new_plaques_global.rds")
channelNames(images) <- rownames(spe_filtered)
# images
# all.equal(names(images), names(masks))

patient_id <- str_extract(names(images), "^[^_]+")
diagnosis <- metadata_filtered$diagnosis[match(patient_id, metadata_filtered$CaseID)]
trem2 <- metadata_filtered$trem2[match(patient_id, metadata_filtered$CaseID)]
# Save the matched sample_id, patient_id, diagnosis, and trem2 information within the elementMetadata slot of
# of the multi-channel images and segmentation masks objects.
mcols(images) <- mcols(masks) <- mcols(new_plaques_global)<- DataFrame(sample_id = names(images),
                                                                                        patient_id = patient_id,
                                                                                        diagnosis = diagnosis,
                                                                                        trem2 = trem2)

## Generate single-cell data from images
# sce <- measureObjects(masks, image = images, img_id = "sample_id")
# saveRDS(sce, "sce.rds")
sce <- readRDS("sce.rds")
### Single-cell visualisation
#### Dimensionality reduction
#library(scater)

#set.seed(220225)
#spe <- runUMAP(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs") 
#spe <- runTSNE(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs") 

#reducedDims(spe)
#spe@int_colData$reducedDims$UMAP
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
  plot <- plotSpatial(spe_filtered[, grepl(patient_id, spe_filtered$sample_id)], 
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
  plot <- plotSpatial(spe_filtered[, grepl(patient_id, spe_filtered$sample_id)], 
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
  plot <- plotSpatial(spe_filtered[, grepl(patient_id, spe_filtered$sample_id)], 
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
  plot <- plotSpatial(spe_filtered[, grepl(patient_id, spe_filtered$sample_id)], 
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
  plot <- plotSpatial(spe_filtered[, grepl(patient_id, spe_filtered$sample_id)], 
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
  plot <- plotSpatial(spe_filtered[, grepl(patient_id, spe_filtered$sample_id)], 
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
  plot <- plotSpatial(spe_filtered[, grepl(patient_id, spe_filtered$sample_id)], 
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
  plot <- plotSpatial(spe_filtered[, grepl(patient_id, spe_filtered$sample_id)], 
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
  plot <- plotSpatial(spe_filtered[, grepl(patient_id, spe_filtered$sample_id)], 
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
  plot <- plotSpatial(spe_filtered[, grepl(patient_id, spe_filtered$sample_id)], 
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
#### Outlining plaques on image
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

set.seed(123)
for(i in seq(1,176,by=4)){
  norm_images <- cytomapper::normalize(images[i:(i+3)])
  new_plaques_global1 <- new_plaques_global[i:(i+3)]
  pdf(file = paste0(names(norm_images)[1],"_outline_plaques.pdf"),
      paper='a4')# The height of the plot in inches
  plotPixels(norm_images, mask = new_plaques_global1,
             object = sce, img_id = "sample_id",
             cell_id = "object_id",
             colour_by = c("Ab","GFAP","Iba1","CD68"),
             bcg = list(Ab = c(0, 45, 5),
                        GFAP = c(0, 10, 1),
                        Iba1 = c(0, 20, 1),
                        CD68 = c(0, 15, 1)),
             outline_by="mask",
             scale=T,thick = TRUE) 
  dev.off()
  
}

## DISPLAYING REPRESENTATIVE IMAGES (6)
# CV -> AD # 19920035_001 11 (4); 19920077_002 11 (5); 19930154_002 11 (9);
set.seed(1234)
for(i in seq(13,16,by=4)){
  norm_images <- cytomapper::normalize(images[i:(i+3)])
  new_plaques_global1 <- new_plaques_global[i:(i+3)]
  pdf(file = paste0(names(norm_images)[1],"_outline_plaques.pdf"),
      paper='a4')# The height of the plot in inches
  plotPixels(norm_images, mask = new_plaques_global1,
             object = sce, img_id = "sample_id",
             cell_id = "object_id",
             colour_by = c("Ab","GFAP","Iba1","CD68"),
             bcg = list(Ab = c(0, 30, 5),
                        GFAP = c(0, 10, 1),
                        Iba1 = c(0, 20, 1),
                        CD68 = c(0, 15, 1)),
             outline_by="mask",
             display = "single",
             scale=T,thick = TRUE) 
  dev.off()
}
for(i in seq(17,20,by=4)){
  norm_images <- cytomapper::normalize(images[i:(i+3)])
  new_plaques_global1 <- new_plaques_global[i:(i+3)]
  pdf(file = paste0(names(norm_images)[1],"_outline_plaques.pdf"),
      paper='a4')# The height of the plot in inches
  plotPixels(norm_images, mask = new_plaques_global1,
             object = sce, img_id = "sample_id",
             cell_id = "object_id",
             colour_by = c("Ab","GFAP","Iba1","CD68"),
             bcg = list(Ab = c(0, 30, 5),
                        GFAP = c(0, 10, 1),
                        Iba1 = c(0, 20, 1),
                        CD68 = c(0, 15, 1)),
             outline_by="mask",
             display = "single",
             scale=T,thick = TRUE) 
  dev.off()
}
for(i in seq(33,36,by=4)){
  norm_images <- cytomapper::normalize(images[i:(i+3)])
  new_plaques_global1 <- new_plaques_global[i:(i+3)]
  pdf(file = paste0(names(norm_images)[1],"_outline_plaques.pdf"),
      paper='a4')# The height of the plot in inches
  plotPixels(norm_images, mask = new_plaques_global1,
             object = sce, img_id = "sample_id",
             cell_id = "object_id",
             colour_by = c("Ab","GFAP","Iba1","CD68"),
             bcg = list(Ab = c(0, 900, 5),
                        GFAP = c(0, 10, 1),
                        Iba1 = c(0, 20, 1),
                        CD68 = c(0, 15, 1)),
             outline_by="mask",
             display = "single",
             scale=T,thick = TRUE) 
  dev.off()
}

# CV -> Control # A297.16_004 2; PDC26_001 2
for(i in seq(117,120,by=4)){
  norm_images <- cytomapper::normalize(images[i:(i+3)])
  new_plaques_global1 <- new_plaques_global[i:(i+3)]
  pdf(file = paste0(names(norm_images)[1],"_outline_plaques.pdf"),
      paper='a4')# The height of the plot in inches
  plotPixels(norm_images, mask = new_plaques_global1,
             object = sce, img_id = "sample_id",
             cell_id = "object_id",
             colour_by = c("Ab","GFAP","Iba1","CD68"),
             bcg = list(Ab = c(0, 30, 5),
                        GFAP = c(0, 10, 1),
                        Iba1 = c(0, 20, 1),
                        CD68 = c(0, 15, 1)),
             outline_by="mask",
             display = "single",
             scale=T,thick = TRUE) 
  dev.off()
}

for(i in seq(173,176,by=4)){
  norm_images <- cytomapper::normalize(images[i:(i+3)])
  new_plaques_global1 <- new_plaques_global[i:(i+3)]
  pdf(file = paste0(names(norm_images)[1],"_outline_plaques.pdf"),
      paper='a4')# The height of the plot in inches
  plotPixels(norm_images, mask = new_plaques_global1,
             object = sce, img_id = "sample_id",
             cell_id = "object_id",
             colour_by = c("Ab","GFAP","Iba1","CD68"),
             bcg = list(Ab = c(0, 30, 5),
                        GFAP = c(0, 10, 1),
                        Iba1 = c(0, 20, 1),
                        CD68 = c(0, 15, 1)),
             outline_by="mask",
             display = "single",
             scale=T,thick = TRUE) 
  dev.off()
}

# R47H -> AD # 19920194_001 12; 19920194_003 12 (6); A166.04_004 11 (22); A220.11_004 11 (27); 
for(i in seq(21,24,by=4)){
  norm_images <- cytomapper::normalize(images[i:(i+3)])
  new_plaques_global1 <- new_plaques_global[i:(i+3)]
  pdf(file = paste0(names(norm_images)[1],"_outline_plaques.pdf"),
      paper='a4')# The height of the plot in inches
  plotPixels(norm_images, mask = new_plaques_global1,
             object = sce, img_id = "sample_id",
             cell_id = "object_id",
             colour_by = c("Ab","GFAP","Iba1","CD68"),
             bcg = list(Ab = c(0, 30, 5),
                        GFAP = c(0, 10, 1),
                        Iba1 = c(0, 20, 1),
                        CD68 = c(0, 15, 1)),
             outline_by="mask",
             display = "single",
             scale=T,thick = TRUE) 
  dev.off()
}
for(i in seq(85,88,by=4)){
  norm_images <- cytomapper::normalize(images[i:(i+3)])
  new_plaques_global1 <- new_plaques_global[i:(i+3)]
  pdf(file = paste0(names(norm_images)[1],"_outline_plaques.pdf"),
      paper='a4')# The height of the plot in inches
  plotPixels(norm_images, mask = new_plaques_global1,
             object = sce, img_id = "sample_id",
             cell_id = "object_id",
             colour_by = c("Ab","GFAP","Iba1","CD68"),
             bcg = list(Ab = c(0, 30, 5),
                        GFAP = c(0, 10, 1),
                        Iba1 = c(0, 20, 1),
                        CD68 = c(0, 15, 1)),
             outline_by="mask",
             display = "single",
             scale=T,thick = TRUE) 
  dev.off()
}
for(i in seq(105,108,by=4)){
  norm_images <- cytomapper::normalize(images[i:(i+3)])
  new_plaques_global1 <- new_plaques_global[i:(i+3)]
  pdf(file = paste0(names(norm_images)[1],"_outline_plaques.pdf"),
      paper='a4')# The height of the plot in inches
  plotPixels(norm_images, mask = new_plaques_global1,
             object = sce, img_id = "sample_id",
             cell_id = "object_id",
             colour_by = c("Ab","GFAP","Iba1","CD68"),
             bcg = list(Ab = c(0, 30, 5),
                        GFAP = c(0, 10, 1),
                        Iba1 = c(0, 20, 1),
                        CD68 = c(0, 15, 1)),
             outline_by="mask",
             display = "single",
             scale=T,thick = TRUE) 
  dev.off()
}

# R47H -> Control # C19.93_001 1
for(i in seq(145,148,by=4)){
  norm_images <- cytomapper::normalize(images[i:(i+3)])
  new_plaques_global1 <- new_plaques_global[i:(i+3)]
  pdf(file = paste0(names(norm_images)[1],"_outline_plaques.pdf"),
      paper='a4')# The height of the plot in inches
  plotPixels(norm_images, mask = new_plaques_global1,
             object = sce, img_id = "sample_id",
             cell_id = "object_id",
             colour_by = c("Ab","GFAP","Iba1","CD68"),
             bcg = list(Ab = c(0, 30, 5),
                        GFAP = c(0, 10, 1),
                        Iba1 = c(0, 20, 1),
                        CD68 = c(0, 15, 1)),
             outline_by="mask",
             display = "single",
             scale=T,thick = TRUE) 
  dev.off()
}

# R62H -> AD # A138.12_001 13 (20)
for(i in seq(77,80,by=4)){
  norm_images <- cytomapper::normalize(images[i:(i+3)])
  new_plaques_global1 <- new_plaques_global[i:(i+3)]
  pdf(file = paste0(names(norm_images)[1],"_outline_plaques.pdf"),
      paper='a4')# The height of the plot in inches
  plotPixels(norm_images, mask = new_plaques_global1,
             object = sce, img_id = "sample_id",
             cell_id = "object_id",
             colour_by = c("Ab","GFAP","Iba1","CD68"),
             bcg = list(Ab = c(0, 40, 5),
                        GFAP = c(0, 10, 1),
                        Iba1 = c(0, 20, 1),
                        CD68 = c(0, 15, 1)),
             outline_by="mask",
             display = "single",
             scale=T,thick = TRUE) 
  dev.off()
}

# R62H -> Control # A319.11_004 7 (31); NP13.012_001 7 (40); DPM16.44_002 5 (38)
for(i in seq(121,124,by=4)){
  norm_images <- cytomapper::normalize(images[i:(i+3)])
  new_plaques_global1 <- new_plaques_global[i:(i+3)]
  pdf(file = paste0(names(norm_images)[1],"_outline_plaques.pdf"),
      paper='a4')# The height of the plot in inches
  plotPixels(norm_images, mask = new_plaques_global1,
             object = sce, img_id = "sample_id",
             cell_id = "object_id",
             colour_by = c("Ab","GFAP","Iba1","CD68"),
             bcg = list(Ab = c(0, 200, 5),
                        GFAP = c(0, 10, 1),
                        Iba1 = c(0, 20, 1),
                        CD68 = c(0, 15, 1)),
             outline_by="mask",
             display = "single",
             scale=T,thick = TRUE) 
  dev.off()
}
for(i in seq(157,160,by=4)){
  norm_images <- cytomapper::normalize(images[i:(i+3)])
  new_plaques_global1 <- new_plaques_global[i:(i+3)]
  pdf(file = paste0(names(norm_images)[1],"_outline_plaques.pdf"),
      paper='a4')# The height of the plot in inches
  plotPixels(norm_images, mask = new_plaques_global1,
             object = sce, img_id = "sample_id",
             cell_id = "object_id",
             colour_by = c("Ab","GFAP","Iba1","CD68"),
             bcg = list(Ab = c(0, 30, 5),
                        GFAP = c(0, 10, 1),
                        Iba1 = c(0, 20, 1),
                        CD68 = c(0, 15, 1)),
             outline_by="mask",
             display = "single",
             scale=T,thick = TRUE) 
  dev.off()
}

# ## Manually adjust the contrast of Ab of specific images
# set.seed(123)
# for(i in seq(161,164,by=4)){
#   norm_images <- cytomapper::normalize(images[i:(i+3)])
#   new_plaques_global1 <- new_plaques_global[i:(i+3)]
#   pdf(file = paste0(names(norm_images)[1],"_outline_plaques.pdf"),
#       paper='a4')# The height of the plot in inches
#   plotPixels(norm_images, mask = new_plaques_global1,
#              object = sce, img_id = "sample_id",
#              cell_id = "object_id",
#              colour_by = c("Ab","GFAP","Iba1","CD68"),
#              bcg = list(Ab = c(0, 300, 5),
#                         GFAP = c(0, 10, 1),
#                         Iba1 = c(0, 20, 1),
#                         CD68 = c(0, 15, 1)),
#              outline_by="mask",
#              scale=T,thick = TRUE) 
#   dev.off()
#   
# }


# count mean cell count:
# by sample_id
mean_num_cells <- colData(spe_filtered) %>%
  as.data.frame() %>%
  group_by(patient_id, cluster_celltype) %>%
  summarise(num_cells = n()) %>%
  group_by(cluster_celltype) %>%
  summarise(mean_num_cells = mean(num_cells))

# Display the mean number of cells
print(mean_num_cells)
# mean_num_cells
# <dbl>
#   1           404.

# by diagnosis: AD
mean_num_cells <- colData(spe_filtered) %>%
  as.data.frame() %>%
  group_by(diagnosis, cluster_celltype) %>%
  filter(diagnosis %in% c("AD")) %>%
  summarise(num_cells=n()) %>%
  group_by(cluster_celltype) %>%
  summarise(mean_num_cells = mean(num_cells))

# Display the mean number of cells for the "AD" group
print(mean_num_cells)
# <dbl>
#   1          49261

# by diagnosis: Control
mean_num_cells <- colData(spe_filtered) %>%
  as.data.frame() %>%
  group_by(diagnosis, cluster_celltype) %>%
  filter(diagnosis %in% c("Control")) %>%
  summarise(num_cells=n()) %>%
  group_by(cluster_celltype) %>%
  summarise(mean_num_cells = mean(num_cells))

# Display the mean number of cells for the "Control" group
print(mean_num_cells)
# <dbl>
#   1          21904

# by AD_TREM2
mean_num_cells <- colData(spe_filtered) %>%
  as.data.frame() %>%
  group_by(groups, cluster_celltype) %>%
  filter(groups %in% c("AD_TREM2")) %>%
  #group_by(groups) %>%
  summarise(num_cells=n()) %>%
  group_by(cluster_celltype) %>%
  summarise(mean_num_cells = mean(num_cells))

# Display the mean number of cells for the "AD_TREM2" group
print(mean_num_cells)
# <dbl>
#   1          19173

# by AD_CV
mean_num_cells <- colData(spe_filtered) %>%
  as.data.frame() %>%
  group_by(groups, cluster_celltype) %>%
  filter(groups %in% c("AD_CV")) %>%
  summarise(num_cells=n()) %>%
  group_by(cluster_celltype) %>%
  summarise(mean_num_cells = mean(num_cells))

# Display the mean number of cells for the "AD_CV" group
print(mean_num_cells)
# <dbl>
#   1          30088

# by CTRL_TREM2
mean_num_cells <- colData(spe_filtered) %>%
  as.data.frame() %>%
  group_by(groups, cluster_celltype) %>%
  filter(groups %in% c("CTRL_TREM2")) %>%
  summarise(num_cells=n()) %>%
  group_by(cluster_celltype) %>%
  summarise(mean_num_cells = mean(num_cells))

# Display the mean number of cells for the "CTRL_TREM2" group
print(mean_num_cells)
# <dbl>
#   1          13152

# by CTRL_CV
mean_num_cells <- colData(spe_filtered) %>%
  as.data.frame() %>%
  group_by(groups, cluster_celltype) %>%
  filter(groups %in% c("CTRL_CV")) %>%
  #group_by(groups) %>%
  summarise(num_cells=n()) %>%
  group_by(cluster_celltype) %>%
  summarise(mean_num_cells = mean(num_cells))

# Display the mean number of cells for the "CTRL_CV" group
print(mean_num_cells)
# <dbl>
#   1           8752

# spe_filtered_AD only
# Braak 3-4 
mean_num_cells <- colData(spe_filtered_AD) %>%
  as.data.frame() %>%
  group_by(BraakGroup, cluster_celltype) %>%
  filter(BraakGroup %in% c("Braak 3-4")) %>%
  #group_by(BraakGroup) %>%
  summarise(num_cells=n()) %>%
  group_by(cluster_celltype) %>%
  summarise(mean_num_cells = mean(num_cells))

# Display the mean number of cells for the group
print(mean_num_cells)
#            <dbl>
# 1           4803


# spe_filtered_AD only
# TREM2 - CV
mean_num_cells <- colData(spe_filtered_AD) %>%
  as.data.frame() %>%
  group_by(trem2, cluster_celltype) %>%
  filter(trem2 %in% c("CV")) %>%
  #group_by(BraakGroup) %>%
  summarise(num_cells=n()) %>%
  group_by(cluster_celltype) %>%
  summarise(mean_num_cells = mean(num_cells))

# Display the mean number of cells for the "CV" group
print(mean_num_cells)

# TREM2 - R47H
mean_num_cells <- colData(spe_filtered_AD) %>%
  as.data.frame() %>%
  group_by(trem2, cluster_celltype) %>%
  filter(trem2 %in% c("R47H")) %>%
  summarise(num_cells=n()) %>%
  group_by(cluster_celltype) %>%
  summarise(mean_num_cells = mean(num_cells))

# Display the mean number of cells for the "R47H" group
print(mean_num_cells)

# TREM2 - R62H
mean_num_cells <- colData(spe_filtered_AD) %>%
  as.data.frame() %>%
  group_by(trem2, cluster_celltype) %>%
  filter(trem2 %in% c("R62H")) %>%
  summarise(num_cells=n()) %>%
  group_by(cluster_celltype) %>%
  summarise(mean_num_cells = mean(num_cells))

# Display the mean number of cells for the "R62H" group
print(mean_num_cells)
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

spe_filtered$OnPlaqueObject=ifelse(rownames(colData(spe_filtered)) %in% all_df$DataFrame,'OnPlaque','OutsidePlaque')
spe_filtered$PlaqueObj=all_df$masks_plaques[match(rownames(colData(spe_filtered)),all_df$DataFrame)]

#sanity check
all(is.na(spe_filtered$PlaqueObj) == (spe_filtered$OnPlaqueObject=='OutsidePlaque'))
spe_filtered$PlaqueObj[rownames(colData(spe_filtered))=="19891053_001_2"]


