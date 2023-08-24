# sce_new<- read_steinbock("~/mtg_dna/deepcell",
#                           return_as="sce_new") 
# 
# rm(sce_new)
# sce_new<- readSCEfromTXT("~/mtg_dna/deepcell/test_txt") 
# 
# assay(sce_new, "exprs") <- asinh(counts(sce_new)/5)
# rowData(sce_new)$channel_id=rowData(sce_new)$channel
# plotSpotHeatmap(sce_new,
#                 spot_id = 'sample_id',
#                 channel_id='channel_id',
#                 assay_type = 'exprs')
# 
# help("plotSpotHeatmap")
# 
# colData(sce_new)
# plotSpotHeatmap(sce_new,
#                 spot_id = 'sample_id',
#                 channel_id='channel_id',
#                 log = FALSE, threshold = 200)
# 
# bc_key <- as.numeric(unique(sce_new$sample_mass))
# bc_key <- bc_key[order(bc_key)]
# 
# sce_new <- assignPrelim(sce_new, bc_key = bc_key)
# sce_new <- estCutoffs(sce_new)
# sce_new <- applyCutoffs(sce_new)
# 
# library(pheatmap)
# cur_table <- table(sce_new$bc_id, sce_new$sample_mass)
# 
# pheatmap(log10(cur_table + 1),
#          cluster_rows = FALSE,
#          cluster_cols = FALSE)
# 

# ### Define interesting channels
# #### Store an additional entry in the rowData slot that specifies the markers of interest. Deselect nuclear markers which were used for cell segmentation, and keep all other biological targets.
# rowData(spe)$use_channel <- !grepl("DNA 1|DNA 2", rownames(spe))
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
# merge_info <- metadata(out)$merge.info
# 
# DataFrame(left = merge_info$left,
#           right = merge_info$right,
#           batch.size = merge_info$batch.size,
#           max_lost_var = rowMax(merge_info$lost.var))
# # Visualisation
# library(scater)
# 
# set.seed(220228)
# spe <- runUMAP(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs") 
# spe <- runUMAP(spe, dimred= "fastMNN", name = "UMAP_mnnCorrected") 
# library(cowplot)
# library(dittoSeq)
# library(viridis)
# 
# # visualize patient id 
# p1 <- dittoDimPlot(spe, var = "patient_id", 
#                    reduction.use = "UMAP", size = 0.2) +
#   ggtitle("Patient ID on UMAP before correction")
# p2 <- dittoDimPlot(spe, var = "patient_id", 
#                    reduction.use = "UMAP_mnnCorrected", size = 0.2) +
#   ggtitle("Patient ID on UMAP after correction")
# 
# plot_grid(p1, p2)

### Define interesting channels
#### Store an additional entry in the rowData slot that specifies the markers of interest. Deselect nuclear markers which were used for cell segmentation, and keep all other biological targets.
rowData(spe_filtered)$use_channel <- !grepl("DNA 1|DNA 2", rownames(spe_filtered))
col=rownames(spe_filtered)[rowData(spe_filtered)$use_channel]
library(batchelor)
sum(rowData(spe_filtered)$use_channel)
set.seed(220228)
out <- fastMNN(spe_filtered, batch = spe_filtered$patient_id,
               auto.merge = TRUE,
               subset.row = rowData(spe_filtered)$use_channel,
               assay.type = "exprs")
out <- readRDS("out.rds")
# Transfer the correction results to the main spe object
test <- reducedDim(out, "corrected")
test2 <- out@assays@data@listData[["reconstructed"]]
reducedDim(spe_filtered, "fastMNN") <- reducedDim(out, "corrected")


merge_info <- metadata(out)$merge.info

DataFrame(left = merge_info$left,
          right = merge_info$right,
          batch.size = merge_info$batch.size,
          max_lost_var = rowMax(merge_info$lost.var))
# saveRDS(out, "out.rds")

# Visualisation
library(scater)

set.seed(220228)
spe_filtered <- runUMAP(spe_filtered, subset_row = rowData(spe_filtered)$use_channel, exprs_values = "exprs") 
spe_filtered <- runUMAP(spe_filtered, dimred= "fastMNN", name = "UMAP_mnnCorrected") 
rowData(spe_filtered)$use_channel <- !grepl("DNA 1|DNA 2|80ArAr|131Xe|127I|134Xe|138Ba|208Pb", rownames(spe_filtered))

spe_filtered <- runUMAP(spe_filtered, subset_row = rowData(spe_filtered)$use_channel, dimred= "fastMNN", name = "UMAP_mnnCorrected2")

saveRDS(spe_filtered, "spe_filtered.rds") #saved on 20.6
library(cowplot)
library(dittoSeq)
library(viridis)

# visualize patient id 
p1 <- dittoDimPlot(spe_filtered, var = "patient_id", 
                   reduction.use = "UMAP", size = 0.2) +
  ggtitle("Patient ID on UMAP before correction")
p2 <- dittoDimPlot(spe_filtered, var = "patient_id", 
                   reduction.use = "UMAP_mnnCorrected", size = 0.2) +
  ggtitle("Patient ID on UMAP after correction")

plot_grid(p1, p2)

markers <- c("Ab", "ApoE", "CD68", "CD163", "GFAP", "HLA-DR", "Iba1", "MAP2", "NfL", "TSPO")

# Before correction
plot_list <- multi_dittoDimPlot(spe_filtered, var = markers, reduction.use = "UMAP", 
                                assay = "exprs", size = 0.2, list.out = TRUE) 

# Modify each plot in the list
plot_list <- lapply(plot_list, function(x) {
  x + scale_color_viridis()
})

# Combine the plots into a grid with a single title
plot_grid(plotlist = plot_list, ncol = 2) +
  plot_annotation(title = "Markers Across All Cells - Before Correction")

# After correction
plot_list <- multi_dittoDimPlot(spe, var = markers, reduction.use = "UMAP_mnnCorrected", 
                                assay = "exprs", size = 0.2, list.out = TRUE) 
plot_list <- lapply(plot_list, function(x) {
  x + scale_color_viridis()
})

plot_grid(plotlist = plot_list, ncol=2) + 
  plot_annotation(title = "Markers Across All Cells - After Correction")
