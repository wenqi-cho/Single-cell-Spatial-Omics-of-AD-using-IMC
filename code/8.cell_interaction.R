library(scales)
library(SpatialExperiment)
library(imcRtools)
spe_filtered_AD <- buildSpatialGraph(spe_filtered_AD, img_id = "sample_id", type = "knn", k = 20)
# spe_filtered_AD <- buildSpatialGraph(spe_filtered_AD, img_id = "sample_id", type = "expansion", threshold = 30)
spe_filtered_AD$cluster_names=spe_filtered_AD$nn_clusters_corrected14

out2 <- testInteractions(spe_filtered_AD, 
                         group_by = "sample_id",
                         label = "cluster_names", 
                         colPairName = "knn_interaction_graph",
                         BPPARAM = SerialParam(RNGseed = 221029))


out2$trem2=spe_filtered_AD$trem2[match(out2$group_by,spe_filtered_AD$sample_id)]

out2%>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


out2%>% as_tibble() %>%
  filter(trem2=='R47H') %>% 
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

out2%>% as_tibble() %>%
  filter(trem2=='R62H') %>% 
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


out2%>% as_tibble() %>%
  filter(trem2=='CV') %>% 
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


R47H=out2%>% as_tibble() %>%
  filter(trem2=='R47H') %>% 
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE),
            trem2='R47H')


R62H=out2%>% as_tibble() %>%
  filter(trem2=='R62H') %>% 
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE),
            trem2='R62H')


CV=out2%>% as_tibble() %>%
  filter(trem2=='CV') %>% 
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE),
            trem2='CV')

all=rbind(R47H,R62H,CV)
library(circlize)
set.seed(123)

duplicate_pairs<-function(df){
  # Create unique identifier for each combination
  df$identifier <- apply(df[, c("from_label", "to_label")], 1, function(x) {
    sorted_pair <- sort(x)
    paste(sorted_pair, collapse = "-")
  })
  
  # Remove duplicated pairs
  df_unique <- df[!duplicated(df$identifier), ]
  return(df_unique)
}


# chordDiagram(CV,
#              annotationTrack =  c("name", "grid"),
#              grid.col=c("#E6A005", "#56B4E9", "#009E73","#F0E442","#0072B2", "#D55E00" ,"#CC79A7", "#666666","#AD7700","#1C91D4","#007756","#D5C711","#005685", "#A04700", "#B14380"),
#              transparency = 0.5)
# 
# 
# chordDiagram(R62H,
#              annotationTrack =  c("name", "grid"),
#              grid.col=c("#E6A005", "#56B4E9", "#009E73","#F0E442","#0072B2", "#D55E00" ,"#CC79A7", "#666666","#AD7700","#1C91D4","#007756","#D5C711","#005685", "#A04700", "#B14380"),
#              transparency = 0.5)
# 
# chordDiagram(R47H,
#              annotationTrack =  c("name", "grid"),
#              grid.col=c("#E6A005", "#56B4E9", "#009E73","#F0E442","#0072B2", "#D55E00" ,"#CC79A7", "#666666","#AD7700","#1C91D4","#007756","#D5C711","#005685", "#A04700", "#B14380"),
#              transparency = 0.5)


### removing self links
all$self=ifelse(all$from_label==all$to_label,'self','normal')

all_normal=all%>%
  filter(self!='self')

CV$self=ifelse(CV$from_label==CV$to_label,'self','normal')
CV$double=paste(CV$from_label,CV$to_label)

CV_normal=CV%>%
  filter(self!='self')

R62H$self=ifelse(R62H$from_label==R62H$to_label,'self','normal')

R62H_normal=R62H%>%
  filter(self!='self')

R47H$self=ifelse(R47H$from_label==R47H$to_label,'self','normal')

R47H_normal=R47H%>%
  filter(self!='self')

library(RColorBrewer)

# print(unique(CV_normal$from_label))
# par(cex = 0.2, mar = c(0, 0, 0, 0))

## CV
new_CV=duplicate_pairs(CV_normal)
cols = colorRamp2(c(min(new_CV$sum_sigval),0,max(new_CV$sum_sigval)),c("blue","white","red"),transparency = 0.5)
chordDiagram(new_CV,
             annotationTrack =  c("name", "grid"),
             grid.col=c("#E6A005", "#56B4E9", "#009E73","#F0E442","#0072B2", "#D55E00" ,"#CC79A7", "#666666","#AD7700","#1C91D4","#007756","#D5C711","#005685", "#A04700", "#B14380"),
             transparency = 0.5,
             col=cols,
             preAllocateTracks = list(track.height = 0.1))
title("CV (AD patients)", line=-4, cex.main = 5)
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
              niceFacing = TRUE, Controlj = c(0, 0.5), col = "black", cex = 4)
  
}, bg.border = NA)

## R62H
new_R62H=duplicate_pairs(R62H_normal)
cols = colorRamp2(c(min(new_R62H$sum_sigval),0,max(new_R62H$sum_sigval)),c("blue","white","red"),transparency = 0.5)
chordDiagram(new_R62H,
             annotationTrack =  c("name", "grid"),
             grid.col=c("#E6A005", "#56B4E9", "#009E73","#F0E442","#0072B2", "#D55E00" ,"#CC79A7", "#666666","#AD7700","#1C91D4","#007756","#D5C711","#005685", "#A04700", "#B14380"),
             transparency = 0.5,
             col=cols,
             preAllocateTracks = list(track.height = 0.1))
title("R62H (AD patients)", line=-4, cex.main = 5)
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
              niceFacing = TRUE, Controlj = c(0, 0.5), col = "black", cex = 4)
  
}, bg.border = NA)

## R47H
new_R47H=duplicate_pairs(R47H_normal)
cols = colorRamp2(c(min(new_R47H$sum_sigval),0,max(new_R47H$sum_sigval)),c("blue","white","red"),transparency = 0.5)
chordDiagram(new_R47H,
             annotationTrack =  c("name", "grid"),
             grid.col=c("#E6A005", "#56B4E9", "#009E73","#F0E442","#0072B2", "#D55E00" ,"#CC79A7", "#666666","#AD7700","#1C91D4","#007756","#D5C711","#005685", "#A04700", "#B14380"),
             transparency = 0.5,
             col=cols,
             preAllocateTracks = list(track.height = 0.1))
title("R47H (AD patients)", line=-4, cex.main = 5)
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
              niceFacing = TRUE, Controlj = c(0, 0.5), col = "black", cex = 4)
  
}, bg.border = NA)
#saveRDS(spe_filtered_AD,'spe_filtered_AD_new.rds')
#####################################################
# using spe_filtered
# For AD and Control only
spe_filtered <- buildSpatialGraph(spe_filtered, img_id = "sample_id", type = "knn", k = 20)

# spe_filtered_AD <- buildSpatialGraph(spe_filtered_AD, img_id = "sample_id", type = "expansion", threshold = 30)

spe_filtered$cluster_names<-as.factor(spe_filtered$nn_clusters_corrected14)
out3 <- testInteractions(spe_filtered, 
                         group_by = "sample_id",
                         label = "cluster_names", 
                         colPairName = "knn_interaction_graph",
                         BPPARAM = SerialParam(RNGseed = 221029)) #221029

out3$diagnosis=spe_filtered$diagnosis[match(out3$group_by,spe_filtered$sample_id)]


out3%>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

AD=out3%>% as_tibble() %>%
  filter(diagnosis=='AD') %>% 
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE),
            diagnosis='AD')
AD = as.data.frame(AD)
Control=out3%>% as_tibble() %>%
  filter(diagnosis=='Control') %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE),
            diagnosis='Control')
Control = as.data.frame(Control)
all=rbind(AD,Control)

set.seed(123)

# chordDiagram(AD,
#              annotationTrack =  c("name", "grid"),
#              grid.col=c(  "yellow", "#C2E699", "springgreen","mediumseagreen","#EFF3FF", "#BDD7E7" ,"#6BAED6", "#2171B5","violetred","salmon","seagreen1","lavender","skyblue", "moccasin", "plum"),
#              transparency = 0.5)
# 
# chordDiagram(Control,
#              annotationTrack =  c("name", "grid"),
#              grid.col=c(  "yellow", "#C2E699", "springgreen","mediumseagreen","#EFF3FF", "#BDD7E7" ,"#6BAED6", "#2171B5","violetred","salmon","seagreen1","lavender","skyblue", "moccasin", "plum"),
#              transparency = 0.5)

### removing self links
all$self=ifelse(all$from_label==all$to_label,'self','normal')

all_normal=all%>%
  filter(self!='self')

AD$self=ifelse(AD$from_label==AD$to_label,'self','normal')
#AD$double=paste(AD$from_label,AD$to_label)
AD_normal=AD%>%
  filter(self!='self')
#AD_normal[,c("from_label", "to_label", "sum_sigval")]

Control$self=ifelse(Control$from_label==Control$to_label,'self','normal')
Control_normal=Control%>%
  filter(self!='self')

library(RColorBrewer)
# brewer.pal(n = 7, name = "YlGn")
# brewer.pal(n = 8, name = "Blues")
print(unique(AD_normal$from_label))
print(unique(Control_normal$from_label))
par(cex = 0.2, mar = c(0, 0, 0, 0))

# # Define sector names and corresponding colors
# sector_names <- c("1: CD163 + HLA-DR + ApoE + Ab + Iba1 + CD68 + pTau + TSPO + GFAP",
#                   "10: pTau",
#                   "11: PLP",
#                   "12: -ApoE -TSPO -NfL -pTau",
#                   "13: -Ab -NfL",
#                   "14: -MAP2 -NfL -Ab -PLP",
#                   "15: NfL + PLP + MAP2 + Ab",
#                   "2: CD68 + Iba1 + TSPO",
#                   "3: pTau + NfL + GFAP + PLP + Ab",
#                   "4: MAP2 + ApoE + NfL + GFAP",
#                   "5: Ab",
#                   "6: GFAP + TSPO",
#                   "7: TSPO + GFAP",
#                   "8: NfL + pTau + PLP",
#                   "9: MAP2"
# )
# 
# sector_colors <- c("yellow", "#C2E699", "springgreen", "mediumseagreen",
#                    "#EFF3FF", "#BDD7E7", "#6BAED6", "#2171B5", "violetred",
#                    "salmon", "seagreen1", "lavender", "skyblue", "moccasin",
#                    "plum"
# )
# 
# # Assign names to colors in grid.col vector
# grid_col <- setNames(sector_colors, sector_names)
# 
# 
# chordDiagram(new_AD,
#              annotationTrack =c("name", "grid"),
#              order=sector_names,
#              grid.col=grid_col,
#              col=cols)
## AD
new_AD=duplicate_pairs(AD_normal)
cols = colorRamp2(c(min(new_AD$sum_sigval),0,max(new_AD$sum_sigval)),c("blue","white","red"),transparency = 0.5)
chordDiagram(new_AD,
             annotationTrack =  c("name", "grid"),
             grid.col=c("#E6A005", "#56B4E9", "#009E73","#F0E442","#0072B2", "#D55E00" ,"#CC79A7", "#666666","#AD7700","#1C91D4","#007756","#D5C711","#005685", "#A04700", "#B14380"),
             transparency = 0.5,
             col=cols,
             preAllocateTracks = list(track.height = 0.1))
title("AD", line=-4, cex.main = 5)
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
              niceFacing = TRUE, Controlj = c(0, 0.5), col = "black", cex = 4)
  
}, bg.border = NA)

## Control
new_Control=duplicate_pairs(Control_normal)
cols = colorRamp2(c(min(new_Control$sum_sigval),0,max(new_Control$sum_sigval)),c("blue","white","red"),transparency = 0.5)
chordDiagram(new_Control,
             annotationTrack =  c("name", "grid"),
             grid.col=c("#E6A005", "#56B4E9", "#009E73","#F0E442","#0072B2", "#D55E00" ,"#CC79A7", "#666666","#AD7700","#1C91D4","#007756","#D5C711","#005685", "#A04700", "#B14380"),
             transparency = 0.5,
             col=cols,
             preAllocateTracks = list(track.height = 0.1))
title("Control", line=-4, cex.main = 5)
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
              niceFacing = TRUE, Controlj = c(0, 0.5), col = "black", cex = 4)
  
}, bg.border = NA)


# Modify chordDiagram function calls
# chordDiagram(new_AD,
#              annotationTrack = c("grid"),
#              order = sector_names_new_AD,
#              grid.col = grid_col,
#              transparency = 0.5,
#              col = cols)
# 
# # Add labels to the chord diagram
# circos.trackText(1:length(sector_names_new_AD),
#                  rep(1, length(sector_names_new_AD)),  # Specify y position as 1 for all labels
#                  sector_names_new_AD,
#                  facing = "inside",
#                  niceFacing = TRUE,
#                  adj = c(0, 0.5),
#                  cex = 0.6)

# circos.text(
#   x = c(0, 1),
#   y = c(0, 1),
#   labels = unique(AD_normal$from_label),
#   sector.index = get.current.sector.index(),
#   facing = "clockwise",
#   niceFacing = TRUE,
#   adj = c(0, 0),
#   cex = 0.6
# )

################################## Cell type
spe_filtered$cell_type <- ifelse(spe_filtered$nn_clusters_corrected14 %in% c(4, 8, 9, 10, 15),"Neuronal cells",
                                ifelse(spe_filtered$nn_clusters_corrected14 %in% c(2),"Activated Microglia",
                                       ifelse(spe_filtered$nn_clusters_corrected14 %in% c(11),"Oligodendrocytes",
                                        ifelse(spe_filtered$nn_clusters_corrected14 %in% c(6, 7),"Astrocytes", "Other cells"))))
# 
# # using spe_filtered
# # For AD and Control only
spe_filtered <- buildSpatialGraph(spe_filtered, img_id = "sample_id", type = "knn", k = 20)

# spe_filtered_AD <- buildSpatialGraph(spe_filtered_AD, img_id = "sample_id", type = "expansion", threshold = 30)
spe_filtered$cell_type<-as.factor(spe_filtered$cell_type)
out3 <- testInteractions(spe_filtered,
                         group_by = "sample_id",
                         label = "cell_type",
                         colPairName = "knn_interaction_graph",
                         BPPARAM = SerialParam(RNGseed = 221029)) #221029

out3$diagnosis=spe_filtered$diagnosis[match(out3$group_by,spe_filtered$sample_id)]


out3%>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

AD=out3%>% as_tibble() %>%
  filter(diagnosis=='AD') %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE),
            diagnosis='AD')
AD = as.data.frame(AD)
Control=out3%>% as_tibble() %>%
  filter(diagnosis=='Control') %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE),
            diagnosis='Control')
Control = as.data.frame(Control)
all=rbind(AD,Control)

set.seed(123)

### removing self links
all$self=ifelse(all$from_label==all$to_label,'self','normal')

all_normal=all%>%
  filter(self!='self')

AD$self=ifelse(AD$from_label==AD$to_label,'self','normal')
#AD$double=paste(AD$from_label,AD$to_label)
AD_normal=AD%>%
  filter(self!='self')
#AD_normal[,c("from_label", "to_label", "sum_sigval")]

Control$self=ifelse(Control$from_label==Control$to_label,'self','normal')
Control_normal=Control%>%
  filter(self!='self')

library(RColorBrewer)
# brewer.pal(n = 7, name = "YlGn")
# brewer.pal(n = 8, name = "Blues")
print(unique(AD_normal$from_label))
print(unique(Control_normal$from_label))
par(cex = 0.2, mar = c(0, 0, 0, 0))

## AD
new_AD=duplicate_pairs(AD_normal)
cols = colorRamp2(c(min(new_AD$sum_sigval),0,max(new_AD$sum_sigval)),c("blue","white","red"),transparency = 0.5)
chordDiagram(new_AD,
             annotationTrack =  c("grid"),
             grid.col=c("#E6A005", "#56B4E9", "#009E73","#F0E442","#0072B2"),
             transparency = 0.5,
             col=cols,
             preAllocateTracks = list(track.height = 0.1))
title("AD", line=-4, cex.main = 5)
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")

  circos.text(mean(xlim), ylim[1], sector.name, facing = "outside",
              niceFacing = TRUE, Controlj = c(0, 0.5), col = "black", cex = 4)

}, bg.border = NA)

## Control
new_Control=duplicate_pairs(Control_normal)
cols = colorRamp2(c(min(new_Control$sum_sigval),0,max(new_Control$sum_sigval)),c("blue","white","red"),transparency = 0.5)
chordDiagram(new_Control,
             annotationTrack =  c("grid"),
             grid.col=c("#E6A005", "#56B4E9", "#009E73","#F0E442","#0072B2"),
             transparency = 0.5,
             col=cols,
             preAllocateTracks = list(track.height = 0.1))
title("Control", line=-4, cex.main = 5)
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")

  circos.text(mean(xlim), ylim[1], sector.name, facing = "outside",
              niceFacing = TRUE, Controlj = c(0, 0.5), col = "black", cex = 4)

}, bg.border = NA)

###### AD only
spe_filtered_AD$cell_type <- ifelse(spe_filtered_AD$nn_clusters_corrected14 %in% c(4, 8, 9, 10, 15),"Neuronal cells",
                                 ifelse(spe_filtered_AD$nn_clusters_corrected14 %in% c(2),"Activated Microglia",
                                        ifelse(spe_filtered_AD$nn_clusters_corrected14 %in% c(11),"Oligodendrocytes",
                                               ifelse(spe_filtered_AD$nn_clusters_corrected14 %in% c(6, 7),"Astrocytes", "Other cells"))))
spe_filtered_AD <- buildSpatialGraph(spe_filtered_AD, img_id = "sample_id", type = "knn", k = 20)
# spe_filtered_AD <- buildSpatialGraph(spe_filtered_AD, img_id = "sample_id", type = "expansion", threshold = 30)
spe_filtered_AD$cell_type<-as.factor(spe_filtered_AD$cell_type)

out2 <- testInteractions(spe_filtered_AD, 
                         group_by = "sample_id",
                         label = "cell_type", 
                         colPairName = "knn_interaction_graph",
                         BPPARAM = SerialParam(RNGseed = 221029))


out2$trem2=spe_filtered_AD$trem2[match(out2$group_by,spe_filtered_AD$sample_id)]

out2%>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


out2%>% as_tibble() %>%
  filter(trem2=='R47H') %>% 
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

out2%>% as_tibble() %>%
  filter(trem2=='R62H') %>% 
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


out2%>% as_tibble() %>%
  filter(trem2=='CV') %>% 
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


R47H=out2%>% as_tibble() %>%
  filter(trem2=='R47H') %>% 
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE),
            trem2='R47H')


R62H=out2%>% as_tibble() %>%
  filter(trem2=='R62H') %>% 
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE),
            trem2='R62H')


CV=out2%>% as_tibble() %>%
  filter(trem2=='CV') %>% 
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE),
            trem2='CV')

all=rbind(R47H,R62H,CV)
library(circlize)
set.seed(123)

duplicate_pairs<-function(df){
  # Create unique identifier for each combination
  df$identifier <- apply(df[, c("from_label", "to_label")], 1, function(x) {
    sorted_pair <- sort(x)
    paste(sorted_pair, collapse = "-")
  })
  
  # Remove duplicated pairs
  df_unique <- df[!duplicated(df$identifier), ]
  return(df_unique)
}


# chordDiagram(CV,
#              annotationTrack =  c("name", "grid"),
#              grid.col=c("#E6A005", "#56B4E9", "#009E73","#F0E442","#0072B2", "#D55E00" ,"#CC79A7", "#666666","#AD7700","#1C91D4","#007756","#D5C711","#005685", "#A04700", "#B14380"),
#              transparency = 0.5)
# 
# 
# chordDiagram(R62H,
#              annotationTrack =  c("name", "grid"),
#              grid.col=c("#E6A005", "#56B4E9", "#009E73","#F0E442","#0072B2", "#D55E00" ,"#CC79A7", "#666666","#AD7700","#1C91D4","#007756","#D5C711","#005685", "#A04700", "#B14380"),
#              transparency = 0.5)
# 
# chordDiagram(R47H,
#              annotationTrack =  c("name", "grid"),
#              grid.col=c("#E6A005", "#56B4E9", "#009E73","#F0E442","#0072B2", "#D55E00" ,"#CC79A7", "#666666","#AD7700","#1C91D4","#007756","#D5C711","#005685", "#A04700", "#B14380"),
#              transparency = 0.5)


### removing self links
all$self=ifelse(all$from_label==all$to_label,'self','normal')

all_normal=all%>%
  filter(self!='self')

CV$self=ifelse(CV$from_label==CV$to_label,'self','normal')
CV$double=paste(CV$from_label,CV$to_label)

CV_normal=CV%>%
  filter(self!='self')

R62H$self=ifelse(R62H$from_label==R62H$to_label,'self','normal')

R62H_normal=R62H%>%
  filter(self!='self')

R47H$self=ifelse(R47H$from_label==R47H$to_label,'self','normal')

R47H_normal=R47H%>%
  filter(self!='self')

library(RColorBrewer)

# print(unique(CV_normal$from_label))
# par(cex = 0.2, mar = c(0, 0, 0, 0))

## CV
new_CV=duplicate_pairs(CV_normal)
cols = colorRamp2(c(min(new_CV$sum_sigval),0,max(new_CV$sum_sigval)),c("blue","white","red"),transparency = 0.5)
chordDiagram(new_CV,
             annotationTrack =  c("grid"),
             grid.col=c("#E6A005", "#56B4E9", "#009E73","#F0E442","#0072B2"),
             transparency = 0.5,
             col=cols,
             preAllocateTracks = list(track.height = 0.1))
title("CV (AD patients)", line=-4, cex.main = 5)
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  circos.text(mean(xlim), ylim[1], sector.name, facing = "outside",
              niceFacing = TRUE, Controlj = c(0, 0.5), col = "black", cex = 4)
  
}, bg.border = NA)

## R62H
new_R62H=duplicate_pairs(R62H_normal)
cols = colorRamp2(c(min(new_R62H$sum_sigval),0,max(new_R62H$sum_sigval)),c("blue","white","red"),transparency = 0.5)
chordDiagram(new_R62H,
             annotationTrack =  c("grid"),
             grid.col=c("#E6A005", "#56B4E9", "#009E73","#F0E442","#0072B2"),
             transparency = 0.5,
             col=cols,
             preAllocateTracks = list(track.height = 0.1))
title("R62H (AD patients)", line=-4, cex.main = 5)
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  circos.text(mean(xlim), ylim[1], sector.name, facing = "outside",
              niceFacing = TRUE, Controlj = c(0, 0.5), col = "black", cex = 4)
  
}, bg.border = NA)

## R47H
new_R47H=duplicate_pairs(R47H_normal)
cols = colorRamp2(c(min(new_R47H$sum_sigval),0,max(new_R47H$sum_sigval)),c("blue","white","red"),transparency = 0.5)
chordDiagram(new_R47H,
             annotationTrack =  c("grid"),
             grid.col=c("#E6A005", "#56B4E9", "#009E73","#F0E442","#0072B2"),
             transparency = 0.5,
             col=cols,
             preAllocateTracks = list(track.height = 0.1))
title("R47H (AD patients)", line=-4, cex.main = 5)
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  circos.text(mean(xlim), ylim[1], sector.name, facing = "outside",
              niceFacing = TRUE, Controlj = c(0, 0.5), col = "black", cex = 4)
  
}, bg.border = NA)

##########################
## neighborhood
out4 <- testInteractions(spe_filtered,
                         group_by = "sample_id",
                         label = "cell_type",
                         colPairName = "neighborhood",
                         BPPARAM = SerialParam(RNGseed = 221029)) #221029

out4$diagnosis=spe_filtered$diagnosis[match(out4$group_by,spe_filtered$sample_id)]


out4%>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

AD=out4%>% as_tibble() %>%
  filter(diagnosis=='AD') %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE),
            diagnosis='AD')
AD = as.data.frame(AD)
Control=out4%>% as_tibble() %>%
  filter(diagnosis=='Control') %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE),
            diagnosis='Control')
Control = as.data.frame(Control)
all=rbind(AD,Control)

set.seed(123)

### removing self links
all$self=ifelse(all$from_label==all$to_label,'self','normal')

all_normal=all%>%
  filter(self!='self')

AD$self=ifelse(AD$from_label==AD$to_label,'self','normal')
#AD$double=paste(AD$from_label,AD$to_label)
AD_normal=AD%>%
  filter(self!='self')
#AD_normal[,c("from_label", "to_label", "sum_sigval")]

Control$self=ifelse(Control$from_label==Control$to_label,'self','normal')
Control_normal=Control%>%
  filter(self!='self')

library(RColorBrewer)
# brewer.pal(n = 7, name = "YlGn")
# brewer.pal(n = 8, name = "Blues")
print(unique(AD_normal$from_label))
print(unique(Control_normal$from_label))
par(cex = 0.2, mar = c(0, 0, 0, 0))

## AD
new_AD=duplicate_pairs(AD_normal)
cols = colorRamp2(c(min(new_AD$sum_sigval),0,max(new_AD$sum_sigval)),c("blue","white","red"),transparency = 0.5)
chordDiagram(new_AD,
             annotationTrack =  c("grid"),
             grid.col=c("#E6A005", "#56B4E9", "#009E73","#F0E442","#0072B2"),
             transparency = 0.5,
             col=cols,
             preAllocateTracks = list(track.height = 0.1))
title("AD", line=-4, cex.main = 5)
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  circos.text(mean(xlim), ylim[1], sector.name, facing = "outside",
              niceFacing = TRUE, Controlj = c(0, 0.5), col = "black", cex = 4)
  
}, bg.border = NA)

## Control
new_Control=duplicate_pairs(Control_normal)
cols = colorRamp2(c(min(new_Control$sum_sigval),0,max(new_Control$sum_sigval)),c("blue","white","red"),transparency = 0.5)
chordDiagram(new_Control,
             annotationTrack =  c("grid"),
             grid.col=c("#E6A005", "#56B4E9", "#009E73","#F0E442","#0072B2"),
             transparency = 0.5,
             col=cols,
             preAllocateTracks = list(track.height = 0.1))
title("Control", line=-4, cex.main = 5)
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  circos.text(mean(xlim), ylim[1], sector.name, facing = "outside",
              niceFacing = TRUE, Controlj = c(0, 0.5), col = "black", cex = 4)
  
}, bg.border = NA)