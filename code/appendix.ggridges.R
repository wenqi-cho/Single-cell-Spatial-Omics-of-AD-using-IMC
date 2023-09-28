library(ggridges)

### SOM
## Plot - naming NAs as 200
metadata=data.frame(colData(spe_filtered))
metadata$distance2=ifelse(metadata$distance=='NA',200,metadata$distance)
# metadata$cluster=metadata$som_clusters_corrected
# metadata$cluster_pg=metadata$pg_clusters_corrected
metadata$cluster=metadata$nn_clusters_corrected14
metadata$distance2=as.numeric(metadata$distance2)

metadata_f=metadata[which(metadata$distance!='NA'),]
metadata_f$distance=as.numeric(metadata_f$distance)

# ggplot(metadata_f, aes(x = distance, y = cluster, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
#   scale_fill_viridis_c(name = "Tail probability", direction = -1) #cluster 3, 4, 6, 8

## Plot - removing NAs and distance = 100
metadata_ff=metadata[which(metadata$distance!='NA' & metadata$distance!="100"),]
metadata_ff$distance=as.numeric(metadata_ff$distance)

metadata_ff$cluster=as.factor(metadata_ff$cluster)
### SOM
# ggplot(metadata_ff, aes(x = distance, y = cluster, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
#   scale_fill_viridis_c(name = "Tail probability", direction = -1) #cluster 3,6, 8
# 
# ggplot(metadata_ff, aes(x = distance, y = trem2, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
#   scale_fill_viridis_c(name = "Tail probability", direction = -1)
# 
# ggplot(metadata_ff, aes(x = distance, y = diagnosis, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
#   scale_fill_viridis_c(name = "Tail probability", direction = -1)
# 
# ###
# metadata_ff1=metadata_ff[which(metadata_ff$cluster=="1"),]
# metadata_ff3=metadata_ff[which(metadata_ff$cluster=="3"),]
# metadata_ff6=metadata_ff[which(metadata_ff$cluster=="6"),]
# metadata_ff8=metadata_ff[which(metadata_ff$cluster=="8"),]
# 
# ggplot(metadata_ff1, aes(x = distance, y = diagnosis, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
#   scale_fill_viridis_c(name = "Tail probability", direction = -1)
# 
# ggplot(metadata_ff3, aes(x = distance, y = diagnosis, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
#   scale_fill_viridis_c(name = "Tail probability", direction = -1)
# 
# ggplot(metadata_ff6, aes(x = distance, y = diagnosis, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
#   scale_fill_viridis_c(name = "Tail probability", direction = -1)
# 
# ggplot(metadata_ff6, aes(x = distance, y = trem2, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
#   scale_fill_viridis_c(name = "Tail probability", direction = -1)
# 
# ggplot(metadata_ff8, aes(x = distance, y = diagnosis, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
#   scale_fill_viridis_c(name = "Tail probability", direction = -1)
# 
# 
# ### Phenograph
# ## Plot - naming NAs as 200
# ggplot(metadata_f, aes(x = distance, y = cluster_pg, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
#   scale_fill_viridis_c(name = "Tail probability", direction = -1) 
# # Removing distance = NA and 100
# ggplot(metadata_ff, aes(x = distance, y = cluster_pg, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
#   scale_fill_viridis_c(name = "Tail probability", direction = -1) #cluster 10
# 
# ###
# metadata_ff_pg10=metadata_ff[which(metadata_ff$cluster_pg=="10"),]
# 
# ggplot(metadata_ff_pg10, aes(x = distance, y = diagnosis, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
#   scale_fill_viridis_c(name = "Tail probability", direction = -1)
# 
# ggplot(metadata_ff_pg10, aes(x = distance, y = trem2, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
#   scale_fill_viridis_c(name = "Tail probability", direction = -1)


### SNN
## Plot - naming NAs as 200
ggplot(metadata_ff, aes(x = distance, y = cluster, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

ggplot(metadata_ff, aes(x = distance, y = cluster_nn, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1) #cluster 2, 5

ggplot(metadata_ff, aes(x = distance, y = trem2, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

ggplot(metadata_ff, aes(x = distance, y = diagnosis, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

###
metadata_ff_nn2=metadata_ff[which(metadata_ff$cluster_nn=="2"),]
metadata_ff_nn5=metadata_ff[which(metadata_ff$cluster_nn=="5"),]
metadata_ff_nn12=metadata_ff[which(metadata_ff$cluster_nn=="12"),]
metadata_ff_nn15=metadata_ff[which(metadata_ff$cluster_nn=="15"),]

ggplot(metadata_ff_nn2, aes(x = distance, y = diagnosis, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

ggplot(metadata_ff_nn2, aes(x = distance, y = trem2, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)



##### Only AD
# spe_filtered_AD=spe_filtered[,colData(spe_filtered)$diagnosis=='AD'] In dirichlet.R
metadata_AD=data.frame(colData(spe_filtered_AD))
metadata_AD$distance2=ifelse(metadata_AD$distance=='NA',200,metadata_AD$distance)

metadata_AD$cluster_nn=metadata_AD$nn_clusters_corrected14
metadata_AD$distance2=as.numeric(metadata_AD$distance2)

metadata_f=metadata_AD[which(metadata_AD$distance!='NA'),]
metadata_f$distance=as.numeric(metadata_f$distance)

## Plot - removing NAs and distance = 100
metadata_ff=metadata_AD[which(metadata_AD$distance!='NA' & metadata_AD$distance!="100"),]
metadata_ff$distance=as.numeric(metadata_ff$distance)

## Plot
ggplot(metadata_f, aes(x = distance, y = cluster_nn, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

ggplot(metadata_ff, aes(x = distance, y = cluster_nn, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1) #cluster 2, 5

ggplot(metadata_ff, aes(x = distance, y = trem2, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

ggplot(metadata_ff, aes(x = distance, y = diagnosis, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

###
metadata_ff_nn2=metadata_ff[which(metadata_ff$cluster_nn=="2"),]
metadata_ff_nn3=metadata_ff[which(metadata_ff$cluster_nn=="3"),]
metadata_ff_nn4=metadata_ff[which(metadata_ff$cluster_nn=="4"),]
metadata_ff_nn5=metadata_ff[which(metadata_ff$cluster_nn=="5"),]
metadata_ff_nn12=metadata_ff[which(metadata_ff$cluster_nn=="12"),]
metadata_ff_nn15=metadata_ff[which(metadata_ff$cluster_nn=="15"),]

ggplot(metadata_ff_nn2, aes(x = distance, y = diagnosis, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

ggplot(metadata_ff_nn2, aes(x = distance, y = trem2, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

ggplot(metadata_ff_nn3, aes(x = distance, y = diagnosis, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

ggplot(metadata_ff_nn3, aes(x = distance, y = trem2, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

ggplot(metadata_ff_nn4, aes(x = distance, y = diagnosis, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

ggplot(metadata_ff_nn4, aes(x = distance, y = trem2, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

ggplot(metadata_ff_nn5, aes(x = distance, y = diagnosis, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

ggplot(metadata_ff_nn5, aes(x = distance, y = trem2, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

ggplot(metadata_ff_nn12, aes(x = distance, y = diagnosis, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

ggplot(metadata_ff_nn12, aes(x = distance, y = trem2, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

ggplot(metadata_ff_nn15, aes(x = distance, y = diagnosis, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

ggplot(metadata_ff_nn15, aes(x = distance, y = trem2, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)
