library(dplyr)
mat <- reducedDim(spe_filtered, "fastMNN")
metadata=data.frame(colData(spe_filtered))
sce <- SingleCellExperiment(assays = list(counts = t(mat)),colData=metadata)
sce$groups=ifelse(sce$diagnosis=='AD'& sce$trem2_all=='TREM2var', 'AD_TREM2',
                  ifelse(sce$diagnosis=='AD'& sce$trem2_all=='CV','AD_CV',
                         ifelse(sce$diagnosis=='Control'& sce$trem2_all=='TREM2var','CTRL_TREM2',
                                ifelse(sce$diagnosis=='Control'& sce$trem2_all=='CV','CTRL_CV',NA))))
sce$distance_group=ifelse(sce$distance==10,"<=10um",
                                 ifelse(sce$distance==20 |sce$distance==30,"10um<cell<=30um",
                                        ifelse(sce$distance==40 |sce$distance==50,"30um<cell<=50um",
                                               ifelse(sce$distance==1,'plaque','>50um'))))
distance_group_order <- c("plaque", "<=10um", "10<cell<=30um", "30<cell<=50um", ">50um")
sce$distance_group <-factor(sce$distance_group, levels=distance_group_order)
# ### SOM
# sce$cluster=sce$som_clusters_corrected
# dirichelet=model_celltype_freqs(sce,
#                                 unique_id_var = "sample_id",
#                                 celltype_var = "cluster",
#                                 dependent_var = "groups",
#                                 ref_class = "CTRL_CV",
#                                 confounding_vars= c('sex','PMD'))
# 
# dirichelet[["dirichlet_plot"]]
# 
# # sce_sub=sce[,sce$cluster %in% c(3,4,6,8)]
# # ## test 
# # dirichelet_t=model_celltype_freqs(sce_sub,
# #                                   unique_id_var = "sample_id",
# #                                   celltype_var = "cluster",
# #                                   dependent_var = "groups",
# #                                   ref_class = "CTRL_CV", 
# #                                   confounding_vars= c('sex','PMD')
# # )


# dirichelet_t[["dirichlet_plot"]]

spe_filtered_AD=spe_filtered[,colData(spe_filtered)$diagnosis=='AD']
mat_AD <- reducedDim(spe_filtered_AD, "fastMNN")
metadata=data.frame(colData(spe_filtered_AD))
sce_AD <- SingleCellExperiment(assays = list(counts = t(mat_AD)),colData=metadata)
sce_AD$distance_group=ifelse(sce_AD$distance==10,"<=10um",
                             ifelse(sce_AD$distance==20 |sce_AD$distance==30,"10um<cell<=30um",
                                    ifelse(sce_AD$distance==40 |sce_AD$distance==50,"30um<cell<=50um",
                                           ifelse(sce_AD$distance==1,'plaque','>50um'))))
# sce_AD$cluster=sce_AD$som_clusters_corrected
# # Without confounding
# dirichelet_AD=model_celltype_freqs(sce_AD,
#                                    unique_id_var = "sample_id",
#                                    celltype_var = "cluster",
#                                    dependent_var = "trem2",
#                                    ref_class = "CV")
# # With confounding
# dirichelet_AD=model_celltype_freqs(sce_AD,
#                                    unique_id_var = "sample_id",
#                                    celltype_var = "cluster",
#                                    dependent_var = "trem2",
#                                    ref_class = "CV",
#                                    confounding_vars= c('sex','PMD'))
# 
# 
# dirichelet_AD[["dirichlet_plot"]]


###########
### Phenograph
# sce$cluster_pg=sce$pg_clusters_corrected
# dirichelet_pg=model_celltype_freqs(sce,
#                                 unique_id_var = "sample_id",
#                                 celltype_var = "cluster_pg",
#                                 dependent_var = "groups",
#                                 ref_class = "CTRL_CV",
#                                 confounding_vars= c('sex','PMD'))
# 
# dirichelet_pg[["dirichlet_plot"]]


# sce_sub=sce[,sce$cluster %in% c(3,4,6,8)]
# ## test 
# dirichelet_t=model_celltype_freqs(sce_sub,
#                                   unique_id_var = "sample_id",
#                                   celltype_var = "cluster",
#                                   dependent_var = "groups",
#                                   ref_class = "CTRL_CV", 
#                                   confounding_vars= c('sex','PMD')
# )
# 
# 
# dirichelet_t[["dirichlet_plot"]]


# sce_AD$cluster_pg=sce_AD$pg_clusters_corrected
# # Without confounding
# dirichelet_AD_pg=model_celltype_freqs(sce_AD,
#                                    unique_id_var = "sample_id",
#                                    celltype_var = "cluster_pg",
#                                    dependent_var = "trem2",
#                                    ref_class = "CV")
# With confounding
# dirichelet_AD_pg=model_celltype_freqs(sce_AD,
#                                    unique_id_var = "sample_id",
#                                    celltype_var = "cluster_pg",
#                                    dependent_var = "trem2",
#                                    ref_class = "CV",
#                                    confounding_vars= c('sex','PMD'))
# 
# 
# dirichelet_AD_pg[["dirichlet_plot"]]

############
### SNN
sce$cluster_nn=sce$nn_clusters_corrected14
sce$cluster_nn=factor(sce$cluster_nn, levels=unique(sce$cluster_nn))
## Without confounding
# dirichelet_nn_noconf=model_celltype_freqs(sce,
#                                    unique_id_var = "sample_id",
#                                    celltype_var = "cluster_nn",
#                                    dependent_var = "groups",
#                                    ref_class = "CTRL_CV")
# dirichelet_nn_noconf[["dirichlet_plot"]]

                                   
dirichelet_nn=model_celltype_freqs(sce,
                                   unique_id_var = "sample_id",
                                   celltype_var = "cluster_nn",
                                   dependent_var = "groups",
                                   ref_class = "CTRL_CV",
                                   confounding_vars= c('age', 'sex','PMD'))

dirichelet_nn[["dirichlet_plot"]]


## distance group
sce$cluster_nn=sce$nn_clusters_corrected14
## With confounding
dirichelet_nn_noconf=model_celltype_freqs(sce,
                                          unique_id_var = "sample_id",
                                          celltype_var = "distance_group",
                                          dependent_var = "trem2",
                                          ref_class = "CV",
                                          confounding_vars= c('age', 'sex','PMD'))
dirichelet_nn_noconf[["dirichlet_plot"]]

dirichelet_nn_noconf=model_celltype_freqs(sce,
                                          unique_id_var = "sample_id",
                                          celltype_var = "distance_group",
                                          dependent_var = "groups",
                                          ref_class = "CTRL_CV",
                                          confounding_vars= c('age', 'sex','PMD'))
dirichelet_nn_noconf[["dirichlet_plot"]]

dirichelet_nn=model_celltype_freqs(sce,
                                   unique_id_var = "sample_id",
                                   celltype_var = "cluster_nn",
                                   dependent_var = "groups",
                                   ref_class = "CTRL_CV",
                                   confounding_vars= c('age', 'sex','PMD'))

dirichelet_nn[["dirichlet_plot"]]

## Astrocytes
sce_cluster6_7 = sce[,colData(sce)$cluster_nn == '6' | colData(sce)$cluster_nn == '7']
### With confounding
dirichelet_cluster6_7=model_celltype_freqs(sce_cluster6_7,
                                              unique_id_var = "sample_id",
                                              celltype_var = "distance_group",
                                              dependent_var = "trem2",
                                              ref_class = "CV",
                                              confounding_vars= c('age','sex','PMD'))

dirichelet_cluster6_7[["dirichlet_plot"]]
### AD PATIENTS ONLY
sce_AD$cluster_nn<-sce_AD$nn_clusters_corrected14
# # Without confounding
dirichelet_AD_nn_noconf=model_celltype_freqs(sce_AD,
                                   unique_id_var = "sample_id",
                                   celltype_var = "cluster_nn",
                                   dependent_var = "trem2",
                                   ref_class = "CV")
dirichelet_AD_nn_noconf[["dirichlet_plot"]]

# With confounding
## cluster 
dirichelet_AD_nn=model_celltype_freqs(sce_AD,
                                      unique_id_var = "sample_id",
                                      celltype_var = "cluster_nn",
                                      dependent_var = "trem2",
                                      ref_class = "CV",
                                      confounding_vars= c('age','sex','PMD'))

dirichelet_AD_nn[["dirichlet_plot"]]

## distance
dirichelet_AD_nn=model_celltype_freqs(sce_AD,
                                      unique_id_var = "sample_id",
                                      celltype_var = "distance",
                                      dependent_var = "trem2",
                                      ref_class = "CV",
                                      confounding_vars= c('age','sex','PMD'))

dirichelet_AD_nn[["dirichlet_plot"]]

## distance_group
dirichelet_AD_nn=model_celltype_freqs(sce_AD,
                                      unique_id_var = "sample_id",
                                      celltype_var = "distance_group",
                                      dependent_var = "trem2",
                                      ref_class = "CV",
                                      confounding_vars= c('age','sex','PMD'))

dirichelet_AD_nn[["dirichlet_plot"]]

## Cluster 1 only
sce_AD_cluster1=sce_AD[,colData(sce_AD)$cluster_nn=='1']
### Without confounding
dirichelet_AD_cluster1=model_celltype_freqs(sce_AD_cluster1,
                                            unique_id_var = "sample_id",
                                            celltype_var = "distance",
                                            dependent_var = "trem2",
                                            ref_class = "CV")

dirichelet_AD_cluster1[["dirichlet_plot"]]

### With confounding
dirichelet_AD_cluster1=model_celltype_freqs(sce_AD_cluster1,
                                            unique_id_var = "sample_id",
                                            celltype_var = "distance",
                                            dependent_var = "trem2",
                                            ref_class = "CV",
                                            confounding_vars= c('age','sex','PMD'))

dirichelet_AD_cluster1[["dirichlet_plot"]]

## Cluster 2 only
sce_AD_cluster2=sce_AD[,colData(sce_AD)$cluster_nn=='2']
### Without confounding
dirichelet_AD_cluster2=model_celltype_freqs(sce_AD_cluster2,
                                            unique_id_var = "sample_id",
                                            celltype_var = "distance",
                                            dependent_var = "trem2",
                                            ref_class = "CV")

dirichelet_AD_cluster2[["dirichlet_plot"]]

### With confounding
dirichelet_AD_cluster2=model_celltype_freqs(sce_AD_cluster2,
                                      unique_id_var = "sample_id",
                                      celltype_var = "distance",
                                      dependent_var = "trem2",
                                      ref_class = "CV",
                                      confounding_vars= c('age','sex','PMD'))

dirichelet_AD_cluster2[["dirichlet_plot"]]

## distance_group
### Without confounding
dirichelet_AD_cluster2=model_celltype_freqs(sce_AD_cluster2,
                                            unique_id_var = "sample_id",
                                            celltype_var = "distance_group",
                                            dependent_var = "trem2",
                                            ref_class = "CV")

dirichelet_AD_cluster2[["dirichlet_plot"]]
### With confounding 
dirichelet_AD_cluster2=model_celltype_freqs(sce_AD_cluster2,
                                            unique_id_var = "sample_id",
                                            celltype_var = "distance_group",
                                            dependent_var = "trem2",
                                            ref_class = "CV",
                                            confounding_vars= c('age','sex','PMD'))

dirichelet_AD_cluster2[["dirichlet_plot"]]

## Cluster 3 only
sce_AD_cluster3=sce_AD[,colData(sce_AD)$cluster_nn=='3']
### Without confounding
dirichelet_AD_cluster3=model_celltype_freqs(sce_AD_cluster3,
                                            unique_id_var = "sample_id",
                                            celltype_var = "distance",
                                            dependent_var = "trem2",
                                            ref_class = "CV")

dirichelet_AD_cluster3[["dirichlet_plot"]]

### With confounding
dirichelet_AD_cluster3=model_celltype_freqs(sce_AD_cluster3,
                                            unique_id_var = "sample_id",
                                            celltype_var = "distance",
                                            dependent_var = "trem2",
                                            ref_class = "CV",
                                            confounding_vars= c('age','sex','PMD'))

dirichelet_AD_cluster3[["dirichlet_plot"]]

## Cluster 5 only
sce_AD_cluster5=sce_AD[,colData(sce_AD)$cluster_nn=='5']

### Without confounding
dirichelet_AD_cluster5=model_celltype_freqs(sce_AD_cluster5,
                                            unique_id_var = "sample_id",
                                            celltype_var = "distance",
                                            dependent_var = "trem2",
                                            ref_class = "CV")

dirichelet_AD_cluster5[["dirichlet_plot"]]
                       
### With confounding
dirichelet_AD_cluster5=model_celltype_freqs(sce_AD_cluster5,
                                            unique_id_var = "sample_id",
                                            celltype_var = "distance",
                                            dependent_var = "trem2",
                                            ref_class = "CV",
                                            confounding_vars= c('age','sex','PMD'))

dirichelet_AD_cluster5[["dirichlet_plot"]]

## Cluster 6 & 7
sce_AD_cluster6_7 = sce_AD[,colData(sce_AD)$cluster_nn == '6' | colData(sce_AD)$cluster_nn == '7']
### With confounding
dirichelet_AD_cluster6_7=model_celltype_freqs(sce_AD_cluster6_7,
                                            unique_id_var = "sample_id",
                                            celltype_var = "distance_group",
                                            dependent_var = "trem2",
                                            ref_class = "CV",
                                            confounding_vars= c('age','sex','PMD'))

dirichelet_AD_cluster6_7[["dirichlet_plot"]]
## Cluster 8 only
sce_AD_cluster8=sce_AD[,colData(sce_AD)$cluster_nn=='8']

dirichelet_AD_cluster8=model_celltype_freqs(sce_AD_cluster8,
                                            unique_id_var = "sample_id",
                                            celltype_var = "distance",
                                            dependent_var = "trem2",
                                            ref_class = "CV",
                                            confounding_vars= c('age','sex','PMD'))

dirichelet_AD_cluster8[["dirichlet_plot"]]
