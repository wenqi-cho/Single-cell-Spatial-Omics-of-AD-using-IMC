library(scProportionTest)
library(forcats)
library(SingleCellExperiment)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(Seurat)
library(SeuratObject)
spe_filtered <- readRDS("spe_filtered1.rds")
spe_filtered$neuronal <- ifelse(spe_filtered$nn_clusters_corrected14 %in% c(4, 8, 9, 10, 15),
                                "Neuronal", "Non-neuronal")

seurat_obj <- as.Seurat(spe_filtered[rowData(spe_filtered)$use_channel], counts = "counts", data = "exprs")
seurat_obj <- AddMetaData(seurat_obj, as.data.frame(colData(spe_filtered)))
seurat_obj$distance_group_num=ifelse(seurat_obj$distance==10,"2",
                                     ifelse(seurat_obj$distance==20 |seurat_obj$distance==30,"3",
                                            ifelse(seurat_obj$distance==40 |seurat_obj$distance==50,"4",
                                                   ifelse(seurat_obj$distance==1,'1','5'))))
seurat_obj$distance_group=ifelse(spe_filtered$distance==10,"<=10um",
                                 ifelse(spe_filtered$distance==20 |spe_filtered$distance==30,"10um<cell<=30um",
                                        ifelse(spe_filtered$distance==40 |spe_filtered$distance==50,"30um<cell<=50um",
                                               ifelse(spe_filtered$distance==1,'plaque','>50um'))))
prop_test <- sc_utils(seurat_obj)

# Permutation Test For Proportions
permutation_test <- function(
  sc_utils_obj,
  cluster_identity = NA,
  sample_1 = NA,
  sample_2 = NA,
  sample_identity = "orig.ident",
  n_permutations = 1000
) {
  
  ## Prepare data.
  meta_data <- copy(sc_utils_obj@meta_data)
  
  meta_data <- meta_data[
    get(sample_identity) %in% c(sample_1, sample_2),
    c(..sample_identity, ..cluster_identity)
  ]
  
  setnames(
    meta_data,
    old = c(sample_identity, cluster_identity),
    new = c("samples", "clusters")
  )
  
  meta_data[, clusters := as.character(clusters)]
  cluster_cases <- unique(meta_data[["clusters"]])
  
  ## Get observed differences in fraction.
  obs_diff <- meta_data[, .(count = .N), by = .(samples, clusters)]
  obs_diff <- obs_diff[
    CJ(samples = samples, clusters = cluster_cases, unique = TRUE),
    on = .(samples, clusters)
  ][
    is.na(count), count := 0
  ][]
  obs_diff[, fraction := count / sum(count), by = samples]
  obs_diff <- dcast(obs_diff, clusters ~ samples, value.var = "fraction")
  obs_diff[, obs_log2FD := log2(get(sample_2)) - log2(get(sample_1))]
  
  ## Permutation test.
  perm_results <- matrix(NA, nrow(obs_diff), n_permutations)
  rownames(perm_results) <- sort(cluster_cases)
  
  for (i in seq_len(n_permutations)) {
    permuted <- copy(meta_data)
    permuted[["samples"]] <- sample(permuted[["samples"]])
    permuted <- permuted[, .(count = .N), by = .(samples, clusters)]
    permuted <- permuted[
      CJ(samples = samples, clusters = cluster_cases, unique = TRUE),
      on = .(samples, clusters)
    ][
      is.na(count), count := 0
    ][]
    permuted[, fraction := count / sum(count), by = samples]
    permuted <- dcast(permuted, clusters ~ samples, value.var = "fraction")
    permuted[, perm_log2FD := log2(get(sample_2)) - log2(get(sample_1))]
    
    perm_results[, i] <- permuted[["perm_log2FD"]]
  }
  
  increased <- rowSums(apply(perm_results, 2, function(x) obs_diff[["obs_log2FD"]] <= x))
  increased <- (increased + 1) / (n_permutations + 1)
  
  decreased <- rowSums(apply(perm_results, 2, function(x) obs_diff[["obs_log2FD"]] >= x))
  decreased <- (decreased + 1) / (n_permutations + 1)
  
  obs_diff[, pval := ifelse(obs_log2FD > 0, increased[.I], decreased[.I])]
  obs_diff[, FDR := p.adjust(pval, "fdr")]
  
  ## Boostrap log2FD CI.
  boot_results <- matrix(NA, nrow(obs_diff), n_permutations)
  rownames(boot_results) <- sort(cluster_cases)
  
  for (i in seq_len(n_permutations)) {
    booted <- copy(meta_data)
    booted[, clusters := sample(clusters, replace = TRUE), by = samples]
    booted <- booted[, .(count = .N), by = .(samples, clusters)]
    booted <- booted[
      CJ(samples = samples, clusters = cluster_cases, unique = TRUE),
      on = .(samples, clusters)
    ][
      is.na(count), count := 0
    ][]
    booted[, fraction := count / sum(count), by = samples]
    booted <- dcast(booted, clusters ~ samples, value.var = "fraction")
    booted[, boot_log2FD := log2(get(sample_2)) - log2(get(sample_1))]
    
    boot_results[, i] <- booted[["boot_log2FD"]]
  }
  
  boot_results[!is.finite(boot_results)] <- NA
  boot_mean <- rowMeans(boot_results, na.rm = TRUE)
  boot_ci <- t(apply(boot_results, 1, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)))
  boot_ci <- as.data.table(boot_ci)
  setnames(boot_ci, old = c(1, 2), new = c("boot_CI_2.5", "boot_CI_97.5"))
  
  obs_diff[, boot_mean_log2FD := boot_mean]
  obs_diff <- cbind(obs_diff, boot_ci)
  
  ## Store results and return object.
  sc_utils_obj@results$permutation <- obs_diff
  return(sc_utils_obj)
}

# Permutation test of neuronal group vs trem2
prop_test1 <- permutation_test(prop_test, cluster_identity = "neuronal",
                               sample_1 = "CV", sample_2 = "R47H",
                               sample_identity = "trem2")

prop_test2 <- permutation_test(prop_test, cluster_identity = "neuronal",
                               sample_1 = "CV", sample_2 = "R62H",
                               sample_identity = "trem2")

##### Function for permutation_plot
permutation_plot2 <- function(
  prop_test1,
  prop_test2,
  FDR_threshold = 0.05,
  log2FD_threshold = log2(1.5),
  order_clusters = TRUE
) {
  
  ## Retrieve results.
  prop_test1@results$permutation$trem2="R47H"
  prop_test2@results$permutation$trem2="R62H"
  names(prop_test1@results$permutation)[3]="VAR"
  names(prop_test2@results$permutation)[3]="VAR"
  
  
  
  plot_data <- rbind(  prop_test1@results$permutation,
                       prop_test2@results$permutation)
  
  ## Mark the significant results.
  plot_data[, significance := ifelse(
    FDR < FDR_threshold & abs(obs_log2FD) > log2FD_threshold,
    paste("FDR <", FDR_threshold, "& \nabs(Log2FD) >", round(log2FD_threshold, 2)),
    "n.s."
  )]
  
  plot_data[, significance := factor(significance, levels = c(
    paste("FDR <", FDR_threshold, "& \nabs(Log2FD) >", round(log2FD_threshold, 2)),
    "n.s."
  ))]
  
  ## Order the clusters by observed log2FD if requested.
  if (order_clusters) {
    plot_data[, clusters := fct_reorder(factor(clusters), desc(obs_log2FD))]
  }
  
  plot_data$VAR
  ## Plot the results.
  p <- ggplot(plot_data, aes(x = clusters, y = obs_log2FD)) +
    geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = significance), size=1.5) +
    theme_bw() +
    geom_hline(yintercept = log2FD_threshold, lty = 2) +
    geom_hline(yintercept = -log2FD_threshold, lty = 2) +
    geom_hline(yintercept = 0) +
    scale_color_manual(values = c("salmon", "darkgrey")) +
    coord_flip()+facet_wrap(~trem2)+
    labs(y = "Observed Log2FD") + 
    theme(  # Increase the font size of axis labels
      axis.text = element_text(size = 22),
      axis.title.x = element_text(size = 22),  # Increase the font size of x-axis label
      axis.title.y = element_blank(),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      strip.text = element_text(size = 20)  # Increase the font size of facet labels
    ) 
  
  return(p)
}

permutation_plot2(prop_test1,
                  prop_test2)



# Permutation test of clusters vs diagnosis
prop_test_diag <- permutation_test(prop_test, cluster_identity = "nn_clusters_corrected14",
                                   sample_1 = "Control", sample_2 = "AD",
                                   sample_identity = "diagnosis")

permutation_plot(prop_test_diag)


# Permutation test of distancegroup vs neuronal
prop_test_diag <- permutation_test(prop_test, cluster_identity = "distance_group",
                                   sample_1 = "Non-neuronal", sample_2 = "Neuronal",
                                   sample_identity = "neuronal")

permutation_plot(prop_test_diag)
##################################################################################################
library(lme4)
MASC <- function(dataset, cluster, contrast, random_effects = NULL, fixed_effects = NULL,
                 verbose = FALSE, save_models = FALSE, save_model_dir = NULL) {
  # Check inputs
  if (is.factor(dataset[[contrast]]) == FALSE) {
    stop("Specified contrast term is not coded as a factor in dataset")
  }
  
  # Generate design matrix from cluster assignments
  cluster <- as.character(cluster)
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)
  
  # Convert cluster assignments to string
  cluster <- as.character(cluster)
  # Prepend design matrix generated from cluster assignments
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)
  # Create output list to hold results
  res <- vector(mode = "list", length = length(unique(cluster)))
  names(res) <- attributes(designmat)$dimnames[[2]]
  
  # Create model formulas
  if (!is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(c(paste0(fixed_effects, collapse = " + "),
                          paste0("(1|", random_effects, ")", collapse = " + ")),
                        collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else if (!is.null(fixed_effects) && is.null(random_effects)) {
    model_rhs <- paste0(fixed_effects, collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
      # For now, do not allow models without mixed effects terms
      stop("No random effects specified")
    }
  } else if (is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0("(1|", random_effects, ")", collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else {
    model_rhs <- "1" # only includes intercept
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
      stop("No random or fixed effects specified")
    }
  }
  
  # Initialize list to store model objects for each cluster
  cluster_models <- vector(mode = "list",
                           length = length(attributes(designmat)$dimnames[[2]]))
  names(cluster_models) <- attributes(designmat)$dimnames[[2]]
  
  # Run nested mixed-effects models for each cluster
  for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
    test_cluster <- attributes(designmat)$dimnames[[2]][i]
    if (verbose == TRUE) {
      message(paste("Creating logistic mixed models for", test_cluster))
    }
    null_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ 1 + "),
                                   model_rhs), collapse = ""))
    full_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ ", contrast, " + "),
                                   model_rhs), collapse = ""))
    # Run null and full mixed-effects models
    null_model <- lme4::glmer(formula = null_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"))
    full_model <- lme4::glmer(formula = full_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"))
    model_lrt <- anova(null_model, full_model)
    # calculate confidence intervals for contrast term beta
    contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
    contrast_ci <- confint.merMod(full_model, method = "Wald",
                                  parm = contrast_lvl2)
    # Save model objects to list
    cluster_models[[i]]$null_model <- null_model
    cluster_models[[i]]$full_model <- full_model
    cluster_models[[i]]$model_lrt <- model_lrt
    cluster_models[[i]]$confint <- contrast_ci
  }
  
  # Organize results into output dataframe
  output <- data.frame(cluster = attributes(designmat)$dimnames[[2]],
                       size = colSums(designmat))
  output$model.pvalue <- sapply(cluster_models, function(x) x$model_lrt[["Pr(>Chisq)"]][2])
  output[[paste(contrast_lvl2, "OR", sep = ".")]] <- sapply(cluster_models, function(x) exp(fixef(x$full)[[contrast_lvl2]]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "2.5 %"]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "97.5 %"]))
  
  # Return MASC results and save models if specified
  if (save_models == TRUE) {
    saveModelObj(cluster_models, save_dir = save_model_dir)
    return(output)
  } else {
    return(output)
  }
}
## MASC
### cluster ~ diagnosis
seu=seurat_obj
test.df <- seu@meta.data %>% as.data.frame() %>%
  mutate(diagnosis = factor(diagnosis, levels = c("Control", "AD")),
         age = as.vector(scale(age)))
test_res <- MASC(data = test.df,
                 cluster = test.df$nn_clusters_corrected14,
                 contrast = "diagnosis",
                 random_effects = "sample_id",
                 fixed_effects = c("sex", "age", "PMD")
)
test_res$cluster_num= as.factor(as.numeric(gsub("cluster", "", test_res$cluster)))

test_res %>% 
  mutate(label = case_when(model.pvalue <= 0.001 ~ "***",
                           model.pvalue <= 0.01 ~ "**",
                           model.pvalue <= 0.05 ~ "*",
                           model.pvalue > 0.05 ~ "ns")) %>%
  ggplot(aes(x = cluster_num, y = diagnosisAD.OR, 
             ymin = diagnosisAD.OR.95pct.ci.lower,
             ymax = diagnosisAD.OR.95pct.ci.upper )) +
  geom_pointrange(aes(col = cluster_num), size = 1.5) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_errorbar(aes(ymin = diagnosisAD.OR.95pct.ci.lower,
                    ymax = diagnosisAD.OR.95pct.ci.upper, 
                    col = cluster_num),
                linewidth = 0.5, cex = 1) + 
  geom_text(aes(x = cluster_num, y = diagnosisAD.OR.95pct.ci.upper + 0.3, 
                label = label), size = 5, angle = 90) +
  coord_flip() +
  xlab("") + 
  ylab("Odds Ratio (95% CI)") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 18),
    axis.ticks.y = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid = element_blank(), 
    legend.position = "none"
  ) 

# Add a distance_group column in seurat_obj. 
# seurat_obj$distance_group=ifelse(seurat_obj$distance==10,"<=10",
#                                  ifelse(seurat_obj$distance==20 |seurat_obj$distance==30,"10<cell<=30",
#                                         ifelse(seurat_obj$distance==40 |seurat_obj$distance==50,"30<cell<=50",
#                                                ifelse(seurat_obj$distance==1,'plaque','>50'))))
# seurat_obj$distance_group_num=ifelse(seurat_obj$distance==10,"2",
#                                      ifelse(seurat_obj$distance==20 |seurat_obj$distance==30,"3",
#                                             ifelse(seurat_obj$distance==40 |seurat_obj$distance==50,"4",
#                                                    ifelse(seurat_obj$distance==1,'1','5'))))
## MASC
### cluster ~ trem2all
test.df <- seu@meta.data %>% as.data.frame() %>%
  mutate(trem2_all = factor(trem2_all, levels = c("CV", "TREM2var")),
         age = as.vector(scale(age)))
test_res <- MASC(data = test.df,
                 cluster = test.df$nn_clusters_corrected14,
                 contrast = "trem2_all",
                 random_effects = "sample_id",
                 fixed_effects = c("sex", "age", "PMD")
)
test_res$cluster_num= as.factor(as.numeric(gsub("cluster", "", test_res$cluster)))

test_res %>% 
  mutate(label = case_when(model.pvalue <= 0.001 ~ "***",
                           model.pvalue <= 0.01 ~ "**",
                           model.pvalue <= 0.05 ~ "*",
                           model.pvalue > 0.05 ~ "ns")) %>%
  ggplot(aes(x = cluster_num, y = trem2_allTREM2var.OR, 
             ymin = trem2_allTREM2var.OR.95pct.ci.lower,
             ymax = trem2_allTREM2var.OR.95pct.ci.upper )) +
  geom_pointrange(aes(col = cluster_num), size = 1.5) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_errorbar(aes(ymin = trem2_allTREM2var.OR.95pct.ci.lower,
                    ymax = trem2_allTREM2var.OR.95pct.ci.upper, 
                    col = cluster_num),
                linewidth = 0.5, cex = 1) + 
  geom_text(aes(x = cluster_num, y = trem2_allTREM2var.OR.95pct.ci.upper + 0.3, 
                label = label), size = 5, angle = 90) +
  coord_flip() +
  xlab("") + 
  ylab("Odds Ratio (95% CI)") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 18),
    axis.ticks.y = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid = element_blank(), 
    legend.position = "none"
  )

## MASC
### distance ~ diagnosis
test.df <- seu@meta.data %>% as.data.frame() %>%
  mutate(diagnosis = factor(diagnosis, levels = c("Control", "AD")),
         age = as.vector(scale(age)))
test_res <- MASC(data = test.df,
                 cluster = test.df$nn_clusters_corrected14,
                 contrast = "diagnosis",
                 random_effects = "sample_id",
                 fixed_effects = c("sex", "age", "PMD")
)
test_res$cluster_num= as.factor(as.numeric(gsub("cluster", "", test_res$cluster)))

test_res %>% 
  mutate(label = case_when(model.pvalue <= 0.001 ~ "***",
                           model.pvalue <= 0.01 ~ "**",
                           model.pvalue <= 0.05 ~ "*",
                           model.pvalue > 0.05 ~ "ns")) %>%
  ggplot(aes(x = cluster_num, y = diagnosisAD.OR, 
             ymin = diagnosisAD.OR.95pct.ci.lower,
             ymax = diagnosisAD.OR.95pct.ci.upper )) +
  geom_pointrange(aes(col = cluster_num), size = 1.5) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_errorbar(aes(ymin = diagnosisAD.OR.95pct.ci.lower,
                    ymax = diagnosisAD.OR.95pct.ci.upper, 
                    col = cluster_num),
                linewidth = 0.5, cex = 1) + 
  geom_text(aes(x = cluster_num, y = diagnosisAD.OR.95pct.ci.upper + 0.3, 
                label = label), size = 5, angle = 90) +
  coord_flip() +
  xlab("") + 
  ylab("Odds Ratio (95% CI)") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 18),
    axis.ticks.y = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid = element_blank(), 
    legend.position = "none"
  )
### MASC
### cluster ~ BraakGroup
seurat_obj_AD <- seurat_obj[, seurat_obj@meta.data$diagnosis == "AD"]
test.df <- seurat_obj_AD@meta.data %>% as.data.frame() %>%
  mutate(BraakGroup = factor(BraakGroup, levels = c("Braak 3-4", "Braak 5-6")),
         age = as.vector(scale(age)))
test_res <- MASC(data = test.df,
                 cluster = test.df$nn_clusters_corrected14,
                 contrast = "BraakGroup",
                 random_effects = "sample_id",
                 fixed_effects = c("sex", "age", "PMD")
)
test_res$cluster_num= as.factor(as.numeric(gsub("cluster", "", test_res$cluster)))

test_res %>% 
  mutate(label = case_when(model.pvalue <= 0.001 ~ "***",
                           model.pvalue <= 0.01 ~ "**",
                           model.pvalue <= 0.05 ~ "*",
                           model.pvalue > 0.05 ~ "ns")) %>%
  ggplot(aes(x = cluster_num, y = `BraakGroupBraak 5-6.OR`, 
             ymin = `BraakGroupBraak 5-6.OR.95pct.ci.lower`,
             ymax = `BraakGroupBraak 5-6.OR.95pct.ci.upper`)) +
  geom_pointrange(aes(col = cluster_num), size = 1.5) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_errorbar(aes(ymin = `BraakGroupBraak 5-6.OR.95pct.ci.lower`,
                    ymax = `BraakGroupBraak 5-6.OR.95pct.ci.upper`, 
                    col = cluster_num),
                linewidth = 0.5, cex = 1) + 
  geom_text(aes(x = cluster_num, y = `BraakGroupBraak 5-6.OR.95pct.ci.upper` + 0.3, 
                label = label), size = 5, angle = 90) +
  coord_flip() +
  xlab("") + 
  ylab("Odds Ratio (95% CI)") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 18),
    axis.ticks.y = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid = element_blank(), 
    legend.position = "none"
  )

### distance ~ BraakGroup
seurat_obj_AD$distance <- ifelse(seurat_obj_AD$distance == "larger than 100", "200", seurat_obj_AD$distance)
test.df <- seurat_obj_AD@meta.data %>% as.data.frame() %>%
  mutate(BraakGroup = factor(BraakGroup, levels = c("Braak 3-4", "Braak 5-6")),
         age = as.vector(scale(age)))
test_res <- MASC(data = test.df,
                 cluster = test.df$distance,
                 contrast = "BraakGroup",
                 random_effects = "sample_id",
                 fixed_effects = c("sex", "age", "PMD")
)
test_res$cluster_num= as.factor(as.numeric(gsub("cluster", "", test_res$cluster)))

test_res %>% 
  mutate(label = case_when(model.pvalue <= 0.001 ~ "***",
                           model.pvalue <= 0.01 ~ "**",
                           model.pvalue <= 0.05 ~ "*",
                           model.pvalue > 0.05 ~ "ns")) %>%
  ggplot(aes(x = cluster_num, y = `BraakGroupBraak 5-6.OR`, 
             ymin = `BraakGroupBraak 5-6.OR.95pct.ci.lower`,
             ymax = `BraakGroupBraak 5-6.OR.95pct.ci.upper`)) +
  geom_pointrange(aes(col = cluster_num), size = 1.5) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_errorbar(aes(ymin = `BraakGroupBraak 5-6.OR.95pct.ci.lower`,
                    ymax = `BraakGroupBraak 5-6.OR.95pct.ci.upper`, 
                    col = cluster_num),
                linewidth = 0.5, cex = 1) + 
  geom_text(aes(x = cluster_num, y = `BraakGroupBraak 5-6.OR.95pct.ci.upper` + 0.3, 
                label = label), size = 5, angle = 90) +
  coord_flip() +
  xlab("") + 
  ylab("Odds Ratio (95% CI)") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 18),
    axis.ticks.y = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid = element_blank(), 
    legend.position = "none"
  )

# distance_group~diagnosis
test.df <- seu@meta.data %>% as.data.frame() %>%
  mutate(diagnosis = factor(diagnosis, levels = c("Control", "AD")),
         age = as.vector(scale(age)))
test_res <- MASC(data = test.df,
                 cluster = test.df$distance_group,
                 contrast = "diagnosis",
                 random_effects = "sample_id",
                 fixed_effects = c("sex", "age", "PMD")
)
test_res$cluster_num= as.factor(as.numeric(gsub("cluster", "", test_res$cluster)))

test_res %>% 
  mutate(label = case_when(model.pvalue <= 0.001 ~ "***",
                           model.pvalue <= 0.01 ~ "**",
                           model.pvalue <= 0.05 ~ "*",
                           model.pvalue > 0.05 ~ "ns")) %>%
  ggplot(aes(x = cluster_num, y = diagnosisAD.OR, 
             ymin = diagnosisAD.OR.95pct.ci.lower,
             ymax = diagnosisAD.OR.95pct.ci.upper )) +
  geom_pointrange(aes(col = cluster_num), size = 1.5) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_errorbar(aes(ymin = diagnosisAD.OR.95pct.ci.lower,
                    ymax = diagnosisAD.OR.95pct.ci.upper, 
                    col = cluster_num),
                linewidth = 0.5, cex = 1) + 
  geom_text(aes(x = cluster_num, y = diagnosisAD.OR.95pct.ci.upper + 0.3, 
                label = label), size = 5, angle = 90) +
  coord_flip() +
  xlab("") + 
  ylab("Odds Ratio (95% CI)") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 18),
    axis.ticks.y = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid = element_blank(), 
    legend.position = "none"
  ) 

#### NEURONAL
seurat_obj$neuronal <- ifelse(seurat_obj$nn_clusters_corrected14 %in% c(3, 8, 9, 15),
                              "Neuronal", "Non-neuronal")
## i change the clusters 
seurat_obj$neuronal <- ifelse(seurat_obj$nn_clusters_corrected14 %in% c(4, 8, 9, 11, 15),
                              "Neuronal", "Non-neuronal")
summary_table <- table(seurat_obj$trem2, seurat_obj$neuronal)


# Calculate proportions within each trem2 group
cv_proportions <- summary_table["CV", ] / sum(summary_table["CV", ])
r47h_proportions <- summary_table["R47H", ] / sum(summary_table["R47H", ])
r62h_proportions <- summary_table["R62H", ] / sum(summary_table["R62H", ])

# Combine proportions into a single data frame
plot_data <- data.frame(trem2 = rep(row.names(summary_table), each = 2),
                        neuronal = rep(c("Neuronal", "Non-neuronal"), times = 3),
                        count = c(summary_table["CV", "Neuronal"], summary_table["CV", "Non-neuronal"],
                                  summary_table["R47H", "Neuronal"], summary_table["R47H", "Non-neuronal"],
                                  summary_table["R62H", "Neuronal"], summary_table["R62H", "Non-neuronal"]),
                        proportions = c(cv_proportions[1], cv_proportions[2],
                                        r47h_proportions[1], r47h_proportions[2],
                                        r62h_proportions[1], r62h_proportions[2]))


# Create the stacked bar plot
ggplot(plot_data, aes(x = trem2, y = proportions, fill = neuronal)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "TREM2", y = "Proportion of cells", fill = "Cell Type") +
  ggtitle("Proportion of Neuronal and Non-neuronal Cells by TREM2")

## NEURONAL - AD ONLY
seurat_obj_AD <- seurat_obj[, seurat_obj@meta.data$diagnosis == "AD"]

summary_table <- table(seurat_obj_AD$trem2, seurat_obj_AD$neuronal)

# Calculate proportions within each trem2 group
cv_proportions <- summary_table["CV", ] / sum(summary_table["CV", ])
r47h_proportions <- summary_table["R47H", ] / sum(summary_table["R47H", ])
r62h_proportions <- summary_table["R62H", ] / sum(summary_table["R62H", ])


# Combine proportions into a single data frame
plot_data <- data.frame(trem2 = rep(row.names(summary_table), each = 2),
                        neuronal = rep(c("Neuronal", "Non-neuronal"), times = 3),
                        count = c(summary_table["CV", "Neuronal"], summary_table["CV", "Non-neuronal"],
                                  summary_table["R47H", "Neuronal"], summary_table["R47H", "Non-neuronal"],
                                  summary_table["R62H", "Neuronal"], summary_table["R62H", "Non-neuronal"]),
                        proportions = c(cv_proportions[1], cv_proportions[2],
                                        r47h_proportions[1], r47h_proportions[2],
                                        r62h_proportions[1], r62h_proportions[2]))


# Create the stacked bar plot
ggplot(plot_data, aes(x = trem2, y = proportions, fill = neuronal)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "TREM2", y = "Proportion of cells", fill = "Cell Type") +
  ggtitle("Proportion of Neuronal and Non-neuronal Cells by TREM2 (AD Only)")



###############################################################################
seu <- seurat_obj
# Convert Seurat object metadata to a data frame
test.df <- as.data.frame(seu@meta.data)

# Make necessary transformations to the data frame
test.df <- test.df %>%
  mutate(neuronal = factor(neuronal, levels = c("Non-neuronal", "Neuronal")),
         age = as.vector(scale(age)))

# Call the MASC function
test_res <- MASC(dataset = test.df,
                 cluster = seu$distance_group_num,
                 contrast = "neuronal",
                 random_effects = "sample_id",
                 fixed_effects = c("sex", "age", "PMD"),
                 verbose = TRUE)
test_res$cluster_num= as.factor(as.numeric(gsub("cluster", "", test_res$cluster)))


# Create a named vector to map numeric values to corresponding labels
labels_map <- c("1" = "plaque",
                "2" = "<=10um",
                "3" = "10um<cell<=30um",
                "4" = "30um<cell<=50um",
                "5" = ">50um")

test_res %>% 
  mutate(label = case_when(model.pvalue <= 0.001 ~ "***",
                           model.pvalue <= 0.01 ~ "**",
                           model.pvalue <= 0.05 ~ "*",
                           model.pvalue > 0.05 ~ "ns")) %>%
  ggplot(aes(x = cluster_num, y = neuronalNeuronal.OR, 
             ymin = neuronalNeuronal.OR.95pct.ci.lower,
             ymax = neuronalNeuronal.OR.95pct.ci.upper )) +
  geom_pointrange(aes(col = cluster_num), size = 1.5) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_errorbar(aes(ymin = neuronalNeuronal.OR.95pct.ci.lower,
                    ymax = neuronalNeuronal.OR.95pct.ci.upper, 
                    col = cluster_num),
                linewidth = 0.5, cex = 1) + 
  geom_text(aes(x = cluster_num, y = neuronalNeuronal.OR.95pct.ci.upper + 0.3, 
                label = label), size = 10, angle = 90) +
  coord_flip() +
  xlab("") + 
  ylab("Odds Ratio (95% CI)") +
  # Set the x-axis labels using scale_x_discrete
  scale_x_discrete(labels = labels_map) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 22),
    axis.title = element_text(color = "black", size = 22),
    axis.ticks.y = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid = element_blank(), 
    legend.position = "none"
  ) 

