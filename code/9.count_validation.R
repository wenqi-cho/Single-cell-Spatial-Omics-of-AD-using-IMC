library(stringr)
manual_counts_Iba1<-read.csv("manual_counts_Iba1.csv")
manual_counts_Iba1$sample_id=paste0(manual_counts_Iba1$sample_id, ' ')
manual_counts_Iba1$sample_id=str_replace_all(manual_counts_Iba1$sample_id,
                                             c(" i " =  '_001', " ii " =  '_002', " iii " =  '_003', " iv " =  '_004'))

library(SpatialExperiment)
# # Retrieve the count data for the Iba1 marker
counts <- assays(spe_filtered)$exprs["Iba1",]
# 
# # Filter sample IDs starting with "19910006_002"
# sample_ids_to_sum <- grep("^962.2_004", names(counts), value = TRUE)
# 
# # Subset counts for the selected sample IDs
# counts_to_sum <- counts[sample_ids_to_sum]
# 
# # Calculate the sum of counts
# total_counts <- sum(counts_to_sum)

library(dplyr)
validation_results <- list()
# Loop through each sample ID in the manual counts data
for (sample_id in manual_counts_Iba1$sample_id) {
  # Filter sample IDs for the current sample_id
  sample_ids_to_sum <- grep(paste0("^", sample_id), names(counts), value = TRUE)
  
  # Subset counts for the selected sample IDs
  counts_to_sum <- counts[sample_ids_to_sum]
  
  # Calculate the sum of counts
  total_counts_spe <- sum(counts_to_sum)
  
  # Find the corresponding row in manual counts data
  manual_row <- manual_counts_Iba1[manual_counts_Iba1$sample_id == sample_id, ]
  
  # Get the DNA count from manual counts data
  dna_count_manual <- manual_row$DNA
  
  # Compare the counts
  validation_results[[sample_id]] <- data.frame(
    sample_id = sample_id,
    dna_count_spe = total_counts_spe,
    dna_count_manual = dna_count_manual
  )
}

# Combine validation results into a data frame
validation_df <- do.call(rbind, validation_results)

# Print the validation results
print(validation_df)

# Prepare a data frame for comparison
comparison_df <- data.frame(
  sample_id = manual_counts_Iba1$sample_id,
  manual_count = manual_counts_Iba1$DNA,
  spe_count = sapply(manual_counts_Iba1$sample_id, function(sample_id) {
    sample_ids_to_sum <- grep(paste0("^", sample_id), names(counts), value = TRUE)
    counts_to_sum <- counts[sample_ids_to_sum]
    sum(counts_to_sum)
  })
)
comparison_df$difference=comparison_df$manual_count-comparison_df$spe_count

# Perform Wilcoxon signed-rank test
wilcoxon_test <- wilcox.test(comparison_df$manual_count, comparison_df$spe_count, paired = TRUE)

# Print the test results
print(wilcoxon_test)

# Create a scatter plot with reference line
ggplot(comparison_df, aes(x = manual_count, y = spe_count)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Manual Counts", y = "SPE Counts") +
  ggtitle("Comparison of Manual and SPE DNA Counts") +
  theme_minimal()

# Prepare a data frame for paired t-test
paired_data <- data.frame(
  manual_count = manual_counts_Iba1$DNA,
  spe_count = sapply(manual_counts_Iba1$sample_id, function(sample_id) {
    sample_ids_to_sum <- grep(paste0("^", sample_id), names(counts), value = TRUE)
    counts_to_sum <- counts[sample_ids_to_sum]
    sum(counts_to_sum)
  })
)

# Perform paired t-test
paired_ttest <- t.test(paired_data$manual_count, paired_data$spe_count, paired = TRUE)

# Print the t-test results
print(paired_ttest)


## Blandr
# Perform Bland-Altman analysis
blandr.stats <- blandr.statistics(comparison_df$manual_count, comparison_df$spe_count, sig.level = 0.95)
mean_difference <- mean(comparison_df$difference)
limits_of_agreement <- c(blandr.stats$lowerLOA, blandr.stats$upperLOA)


# Print the Bland-Altman statistics
print(blandr.stats)

# Create Bland-Altman plot
# blandr.plot <- blandr.draw(comparison_df$manual_count, comparison_df$spe_count, ciDisplay = FALSE, ciShading = FALSE) +
#   geom_abline(intercept = mean_difference, slope = 0, linetype = "dashed") +
#   ggtitle("Bland-Altman Plot: Manual vs. Steinbock DNA Counts") +
#   annotate("text", x = 150, y = 50, label = paste("Mean Difference:", round(mean_difference, 2)), color = "red")


# Load the irr package
library(irr)

# Calculate the ICC
icc_result <- icc(comparison_df[, c("manual_count","spe_count")], model = "twoway", type = "agreement")

# Extract the ICC value
icc_value <- icc_result$value[1]
mean_difference <- mean(comparison_df$difference)
# Calculate Spearman's rank correlation
spearman_corr <- cor.test(comparison_df$manual_count, comparison_df$spe_count, method = "spearman")
correlation_coefficient <- spearman_corr$estimate

# Calculate Pearson's correlation coefficient
pearson_corr <- cor( comparison_df$manual_count, comparison_df$spe_count,method = "pearson")

# Create the Bland-Altman plot
blandr.plot <- blandr.draw(comparison_df$manual_count, comparison_df$spe_count, ciDisplay = FALSE, ciShading = FALSE) +
  ggtitle("Bland-Altman Plot: Steinbock vs. Manual DNA Counts") +
  annotate("text", x = 150, y = 50, label = paste("Mean Difference:", round(mean_difference, 2)), color = "red") +
  annotate("text", x = Inf, y = Inf, label = paste0("ICC =", round(icc_value, 2)),
           hjust = 1, vjust = 11, size = 3, color = "blue") +
  annotate("text", x = Inf, y = Inf, label = paste0("Spearman's rho =", round(correlation_coefficient, 2)),
           hjust = 1, vjust = 13, size = 3, color = "blue") +
  annotate("text", x = Inf, y = Inf, label = paste0("Pearson's corr.=", round(pearson_corr, 2)),
           hjust = 1, vjust = 15, size = 3, color = "blue") +
  annotate("text", x = Inf, y = -Inf, label = paste0("Limits of Agreement = [", round(limits_of_agreement[1], 2), ", ", round(limits_of_agreement[2], 2), "]"),
           hjust = 1, vjust = 0, size = 3, color = "black") +
  geom_hline(yintercept = mean_difference, linetype = "dashed", color = "red") +
  #geom_hline(yintercept = limits_of_agreement[1], linetype = "dotted", color = "purple") +
  #geom_hline(yintercept = limits_of_agreement[2], linetype = "dotted", color = "purple") +
  theme_minimal()

# Display the Bland-Altman plot with ICC value
print(blandr.plot)
# Perform paired t-test
ttest_result <- t.test(df_comparison$num_cells, df_comparison$Count, paired = TRUE)


###############################################################################
library(readr)
comarker_full <- read_csv("CoMarker_Colocalisation_analysis_report2.csv")

spe_filtered_AD_counts= spe_filtered_AD[,spe_filtered_AD$area>=30 ]
automatic_counts=data.frame(colData(spe_filtered_AD_counts))

df_aut=automatic_counts%>%
  group_by(sample_id)%>%
  summarise(num_cells = sum(nn_clusters_corrected14 %in% c(1, 2)))

df_aut=data.frame(df_aut)
rownames(df_aut)=df_aut$sample_id
hist(df_aut$num_cells)


comarker_full=comarker_full[which(comarker_full$Slice=='Iba1 Cell'),]
comarker_full$replicate=paste0(comarker_full$replicate, ' ')
comarker_full$replicate=str_replace_all(comarker_full$replicate,c(" i " =  '_001', " ii " =  '_002', " iii " =  '_003', " iv " =  '_004'))

df_com=comarker_full[,c("replicate","Count")]
df_com$replicate[df_com$replicate %in% automatic_counts$sample_id]
df_com=data.frame(df_com)
rownames(df_com)=df_com$replicate


df_merge=merge(df_com,df_aut,by="row.names")
df_comparison=df_merge[,c('replicate','Count','num_cells')]
# gg_df=melt(df_comparison)
# gg_df$replicate

# 
# ggplot(data = gg_df, aes(x=replicate, y = value, fill = variable)) +
#   geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75)+
#   coord_flip()

df_comparison$difference=df_comparison$num_cells-df_comparison$Count
df_comparison$perc_change=df_comparison$num_cells/df_comparison$Count
df_comparison$id=gsub("\\_.*","",df_comparison$replicate)
library(forcats)
df_comparison %>%
  filter(!is.na(difference)) %>%
  mutate(replicate = fct_reorder(replicate, desc(perc_change))) %>%
  ggplot(., aes(x = replicate, y = perc_change, fill = factor(id))) +
  geom_bar(stat = "identity", show.legend = F)+ 
  geom_text(aes(label=round(perc_change,2)), position=position_dodge(width=0.9), vjust=-0.25,size=1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0,size = 5))+ylab("Ratio steinbock/comarker")


df_comparison$replicate
df_comparison <- na.omit(df_comparison)

library(blandr)

#blandr.statistics(df_comparison$num_cells, df_comparison$Count, sig.level=0.95 )
blandr_result <- blandr.statistics( df_comparison$num_cells,df_comparison$Count, sig.level = 0.95)
mean_difference <- mean(df_comparison$difference)
limits_of_agreement <- c(blandr_result$lowerLOA, blandr_result$upperLOA)
# Load the irr package
library(irr)

# Calculate the ICC
icc_result <- icc(df_comparison[, c("num_cells","Count")], model = "twoway", type = "agreement")

# Extract the ICC value
icc_value <- icc_result$value[1]
mean_difference <- mean(df_comparison$difference)
# Calculate Spearman's rank correlation
spearman_corr <- cor.test(df_comparison$num_cells, df_comparison$Count, method = "spearman")
correlation_coefficient <- spearman_corr$estimate

# Calculate Pearson's correlation coefficient
pearson_corr <- cor( df_comparison$num_cells, df_comparison$Count,method = "pearson")

# Create the Bland-Altman plot
blandr.plot <- blandr.draw(df_comparison$num_cells,df_comparison$Count, ciDisplay = FALSE, ciShading = FALSE) +
  ggtitle("Bland-Altman Plot: Steinbock vs. Comarker DNA Counts") +
  annotate("text", x = 100, y = 50, label = paste("Mean Difference:", round(mean_difference, 2)), color = "red") +
  annotate("text", x = Inf, y = Inf, label = paste0("ICC =", round(icc_value, 2)),
           hjust = 1, vjust = 11, size = 8, color = "blue") +
  annotate("text", x = Inf, y = Inf, label = paste0("Spearman's rho =", round(correlation_coefficient, 2)),
           hjust = 1, vjust = 13, size = 8, color = "blue") +
  annotate("text", x = Inf, y = Inf, label = paste0("Pearson's corr.=", round(pearson_corr, 2)),
           hjust = 1, vjust = 15, size = 8, color = "blue") +
  annotate("text", x = Inf, y = -Inf, label = paste0("Limits of Agreement = [", round(limits_of_agreement[1], 2), ", ", round(limits_of_agreement[2], 2), "]"),
           hjust = 1, vjust = 0, size = 8, color = "black") +
  geom_hline(yintercept = mean_difference, linetype = "dashed", color = "red") +
  #geom_hline(yintercept = limits_of_agreement[1], linetype = "dotted", color = "purple") +
  #geom_hline(yintercept = limits_of_agreement[2], linetype = "dotted", color = "purple") +
  theme_minimal()

# Display the Bland-Altman plot with ICC value
print(blandr.plot)
# Perform paired t-test
ttest_result <- t.test(df_comparison$num_cells, df_comparison$Count, paired = TRUE)



# pdf(file = "Iba1_signal.pdf",    # The directory you want to save the file in
#     width = 4, # The width of the plot in inches
#     height = 4) # The height of the plot in inches
# 
# for ( i in df_comparison$replicate){
#   cur_images <- saved_images[names(saved_images) %in%  i]
#   cur_masks <- masks[names(masks) %in% i]
#   
#   plot=plotPixels(image = cur_images,
#                   object = spe_filtered_AD, 
#                   cell_id = "ObjectNumber", 
#                   img_id = "sample_id",
#                   colour_by = "Iba1",
#                   bcg = list(`Iba1` = c(0, 20, 1)),
#                   mask = cur_masks,
#                   outline_by = "Iba1_cell",
#                   colour=list(Iba1_cell = c(`Iba1+`= "red",
#                                             other="purple")),
#                   thick = TRUE)
#   print(plot)
# }
# dev.off()
## Outliers
# Create box plots
boxplot(df_comparison$num_cells, df_comparison$Count, names = c("Steinbock", "Comarker"))
# Calculate Z-scores
z_scores_steinbock <- (df_comparison$num_cells - mean(df_comparison$num_cells)) / sd(df_comparison$num_cells)
z_scores_comarker <- (df_comparison$Count - mean(df_comparison$Count)) / sd(df_comparison$Count)
# Create scatter plot
plot(df_comparison$num_cells, df_comparison$Count, xlab = "Steinbock", ylab = "Comarker")
# Create Q-Q plot
qqnorm(df_comparison$num_cells - df_comparison$Count)
qqline(df_comparison$num_cells - df_comparison$Count)
# Perform linear regression
linear_model <- lm(num_cells ~ Count, data = df_comparison)
# Steinbock = 11.341 + 1.006*Count 1.006 suggest a + r/s. for each unit increase in the count var, num_cells increases by1.006
# Create a scatter plot with regression line
ggplot(df_comparison, aes(x = Count, y = num_cells)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Linear Regression: Steinbock vs. Comarker",
       x = "Comarker Counts", y = "Steinbock Counts") +
  theme_minimal()
################################################################
## MANUAL
spe_filtered_AD_counts= spe_filtered_AD[,spe_filtered_AD$area>=30 ]
automatic_counts=data.frame(colData(spe_filtered_AD_counts))

df_aut=automatic_counts%>%
  group_by(sample_id)%>%
  summarise(num_cells = sum(nn_clusters_corrected14 %in% c(2)))

df_aut=data.frame(df_aut)
rownames(df_aut)=df_aut$sample_id
hist(df_aut$num_cells)

manual_count <- read.csv("manualCount_mtg_metadata2.csv")
manual_count <- manual_count[,1:5]
manual_count <- manual_count[manual_count$CaseID_imageID != "replicate", ]

manual_count$sample_id=paste0(manual_count$CaseID_imageID, ' ')
manual_count$sample_id=str_replace_all(manual_count$sample_id,
                                             c(" i " =  '_001', " ii " =  '_002', " iii " =  '_003', " iv " =  '_004'))

df_manual=manual_count[,c("sample_id","total.microglia")]
# df_manual$sample_id[df_manual$sample_id %in% automatic_counts$sample_id]
df_manual = df_manual[df_manual$sample_id %in% automatic_counts$sample_id, ]
df_manual=data.frame(df_manual)
rownames(df_manual)=df_manual$sample_id


df_merge=merge(df_manual,df_aut,by="row.names")
df_comparison=df_merge[,c('sample_id.x','total.microglia','num_cells')]
colnames(df_comparison) <- c('sample_id', 'manual_count', 'auto_count')

#scatter plot
ggplot(df_comparison, aes(x = manual_count, y = auto_count)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Manual (N)", y = "Steinbock (N)") +
  ggtitle("Comparison of Microglia - Steinbock vs Manual") +
  theme_minimal()
df_comparison$difference=df_comparison$auto_count-df_comparison$manual_count
##
#blandr.statistics(df_comparison$num_cells, df_comparison$Count, sig.level=0.95 )
blandr_result <- blandr.statistics( df_comparison$manual_count,df_comparison$auto_count, sig.level = 0.95)
mean_difference <- mean(df_comparison$difference)
limits_of_agreement <- c(blandr_result$lowerLOA, blandr_result$upperLOA)
# Load the irr package
library(irr)

# Calculate the ICC
icc_result <- icc(df_comparison[, c("auto_count","manual_count")], model = "twoway", type = "agreement")

# Extract the ICC value
icc_value <- icc_result$value[1]
mean_difference <- mean(df_comparison$difference)
# Calculate Spearman's rank correlation
spearman_corr <- cor.test(df_comparison$manual_count,df_comparison$auto_count, method = "spearman")
correlation_coefficient <- spearman_corr$estimate

# Calculate Pearson's correlation coefficient
pearson_corr <- cor(df_comparison$manual_count,df_comparison$auto_count,method = "pearson")

# Create the Bland-Altman plot
blandr.plot <- blandr.draw(df_comparison$auto_count,df_comparison$manual_count, ciDisplay = FALSE, ciShading = FALSE) +
  ggtitle("Bland-Altman Plot: Steinbock vs. Manual Microglia Count") +
  annotate("text", x = 10, y = 40, label = paste("Mean Difference:", round(mean_difference, 2)), size = 4.5, color = "red") +
  annotate("text", x = Inf, y = Inf, label = paste0("ICC =", round(icc_value, 2)),
           hjust = 1, vjust = 11, size = 4, color = "blue") +
  annotate("text", x = Inf, y = Inf, label = paste0("Spearman's rho =", round(correlation_coefficient, 2)),
           hjust = 1, vjust = 13, size = 4, color = "blue") +
  annotate("text", x = Inf, y = Inf, label = paste0("Pearson's corr.=", round(pearson_corr, 2)),
           hjust = 1, vjust = 15, size = 4, color = "blue") +
  annotate("text", x = Inf, y = -Inf, label = paste0("Limits of Agreement = [", round(limits_of_agreement[1], 2), ", ", round(limits_of_agreement[2], 2), "]"),
           hjust = 1, vjust = 0, size = 4, color = "black") +
  geom_hline(yintercept = mean_difference, linetype = "dashed", color = "red") +
  #geom_hline(yintercept = limits_of_agreement[1], linetype = "dotted", color = "purple") +
  #geom_hline(yintercept = limits_of_agreement[2], linetype = "dotted", color = "purple") +
  theme(axis.text = element_text(size = 4),axis.text.x = element_text(size = 12),  # Adjust x-axis tick label size
        axis.text.y = element_text(size = 12))

# Display the Bland-Altman plot with ICC value
print(blandr.plot)
# Perform paired t-test
ttest_result <- t.test(df_comparison$auto_count,df_comparison$manual_count, paired = TRUE)
# data:  df_comparison$auto_count and df_comparison$manual_count
# t = 9.6386, df = 59, p-value = 1.002e-13
# alternative hypothesis: true mean difference is not equal to 0
# 95 percent confidence interval:
#   15.17443 23.12557
################################################################
