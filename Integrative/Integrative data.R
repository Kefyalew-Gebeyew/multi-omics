# Load necessary libraries
library(mixOmics)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pROC)

# Set random seed
set.seed(12345)

## Improved preprocessing function with stricter filtering
preprocess_and_filter <- function(data, name, is_microbiome = FALSE) {
  cat("\nProcessing", name, "\n")
  
  if(!is_microbiome) {
    # For FPKM data
    # 1. More stringent expression filtering
    mean_expression <- colMeans(data)
    keep <- mean_expression > quantile(mean_expression, 0.75)  # Top 25% expressed genes
    data <- data[, keep]
    
    # 2. More stringent variance filtering
    var_genes <- apply(data, 2, var)
    keep_var <- var_genes > quantile(var_genes, 0.75)  # Top 25% variable genes
    data <- data[, keep_var]
    
    # 3. Log transformation with variance stabilization
    data <- log2(data + 1)
    
    # 4. Remove outliers based on MAD
    mad_scores <- apply(data, 1, mad)
    keep_samples <- mad_scores < (median(mad_scores) + 3*mad(mad_scores))
    data <- data[keep_samples,]
    
  } else {
    # For microbiome data
    # 1. Stricter prevalence filtering
    keep <- colMeans(data > 0) >= 0.75  # Present in 75% of samples
    data <- data[, keep]
    
    # 2. Stricter abundance filtering
    rel_abundance <- sweep(data, 1, rowSums(data), '/')
    abundance_filter <- colMeans(rel_abundance) > 0.1/100  # More than 0.1%
    data <- data[, abundance_filter]
    
    # 3. Remove outliers
    clr_data <- mixOmics::logratio.transfo(data + 1, logratio = 'CLR')
    pca_temp <- prcomp(clr_data)
    scores <- pca_temp$x[,1:2]
    dist_center <- sqrt(rowSums(scores^2))
    keep_samples <- dist_center < (median(dist_center) + 3*mad(dist_center))
    data <- data[keep_samples,]
    
    # 4. Final CLR transformation
    data <- mixOmics::logratio.transfo(data + 1, logratio = 'CLR')
  }
  
  # Scale the data
  data_scaled <- scale(data)
  
  # Print dimensions after filtering
  cat("\nDimensions after filtering:", dim(data_scaled)[1], "x", dim(data_scaled)[2], "\n")
  
  return(data_scaled)
}

# Read and process data
fpkm_data <- read.csv("Colon_FPKM_cleaned.csv", row.names = 1)
genus_abundance <- read.csv("Colon_abundance.csv", row.names = 1)

# Create metadata
metadata <- data.frame(
  Sample = rownames(t(fpkm_data)),
  Group = ifelse(grepl("^E", rownames(t(fpkm_data))), "E", "L")
)

# Process datasets
fpkm_processed <- preprocess_and_filter(t(fpkm_data), "FPKM", FALSE)
genus_processed <- preprocess_and_filter(genus_abundance, "Microbiome", TRUE)

# Update metadata after potential outlier removal
common_samples <- intersect(rownames(fpkm_processed), rownames(genus_processed))
fpkm_processed <- fpkm_processed[common_samples,]
genus_processed <- genus_processed[common_samples,]
metadata <- metadata[metadata$Sample %in% common_samples,]

# Prepare data for DIABLO
X <- list(
  Transcriptome = fpkm_processed,
  Microbiome = genus_processed
)
Y <- metadata$Group

# Modified design matrix with stronger block connections
design <- matrix(c(1, 0.5,
                   0.5, 1), ncol = 2)

# Run DIABLO with modified parameters
diablo_result <- block.splsda(X, Y, 
                              design = design,
                              ncomp = 3,  # Reduced number of components
                              scale = TRUE,
                              near.zero.var = TRUE)

# Visualize results with improved settings
png("DIABLO_Individual_Plot_Integrated_Colon.png", 
    width = 1000, height = 800, 
    res = 120)
plotIndiv(diablo_result, 
          group = Y, 
          legend = TRUE, 
          title = "Integrated DIABLO Analysis - Colon",
          ellipse = TRUE,  # Add confidence ellipses
          style = "graphics",
          conf.level = 0.95)
dev.off()



#####________to create diagnostic plots for multi-block transcriptome and microbiome data analysis_____########


# Create diagnostic plots
# 1. Create the DIABLO diagnostic plot showing correlations between blocks
png("DIABLO_Diagnostic_Plot_Colon.png", 
    width = 1000, height = 1000, 
    res = 120)
plotDiablo(diablo_result, ncomp = 1)
dev.off()


# 2. Evaluate model performance with cross-validation
perf.diablo <- perf(diablo_result, 
                    validation = "Mfold",
                    folds = 5,
                    nrepeat = 10)

# Plot performance
png("DIABLO_Performance_Plot_Colon.png", 
    width = 1000, height = 800, 
    res = 120)
plot(perf.diablo)
dev.off()


#######________________________________________end________________________________###############



#  Feature Selection and Plotting
get_top_features <- function(diablo_result, dataset, is_microbiome = FALSE, top_percentage = NULL, top_n = NULL) {
  loadings <- diablo_result$loadings[[dataset]][, 1]
  importance <- abs(loadings)
  top_features <- data.frame(
    Feature = names(importance),
    Importance = importance,
    Direction = ifelse(loadings > 0, "Positive", "Negative")
  )
  top_features <- top_features %>% arrange(desc(Importance))
  
  if (is_microbiome) {
    top_features <- head(top_features, top_n)
  } else if (!is.null(top_percentage)) {
    n_top <- ceiling(nrow(top_features) * top_percentage)
    top_features <- head(top_features, n_top)
  } else {
    top_features <- head(top_features, top_n)
  }
  
  return(top_features)
}

plot_top_features <- function(top_features, dataset) {
  top_features$Group <- ifelse(top_features$Direction == "Positive", "L", "E")
  p <- ggplot(top_features, aes(x = reorder(Feature, Importance), y = Importance, fill = Group)) +
    geom_bar(stat = "identity", width = 0.1) +
    geom_point(aes(y = Importance, color = Group), size = 3, shape = 16) +
    coord_flip() +
    labs(title = paste("Top Features -", dataset, "- Colon"), x = "Feature", y = "Importance Score") +
    scale_fill_manual(values = c("E" = "blue", "L" = "orange")) +
    scale_color_manual(values = c("E" = "blue", "L" = "orange")) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12))
  
  ggsave(filename = paste0("Top_Features_", dataset, "_Colon.png"), plot = p, width = 6, height = 8, bg = "white", dpi=300)
}

top_features_list <- lapply(names(X), function(dataset) {
  if (dataset == "Microbiome") {
    selected_features <- get_top_features(diablo_result, dataset, is_microbiome = TRUE, top_n = 20)
  } else if (dataset == "Transcriptome") {
    selected_features <- get_top_features(diablo_result, dataset, is_microbiome = FALSE, top_percentage = 0.10)
  }
  plot_top_features(selected_features, dataset)
  return(selected_features)
})
names(top_features_list) <- names(X)

# Save results
saveRDS(diablo_result, "diablo_result_integrated_colon.rds")
write.csv(auc_df, "AUC_Values_Integrated_Colon.csv", row.names = FALSE)

for (dataset in names(top_features_list)) {
  write.csv(top_features_list[[dataset]], paste0("Top_Features_", dataset, "_Colon.csv"), row.names = FALSE)
}




# Generate correlation circle plot
top_genes <- read.csv("Top_Features_Transcriptome_Colon.csv", stringsAsFactors = FALSE)$Feature
top_otus <- read.csv("Top_Features_Microbiome_Colon.csv", stringsAsFactors = FALSE)$Feature

# Subset DIABLO result for each block to include only top features
diablo_result$X$Transcriptome <- diablo_result$X$Transcriptome[, colnames(diablo_result$X$Transcriptome) %in% top_genes]
diablo_result$X$Microbiome <- diablo_result$X$Microbiome[, colnames(diablo_result$X$Microbiome) %in% top_otus]

# Extract data from plotVar() for ggplot customization
var_data <- plotVar(diablo_result, comp = c(1, 2), cutoff = 0.5, var.names = TRUE, legend = TRUE, plot = FALSE)

# Create a data frame from var_data
var_df <- data.frame(
  Comp1 = var_data$x,
  Comp2 = var_data$y,
  Block = var_data$Block,
  Feature = var_data$names
)

# Assign custom shapes and colors for each block
shape_mapping <- c("Microbiome" = 23, "Transcriptome" = 24)
color_mapping <- c("Microbiome" = "blue", "Transcriptome" = "brown2")

# Create customized correlation circle plot using ggplot2
p <- ggplot(var_df, aes(x = Comp1, y = Comp2)) +
  geom_path(data = data.frame(x = cos(seq(0, 2 * pi, length.out = 100)), y = sin(seq(0, 2 * pi, length.out = 100))), aes(x=x, y=y), color="darkgrey", linetype="solid") +
  geom_path(data = data.frame(x = cos(seq(0, 2 * pi, length.out = 100)) * sqrt(0.5), y = sin(seq(0, 2 * pi, length.out = 100)) * sqrt(0.5)), aes(x=x, y=y), color="darkgrey", linetype="solid") +
  geom_hline(yintercept=0, color="darkgrey", linetype="dotted") +
  geom_vline(xintercept=0, color="darkgrey", linetype="dotted") +
  geom_point(aes(shape = Block, color = Block), size = 3, fill=NA) +
  geom_text_repel(aes(label = Feature, color = Block), size=3) +
  scale_shape_manual(values = shape_mapping) +
  scale_color_manual(values = color_mapping) +
  labs(title = "Correlation Circle Plot - Top Features - Colon", x = "Component 1", y = "Component 2") +
  theme_minimal(base_size = 14) +
  theme(panel.border=element_rect(color="black", fill=NA), legend.position="right")

# Save the plot as a PNG file
ggsave(filename = "Improved_Correlation_Circle_Plot_Top_Features_Colon.png", plot = p, width = 7, height = 5.5, bg="white")

# Export correlation results presented on circle plots as CSV file
write.csv(var_df, "Correlation_Results_Top_Features_Colon.csv", row.names=FALSE)



############______________________ROC curve___________________________________#######

# Load required libraries
library(mixOmics)
library(ggplot2)
library(pROC)
library(gridExtra)

# Function to create ROC curve plot using ggplot2
create_roc_plot <- function(diablo_result, dataset_name) {
  # Get scores and class labels for component 1 only
  scores <- diablo_result$variates[[dataset_name]][,1]
  classes <- factor(diablo_result$Y)
  
  # Calculate ROC curve
  roc_obj <- roc(classes, scores, levels = c("E", "L"))
  
  # Create data frame for plotting
  roc_data <- data.frame(
    Specificity = (1 - roc_obj$specificities) * 100,  # Convert to percentage
    Sensitivity = roc_obj$sensitivities * 100         # Convert to percentage
  )
  
  # Create plot
  p <- ggplot(roc_data, aes(x = Specificity, y = Sensitivity)) +
    geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black") +
    geom_path(color = "#FF9999", linewidth = 1) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = paste0("Block: ", dataset_name, ", Using Comp(s): 1"),
      x = "100 - Specificity (%)",
      y = "Sensitivity (%)"
    ) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    annotate("text", x = 75, y = 25,
             label = sprintf("E vs L: %.4f", auc(roc_obj)),
             size = 5)
  
  # Save plot
  ggsave(
    filename = paste0("ROC-1-", tolower(dataset_name), "-colon.png"),
    plot = p,
    width = 6,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  return(p)
}




###___________perform the LOO cross-validation________ optional step________####


library(mixOmics)

# Run DIABLO with correct number of components
diablo_result <- block.splsda(X, Y, 
                              ncomp = 3,  # Specify number of components
                              design = design,
                              scale = TRUE,
                              near.zero.var = TRUE)

# Perform leave-one-out cross-validation
loo_result <- perf(diablo_result, validation = 'loo', nrepeat = 1)

# Rest of your error rate calculation code remains the same...

# Enhanced visualization
library(ggplot2)

# Create separate plots for BER and group-specific error rates
plot_error_rates <- function(error_data) {
  # Plot for BER
  p1 <- ggplot(error_data, aes(x = Component, y = BER, color = Dataset)) +
    geom_point(size = 3) +
    geom_line(aes(group = Dataset)) +
    facet_wrap(~ Distance) +
    theme_minimal() +
    labs(title = "Balanced Error Rate (BER) by Component",
         y = "BER") +
    ylim(0, 1)
  
  # Plot for group-specific errors
  plot_data <- reshape2::melt(error_data[, c("Dataset", "Distance", "Component", "E", "L")],
                              id.vars = c("Dataset", "Distance", "Component"))
  
  p2 <- ggplot(plot_data, aes(x = Component, y = value, color = variable, shape = Dataset)) +
    geom_point(size = 3) +
    geom_line(aes(group = interaction(Dataset, variable))) +
    facet_wrap(~ Distance) +
    theme_minimal() +
    labs(title = "Group-Specific Error Rates",
         y = "Error Rate",
         color = "Group") +
    ylim(0, 1)
  
  return(list(ber_plot = p1, group_plot = p2))
}

# Create and save plots
plots <- plot_error_rates(all_error_rates)
ggsave("BER_Plot_Colon.png", plots$ber_plot, width = 10, height = 6, bg = "white")
ggsave("Group_Error_Rates_Plot_Colon.png", plots$group_plot, width = 10, height = 6, bg = "white")

#  Export Error Rates
# Format error rates for export
error_rates_export <- data.frame(
  Dataset = all_error_rates$Dataset,
  Component = all_error_rates$Component,
  Distance_Method = all_error_rates$Distance,
  Error_Rate_E = all_error_rates$E,
  Error_Rate_L = all_error_rates$L,
  BER = all_error_rates$BER
)

# Export error rates
write.csv(error_rates_export, 
          "Error_Rates_Summary_Colon.csv", 
          row.names = FALSE)

