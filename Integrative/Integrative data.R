#############################################################
## DIABLO Analysis – Minimal Script (Required Outputs)
#############################################################

## Load libraries ----
library(mixOmics)   # v6.24 or later
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pROC)
library(gridExtra)
set.seed(12345)


## -------- 1. Pre-processing helper -----------------------------------------
preprocess_and_filter <- function(data, name, is_microbiome = FALSE) {
  message("\nProcessing ", name)
  
  if (!is_microbiome) {                                # Transcriptome
    mean_expression <- colMeans(data)
    keep <- mean_expression > quantile(mean_expression, 0.75)
    data <- data[, keep]
    
    var_genes <- apply(data, 2, var)
    keep_var <- var_genes > quantile(var_genes, 0.75)
    data <- data[, keep_var]
    
    data <- log2(data + 1)
    
    mad_scores <- apply(data, 1, mad)
    keep_samples <- mad_scores < (median(mad_scores) + 3 * mad(mad_scores))
    data <- data[keep_samples, ]
  } else {                                             # Microbiome
    keep <- colMeans(data > 0) >= 0.75
    data <- data[, keep]
    
    rel_abundance <- sweep(data, 1, rowSums(data), "/")
    abundance_filter <- colMeans(rel_abundance) > 0.1 / 100
    data <- data[, abundance_filter]
    
    clr_data <- mixOmics::logratio.transfo(data + 1, logratio = "CLR")
    pca_temp <- prcomp(clr_data)
    scores <- pca_temp$x[, 1:2]
    dist_center <- sqrt(rowSums(scores^2))
    keep_samples <- dist_center < (median(dist_center) + 3 * mad(dist_center))
    data <- data[keep_samples, ]
    
    data <- mixOmics::logratio.transfo(data + 1, logratio = "CLR")
  }
  
  data_scaled <- scale(data)
  message("Dimensions after filtering: ", nrow(data_scaled), " × ", ncol(data_scaled))
  data_scaled
}

## -------- 2. Data import ----------------------------------------------------
fpkm_data       <- read.csv("Colon_FPKM_cleaned.csv", row.names = 1)
genus_abundance <- read.csv("Colon_abundance.csv",   row.names = 1)

## Create metadata (E vs L)
metadata <- data.frame(
  Sample = colnames(fpkm_data),
  Group  = ifelse(grepl("^E", colnames(fpkm_data)), "E", "L"),
  row.names = colnames(fpkm_data)
)

## -------- 3. Pre-processing -------------------------------------------------
fpkm_proc  <- preprocess_and_filter(t(fpkm_data),       "FPKM")
genus_proc <- preprocess_and_filter(genus_abundance,    "Microbiome", TRUE)

## Align sample IDs after outlier removal
common_ids   <- intersect(rownames(fpkm_proc), rownames(genus_proc))
fpkm_proc    <- fpkm_proc[common_ids, ]
genus_proc   <- genus_proc[common_ids, ]
Y            <- metadata[common_ids, "Group", drop = TRUE]

## -------- 4. DIABLO model ---------------------------------------------------
X      <- list(Transcriptome = fpkm_proc, Microbiome = genus_proc)
design <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)

diablo_res <- block.splsda(
  X      = X,
  Y      = Y,
  design = design,
  ncomp  = 3,
  scale  = TRUE,
  near.zero.var = TRUE
)

## -------- 5. Plots & diagnostics -------------------------------------------
## 5a. Individual plot
png("DIABLO_Individual_Plot_Integrated_Colon.png", width = 1000, height = 800, res = 120)
plotIndiv(diablo_res, group = Y, legend = TRUE,
          title = "Integrated DIABLO Analysis – Colon",
          ellipse = TRUE, style = "graphics", conf.level = 0.95)
dev.off()

## 5b. DIABLO diagnostic plot
png("DIABLO_Diagnostic_Plot_Colon.png", width = 1000, height = 1000, res = 120)
plotDiablo(diablo_res, ncomp = 1)
dev.off()

## 5c. Global correlation circle (all blocks)
png("DIABLO_Correlation_Circle_Colon.png", width = 800, height = 800, res = 120)
plotVar(diablo_res, comp = c(1, 2), var.names = TRUE, legend = TRUE, pch = c(16, 16))
dev.off()


## 5D. ROC curve
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

# Assuming diablo_result is already available from previous steps
# Generate ROC curves for both datasets
roc_transcriptome <- create_roc_plot(diablo_res, "Transcriptome")
roc_microbiome <- create_roc_plot(diablo_res, "Microbiome")

# Combine plots side by side
combined_plot <- grid.arrange(roc_transcriptome, roc_microbiome, ncol = 2)

# Save the combined plot
ggsave("Combined_ROC_Curves_Colon.png", combined_plot, 
       width = 10, height = 6, dpi = 300, bg = "white")

## -------- 6. Top-feature extraction ----------------------------------------
get_top_features <- function(model, block, is_microbiome = FALSE,
                             top_percentage = NULL, top_n = NULL) {
  loadings   <- model$loadings[[block]][, 1]
  importance <- abs(loadings)
  feats  <- data.frame(
    Feature   = names(importance),
    Importance = importance,
    Direction  = ifelse(loadings > 0, "Positive", "Negative")
  ) %>% arrange(desc(Importance))
  
  if (is_microbiome) feats <- head(feats, top_n)
  else if (!is.null(top_percentage)) {
    n_top <- ceiling(nrow(feats) * top_percentage)
    feats <- head(feats, n_top)
  } else feats <- head(feats, top_n)
  feats
}

plot_top_features <- function(df, block) {
  df$Group <- ifelse(df$Direction == "Positive", "L", "E")
  ggplot(df, aes(reorder(Feature, Importance), Importance, fill = Group)) +
    geom_bar(stat = "identity", width = 0.1) +
    geom_point(aes(color = Group), size = 3, shape = 16) +
    coord_flip() +
    labs(title = paste("Top Features –", block, "– Colon"),
         x = "Feature", y = "Importance Score") +
    scale_fill_manual(values = c(E = "blue", L = "orange")) +
    scale_color_manual(values = c(E = "blue", L = "orange")) +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank())
}

## Extract & save
top_micro <- get_top_features(diablo_res, "Microbiome",  TRUE,  top_n        = 20)
top_trans <- get_top_features(diablo_res, "Transcriptome",FALSE, top_percentage = 0.10)

write.csv(top_micro,  "Top_Features_Microbiome_Colon.csv",   row.names = FALSE)
write.csv(top_trans,  "Top_Features_Transcriptome_Colon.csv",row.names = FALSE)

ggsave("Top_Features_Microbiome_Colon.png",
       plot_top_features(top_micro, "Microbiome"),
       width = 6, height = 8, dpi = 300, bg = "white")

ggsave("Top_Features_Transcriptome_Colon.png",
       plot_top_features(top_trans, "Transcriptome"),
       width = 6, height = 8, dpi = 300, bg = "white")

## -------- 7. Correlation circle – selected features ------------------------
# Restrict model data matrices to top features for clarity
diablo_res$X$Transcriptome <- diablo_res$X$Transcriptome[, colnames(diablo_res$X$Transcriptome) %in% top_trans$Feature]
diablo_res$X$Microbiome    <- diablo_res$X$Microbiome[,    colnames(diablo_res$X$Microbiome)    %in% top_micro$Feature]

var_dat <- plotVar(diablo_res, comp = c(1, 2), cutoff = 0.5,
                   var.names = TRUE, legend = TRUE, plot = FALSE)

var_df <- data.frame(
  Comp1   = var_dat$x,
  Comp2   = var_dat$y,
  Block   = var_dat$Block,
  Feature = var_dat$names
)
write.csv(var_df, "Correlation_Results_Top_Features_Colon.csv", row.names = FALSE)

## Highlight only designated features
roi_features <- c("COL3A1","SPARC","HMGA1","HSP90AB1","S100A11",
                  "COL1A1","COL1A2","OGDH","LENG8","Roseburia","Prevotella")
var_df$Label <- ifelse(var_df$Feature %in% roi_features, var_df$Feature, "")

shape_map  <- c(Microbiome = 23, Transcriptome = 24)
color_map  <- c(Microbiome = "blue", Transcriptome = "aquamarine4")

corr_plot <- ggplot(var_df, aes(Comp1, Comp2)) +
  # outer & inner circles
  geom_path(data = data.frame(x = cos(seq(0, 2*pi, length.out = 100)),
                              y = sin(seq(0, 2*pi, length.out = 100))),
            aes(x, y), colour = "darkgrey") +
  geom_path(data = data.frame(x = cos(seq(0, 2*pi, length.out = 100))*sqrt(0.5),
                              y = sin(seq(0, 2*pi, length.out = 100))*sqrt(0.5)),
            aes(x, y), colour = "darkgrey") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "darkgrey") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "darkgrey") +
  geom_point(aes(shape = Block, colour = Block), size = 3, fill = NA) +
  geom_text_repel(data = subset(var_df, Label != ""),
                  aes(label = Label, colour = Block), size = 3.5) +
  scale_shape_manual(values = shape_map) +
  scale_colour_manual(values = color_map) +
  coord_fixed() +
  labs(title = "Correlation Circle Plot – Selected Features – Colon",
       x = "Component 1", y = "Component 2") +
  theme_minimal(base_size = 14) +
  theme(panel.border = element_rect(colour = "black", fill = NA))

ggsave("Correlation_Circle_Plot2_Selected_Features_Colon.png",
       corr_plot, width = 8, height = 6, dpi = 300, bg = "white")

## -------- 8. Save full model (optional) ------------------------------------
saveRDS(diablo_res, "diablo_result_integrated_colon.rds")



#________calculate the CLR-transformed abundance for Roseburia and Prevotella__________####

# Load required library
library(mixOmics)

# Assuming your data is in a CSV file with samples as rows and genera as columns
genus_abundance <- read.csv("Colon_abundance.csv", row.names = 1)

# Filter for samples in E and L groups
E_samples <- genus_abundance[grep("^E", rownames(genus_abundance)), ]
L_samples <- genus_abundance[grep("^L", rownames(genus_abundance)), ]

# Add pseudocount of 1 to handle zeros
genus_abundance_plus1 <- genus_abundance + 1

# Apply CLR transformation
clr_transformed <- mixOmics::logratio.transfo(genus_abundance_plus1, 
                                              logratio = 'CLR')

# Extract Roseburia and Prevotella abundances
roseburia_clr <- clr_transformed[, "Roseburia"]
prevotella_clr <- clr_transformed[, "Prevotella"]


# Prepare data for plotting
plot_data <- data.frame(
  Abundance = c(roseburia_clr, prevotella_clr),
  Group = rep(c(rep("E", sum(grepl("^E", rownames(clr_transformed)))),
                rep("L", sum(grepl("^L", rownames(clr_transformed))))), 2),
  Genus = rep(c("Roseburia", "Prevotella"), each = nrow(clr_transformed))
)


# Assuming your data is already CLR transformed
# Create a data frame for plotting with proper structure
plot_data <- data.frame(
  Abundance = c(roseburia_clr, prevotella_clr),
  Group = rep(c(rep("E", sum(grepl("^E", rownames(clr_transformed)))),
                rep("L", sum(grepl("^L", rownames(clr_transformed))))), 2),
  Genus = rep(c("Roseburia", "Prevotella"), each = nrow(clr_transformed))
)



# Create improved plot with specified colors
ggplot(plot_data, aes(x = Group, y = Abundance, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_point(aes(color = Group, shape = Group), 
             position = position_jitter(width = 0.2),
             size = 3) +
  scale_shape_manual(values = c("E" = 16, "L" = 17)) +
  scale_fill_manual(values = c("E" = "#1F77B4", "L" = "#FF7F0E")) +
  scale_color_manual(values = c("E" = "#1F77B4", "L" = "#FF7F0E")) +
  facet_wrap(~Genus, scales = "free_y") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "plain"),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    legend.position = "none"
  ) +
  labs(
    y = "CLR transformed abundance",
    x = "Group"
  )

# Save the plot
ggsave("CLR_Boxplots_Colored.png", 
       width = 6, 
       height = 4, 
       dpi = 300)



#_________calculate the CLR-transformed expression for selected genes__________####


# Load necessary libraries
library(ggplot2)
library(reshape2)

# Assuming your data is in a CSV file with samples as rows and genes as columns
gene_expression <- read.csv("Colon_FPKM_targeted.csv", row.names = 1)

# Transpose the data
gene_expression_t <- t(gene_expression)

# Convert to data frame
gene_expression_df <- as.data.frame(gene_expression_t)

# Apply log2(FPKM + 1) transformation
gene_expression_plus1 <- gene_expression_df + 1
log2_transformed <- log2(gene_expression_plus1)

# Melt the data frame to long format
long_data <- melt(log2_transformed, id.vars = NULL, variable.name = "Gene", value.name = "Expression")

# Assuming the row names of the transposed data correspond to different samples and are related to "E" or "L"
long_data$Sample <- rownames(log2_transformed)
long_data$Group <- ifelse(grepl("^E", long_data$Sample), "E", "L")

# Create the improved plot
ggplot(long_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_point(aes(color = Group, shape = Group), 
             position = position_jitter(width = 0.2),
             size = 3) +
  scale_shape_manual(values = c("E" = 16, "L" = 17)) +
  scale_fill_manual(values = c("E" = "#1F77B4", "L" = "#FF7F0E")) +
  scale_color_manual(values = c("E" = "#1F77B4", "L" = "#FF7F0E")) +
  facet_wrap(~Gene, scales = "free_y", ncol = 4) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    y = "Log2 transformed expression",
    x = "Group",
    title = "Log2-transformed Gene Expression in E and L Groups"
  )


# Save the plot
ggsave("Log2_Boxplots_Targeted_Genes.png", 
       width = 12, 
       height = 10, 
       dpi = 300)
