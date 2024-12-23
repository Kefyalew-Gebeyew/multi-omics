
###---------------------------Separate age time point for each tissue------------------########---------------------######


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(vegan)
library(RColorBrewer)
library(grid)

# Load the long format data
data <- read.csv('OTU_new.csv')

# Ensure Count is numeric
data$Count <- as.numeric(data$Count)

# Exclude rows with NA values
data <- na.omit(data)

# Adding a metadata frame, if you have additional metadata such as Age or Tissue
metadata <- read.csv('SampleID.csv')

# Exclude rows with NA values from metadata as well
metadata <- na.omit(metadata)

# Merge the OTU data with metadata
merged_data <- merge(data, metadata, by = "SampleID")

# Create the species table
species_matrix <- xtabs(Count ~ SampleID + Taxonomy, data = merged_data)

# Calculate Shannon and Simpson diversity indices
shannon_diversity <- diversity(species_matrix, index = "shannon")
simpson_diversity <- diversity(species_matrix, index = "simpson")

# Add diversity indices to your data
merged_data$Shannon <- shannon_diversity[match(merged_data$SampleID, names(shannon_diversity))]
merged_data$Simpson <- simpson_diversity[match(merged_data$SampleID, names(simpson_diversity))]

# Exclude rows with NA values after adding diversity indices
merged_data <- na.omit(merged_data)

# Perform Kruskal-Wallis test for Shannon diversity
kw_shannon <- kruskal.test(Shannon ~ interaction(Age, Tissue), data = merged_data)
print(kw_shannon)

# Perform Dunn's post hoc test if Kruskal-Wallis is significant for Shannon diversity
dunn_shannon_results <- NULL
if (kw_shannon$p.value < 0.05) {
  group <- interaction(merged_data$Age, merged_data$Tissue)
  dunn_shannon_results <- dunn.test(x = merged_data$Shannon, g = group, method="bonferroni")
  print(dunn_shannon_results)
}

# Perform Kruskal-Wallis test for Simpson diversity
kw_simpson <- kruskal.test(Simpson ~ interaction(Age, Tissue), data = merged_data)
print(kw_simpson)

# Perform Dunn's post hoc test if Kruskal-Wallis is significant for Simpson diversity
dunn_simpson_results <- NULL
if (kw_simpson$p.value < 0.05) {
  group <- interaction(merged_data$Age, merged_data$Tissue)
  dunn_simpson_results <- dunn.test(x = merged_data$Simpson, g = group, method="bonferroni")
  print(dunn_simpson_results)
}

# Define a color palette for age time points
color_palette <- c("firebrick2", "chartreuse1", "hotpink4", "maroon2", "aquamarine4", "cyan1")

# Define shapes for each tissue
shape_palette <- c("Rumen" = 16, "Jejunum" = 17, "Colon" = 15)  # Circle, triangle, square

# Boxplot for Shannon Diversity by Age and Tissue with adjusted outlier shape and size
shannon_plot <- ggplot(merged_data, aes(x = as.factor(Age), y = Shannon, color = as.factor(Age), shape = Tissue)) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.6, lwd = 0.5, position = position_dodge(width = 0.75)) +  # Remove default outliers and fill color, adjust cap size, reduce gap
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 2, aes(color = as.factor(Age))) +  # Add jittered points with color
  labs(title = "Shannon Diversity by Age and Tissue", x = "Age", y = "Shannon Index") +
  theme_minimal() +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +  # Different shapes for each tissue
  facet_wrap(~ Tissue, scales = "free_x", labeller = as_labeller(c("Colon" = "Colon", "Jejunum" = "Jejunum", "Rumen" = "Rumen"))) +  # Separate panels for each tissue with labels
  theme(strip.text = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border around each panel
        axis.line = element_line(colour = "black"),
        strip.background = element_rect(fill = "grey90", color = "black"),
        panel.spacing = unit(0.5, "lines"))  # Adjust spacing between panels

# Boxplot for Simpson Diversity by Age and Tissue with adjusted outlier shape and size
simpson_plot <- ggplot(merged_data, aes(x = as.factor(Age), y = Simpson, color = as.factor(Age), shape = Tissue)) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.6, lwd = 0.5, position = position_dodge(width = 0.75)) +  # Remove default outliers and fill color, adjust cap size, reduce gap
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 2, aes(color = as.factor(Age))) +  # Add jittered points with color
  labs(title = "Simpson Diversity by Age and Tissue", x = "Age", y = "Simpson Index") +
  theme_minimal() +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +  # Different shapes for each tissue
  facet_wrap(~ Tissue, scales = "free_x", labeller = as_labeller(c("Colon" = "Colon", "Jejunum" = "Jejunum", "Rumen" = "Rumen"))) +  # Separate panels for each tissue with labels
  theme(strip.text = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border around each panel
        axis.line = element_line(colour = "black"),
        strip.background = element_rect(fill = "grey90", color = "black"),
        panel.spacing = unit(0.5, "lines"))  # Adjust spacing between panels

# Save the plots to files
ggsave("shannon_plot.png", plot = shannon_plot, width = 12, height = 6, dpi = 300, bg = "white")



###---------------------------PCoA and PCA------------------########---------------------######

# Load necessary libraries
install.packages(c("vegan", "ggplot2", "ape"))
library(vegan)
library(ggplot2)
library(ape)

# Load the data
otu_data <- read.csv('OTU_new.csv')
metadata <- read.csv('SampleID.csv')

# Ensure Count is numeric and exclude rows with NA values
otu_data$Count <- as.numeric(otu_data$Count)
otu_data <- na.omit(otu_data)
metadata <- na.omit(metadata)

# Merge the OTU data with metadata
merged_data <- merge(otu_data, metadata, by = "SampleID")

# Create species matrix
species_matrix <- xtabs(Count ~ SampleID + Taxonomy, data = merged_data)

# Calculate Bray-Curtis distance matrix
bray_curtis_dist <- vegdist(species_matrix, method = "bray")

# Perform PCoA
pcoa_result <- cmdscale(bray_curtis_dist, eig = TRUE, k = 2)

# Calculate the percentage of variance explained for PCoA
pcoa_var_explained <- pcoa_result$eig[1:2] / sum(pcoa_result$eig) * 100

# Extract coordinates
pcoa_coords <- as.data.frame(pcoa_result$points)
colnames(pcoa_coords) <- c("PCoA1", "PCoA2")
pcoa_coords$SampleID <- rownames(pcoa_coords)

# Merge with metadata
pcoa_data <- merge(pcoa_coords, metadata, by = "SampleID")

# Convert Age to a factor
pcoa_data$Age <- as.factor(pcoa_data$Age)

# Define a color palette for age time points
color_palette <- c("firebrick2", "chartreuse1", "hotpink4", "maroon2", "aquamarine4", "cyan1")

# Define shapes for each tissue
shape_palette <- c("Rumen" = 16, "Jejunum" = 17, "Colon" = 15)

# Plot PCoA
pcoa_plot <- ggplot(pcoa_data, aes(x = PCoA1, y = PCoA2, color = Age, shape = Tissue)) +
  geom_point(size = 3) +
  labs(title = "PCoA of Microbial Communities", 
       x = paste0("PCoA1: ", round(pcoa_var_explained[1], 1), "%"), 
       y = paste0("PCoA2: ", round(pcoa_var_explained[2], 1), "%")) +
  theme_minimal() +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  theme(panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA))

# Save PCoA plot
ggsave("pcoa_plot.png", plot = pcoa_plot, width = 12, height = 6, dpi = 300, bg = "white")

# Perform PCA
pca_result <- prcomp(species_matrix, scale. = TRUE)

# Calculate the percentage of variance explained for PCA
pca_var_explained <- summary(pca_result)$importance[2, 1:2] * 100

# Extract PCA coordinates
pca_coords <- as.data.frame(pca_result$x)
pca_coords$SampleID <- rownames(pca_coords)

# Merge with metadata
pca_data <- merge(pca_coords, metadata, by = "SampleID")

# Convert Age to a factor
pca_data$Age <- as.factor(pca_data$Age)

# Plot PCA
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Age, shape = Tissue)) +
  geom_point(size = 3) +
  labs(title = "PCA of Microbial Communities", 
       x = paste0("PC1: ", round(pca_var_explained[1], 1), "%"), 
       y = paste0("PC2: ", round(pca_var_explained[2], 1), "%")) +
  theme_minimal() +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  theme(panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA))

# Save PCA plot
ggsave("pca_plot.png", plot = pca_plot, width = 12, height = 6, dpi = 300, bg = "white")

# Perform PERMANOVA
permanova_result <- adonis2(bray_curtis_dist ~ Age * Tissue, data = metadata, permutations = 999)

# Extract the ANOVA table from the PERMANOVA result
permanova_table <- as.data.frame(permanova_result$aov.tab)

# Save PERMANOVA results to CSV
write.csv(permanova_table, "permanova_results.csv")

ggsave("simpson_plot.png", plot = simpson_plot, width = 12, height = 6, dpi = 300, bg = "white")