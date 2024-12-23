
######--------Age associated bacteria in a Separate tissue------------------## for shared bacteria------##


# Load necessary libraries
library(data.table)
library(Maaslin2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(ggrepel)

# Load the data
otu_data <- fread('OTU_kefe.csv')
metadata <- fread('SampleID.csv')
taxa <- fread('Taxa.csv')

# Ensure the column names are consistent
colnames(metadata)[1] <- "SampleID"
colnames(taxa)[1] <- "OTU_ID"

# Reshape OTU data to long format
otu_long <- melt(otu_data, id.vars = "SampleID", variable.name = "OTU_ID", value.name = "Count")
otu_long <- merge(otu_long, taxa, by = "OTU_ID")
otu_long$Phylum_Genus <- paste(otu_long$Phylum, otu_long$Genus, sep = "_")
otu_long <- otu_long[!is.na(otu_long$Phylum) & !is.na(otu_long$Genus), ]

# Merge OTU Data and Metadata
full_data <- merge(otu_long, metadata, by = "SampleID")
full_data$Tissue <- factor(full_data$Tissue, levels = c("Colon", "Jejunum", "Rumen"))

# Convert full_data to data.table
setDT(full_data)

# CSS Normalization
total_counts <- full_data[, .(TotalCount = sum(Count)), by = "SampleID"]
full_data <- merge(full_data, total_counts, by = "SampleID", suffixes = c("", ".total"))
median_total_counts <- median(full_data$TotalCount)
full_data[, NormalizedCount := (Count / TotalCount) * median_total_counts]

# Convert Data to Wide Format for MaAsLin2
wide_data <- dcast(full_data, SampleID + Tissue + Age ~ Phylum_Genus, value.var = "NormalizedCount", fun.aggregate = sum, fill = 0)

# Ensure wide_data is a data.table
setDT(wide_data)

# Create input_metadata and input_data
input_metadata <- wide_data[, .(SampleID, Tissue, Age)]
input_data <- wide_data[, !c("SampleID", "Tissue", "Age"), with = FALSE]
rownames(input_data) <- input_metadata$SampleID
rownames(input_metadata) <- input_metadata$SampleID
input_metadata$Tissue <- factor(input_metadata$Tissue, levels = c("Colon", "Jejunum", "Rumen"))

# Run MaAsLin2 analysis for each tissue type
results_list <- list()
tissues <- unique(input_metadata$Tissue)

for (tissue in tissues) {
  tissue_data <- input_data[input_metadata$Tissue == tissue, ]
  tissue_metadata <- input_metadata[input_metadata$Tissue == tissue, ]
  
  # Ensure row names are correctly set
  rownames(tissue_data) <- tissue_metadata$SampleID
  rownames(tissue_metadata) <- tissue_metadata$SampleID
  
  fit <- Maaslin2(
    input_data = tissue_data, 
    input_metadata = tissue_metadata, 
    output = paste0("Maaslin2_results_", tissue),
    fixed_effects = c("Age"),
    normalization = "NONE",
    transform = "LOG",
    max_significance = 0.05,
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )
  
  results_list[[tissue]] <- fit$results %>%
    filter(qval < 0.05) %>%
    mutate(Tissue = tissue)
}

# Combine results from all tissues
combined_results <- bind_rows(results_list)

# Filter for top 15 bacteria in each tissue
top_bacteria <- combined_results %>%
  group_by(Tissue) %>%
  top_n(15, abs(coef)) %>%
  ungroup()

# Reshape data for bubble plot
bubble_data <- top_bacteria %>%
  select(feature, Tissue, coef)

# Adjust the color scheme to match the previous example (Colon = green, Jejunum = orange, Rumen = purple)
color_scheme <- c("Colon" = "#66C2A5", "Jejunum" = "#FC8D62", "Rumen" = "#8DA0CB")

# Create bubble plot with ggplot2, matching previous visualization style
bubble_plot <- ggplot(bubble_data, aes(x = Tissue, y = feature, size = abs(coef), color = Tissue)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(3, 9), name = "Coefficient (abs)") +
  scale_color_manual(values = color_scheme) + # Use the custom color scheme
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dashed", color = "black", size = 0.5) +  # Add vertical lines between tissues
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border
    panel.background = element_rect(fill = "white"),  # White background
    plot.background = element_rect(fill = "white")  # Ensure full plot background is white
  ) +
  labs(title = "Top 15 Bacteria Associations in Each Tissue",
       x = "Tissue", y = "Phylum_Genus")

# Save the plot as a PNG file with a white background
ggsave("bubble_plot_with_vline.png", plot = bubble_plot, width = 6, height = 8, dpi = 300, bg = "white")

# Export significant results for each tissue type
for (tissue in tissues) {
  tissue_results <- combined_results %>%
    filter(Tissue == tissue)
  
  write.csv(tissue_results, paste0("age_associated_genera_", tissue, ".csv"), row.names = FALSE)
}

######--------visualize genera shows significantly changes over time in at least in a tissue------------------## for shared bacteria------##

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(reshape2)

# Assuming `full_data` is already prepared, we will aggregate it by Age, Tissue, and Genera (Phylum_Genus)
# Summarize the data by Age, Tissue, and Genera
age_abundance_data <- full_data %>%
  group_by(Tissue, Phylum_Genus) %>%
  summarize(Total_Abundance = sum(NormalizedCount, na.rm = TRUE)) %>%
  ungroup()

# Select the top 30 genera for each tissue based on total abundance
top_30_genera <- age_abundance_data %>%
  group_by(Tissue) %>%
  top_n(30, Total_Abundance) %>%
  pull(Phylum_Genus)

# Filter the original dataset to include only the top 30 genera in each tissue
filtered_age_abundance_data <- full_data %>%
  filter(Phylum_Genus %in% top_30_genera) %>%
  group_by(Age, Phylum_Genus, Tissue) %>%
  summarize(Mean_Abundance = mean(NormalizedCount, na.rm = TRUE)) %>%
  ungroup()

# Adjust the color scheme to match the previous example (Colon = green, Jejunum = orange, Rumen = purple)
color_scheme <- c("Colon" = "#66C2A5", "Jejunum" = "#FC8D62", "Rumen" = "#8DA0CB")

# Create the bubble plot with ggplot2 for the top 30 genera in each tissue
bubble_plot_age <- ggplot(filtered_age_abundance_data, aes(x = Age, y = Phylum_Genus, size = Mean_Abundance, color = Tissue)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(2, 12), name = "Mean Abundance") +
  scale_color_manual(values = color_scheme) +  # Use the custom color scheme for tissues
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border
    panel.background = element_rect(fill = "white"),  # White background
    plot.background = element_rect(fill = "white")  # Ensure full plot background is white
  ) +
  labs(title = "Top 30 Abundant Genera Across Age Time Points", 
       x = "Age (Months)", 
       y = "Genus") +
  facet_wrap(~ Tissue)  # Create separate panels for each tissue

# Save the plot as a PNG file with a white background
ggsave("bubble_plot_top30_age_timepoints.png", plot = bubble_plot_age, width = 12, height = 10, dpi = 300, bg = "white")