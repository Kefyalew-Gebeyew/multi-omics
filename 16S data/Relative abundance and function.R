
# Load required libraries
library(streamgraph)
library(reshape2)
library(RColorBrewer)
library(htmlwidgets)
library(webshot)
library(magick)
library(grid)
library(ggplot2)
library(png)


# Load the data
data <- read.csv("/Users/kefe/Desktop/Data/Bacterial genera/Stacked bar/Phylum.csv")

# Melt the data to long format
long_data <- melt(data, id.vars = c("sample_ID", "Tissue", "Age"),
                  variable.name = "Phylum", value.name = "Abundance")

# Create a sequence variable for the x-axis
long_data$Sequence <- as.numeric(as.factor(long_data$sample_ID))

# Generate color palettes
palette1 <- brewer.pal(12, "Paired")
palette2 <- brewer.pal(9, "Set1")
palette3 <- brewer.pal(8, "Set3")
additional_colors <- c("#A52A2A", "#5F9EA0", "#D2691E", "#FF7F50", "#6495ED", "#FFF8DC", "#DC143C", "#00FFFF", "#00008B")
combined_palette <- c(palette1, palette2, palette3, additional_colors)

# Define actual phyla names
actual_phyla_names <- unique(long_data$Phylum)

# Verify that we have enough colors for all phyla
if (length(actual_phyla_names) > length(combined_palette)) {
  stop("Not enough colors in the combined palette for all phyla.")
}

# Assign the real names to the colors
phyla_colors <- setNames(combined_palette[1:length(actual_phyla_names)], actual_phyla_names)

# Create the stacked area plot
color_bar_plot <- ggplot(long_data, aes(x = Sequence, y = Abundance, fill = Phylum)) +
  geom_area(position = "fill") +
  scale_fill_manual(values = phyla_colors) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Phylum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~ Tissue, scales = "free_x") +
  scale_x_continuous(expand = c(0, 0))

# Save the plot
ggsave("Phylum_stacked_area.png", width = 10, height = 6, bg = "White")



###---------------------------------Functional analysis and visualization---------------------------------------#######

# Load your data
data <- read.csv("/Users/kefe/Desktop/Data/Bacterial genera/Stacked bar/cog_Level2_Function.csv")

# Check the column names
colnames(data)

# Extract tissue information from the first row
tissue_info <- as.character(data[1, -1])  # Exclude the first column (SampleID) and convert to character
data <- data[-1, ]  # Remove the first row

# Reshape the data to include tissue information in a separate column
data_long <- melt(data, id.vars = "SampleID", variable.name = "Sample", value.name = "Abundance")
data_long$Tissue <- rep(tissue_info, each = nrow(data))

# Convert Abundance to numeric
data_long$Abundance <- as.numeric(data_long$Abundance)

# Create a sequence variable for the x-axis
data_long$Sequence <- as.numeric(as.factor(data_long$Sample))

# Check the structure of the long_data dataframe
str(data_long)

# Generate color palettes
palette1 <- brewer.pal(12, "Paired")
palette2 <- brewer.pal(9, "Set1")
palette3 <- brewer.pal(8, "Set3")
additional_colors <- c("#A52A2A", "#5F9EA0", "#D2691E", "#FF7F50", "#6495ED", "#FFF8DC", "#DC143C", "#00FFFF", "#00008B")
more_colors <- colorRampPalette(brewer.pal(9, "Set1"))(100) # Generate 100 more colors
combined_palette <- c(palette1, palette2, palette3, additional_colors, more_colors)

# Define actual COG names
actual_cog_names <- unique(data_long$SampleID)

# Verify that we have enough colors for all COG categories
if (length(actual_cog_names) > length(combined_palette)) {
  stop("Not enough colors in the combined palette for all COG categories.")
}

# Assign the real names to the colors
cog_colors <- setNames(combined_palette[1:length(actual_cog_names)], actual_cog_names)

# Print the colors (optional)
print(cog_colors)

# Create the stacked area chart
ggplot(data_long, aes(x = Sequence, y = Abundance, fill = SampleID)) +
  geom_area(position = "fill") +
  scale_fill_manual(values = cog_colors) +
  labs(x = "Sample", y = "Relative Abundance", fill = "PicrustFunction") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~ Tissue, scales = "free_x") +  # Add facet_wrap to separate by tissue and allow free x-axis scaling
  scale_x_continuous(expand = c(0, 0))  # Remove spaces on the x-axis

# Save the plot
ggsave("cog_stacked_area_chart_with_tissue.png", width = 10, height = 6, bg = "White")