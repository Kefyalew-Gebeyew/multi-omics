

###-------------------------enterotypes transition pie chart across age time point regardless of tissue---------------

# Load necessary libraries
library(phyloseq)
library(metagenomeSeq)
library(DirichletMultinomial)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(viridis)

# Load OTU Data
otu_path <- "/Users/kefe/Desktop/Data/Bacterial genera/DMN_Genera/OTU_kefe.csv"
otu_table <- read.csv(otu_path, check.names = FALSE)
rownames(otu_table) <- otu_table[, 1]
otu_table <- otu_table[, -1]
otu_matrix <- as.matrix(otu_table)

# Load Metadata
metadata_path <- "/Users/kefe/Desktop/Data/Bacterial genera/DMN_Genera/SampleID.csv"
metadata <- read.csv(metadata_path)
metadata$SampleID <- as.character(metadata$SampleID)
metadata$Tissue <- as.factor(metadata$Tissue)
metadata$Age <- as.numeric(as.character(metadata$Age))

# Align OTU Matrix with Metadata
common_samples <- intersect(rownames(otu_matrix), metadata$SampleID)
otu_matrix <- otu_matrix[common_samples, ]
metadata <- metadata[metadata$SampleID %in% common_samples, ]
sample_data <- sample_data(metadata)
rownames(sample_data) <- metadata$SampleID

# Load Taxa Data
taxa_path <- "/Users/kefe/Desktop/Data/Bacterial genera/DMN_Genera/Taxa.csv"
taxa_data <- read.csv(taxa_path, check.names = FALSE)
rownames(taxa_data) <- taxa_data$OTU_ID
taxa_data$OTU_ID <- NULL
taxa_data[is.na(taxa_data)] <- "Unclassified"
taxa_matrix <- as.matrix(taxa_data)
tax_table <- tax_table(taxa_matrix)

# Create phyloseq object
physeq <- phyloseq(otu_table(otu_matrix, taxa_are_rows = FALSE), sample_data, tax_table)

# Normalize Using CSS Method
physeq_metaseq <- phyloseq_to_metagenomeSeq(physeq)
physeq_metaseq_norm <- cumNorm(physeq_metaseq, p = cumNormStat(physeq_metaseq))
CSS_phyloseq <- MRcounts(physeq_metaseq_norm, norm = TRUE)
CSS_phyloseq_ps <- merge_phyloseq(phyloseq(otu_table(CSS_phyloseq, taxa_are_rows = TRUE), sample_data, tax_table))

# Filter Low Abundant Taxa and Agglomerate Data at Genus Level
abundance_threshold <- 0.005
taxa_sums_before <- taxa_sums(CSS_phyloseq_ps)
taxa_to_keep <- taxa_sums_before > (abundance_threshold * sum(taxa_sums_before))
pseq_taxa_fil <- prune_taxa(taxa_to_keep, CSS_phyloseq_ps)
pseq_genus <- tax_glom(pseq_taxa_fil, "Genus")

# Ensure Genus names are retained and unique
glomTax <- tax_table(pseq_genus)[, "Genus"]
unique_genus <- make.unique(as.character(glomTax))
tax_table(pseq_genus)[, "Genus"] <- unique_genus

# Prepare Input Data for DMM Analysis
glomOTU <- otu_table(pseq_genus)
glomTable <- as.data.frame(glomOTU)
glomTable$Genus <- unique_genus

# Set rownames to unique genus names
rownames(glomTable) <- glomTable$Genus
glomTable <- glomTable[, -which(names(glomTable) == "Genus")]

# Transpose for DMM analysis
dmm_input_genus <- t(glomTable)

# Perform multiple runs to get a stable estimate of the number of components
set.seed(123)  # For reproducibility
num_runs <- 5
optimal_components_list <- numeric(num_runs)

for (i in 1:num_runs) {
  model_fits <- sapply(2:10, function(k) {
    model <- dmn(dmm_input_genus, k = k, verbose = TRUE)
    laplace(model)
  })
  model_fit_df <- data.frame(Number_of_Components = 2:10, Model_Fit = model_fits)
  optimal_components_list[i] <- model_fit_df$Number_of_Components[which.min(model_fit_df$Model_Fit)]
}

# Average the results to get a stable estimate
optimal_components <- round(mean(optimal_components_list))
cat("Stable estimate of the optimal number of Dirichlet components:", optimal_components, "\n")

# Fit DMM Model with Optimal Number of Components
best_model <- dmn(dmm_input_genus, k = optimal_components, verbose = TRUE)

# Assign Enterotypes
prob_matrix <- predict(best_model, newdata = dmm_input_genus, type = "prob")
enterotype_assignments <- max.col(prob_matrix, ties.method = "first")
metadata <- metadata[metadata$SampleID %in% rownames(dmm_input_genus), ]
if (length(enterotype_assignments) != nrow(metadata)) {
  enterotype_assignments <- enterotype_assignments[1:nrow(metadata)]
}
metadata$Enterotype <- factor(enterotype_assignments, levels = 1:optimal_components)

# Print enterotype assignments for debugging
print(table(metadata$Enterotype))

# Save processed data for consistency in both plots
write.csv(metadata, "processed_metadata_combined.csv", row.names = FALSE)
saveRDS(pseq_genus, "processed_pseq_genus_combined.rds")

# Prepare data for pie chart visualization
pie_data <- metadata %>%
  group_by(Age, Enterotype) %>%
  summarise(Count = n()) %>%
  spread(Enterotype, Count, fill = 0)

# Convert to long format for ggplot
pie_data_long <- pie_data %>%
  gather(key = "Enterotype", value = "Count", -Age)

# Adjust the plot direction to be horizontal and the legend on the right side
ggplot(pie_data_long, aes(x = "", y = Count, fill = Enterotype)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ Age, nrow = 1, strip.position = "bottom") +  # Position the age labels at the bottom
  scale_fill_viridis_d() +
  theme_void() +
  theme(
    legend.position = "right",  # Place legend on the right side
    strip.text = element_text(size = 8, face = "bold", color = "black"),
    strip.background = element_rect(fill = "lightgrey", color = "black"),  # Keep strip background
    panel.spacing = unit(1, "lines"),  # Space between panels
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.placement = "outside",  # Place strips outside the plot
    plot.title = element_text(hjust = 0.5, size = 14, margin = margin(t = 10, b = 10)),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)  # Adjust the plot margins
  ) +
  labs(title = "Enterotype Distribution Across Ages", fill = "Enterotype")

# Save the plot with adjusted dimensions to reduce white space
ggsave("horizontal_enterotype_pie_charts_with_labels_bottom.png", bg = "white", width = 14, height = 6)





####--------------------------Sanky plot for all age time point for jejunum------------------------####

# Load necessary libraries
library(dplyr)
library(networkD3)
library(webshot)

# Load the processed metadata
metadata <- read.csv("processed_metadata.csv")

# Extract the numeric part of SampleID to serve as SubjectID
metadata$SubjectID <- sub("^[A-Z]+|[A-Z]+$", "", metadata$SampleID)

# Define age time points
metadata$AgeTimePoint <- as.factor(metadata$Age)

# Prepare data for all tissues combined
cat("Processing all tissues combined\n")

transitions <- data.frame(Source = character(), Target = character(), Value = integer())
age_time_points <- c("3", "4", "5", "6", "7", "8")

for (i in 1:(length(age_time_points) - 1)) {
  age_time_point_1 <- age_time_points[i]
  age_time_point_2 <- age_time_points[i + 1]
  
  for (enterotype_1 in unique(metadata$Enterotype)) {
    for (enterotype_2 in unique(metadata$Enterotype)) {
      filtered_data_1 <- metadata %>% 
        filter(AgeTimePoint == age_time_point_1 & Enterotype == enterotype_1) %>% 
        distinct(SubjectID, .keep_all = TRUE) # Ensure unique SubjectID
      
      filtered_data_2 <- metadata %>% 
        filter(AgeTimePoint == age_time_point_2 & Enterotype == enterotype_2) %>% 
        distinct(SubjectID, .keep_all = TRUE) # Ensure unique SubjectID
      
      if (nrow(filtered_data_1) > 0 & nrow(filtered_data_2) > 0) {
        count <- nrow(filtered_data_1 %>% inner_join(filtered_data_2, by = "SubjectID"))
        
        if (count > 0) {
          transitions <- rbind(transitions, data.frame(Source = paste(enterotype_1, age_time_point_1, sep = "_"),
                                                       Target = paste(enterotype_2, age_time_point_2, sep = "_"),
                                                       Value = count))
        }
      }
    }
  }
}

# Create Nodes Data Frame
nodes <- data.frame(name = unique(c(transitions$Source, transitions$Target)))

# Convert Source and Target to Indices
transitions$source <- match(transitions$Source, nodes$name) - 1
transitions$target <- match(transitions$Target, nodes$name) - 1

# Define a color scale for the enterotypes
color_scale <- 'd3.scaleOrdinal().domain(["1_3", "2_3", "3_3", "1_4", "2_4", "3_4",
                                          "1_5", "2_5", "3_5", "1_6", "2_6", "3_6",
                                          "1_7", "2_7", "3_7", "1_8", "2_8", "3_8"])
                        .range(["#1f77b4", "#ff7f0e", "#2ca02c", "#1f77b4", "#ff7f0e", "#2ca02c",
                                "#1f77b4", "#ff7f0e", "#2ca02c", "#1f77b4", "#ff7f0e", "#2ca02c",
                                "#1f77b4", "#ff7f0e", "#2ca02c"])'

# Create the Sankey Diagram
sankey <- sankeyNetwork(Links = transitions, Nodes = nodes, Source = "source", Target = "target", 
                        Value = "Value", NodeID = "name", units = "Subjects", fontSize = 12, nodeWidth = 30,
                        colourScale = color_scale, sinksRight = FALSE, nodePadding = 10)

# Customizing the HTML to hide the labels using CSS
# Save the Sankey diagram as an HTML file
html_file <- "sankey_diagram_all_tissues.html"
saveNetwork(sankey, html_file)

# Read the HTML content
html_content <- readLines(html_file)

# Add custom CSS to hide node labels
css <- '<style> .node text { display: none; } </style>'
html_content <- gsub("</head>", paste0(css, "</head>"), html_content)

# Write back the modified HTML content
writeLines(html_content, html_file)

# Convert the HTML to PNG
png_file <- "sankey_diagram_all_tissues.png"
webshot(html_file, file = png_file, vwidth = 1000, vheight = 800)





####--------------------------Markov Chain plot across age time for each tissue------------------------####

# Load necessary libraries
library(dplyr)
library(markovchain)
library(igraph)

# Read the CSV file
metadata <- read.csv("processed_metadata.csv", stringsAsFactors = FALSE)

# Check the structure of the dataframe
str(metadata)

# Convert Enterotype to factor if not already
metadata$Enterotype <- as.factor(metadata$Enterotype)

# Convert Age to factor if needed
metadata$Age <- as.factor(metadata$Age)

# Define age time points
metadata$AgeTimePoint <- as.factor(metadata$Age)

# Define tissues
tissues <- unique(metadata$Tissue)

# Iterate over each tissue
for (tissue in tissues) {
  
  # Subset data for the current tissue
  tissue_data <- metadata %>% filter(Tissue == tissue)
  
  # Create a data frame to store transitions
  transitions <- data.frame(Source = character(), Target = character(), Value = integer(), stringsAsFactors = FALSE)
  
  # Calculate transitions for each pair of consecutive age time points
  age_time_points <- c("3", "4", "5", "6", "7", "8")
  for (i in 1:(length(age_time_points) - 1)) {
    age_time_point_1 <- age_time_points[i]
    age_time_point_2 <- age_time_points[i + 1]
    
    for (enterotype_1 in unique(tissue_data$Enterotype)) {
      for (enterotype_2 in unique(tissue_data$Enterotype)) {
        count <- nrow(tissue_data %>% filter(AgeTimePoint == age_time_point_1 & Enterotype == enterotype_1)) *
          nrow(tissue_data %>% filter(AgeTimePoint == age_time_point_2 & Enterotype == enterotype_2))
        if (count > 0) {
          transitions <- rbind(transitions, data.frame(Source = as.character(enterotype_1),
                                                       Target = as.character(enterotype_2),
                                                       Value = count))
        }
      }
    }
  }
  
  # Verify the structure of the transitions dataframe
  str(transitions)
  
  if (nrow(transitions) == 0) {
    next
  }
  
  # Create Transition Matrix
  CSTs <- unique(c(transitions$Source, transitions$Target))
  nstates <- length(CSTs)
  trans_matrix <- matrix(0, nrow = nstates, ncol = nstates)
  rownames(trans_matrix) <- CSTs
  colnames(trans_matrix) <- CSTs
  
  for (i in 1:nrow(transitions)) {
    trans_matrix[transitions$Source[i], transitions$Target[i]] <- transitions$Value[i]
  }
  
  # Normalize the transition matrix
  row_sums <- rowSums(trans_matrix, na.rm = TRUE)
  trans_matrix <- trans_matrix / row_sums
  trans_matrix[is.na(trans_matrix)] <- 0
  
  # Filter transitions greater than 0.2 and normalize again to ensure rows sum to 1
  trans_matrix <- ifelse(trans_matrix > 0.2, trans_matrix, 0)
  row_sums <- rowSums(trans_matrix, na.rm = TRUE)
  trans_matrix <- trans_matrix / row_sums
  trans_matrix[is.na(trans_matrix)] <- 0
  
  # Ensure states are character vectors
  CSTs <- as.character(CSTs)
  
  # Create Markov chain object
  mc <- new("markovchain", states = CSTs, transitionMatrix = trans_matrix, name = paste("Enterotype Transitions -", tissue))
  
  # Verify the Markov chain object
  print(mc)
  
  # Create igraph object for plotting
  graph <- graph_from_adjacency_matrix(trans_matrix, mode = "directed", weighted = TRUE)
  
  # Customize plot parameters
  # Calculate vertex sizes
  tissue_data$Enterotype <- as.character(tissue_data$Enterotype)
  vert.sz <- sapply(states(mc), function(x) nrow(tissue_data[tissue_data$Enterotype == x, , drop = FALSE]))
  
  # Set vertex sizes proportional to the log of the number of subjects, scaled
  V(graph)$size <- log(vert.sz + 1) * 15
  
  # Set vertex colors
  V(graph)$color <- c("#CDBE6A", "#89CFBE", "#86A4CF") # Adjusted to match the number of enterotypes
  
  # Set vertex labels to E1, E2, E3
  V(graph)$label <- paste("E", V(graph)$name, sep = "")
  
  # Set edge widths proportional to transition probabilities, with increased thickness
  E(graph)$width <- E(graph)$weight * 10
  
  # Set edge labels to transition probabilities, rounded to two decimal places
  E(graph)$label <- round(E(graph)$weight, 2)
  
  # Set vertex label color
  V(graph)$label.color <- "white"
  
  # Define plot dimensions and margin
  plot_width <- 990
  plot_height <- 800
  margin <- 50
  
  # Plot the Markov chain and save as PNG
  png(filename = paste0("/Users/kefe/Desktop/Data/Bacterial genera/DMN_Genera/Sanky plot/Markov_Chain_Plot_", tissue, ".png"), width = plot_width, height = plot_height)
  par(mar = c(5, 5, 5, 5))
  
  # Increase edge.arrow.size to make arrow heads more prominent
  plot(graph, 
       edge.arrow.size = 1.5,  # Increase the arrow size
       edge.curved = 0.2, 
       vertex.label.font = 2, 
       vertex.label.cex = 1.5, 
       vertex.frame.color = NA,
       layout = layout_in_circle)
  
  # Add legend with adjusted line widths
  legend("topright", 
         legend = c("0.2", "0.4", "0.6", "0.8"), 
         title = "Transition rate", 
         col = "gray", 
         lwd = c(4, 8, 12, 16),
         bty = "n")
  
  dev.off()
  
  # Save transition matrix as CSV
  write.csv(trans_matrix, file = paste0("/Users/kefe/Desktop/Data/Bacterial genera/DMN_Genera/Sanky plot/Transition_Matrix_", tissue, ".csv"), row.names = TRUE)
}
