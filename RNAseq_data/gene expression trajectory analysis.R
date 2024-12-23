# Load required libraries
library(Mfuzz)
library(reshape2)
library(ggplot2)
library(circlize)
library(clusterProfiler)
library(AnnotationDbi)
library(ggthemes)

# Set working directory
setwd("/Users/vivi/Desktop/all")

# Data Preprocessing
#-------------------

# Process gene expression data for four tissues: colon, jejunum, liver, and rumen
process_tissue_data <- function(data) {
  data <- data.frame(data)
  data[, 37] <- apply(data[, 1:36], 1, mean) # Calculate mean expression across samples
  data <- data[order(data[, 37], decreasing = TRUE), ] # Sort by mean expression
  return(data)
}

assay_dds_colon_2 <- process_tissue_data(assay_dds_colon)
assay_dds_jejunum_2 <- process_tissue_data(assay_dds_jejunum)
assay_dds_liver_2 <- process_tissue_data(assay_dds_liver)
assay_dds_rumen_2 <- process_tissue_data(assay_dds_rumen)

# Identify common genes across all tissues (top 15,000 genes by expression)
loess_genes <- Reduce(intersect, list(
  rownames(assay_dds_colon_2[1:15000,]),
  rownames(assay_dds_jejunum_2[1:15000,]),
  rownames(assay_dds_liver_2[1:15000,]),
  rownames(assay_dds_rumen_2[1:15000,])
))

# Normalize and scale data for all tissues
assay_dds_all_4 <- data.frame(scale(data.frame(assay_dds_all)[loess_genes,]))

# Fit LOESS regression to each gene using median expression per age group
col_num <- 0
for (i in 1:24) {
  assay_dds_all_4[, 144 + i] <- apply(assay_dds_all_4[, (1 + col_num):(i * 6)], 1, median)
  col_num <- col_num + 6
}

# Estimate whole-organism trajectory per gene by averaging trajectories across tissues
gene_mean <- as.data.frame(matrix(0, nrow = length(loess_genes), ncol = 6))
colnames(gene_mean) <- c("m3", "m4", "m5", "m6", "m7", "m8")
rownames(gene_mean) <- loess_genes

for (i in loess_genes) {
  gene_mean[i, ] <- apply(
    assay_dds_all_4[i, c(145:150)], 
    MARGIN = 2,
    FUN = mean
  )
}

# Clustering Analysis with Mfuzz
#-------------------------------

df_ddgs_3a <- as.matrix(gene_mean)
eset_ddgs <- ExpressionSet(df_ddgs_3a)

# Filter missing values and standardize the data
eset_ddgs <- filter.NA(eset_ddgs, thres = 0.25)
eset_ddgs <- standardise(eset_ddgs)

# Perform Mfuzz clustering
set.seed(2022)
c_ddgs <- 9 # Number of clusters
m <- mestimate(eset_ddgs)
cl_ddgs <- mfuzz(eset_ddgs, c = c_ddgs, m = m)

# Plot cluster overlap
O_ddgs <- overlap(cl_ddgs)
par(mfrow = c(1, 1))
overlap.plot(cl_ddgs, over = O_ddgs, thres = 0.05)

# Visualize clusters with Mfuzz plots
mfuzz.plot2(
  eset_ddgs,
  cl = cl_ddgs,
  mfrow = c(3,3),
  centre = TRUE,
  time.labels = c("3m", "4m", "5m", "6m", "7m", "8m"),
  x11 = FALSE
)

# Analyze Clustering Results
#---------------------------

# Extract cluster sizes and membership information
cluster_size <- cl_ddgs$size
names(cluster_size) <- seq_len(c_ddgs)

raw_cluster_anno <- cbind(df_ddgs_3a, cluster = cl_ddgs$cluster)

# Prepare normalized and membership matrices for visualization
norm_cluster_anno <- cbind(eset_ddgs@assayData$exprs, cluster = cl_ddgs$cluster)

membership_info <- cl_ddgs$membership %>%
  as.data.frame() %>%
  mutate(gene = rownames(.))

dnorm <- norm_cluster_anno %>%
  as.data.frame() %>%
  mutate(gene = rownames(.))

final_res <- merge(dnorm, membership_info, by = 'gene')

df_long <- melt(
  final_res,
  id.vars = c('cluster', 'gene', 'membership'),
  variable.name = 'time_point',
  value.name = 'norm_value'
)

df_long$cluster_name <- paste('Cluster', df_long$cluster)

# Visualize gene expression changes across clusters using ggplot2
pdf('mfuzz.pdf', width = 10)
ggplot(df_long, aes(x = time_point, y = norm_value)) +
  geom_line(aes(color = membership, group = gene), size = 0.5) +
  theme_classic() +
  facet_wrap(~cluster_name) +
  ylab('Expression Changes') +
  xlab('Age [months]')
dev.off()

# Functional Enrichment Analysis (GO and KEGG Pathways)
#------------------------------------------------------

sheep_db <- loadDb("/Users/vivi/Desktop/colon/sheep.sqlite")

for (i in seq_len(c_ddgs)) {
  
  genes_in_cluster <- raw_cluster_anno[raw_cluster_anno$cluster == i, ]
  
  # Perform GO enrichment analysis
  enrich_go_results <- enrichGO(
    gene = genes_in_cluster,
    OrgDb = sheep_db,
    keyType = "ENSEMBL",
    ont = "ALL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  
  write.table(as.data.frame(enrich_go_results), paste0("cluster_", i,".txt"))
}

####______________calculating and visualizing the amplitude and variability index____________#####

library(ggplot2)
library(ggrepel)
library(ggthemes)

setwd("/Users/vivi/Desktop/all")

organ_cluster <- data.frame(organ = rep(NA, c_ddgs*4),
                            cluster = rep(NA, c_ddgs*4),
                            Amp = rep(NA, c_ddgs*4),
                            Var = rep(NA, c_ddgs*4))
for (i in 1:c_ddgs) {
  genes <- unique(df_2[which(df_2$cluster == i), "gene"])
  gene_median <- assay_dds_all_4[genes, 145:150]
  df_ddgs_3 <- gene_median
  df_ddgs_3a <- as.matrix(df_ddgs_3)
  eset_ddgs <- ExpressionSet(df_ddgs_3a)
  eset_ddgs <- filter.NA(eset_ddgs, thres = 0.25)
  eset_ddgs <- standardise(eset_ddgs)
  
  set.seed(2022)
  m <- mestimate(eset_ddgs)
  organ_membership <- data.frame(t(membership(as.matrix(eset_ddgs@assayData[["exprs"]]),
                                              clusters = cl_ddgs[[1]], m = m)))
  organ_cluster[i, "Var"] <- apply(organ_membership[i,], 1, mean)
  organ_cluster[i, "Amp"] <- sd(abs(apply(data.frame(t(eset_ddgs@assayData[["exprs"]])),1,mean)))
  organ_cluster[i, "organ"] <- "Colon"
  organ_cluster[i, "cluster"] <- i
}
for (i in 1:c_ddgs) {
  genes <- unique(df_2[which(df_2$cluster == i), "gene"])
  gene_median <- assay_dds_all_4[genes, 151:156]
  df_ddgs_3 <- gene_median
  df_ddgs_3a <- as.matrix(df_ddgs_3)
  eset_ddgs <- ExpressionSet(df_ddgs_3a)
  eset_ddgs <- filter.NA(eset_ddgs, thres = 0.25)
  eset_ddgs <- standardise(eset_ddgs)
  
  set.seed(2022)
  m <- mestimate(eset_ddgs)
  organ_membership <- data.frame(t(membership(as.matrix(eset_ddgs@assayData[["exprs"]]),
                                              clusters = cl_ddgs[[1]], m = m)))
  organ_cluster[i+9, "Var"] <- apply(organ_membership[i,], 1, mean)
  organ_cluster[i+9, "Amp"] <- sd(abs(apply(data.frame(t(eset_ddgs@assayData[["exprs"]])),1,mean)))
  organ_cluster[i+9, "organ"] <- "Jejunum"
  organ_cluster[i+9, "cluster"] <- i
  
}
for (i in 1:c_ddgs) {
  genes <- unique(df_2[which(df_2$cluster == i), "gene"])
  gene_median <- assay_dds_all_4[genes, 157:162]
  df_ddgs_3 <- gene_median
  df_ddgs_3a <- as.matrix(df_ddgs_3)
  eset_ddgs <- ExpressionSet(df_ddgs_3a)
  eset_ddgs <- filter.NA(eset_ddgs, thres = 0.25)
  eset_ddgs <- standardise(eset_ddgs)
  
  set.seed(2022)
  m <- mestimate(eset_ddgs)
  organ_membership <- data.frame(t(membership(as.matrix(eset_ddgs@assayData[["exprs"]]),
                                              clusters = cl_ddgs[[1]], m = m)))
  organ_cluster[i+18, "Var"] <- apply(organ_membership[i,], 1, mean)
  organ_cluster[i+18, "Amp"] <- sd(abs(apply(data.frame(t(eset_ddgs@assayData[["exprs"]])),1,mean)))
  organ_cluster[i+18, "organ"] <- "Liver"
  organ_cluster[i+18, "cluster"] <- i
  
}
for (i in 1:c_ddgs) {
  genes <- unique(df_2[which(df_2$cluster == i), "gene"])
  gene_median <- assay_dds_all_4[genes, 163:168]
  df_ddgs_3 <- gene_median
  df_ddgs_3a <- as.matrix(df_ddgs_3)
  eset_ddgs <- ExpressionSet(df_ddgs_3a)
  eset_ddgs <- filter.NA(eset_ddgs, thres = 0.25)
  eset_ddgs <- standardise(eset_ddgs)
  
  set.seed(2022)
  m <- mestimate(eset_ddgs)
  organ_membership <- data.frame(t(membership(as.matrix(eset_ddgs@assayData[["exprs"]]),
                                              clusters = cl_ddgs[[1]], m = m)))
  organ_cluster[i+27, "Var"] <- apply(organ_membership[i,], 1, mean)
  organ_cluster[i+27, "Amp"] <- sd(abs(apply(data.frame(t(eset_ddgs@assayData[["exprs"]])),1,mean)))
  organ_cluster[i+27, "organ"] <- "Rumen"
  organ_cluster[i+27, "cluster"] <- i
  
}

c_data <- data.frame(cl_ddgs[[1]])
c_data$cluster <- c(8,3,6,7,9,2,1,4,5)
amp_via <- data.frame(Amplitude = rep(NA, 9), Variability = rep(NA, 9))
for (i in 1:c_ddgs) {
  sd_data <- sd(abs(c_data[which(c_data$cluster == i),1:6]))
  mean_mem <- mean(organ_cluster[which(organ_cluster$cluster == i),"Var"])
  amp_via[i, "Amplitude"] <- sd_data
  amp_via[i, "Variability"] <- mean_mem
}
amp_via$name <- paste("cluster", 1:9, sep = " ")

##########################################################################################
pdf("amp_var.pdf")
ggplot(data = amp_via, aes(x = Variability, y = Amplitude)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = name),
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.5, "lines"),
                  segment.color = 'grey50',
                  size = 7) +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        axis.ticks = element_line(size = 1))
dev.off()
##########################################################################################
organ_cluster$cluster_name <- paste("cluster", organ_cluster$cluster, sep = " ")
pdf("organ_cluster.pdf", width = 12, height = 5)
ggplot(data = organ_cluster, aes(x = Var, y = Amp, color = organ, group = cluster)) +
  geom_point(size = 3.5) + 
  theme_few() +
  facet_wrap(~cluster_name ,ncol = 5, scales = 'free') +
  ylab('Amplitude') + xlab('Variability') +
  theme(text = element_text(size = 12, face = "bold"),
        legend.position = c(0.9, 0.2),
        title = element_text(size = 15, face = "bold"),
        strip.text.x = element_text(size = 15, face = "bold"),
        plot.margin = margin(10,10,10,10, "pt"))
dev.off()
##########################################################################################

