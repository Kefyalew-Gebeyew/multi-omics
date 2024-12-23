
# Retrieve-counts

setwd("/Users/vivi/Desktop/colon")

sampleinfo <- read.table("/Users/vivi/Desktop/各种数据/sampleinfo.csv",sep = ",",header = T)
samplecolon <- sampleinfo[c(1:36),]
group <- paste(samplecolon$Group,samplecolon$Timepoint,sep = ".")
# Create the two variables
group <- as.character(group)
type <- sapply(strsplit(group, ".", fixed=T), function(x) x[1])
status <- sapply(strsplit(group, ".", fixed=T), function(x) x[2])
# Specify a design matrix with an intercept term
design <- model.matrix(~0+ status)

# retrieve counts
count_all <- read.table("/Users/vivi/Desktop/各种数据/oar_v31.trimmomatic.STAR.featurecounts.count.merge.txt",
                        sep = "\t",header = T)
rownames(count_all) <- count_all$X
count_all <- count_all[,-1]
count_colon <- count_all[,c(1:36)]
group <- factor(group)


#############__________________Normalization____________________________####

# Load required libraries
library(RUVSeq)
library(RColorBrewer)
library(DESeq2)

# Set working directory
setwd("/Users/vivi/Desktop/colon")

# Data Preprocessing
#-------------------

# Filter genes with at least 2 samples having counts > 5
filter_colon <- apply(count_colon, 1, function(x) length(x[x > 5]) >= 2)
filtered_colon <- count_colon[filter_colon,]

# Separate genes and spike-ins
genes_colon <- rownames(filtered_colon)[grep("^ENS", rownames(filtered_colon))]
spikes_colon <- rownames(filtered_colon)[grep("^ERCC", rownames(filtered_colon))]

# Create SeqExpressionSet object
x <- factor(group)
set <- newSeqExpressionSet(as.matrix(filtered_colon), 
                           phenoData = data.frame(x, row.names = colnames(filtered_colon)))

# Visualization before normalization
colors <- brewer.pal(6, "Set2")
pdf('before_normalization_RLE.pdf')
plotRLE(set, outline = FALSE, ylim = c(-4, 4), col = colors[x])
dev.off()

pdf('before_normalization_RUV_PCA.pdf')
plotPCA(set, col = colors[x], cex = 1.2)
dev.off()

# Normalization
#--------------

# Between-lane normalization
set <- betweenLaneNormalization(set, which = "upper")

# Prepare design matrix and perform differential expression analysis
design <- model.matrix(~x, data = pData(set))
y <- DGEList(counts = counts(set), group = x)
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type = "deviance")

# RUVr normalization
set4 <- RUVr(set, genes_colon, k = 1, res)

# Visualization after normalization
pdf('after_normalization_RLE.pdf')
plotRLE(set4, outline = FALSE, ylim = c(-4, 4), col = colors[x])
dev.off()

pdf('after_normalization_RUV_PCA.pdf')
plotPCA(set4, col = colors[x], cex = 1.2)
dev.off()

set4_colon <- normCounts(set4)

# Differential Expression Analysis
#---------------------------------

# Perform DESeq2 analysis
degobj_colon <- DESeq(dds_colon, betaPrior = TRUE, quiet = TRUE)

# Function to perform pairwise comparisons
perform_de_analysis <- function(condition1, condition2) {
  res <- results(degobj_colon, contrast = c("x", condition1, condition2), pAdjustMethod = 'BH')
  base1 <- get_base_mean(condition1)
  base2 <- get_base_mean(condition2)
  res_combined <- cbind(base1, base2, as.data.frame(res))
  res_combined <- cbind(ID = rownames(res_combined), res_combined)
  res_combined$baseMean <- rowMeans(cbind(base1, base2))
  res_combined$padj[is.na(res_combined$padj)] <- 1
  res_combined <- res_combined[order(res_combined$padj),]
  return(res_combined)
}

# Function to get base mean for a condition
get_base_mean <- function(condition) {
  base <- counts(degobj_colon, normalized = TRUE)[, colData(degobj_colon)$x == condition]
  if (is.vector(base)) {
    base_mean <- as.data.frame(base)
  } else {
    base_mean <- as.data.frame(rowMeans(base))
  }
  colnames(base_mean) <- condition
  return(base_mean)
}

# Perform pairwise comparisons
comparisons <- list(
  c("colon.3", "colon.4"),
  c("colon.3", "colon.5"),
  c("colon.3", "colon.6"),
  c("colon.3", "colon.7"),
  c("colon.3", "colon.8")
)

results_list <- lapply(comparisons, function(comp) perform_de_analysis(comp[1], comp[2]))

# Extract differentially expressed genes
extract_de_genes <- function(res, condition1, condition2) {
  de_genes <- subset(res, res$padj < 0.05, 
                     select = c('ID', condition1, condition2, 'log2FoldChange', 'padj'))
  up_regulated <- subset(de_genes, de_genes$log2FoldChange >= 1)
  down_regulated <- subset(de_genes, de_genes$log2FoldChange <= -1)
  return(list(de_genes = de_genes, up = up_regulated, down = down_regulated))
}

de_results <- lapply(seq_along(comparisons), function(i) {
  extract_de_genes(results_list[[i]], comparisons[[i]][1], comparisons[[i]][2])
})

# Write results to files
for (i in seq_along(comparisons)) {
  write.table(de_results[[i]]$de_genes, 
              file = paste0('res_colon_de_', comparisons[[i]][1], '_', comparisons[[i]][2], '.txt'), 
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# Combine all differentially expressed gene IDs
all_de_ids <- Reduce(union, lapply(de_results, function(x) c(x$up$ID, x$down$ID)))



####_______DEG relative to adjacent month_______####

setwd("/Users/vivi/Desktop/colon")

#################################################################
#resultsNames(degobj)
res_colonA_others <- DESeq2::results(degobj, 
                                     contrast = c(0, 0, 1, -1/5, -1/5, -1/5, -1/5, -1/5), 
                                     pAdjustMethod = 'BH')

res_colonA_others$padj[is.na(res_colonA_others$padj)] <- 1
res_colonA_others <- res_colonA_others[order(res_colonA_others$padj),]
res_colon_de_A_others <- subset(res_colonA_others, res_colonA_others$padj<0.05)
res_colon_de_A_others$ID <- rownames(res_colon_de_A_others)

res_colon_de_upa_others <- subset(res_colon_de_A_others, res_colon_de_A_others$log2FoldChange>=1)
res_colon_de_dwa_others <- subset(res_colon_de_A_others, res_colon_de_A_others$log2FoldChange<=(-1)*1)

#######################################

res_colonB_others <- DESeq2::results(degobj, 
                                     contrast = c(0, 0, -1/5, 1, -1/5, -1/5, -1/5, -1/5), 
                                     pAdjustMethod = 'BH')

res_colonB_others$padj[is.na(res_colonB_others$padj)] <- 1
res_colonB_others <- res_colonB_others[order(res_colonB_others$padj),]
res_colon_de_B_others <- subset(res_colonB_others, res_colonB_others$padj<0.05)
res_colon_de_B_others$ID <- rownames(res_colon_de_B_others)
# foldchange > 0.5
res_colon_de_upb_others <- subset(res_colon_de_B_others, res_colon_de_B_others$log2FoldChange>=1)
res_colon_de_dwb_others <- subset(res_colon_de_B_others, res_colon_de_B_others$log2FoldChange<=(-1)*1)

#######################################
res_colonC_others <- DESeq2::results(degobj, 
                                     contrast = c(0, 0, -1/5, -1/5, 1, -1/5, -1/5, -1/5), 
                                     pAdjustMethod = 'BH')

res_colonC_others$padj[is.na(res_colonC_others$padj)] <- 1
res_colonC_others <- res_colonC_others[order(res_colonC_others$padj),]
res_colon_de_C_others <- subset(res_colonC_others, res_colonC_others$padj<0.05)
res_colon_de_C_others$ID <- rownames(res_colon_de_C_others)
# foldchange > 0.5
res_colon_de_upc_others <- subset(res_colon_de_C_others, res_colon_de_C_others$log2FoldChange>=1)
res_colon_de_dwc_others <- subset(res_colon_de_C_others, res_colon_de_C_others$log2FoldChange<=(-1)*1)

#######################################
res_colonD_others <- DESeq2::results(degobj, 
                                     contrast = c(0, 0, -1/5, -1/5, -1/5, 1, -1/5, -1/5), 
                                     pAdjustMethod = 'BH')

res_colonD_others$padj[is.na(res_colonD_others$padj)] <- 1
res_colonD_others <- res_colonD_others[order(res_colonD_others$padj),]
res_colon_de_D_others <- subset(res_colonD_others, res_colonD_others$padj<0.05)
res_colon_de_D_others$ID <- rownames(res_colon_de_D_others)
# foldchange > 0.5
res_colon_de_upd_others <- subset(res_colon_de_D_others, res_colon_de_D_others$log2FoldChange>=1)
res_colon_de_dwd_others <- subset(res_colon_de_D_others, res_colon_de_D_others$log2FoldChange<=(-1)*1)

#######################################
#######################################
#######################################
res_colonE_others <- DESeq2::results(degobj, 
                                     contrast = c(0, 0, -1/5, -1/5, -1/5, -1/5, 1, -1/5), 
                                     pAdjustMethod = 'BH')

res_colonE_others$padj[is.na(res_colonE_others$padj)] <- 1
res_colonE_others <- res_colonE_others[order(res_colonE_others$padj),]
res_colon_de_E_others <- subset(res_colonE_others, res_colonE_others$padj<0.05)
res_colon_de_E_others$ID <- rownames(res_colon_de_E_others)
# foldchange > 0.5
res_colon_de_upe_others <- subset(res_colon_de_E_others, res_colon_de_E_others$log2FoldChange>=1)
res_colon_de_dwe_others <- subset(res_colon_de_E_others, res_colon_de_E_others$log2FoldChange<=(-1)*1)

#######################################
res_colonF_others <- DESeq2::results(degobj, 
                                     contrast = c(0, 0, -1/5, -1/5, -1/5, -1/5, -1/5, 1), 
                                     pAdjustMethod = 'BH')

res_colonF_others$padj[is.na(res_colonF_others$padj)] <- 1
res_colonF_others <- res_colonF_others[order(res_colonF_others$padj),]
res_colon_de_F_others <- subset(res_colonF_others, res_colonF_others$padj<0.05)
res_colon_de_F_others$ID <- rownames(res_colon_de_F_others)
# foldchange > 0.5
res_colon_de_upf_others <- subset(res_colon_de_F_others, res_colon_de_F_others$log2FoldChange>=1)
res_colon_de_dwf_others <- subset(res_colon_de_F_others, res_colon_de_F_others$log2FoldChange<=(-1)*1)

write.table(as.data.frame(res_colon_de_upa_others), file='res_colon_de_upa_others.txt', sep="\t", quote=F, row.names=F)
write.table(as.data.frame(res_colon_de_upb_others), file='res_colon_de_upb_others.txt', sep="\t", quote=F, row.names=F)
write.table(as.data.frame(res_colon_de_upc_others), file='res_colon_de_upc_others.txt', sep="\t", quote=F, row.names=F)
write.table(as.data.frame(res_colon_de_upd_others), file='res_colon_de_upd_others.txt', sep="\t", quote=F, row.names=F)
write.table(as.data.frame(res_colon_de_upe_others), file='res_colon_de_upe_others.txt', sep="\t", quote=F, row.names=F)
write.table(as.data.frame(res_colon_de_upf_others), file='res_colon_de_upf_others.txt', sep="\t", quote=F, row.names=F)

write.table(as.data.frame(res_colon_de_A_others), file='res_colon_de_A_others.txt', sep="\t", quote=F, row.names=F)
write.table(as.data.frame(res_colon_de_B_others), file='res_colon_de_B_others.txt', sep="\t", quote=F, row.names=F)
write.table(as.data.frame(res_colon_de_C_others), file='res_colon_de_C_others.txt', sep="\t", quote=F, row.names=F)
write.table(as.data.frame(res_colon_de_D_others), file='res_colon_de_D_others.txt', sep="\t", quote=F, row.names=F)
write.table(as.data.frame(res_colon_de_E_others), file='res_colon_de_E_others.txt', sep="\t", quote=F, row.names=F)
write.table(as.data.frame(res_colon_de_F_others), file='res_colon_de_F_others.txt', sep="\t", quote=F, row.names=F)



####_________generating expression patterns of age-correlated genes in four tissues__________###

# Load required libraries
library(ggplot2)

# Set working directory
setwd("/Users/vivi/Desktop/colon")

# Define function to create data frames for each tissue
create_tissue_df <- function(tissue, up_genes, down_genes) {
  data.frame(
    x_month = x_month,
    y_genes_up = up_genes,
    y_genes_down = -down_genes,
    tissue = tissue
  )
}

# Define month comparisons
x_month <- c("M4:3", "M5:3", "M6:3", "M7:3", "M8:3")

# Create data frames for consecutive month comparisons
consecutive_data <- rbind(
  create_tissue_df("Colon", 
                   c(nrow(res_colon_de_upab), nrow(res_colon_de_upac), nrow(res_colon_de_upad), nrow(res_colon_de_upae), nrow(res_colon_de_upaf)),
                   c(nrow(res_colon_de_dwab), nrow(res_colon_de_dwac), nrow(res_colon_de_dwad), nrow(res_colon_de_dwae), nrow(res_colon_de_dwaf))),
  create_tissue_df("Rumen", 
                   c(nrow(res_rumen_de_upab), nrow(res_rumen_de_upac), nrow(res_rumen_de_upad), nrow(res_rumen_de_upae), nrow(res_rumen_de_upaf)),
                   c(nrow(res_rumen_de_dwab), nrow(res_rumen_de_dwac), nrow(res_rumen_de_dwad), nrow(res_rumen_de_dwae), nrow(res_rumen_de_dwaf))),
  create_tissue_df("Liver", 
                   c(nrow(res_liver_de_upab), nrow(res_liver_de_upac), nrow(res_liver_de_upad), nrow(res_liver_de_upae), nrow(res_liver_de_upaf)),
                   c(nrow(res_liver_de_dwab), nrow(res_liver_de_dwac), nrow(res_liver_de_dwad), nrow(res_liver_de_dwae), nrow(res_liver_de_dwaf))),
  create_tissue_df("Jejunum", 
                   c(nrow(res_jejunum_de_upab), nrow(res_jejunum_de_upac), nrow(res_jejunum_de_upad), nrow(res_jejunum_de_upae), nrow(res_jejunum_de_upaf)),
                   c(nrow(res_jejunum_de_dwab), nrow(res_jejunum_de_dwac), nrow(res_jejunum_de_dwad), nrow(res_jejunum_de_dwae), nrow(res_jejunum_de_dwaf)))
)

# Plot consecutive month comparisons
pdf('consecutive_up_down_line_chart.pdf', width = 8, height = 6)
ggplot(consecutive_data, aes(x = x_month, group = tissue)) +
  geom_line(aes(y = y_genes_up, color = tissue), size = 1) +
  geom_point(aes(y = y_genes_up, color = tissue)) +
  geom_line(aes(y = y_genes_down, color = tissue), linetype = "dashed", size = 1) +
  geom_point(aes(y = y_genes_down, color = tissue), shape = 1) +
  scale_color_manual(values = c("#ff595e", "#ffca3a", "#8ac926", "#1982c4")) +
  theme_classic() +
  labs(x = "Month Comparison", y = "Number of Differentially Expressed Genes", 
       title = "Age-correlated Gene Expression Patterns (Consecutive Months)") +
  theme(legend.position = "right")
dev.off()

# Define month comparisons for adjacent months
x_month <- c("M4:3", "M5:4", "M6:5", "M7:6", "M8:7")

# Create data frames for adjacent month comparisons
adjacent_data <- rbind(
  create_tissue_df("Colon", 
                   c(nrow(res_colon_de_upab), nrow(res_colon_de_upbc), nrow(res_colon_de_upcd), nrow(res_colon_de_upde), nrow(res_colon_de_upef)),
                   c(nrow(res_colon_de_dwab), nrow(res_colon_de_dwbc), nrow(res_colon_de_dwcd), nrow(res_colon_de_dwde), nrow(res_colon_de_dwef))),
  create_tissue_df("Rumen", 
                   c(nrow(res_rumen_de_upab), nrow(res_rumen_de_upbc), nrow(res_rumen_de_upcd), nrow(res_rumen_de_upde), nrow(res_rumen_de_upef)),
                   c(nrow(res_rumen_de_dwab), nrow(res_rumen_de_dwbc), nrow(res_rumen_de_dwcd), nrow(res_rumen_de_dwde), nrow(res_rumen_de_dwef))),
  create_tissue_df("Liver", 
                   c(nrow(res_liver_de_upab), nrow(res_liver_de_upbc), nrow(res_liver_de_upcd), nrow(res_liver_de_upde), nrow(res_liver_de_upef)),
                   c(nrow(res_liver_de_dwab), nrow(res_liver_de_dwbc), nrow(res_liver_de_dwcd), nrow(res_liver_de_dwde), nrow(res_liver_de_dwef))),
  create_tissue_df("Jejunum", 
                   c(nrow(res_jejunum_de_upab), nrow(res_jejunum_de_upbc), nrow(res_jejunum_de_upcd), nrow(res_jejunum_de_upde), nrow(res_jejunum_de_upef)),
                   c(nrow(res_jejunum_de_dwab), nrow(res_jejunum_de_dwbc), nrow(res_jejunum_de_dwcd), nrow(res_jejunum_de_dwde), nrow(res_jejunum_de_dwef)))
)

# Plot adjacent month comparisons
pdf('adjacent_up_down_line_chart.pdf', width = 8, height = 6)
ggplot(adjacent_data, aes(x = x_month, group = tissue)) +
  geom_line(aes(y = y_genes_up, color = tissue), size = 1) +
  geom_point(aes(y = y_genes_up, color = tissue)) +
  geom_line(aes(y = y_genes_down, color = tissue), linetype = "dashed", size = 1) +
  geom_point(aes(y = y_genes_down, color = tissue), shape = 1) +
  scale_color_manual(values = c("#ff595e", "#ffca3a", "#8ac926", "#1982c4")) +
  theme_classic() +
  labs(x = "Month Comparison", y = "Number of Differentially Expressed Genes", 
       title = "Age-correlated Gene Expression Patterns (Adjacent Months)") +
  theme(legend.position = "right")
dev.off()

