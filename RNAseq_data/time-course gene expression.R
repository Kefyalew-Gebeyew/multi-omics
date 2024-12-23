####_____time-course gene expression data analysis using maSigPro package_____####

library(maSigPro)

setwd("/Users/vivi/Desktop/colon")

sample.group <- read.table("colon_group_masigpro.csv",sep = ',',header = T)
rownames(sample.group) <- sample.group$X
sample.group <- sample.group[,-1]

design <- make.design.matrix(
  sample.group,
  degree = length(unique(sample.group$Time)) - 1)

# Removing rows with all zeros:
count_colon.norm_2 <- assay_dds_colon
count_colon.norm_2$sd <- apply(assay_dds_colon,1,sd)
count_colon.norm_2$mean <- apply(assay_dds_colon,1,mean)
rsds <- count_colon.norm_2$sd/count_colon.norm_2$mean
rsds0.2 <- names(which(rsds > 0.2))

count_colon.norm <- assay_dds_colon[rsds0.2,]
sumatot <- apply(count_colon.norm, 1, sum)
counts0 <- which(sumatot == 0)
if (length(counts0) > 0) { count_colon.norm <- count_colon.norm[-counts0,] }

# regression
fit <- p.vector(
  count_colon.norm,
  design,
  counts = T,
  #theta = theta,
  Q = 0.05,
  MT.adjust = "BH",
  min.obs = 20
)
#dim(fit$SELEC)
#save(fit,file = "fit.RData")

# step regresion
tstep <- T.fit(
  fit,
  step.method = "backward",
  alfa = 0.05)
#save(tstep,file = "tstep.RData")

# get siggenes; default rsq = 0.7
sigs <- get.siggenes(
  tstep,
  rsq = 0.7,
  vars = "groups"
)

pdf('masigpro-clusters_of_rsp_0.7.pdf')
see.genes(sigs$sig.genes$colon,
          show.fit = T, 
          dis =design$dis, 
          #cluster.method="hclust",
          cluster.data = 1, 
          k = 6,
          #k.mclust = T,
          newX11 = F)
dev.off()


#####__________ perform soft clustering on time-course gene expression data using Mfuzz package___######

# Load required libraries
library(Mfuzz)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(circlize)

# Set working directory
setwd("/Users/vivi/Desktop/colon")

# Data Preparation
#-----------------

# Extract significant genes from previous analysis
symbols <- sigs$summary$`summary[apply(summary, 1, not.all.empty), ]`

# Calculate mean expression values for each group
df_ddgs_1 <- t(assay_dds_colon[symbols,])
df_ddgs_2 <- aggregate(df_ddgs_1, by = list(group), mean, na.rm = TRUE)
rownames(df_ddgs_2) <- df_ddgs_2[,1]
df_ddgs_3 <- data.frame(t(df_ddgs_2[,-1]))
df_ddgs_3a <- as.matrix(df_ddgs_3)

# Create ExpressionSet object and preprocess data
eset_ddgs <- ExpressionSet(df_ddgs_3a)
eset_ddgs <- filter.NA(eset_ddgs, thres = 0.25)
eset_ddgs <- fill.NA(eset_ddgs, mode = 'mean')
eset_ddgs <- standardise(eset_ddgs)

# Mfuzz Clustering
#-----------------

# Set number of clusters and seed for reproducibility
c_ddgs <- 6
set.seed(2022)

# Estimate fuzzification parameter and perform clustering
m <- mestimate(eset_ddgs)
cl_ddgs <- mfuzz(eset_ddgs, c = c_ddgs, m = m)

# Analyze cluster overlap
O_ddgs <- overlap(cl_ddgs)

# Visualization
#--------------

# Plot cluster overlap
par(mfrow = c(1,1))
overlap.plot(cl_ddgs, over = O_ddgs, thres = 0.05)

# Create Mfuzz plots
color <- colorRampPalette(rev(c("#E02401", "Yellow","#0F52BA", "#3E7C17")))(1000)
mfuzz.plot2(eset_ddgs,
            cl_ddgs,
            mfrow = c(2,3),
            centre = TRUE,
            time.labels = unique(group),
            x11 = FALSE,
            colo = color)

# Analysis of Clustering Results
#-------------------------------

# Display cluster sizes
cluster_size <- cl_ddgs$size
names(cluster_size) <- 1:c_ddgs
print(cluster_size)

# Combine clustering results with original and normalized data
raw_cluster_anno <- cbind(df_ddgs_3, cluster = cl_ddgs$cluster)
norm_cluster_anno <- cbind(eset_ddgs@assayData$exprs, cluster = cl_ddgs$cluster)

# Save results
write.csv(raw_cluster_anno, file = 'mean_mfuzz_cluster_anno.csv', row.names = TRUE)
write.csv(norm_cluster_anno, file = 'mean_mfuzz_norm_cluster_anno.csv', row.names = TRUE)

# Custom ggplot2 Visualization
#-----------------------------

# Prepare data for plotting
mem <- cbind(cl_ddgs$membership, cluster2 = cl_ddgs$cluster) %>% 
  as.data.frame() %>% 
  mutate(gene = rownames(.))

membership_info <- lapply(1:6, function(x){
  ms <- mem %>% filter(cluster2 == x)
  data.frame(membership = ms[[x]], gene = ms$gene, cluster2 = ms$cluster2)
}) %>% do.call('rbind', .)

dnorm <- cbind(eset_ddgs@assayData$exprs, cluster = cl_ddgs$cluster) %>% 
  as.data.frame() %>% 
  mutate(gene = rownames(.))

final_res <- merge(dnorm, membership_info, by = 'gene')
final_res <- final_res[,-10]

# Convert to long format
df <- melt(final_res, id.vars = c('cluster', 'gene', 'membership'),
           variable.name = 'cell_type', value.name = 'norm_value')
df$cluster_name <- paste('cluster ', df$cluster, sep = '')
df$cluster_name <- factor(df$cluster_name, levels = paste('cluster ', 1:6, sep = ''))

# Generate plot
pdf('mfuzz.pdf', width = 10)
ggplot(df, aes(x = cell_type, y = norm_value)) +
  geom_line(aes(color = membership, group = gene), size = 1) +
  geom_line(stat = "summary", fun = "median", colour = "brown", size = 2, aes(group = 1)) +
  theme_bw(base_size = 18) +
  ylab('Expression Changes') + xlab('') +
  theme(axis.ticks.length = unit(0.25, 'cm'),
        axis.text.x = element_text(angle = 45, hjust = 1, color = 'black')) +
  facet_wrap(~cluster_name, ncol = 3, scales = 'free')
dev.off()

###________________ Functional Enrichment Analysis______________________________###

for (i in 1:c_ddgs) {
  genes <- rownames(raw_cluster_anno[which(raw_cluster_anno$cluster == i),])
  
  # GO Enrichment
  enrich.go <- enrichGO(gene = genes,
                        OrgDb = sheep,
                        keyType = "ENSEMBL",
                        ont = 'ALL',
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.01,
                        readable = FALSE)
  
  BP <- enrich.go[enrich.go$ONTOLOGY == 'BP', ]
  CC <- enrich.go[enrich.go$ONTOLOGY == 'CC', ]
  MF <- enrich.go[enrich.go$ONTOLOGY == 'MF', ]
  
  write.table(as.data.frame(BP), paste('mfuzz', i, '_go.BP.txt', sep = ""), sep = '\t', row.names = FALSE, quote = FALSE)
  write.table(as.data.frame(CC), paste('mfuzz', i, '_go.CC.txt', sep = ""), sep = '\t', row.names = FALSE, quote = FALSE)
  write.table(as.data.frame(MF), paste('mfuzz', i, '_go.MF.txt', sep = ""), sep = '\t', row.names = FALSE, quote = FALSE)
  
  # KEGG Enrichment
  L_ids_kegg <- data.frame(na.omit(ls_31[genes,]$entrezgene_id))
  colnames(L_ids_kegg) <- 'entrez_id'
  
  kegg <- enrichKEGG(gene = L_ids_kegg$entrez_id,
                     keyType = 'kegg',
                     organism = "oas",
                     pAdjustMethod = 'BH',
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.01)
  
  write.table(kegg, paste('mfuzz', i, '_go.kegg.txt', sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
}
