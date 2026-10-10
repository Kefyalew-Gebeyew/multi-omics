
###---------------------------Separate age time point for each tissue------------------########---------------------######

# ---- Step 0. Configuration --------------------------------------------------------
otu_file      <- "OTU_kefe.csv"
meta_file     <- "SampleID.csv"
out_dir       <- "output"
do_rarefy     <- FALSE        # TRUE = rarefy to a common depth before analysis
rarefy_depth  <- NULL         # NULL = minimum sample depth; or a fixed number
alpha_level   <- 0.05         # Kruskal-Wallis threshold for running Dunn's test
seed          <- 123
tissue_levels <- c("Colon", "Jejunum", "Rumen")
age_cols      <- c("firebrick2", "chartreuse1", "hotpink4",
                   "maroon2", "aquamarine4", "cyan1")

# ---- 1. Packages -------------------------------------------------------------
pkgs <- c("dplyr", "vegan", "ggplot2", "dunn.test")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
invisible(lapply(pkgs, library, character.only = TRUE))

dir.create(out_dir, showWarnings = FALSE)
set.seed(seed)

# ---- 2. Load and validate data -----------------------------------------------
otu  <- read.csv(otu_file, row.names = 1, check.names = FALSE)
meta <- read.csv(meta_file)

stopifnot(
  "Metadata must contain SampleID, Age and Tissue" =
    all(c("SampleID", "Age", "Tissue") %in% names(meta)),
  "Some metadata samples are missing from the OTU table" =
    all(meta$SampleID %in% rownames(otu)),
  "OTU table contains missing values" = !anyNA(otu)
)
otu <- otu[meta$SampleID, ]   # align row order with metadata

# ---- 3. Optional rarefaction -------------------------------------------------
if (do_rarefy) {
  depth <- if (is.null(rarefy_depth)) min(rowSums(otu)) else rarefy_depth
  message("Rarefying to ", depth, " reads per sample")
  otu <- vegan::rrarefy(otu, sample = depth)
}

# ---- 4. Alpha diversity (one row per sample) ---------------------------------
alpha <- meta %>%
  mutate(
    Shannon  = as.numeric(diversity(otu, "shannon")),
    Simpson  = as.numeric(diversity(otu, "simpson")),   # 1 - D
    Richness = as.numeric(specnumber(otu)),
    Depth    = as.numeric(rowSums(otu)),
    Age      = as.numeric(Age),
    Tissue   = factor(Tissue, levels = tissue_levels)
  )
write.csv(alpha, file.path(out_dir, "alpha_diversity_per_sample.csv"),
          row.names = FALSE)

# ---- 5. Sequencing depth check -----------------------------------------------
depth_test <- cor.test(alpha$Depth, alpha$Shannon,
                       method = "spearman", exact = FALSE)
message(sprintf("Depth vs Shannon: Spearman rho = %.2f, p = %.2g",
                depth_test$estimate, depth_test$p.value))

# ---- 6. Statistics within each tissue ----------------------------------------
run_tests <- function(data, metric, tissue) {
  d  <- filter(data, Tissue == tissue)
  kw <- kruskal.test(d[[metric]] ~ factor(d$Age))
  kw_row <- data.frame(Metric = metric, Tissue = tissue,
                       chi_sq = unname(kw$statistic),
                       df     = unname(kw$parameter),
                       p      = kw$p.value)
  dunn_rows <- NULL
  if (kw$p.value < alpha_level) {
    invisible(capture.output(
      res <- dunn.test(d[[metric]], d$Age, method = "bh",
                       kw = FALSE, table = FALSE)))
    dunn_rows <- data.frame(Metric = metric, Tissue = tissue,
                            Comparison = res$comparisons,
                            Z = res$Z, p_adj_BH = res$P.adjusted)
  }
  list(kw = kw_row, dunn = dunn_rows)
}

combos  <- expand.grid(metric = c("Shannon", "Simpson"),
                       tissue = tissue_levels, stringsAsFactors = FALSE)
results <- Map(function(m, t) run_tests(alpha, m, t), combos$metric, combos$tissue)

kw_table   <- do.call(rbind, lapply(results, `[[`, "kw"))
dunn_table <- do.call(rbind, lapply(results, `[[`, "dunn"))
print(kw_table, row.names = FALSE)
write.csv(kw_table, file.path(out_dir, "kruskal_results.csv"), row.names = FALSE)
if (!is.null(dunn_table)) {
  print(dunn_table, row.names = FALSE)
  write.csv(dunn_table, file.path(out_dir, "dunn_results.csv"), row.names = FALSE)
}

# ---- 7. Plots ----------------------------------------------------------------
plot_alpha <- function(data, metric, ylim = NULL, breaks = waiver()) {
  ggplot(data, aes(factor(Age), .data[[metric]], color = factor(Age))) +
    geom_boxplot(outlier.shape = NA, fill = NA, width = 0.6, linewidth = 0.5) +
    geom_jitter(width = 0.15, size = 1.5) +
    facet_wrap(~Tissue) +
    scale_y_continuous(limits = ylim, breaks = breaks) +
    scale_color_manual(values = age_cols) +
    labs(title = paste(metric, "Diversity by Age and Tissue"),
         x = "Age", y = paste(metric, "Index"), color = "Age") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1),
          strip.background = element_rect(fill = "grey90", colour = "black"),
          strip.text       = element_text(size = 12, face = "bold"),
          panel.spacing    = unit(0.5, "lines"))
}

ggsave(file.path(out_dir, "shannon_plot.png"),
       plot_alpha(alpha, "Shannon", ylim = c(1, 7), breaks = 1:7),
       width = 12, height = 6, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "simpson_plot.png"),
       plot_alpha(alpha, "Simpson", ylim = c(0, 1), breaks = seq(0, 1, 0.2)),
       width = 12, height = 6, dpi = 300, bg = "white")

# ---- 8. Reproducibility ------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))








###---------------------------PCoA and PCA------------------########---------------------######

# ---- Step 0. Configuration --------------------------------------------------------
otu_file      <- "OTU_kefe.csv"
meta_file     <- "SampleID.csv"
out_dir       <- "output"
do_rarefy     <- FALSE      # TRUE = rarefy to a common depth before analysis
rarefy_depth  <- NULL       # NULL = minimum sample depth
n_perm        <- 999
seed          <- 123
tissue_levels <- c("Colon", "Jejunum", "Rumen")
color_palette <- c("firebrick2", "chartreuse1", "hotpink4",
                   "maroon2", "aquamarine4", "cyan1")
shape_palette <- c("Rumen" = 16, "Jejunum" = 17, "Colon" = 15)

# ---- 1. Packages -------------------------------------------------------------
pkgs <- c("vegan", "ggplot2")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
invisible(lapply(pkgs, library, character.only = TRUE))
dir.create(out_dir, showWarnings = FALSE)
set.seed(seed)

# ---- 2. Load and validate data -----------------------------------------------
otu  <- read.csv(otu_file, row.names = 1, check.names = FALSE)
meta <- read.csv(meta_file)
stopifnot(
  all(c("SampleID", "Age", "Tissue") %in% names(meta)),
  all(meta$SampleID %in% rownames(otu)),
  !anyNA(otu)
)
otu  <- otu[meta$SampleID, ]                 # same order as metadata
meta$Age    <- factor(meta$Age)
meta$Tissue <- factor(meta$Tissue, levels = tissue_levels)
rownames(meta) <- meta$SampleID

# ---- 3. Optional rarefaction -------------------------------------------------
if (do_rarefy) {
  depth <- if (is.null(rarefy_depth)) min(rowSums(otu)) else rarefy_depth
  message("Rarefying to ", depth, " reads per sample")
  otu <- rrarefy(otu, sample = depth)
}

# ---- 4. Distances and ordinations --------------------------------------------
# Bray-Curtis on relative abundances (removes library-size effects)
rel_abund <- decostand(otu, method = "total")
bray      <- vegdist(rel_abund, method = "bray")

# PCoA
pcoa     <- cmdscale(bray, eig = TRUE, k = 2)
pos_eig  <- pcoa$eig[pcoa$eig > 0]
pcoa_var <- pcoa$eig[1:2] / sum(pos_eig) * 100
pcoa_df  <- data.frame(SampleID = rownames(pcoa$points),
                       PCoA1 = pcoa$points[, 1], PCoA2 = pcoa$points[, 2])
pcoa_df  <- merge(pcoa_df, meta, by = "SampleID")

# PCA on Hellinger-transformed data (OTUs with zero variance removed)
hell     <- decostand(otu, method = "hellinger")
hell     <- hell[, apply(hell, 2, var) > 0]
pca      <- prcomp(hell, scale. = FALSE)
pca_var  <- summary(pca)$importance[2, 1:2] * 100
pca_df   <- data.frame(SampleID = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
pca_df   <- merge(pca_df, meta, by = "SampleID")

write.csv(pcoa_df, file.path(out_dir, "pcoa_coordinates.csv"), row.names = FALSE)
write.csv(pca_df,  file.path(out_dir, "pca_coordinates.csv"),  row.names = FALSE)

# ---- 5. Plots ----------------------------------------------------------------
plot_ordination <- function(df, x, y, xlab, ylab, title) {
  ggplot(df, aes(.data[[x]], .data[[y]], color = Age, shape = Tissue)) +
    geom_point(size = 3) +
    labs(title = title, x = xlab, y = ylab) +
    scale_color_manual(values = color_palette) +
    scale_shape_manual(values = shape_palette) +
    theme_bw() +
    theme(panel.grid   = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
}

pcoa_plot <- plot_ordination(pcoa_df, "PCoA1", "PCoA2",
                             sprintf("PCoA1 (%.1f%%)", pcoa_var[1]),
                             sprintf("PCoA2 (%.1f%%)", pcoa_var[2]),
                             "PCoA of Microbial Communities (Bray-Curtis)")
pca_plot  <- plot_ordination(pca_df, "PC1", "PC2",
                             sprintf("PC1 (%.1f%%)", pca_var[1]),
                             sprintf("PC2 (%.1f%%)", pca_var[2]),
                             "PCA of Microbial Communities (Hellinger)")

ggsave(file.path(out_dir, "pcoa_plot.png"), pcoa_plot, width = 9, height = 6, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "pca_plot.png"),  pca_plot,  width = 9, height = 6, dpi = 300, bg = "white")

# ---- 6. PERMANOVA and dispersion tests ---------------------------------------
perm <- adonis2(bray ~ Age * Tissue, data = meta, permutations = n_perm, by = "terms")
print(perm)
write.csv(as.data.frame(perm), file.path(out_dir, "permanova_results.csv"))

# PERMANOVA is sensitive to unequal dispersion, so test it
disp_tissue <- betadisper(bray, meta$Tissue)
disp_age    <- betadisper(bray, meta$Age)
disp_tab <- rbind(
  data.frame(Factor = "Tissue", as.data.frame(permutest(disp_tissue, permutations = n_perm)$tab)[1, ]),
  data.frame(Factor = "Age",    as.data.frame(permutest(disp_age,    permutations = n_perm)$tab)[1, ]))
print(disp_tab)
write.csv(disp_tab, file.path(out_dir, "betadisper_results.csv"), row.names = FALSE)

# Age effect within each tissue
within_tissue <- do.call(rbind, lapply(tissue_levels, function(t) {
  idx <- meta$Tissue == t
  res <- adonis2(as.dist(as.matrix(bray)[idx, idx]) ~ Age,
                 data = meta[idx, ], permutations = n_perm)
  data.frame(Tissue = t, R2 = res$R2[1], F = res$F[1], p = res$`Pr(>F)`[1])
}))
within_tissue$p_adj_BH <- p.adjust(within_tissue$p, method = "BH")
print(within_tissue)
write.csv(within_tissue, file.path(out_dir, "permanova_within_tissue.csv"), row.names = FALSE)

# ---- 7. Reproducibility ------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
