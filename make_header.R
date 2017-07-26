# Childhood Cancer Data Lab 2017
# make_header.R
#
# Usage:
#     Rscript make_header.R
#
# Output:
#     Generate gene expression heatmap header from TARGET data

library(ggplot2)
library(dplyr)
library(magick)

set.seed(123)

# Load data
data_df <- readr::read_tsv('target_RSEM_Hugo_norm_count')
phenotype_df <- readr::read_tsv('TARGET_phenotype')

mad_genes <- apply(data_df[, 2:ncol(data_df)], 1, mad)
names(mad_genes) <- data_df$sample
mad_genes <- sort(mad_genes, decreasing = TRUE)

# Take top genes and random samples
num_genes <- 500
num_samples <- 100
use_genes <- names(mad_genes)[1:num_genes]
random_samples <- sample(2:ncol(data_df), num_samples)

subset_data_df <- data_df[data_df$sample %in% use_genes, c(1, random_samples)]

# Find order of hiearchical clustering of samples and genes
gene_order <- hclust(dist(subset_data_df[, -1], method = 'euclidean'))$order
sample_order <- hclust(dist(t(subset_data_df[, -1]), method = 'euclidean'))$order

gene_order <- subset_data_df$sample[gene_order]
sample_order <- colnames(subset_data_df)[2:ncol(subset_data_df)][sample_order]

subset_data_long_df <- reshape2::melt(subset_data_df)
colnames(subset_data_long_df) <- c('gene', 'sample', 'value')

subset_data_long_df$gene <- factor(subset_data_long_df$gene,
                                   levels = rev(gene_order))
subset_data_long_df$sample <- factor(subset_data_long_df$sample,
                                     levels = sample_order)

ggplot(subset_data_long_df, aes(x = sample, y = gene)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = 'white', high = '#71B9F6') +
  xlab("") + ylab("") + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") + 
  ggsave('real_data_heatmap.png', width = 8, heigh = 3)

# What is the phenotype data
phenotype_sub_df <- phenotype_df %>% dplyr::filter(sample_id %in% sample_order)

table(phenotype_sub_df$`_primary_disease`)

# Combine background heatmap with header base
header_raw <- magick::image_read('cures_transparency.png')
header <- header_raw %>% magick::image_scale("1000")
background <- magick::image_read('real_data_heatmap.png')
composite <- magick::image_composite(background, header, offset = "+1400+390")

magick::image_write(composite, "composite_header.png")
