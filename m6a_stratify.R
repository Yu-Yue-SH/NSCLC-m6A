rm(list = ls())


# load packages
library(Seurat)
library(ggplot2)
# library(dplyr)
# library(tidyverse)
# library(scCustomize)
library(wesanderson)


# load data
NSCLC <- readRDS('outputs/NSCLC_new_umap.rds')


# m6a
# 1	Jiang, X. et al. The role of m6A modification in the biological functions
# and diseases. Signal Transduct Target Ther 6, 74 (2021).
# https://doi.org/10.1038/s41392-020-00450-x
# *2	Fang, Z. et al. Role of m6A writers, erasers and readers in cancer. Exp
# Hematol Oncol 11, 45 (2022). https://doi.org/10.1186/s40164-022-00298-7
# 3	Zaccara, S., Ries, R. J. & Jaffrey, S. R. Reading, writing and erasing mRNA
# methylation. Nat Rev Mol Cell Biol 20, 608-624 (2019).
# https://doi.org/10.1038/s41580-019-0168-5
writers <- c('METTL3', 'METTL5', 'METTL14', 'METTL16', 'WTAP', 'RBM15',
             'RBM15B', 'ZC3H13', 'VIRMA', 'CBLL1')
# ZNF217 not found, VIRMA == KIAA1429, CBLL1 == HAKAI
erasers <- c('FTO', 'ALKBH5')
readers <- c('YTHDF1', 'YTHDF2', 'YTHDF3', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3',
             'YTHDC1', 'YTHDC2', 'ELAVL1', 'HNRNPC', 'RBMX', 'HNRNPA2B1')
# RBMX==HNRNPG
# all HNRNPs, while only important genes are selected
# rownames(NSCLC)[grep('^HNRNP', rownames(NSCLC))]
# EIF3s are not included
m6a <- c(writers, erasers, readers)


# subset cancer samples and cancer cells ####
cancer_samples <- subset(
  NSCLC,
  platform %in% c('10x', 'Singleron') &
    origin %in% c('tumor_matastasis', 'tumor_primary')
)
cancer_cells <- subset(cancer_samples, ann_fine == 'Cancer cells')

# saveRDS(cancer_samples, 'outputs/NSCLC_cancer_samples.rds')
# saveRDS(cancer_cells, 'outputs/NSCLC_cancer_cells.rds')
cancer_samples <- readRDS('outputs/NSCLC_cancer_samples.rds')
cancer_cells <- readRDS('outputs/NSCLC_cancer_cells.rds')


# extract data ####
data_cancer_cells <- as.data.frame(GetAssayData(cancer_cells, layer = 'data'))
data_cancer_cells <- data_cancer_cells[m6a, ]
data_cancer_cells <- t(data_cancer_cells)
data_cancer_cells <- cbind(data_cancer_cells, cancer_cells@meta.data)

# saveRDS(data_cancer_cells, 'outputs/data_cancer_cells.rds')
data_cancer_cells <- readRDS('outputs/data_cancer_cells.rds')


# stratify m6a levels ####
for (i in m6a) {
  data_cancer_cells[paste0(i, '_level')] <- NA
  data_cancer_cells[paste0(i, '_level')][data_cancer_cells[[i]] == 0, ] <-
    'Not detected'
  for (j in c('10x', 'Singleron')) {
    cells <- data_cancer_cells[[i]] != 0 & data_cancer_cells[['platform']] == j
    data_cancer_cells[cells, paste0(i, '_level')] <- as.character(
      cut(
        data_cancer_cells[cells, i],
        breaks = quantile(data_cancer_cells[cells, i], c(0, 0.4, 0.6, 1)),
        labels = c('Low', 'Medium', 'High'),
        include.lowest = TRUE
      )
    )
  }
}

# saveRDS(data_cancer_cells, 'outputs/data_cancer_cells_stratified.rds')
data_cancer_cells <- readRDS('outputs/data_cancer_cells_stratified.rds')

# add to seurat object
cancer_samples[['new_broad_type']] <- cancer_samples[['ann_coarse']]
cancer_samples@meta.data[['new_broad_type']] <-
  as.character(cancer_samples@meta.data[['new_broad_type']])
cancer_samples[['new_broad_type']][cancer_samples[['ann_fine']] ==
                                          'Cancer cells'] <- 'Cancer cells'

for (i in m6a) {
  cancer_cells[[paste0(i, '_level')]] <- NA
  cancer_cells[[paste0(i, '_level')]]<- data_cancer_cells[[paste0(i, '_level')]]
  cancer_samples[[paste0(i, '_level')]] <- cancer_samples[['new_broad_type']]
  cancer_samples[[paste0(i, '_level')]] <- cancer_cells[[paste0(i, '_level')]]
}

# saveRDS(cancer_samples, 'outputs/NSCLC_cancer_samples_stratified.rds')
# saveRDS(cancer_cells, 'outputs/NSCLC_cancer_cells_stratified.rds')
cancer_samples <- readRDS('outputs/NSCLC_cancer_samples_stratified.rds')
cancer_cells <- readRDS('outputs/NSCLC_cancer_cells_stratified.rds')

# rerun umap
cancer_samples <- RunUMAP(
  cancer_samples,
  reduction = 'scANVI',
  dims = 1:10,
  reduction.name = 'cancer_sample_umap'
)
cancer_cells <- RunUMAP(
  cancer_cells,
  reduction = 'scANVI',
  dims = 1:10,
  reduction.name = 'cancer_cell_umap'
)

# saveRDS(cancer_samples, 'outputs/NSCLC_cancer_samples_stratified_umap.rds')
# saveRDS(cancer_cells, 'outputs/NSCLC_cancer_cells_stratified_umap.rds')
cancer_samples <-
  readRDS('outputs/NSCLC_cancer_samples_stratified_umap.rds')
cancer_cells <- readRDS('outputs/NSCLC_cancer_cells_stratified_umap.rds')


# visualize ####
for (i in m6a) {
  p <- DimPlot(
    cancer_cells,
    group.by = paste0(i, '_level'),
    reduction = 'cancer_cell_umap',
    cols = c('#F21A00', '#3B9AB2', '#EBCC2A', 'grey'),
    raster = FALSE,
    shuffle = TRUE
  ) +
    labs(title = i, x = 'UMAP_1', y = 'UMAP_2')

  ggsave(
    paste0('cancer_cells_stratified_umap_', i, '.pdf'),
    p,
    path = 'outputs/figures/fig1/cancer_cells_stratified_umap/',
    width = 5,
    height = 5
  )
}

# barplot
for (i in m6a) {
  data_cancer_cells[[paste0(i, '_level')]] <- factor(
    data_cancer_cells[[paste0(i, '_level')]],
    levels = c('High', 'Medium', 'Low', 'Not detected')
  )

  p <- ggplot(
    data_cancer_cells,
    aes(x = platform, fill = .data[[paste0(i, '_level')]])
  ) +
    geom_bar(position = 'fill', alpha = 0.5) +
    labs(
      title = paste0(i, ' expression in different platforms'),
      x = NULL,
      y = 'Percentage',
      fill = i
    ) +
    scale_fill_manual(values = c('#F21A00', '#EBCC2A', '#3B9AB2', 'grey')) +
    theme_classic() +
    scale_y_continuous(labels = scales::percent)

  ggsave(
    paste0('m6a_expression_in_platforms_', i, '.pdf'),
    p,
    path = 'outputs/figures/fig1/cancer_cells_m6a_barplot_platform',
    width = 5,
    height = 5
  )

  p <- ggplot(
    data_cancer_cells,
    aes(x = dataset, fill = .data[[paste0(i, '_level')]])
  ) +
    geom_bar(position = 'fill', alpha = 0.5) +
    labs(
      title = paste0(i, ' expression in different datasets'),
      x = NULL,
      y = 'Percentage',
      fill = i
    ) +
    scale_fill_manual(values = c('#F21A00', '#EBCC2A', '#3B9AB2', 'grey')) +
    theme_classic() +
    scale_y_continuous(labels = scales::percent) +
  coord_flip()

  ggsave(
    paste0('m6a_expression_in_datasets_', i, '.pdf'),
    p,
    path = 'outputs/figures/fig1/cancer_cells_m6a_barplot_dataset',
    width = 8,
    height = 5
  )
}

# vlnplot
for (i in m6a) {
  p <-ggplot(
    data_cancer_cells[data_cancer_cells[[i]] != 0, ],
    aes(x = platform, y = .data[[i]])
  ) +
    geom_violin(
      fill = '#F21A00',
      width = 0.6,
      alpha = 0.2,
      aes(linetype = NA)
    ) +
    geom_boxplot(width = 0.1, fill = '#F21A00', alpha = 0.5) +
    labs(
      title = paste0(i, ' expression in different platforms'),
      x = NULL,
      y = 'Log normalized data'
    ) +
    theme_classic() +
    theme(legend.position = 'none')

  ggsave(
    paste0('m6a_expression_in_platforms_', i, '.pdf'),
    p,
    path = 'outputs/figures/fig1/cancer_cells_m6a_vlnplot_platform',
    width = 5,
    height = 5
  )

  p <-ggplot(
    data_cancer_cells[data_cancer_cells[[i]] != 0, ],
    aes(x = dataset, y = .data[[i]])
  ) +
    geom_violin(
      fill = '#F21A00',
      width = 1,
      alpha = 0.2,
      aes(linetype = NA)
    ) +
    geom_boxplot(width = 0.1, fill = '#F21A00', alpha = 0.5) +
    labs(
      title = paste0(i, ' expression in different datasets'),
      x = NULL,
      y = 'Log normalized data'
    ) +
    theme_classic() +
    theme(legend.position = 'none') +
    coord_flip()

  ggsave(
    paste0('m6a_expression_in_datasets_', i, '.pdf'),
    p,
    path = 'outputs/figures/fig1/cancer_cells_m6a_vlnplot_dataset',
    width = 8,
    height = 5
  )
}
