rm(list = ls())


# load packages
library(Seurat)
library(ggplot2)


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


# function: stratify_m6a_levels_in_platforms ####
# output: a seurat object of the samples with m6a levels stratified, saved
stratify_m6a_levels_in_platforms <- function (
  cancer_samples,
  platforms = c('10x', 'Singleron'),
  thresholds = c(0, 0.4, 0.6, 1)
) {
  # extract data
  cancer_cells <- subset(cancer_samples, ann_fine == 'Cancer cells')

  data_cancer_cells <- as.data.frame(GetAssayData(cancer_cells, layer = 'data'))
  data_cancer_cells <- data_cancer_cells[m6a, ]
  data_cancer_cells <- t(data_cancer_cells)
  data_cancer_cells <- cbind(data_cancer_cells, cancer_cells@meta.data)

  # stratify m6a levels
  for (i in m6a) {
    data_cancer_cells[paste0(i, '_level')] <- NA
    data_cancer_cells[paste0(i, '_level')][data_cancer_cells[[i]] == 0, ] <-
      'Not detected'
    for (j in platforms) {
      cells <- data_cancer_cells[[i]] != 0 &
        data_cancer_cells[['platform']] == j
      data_cancer_cells[cells, paste0(i, '_level')] <- as.character(
        cut(
          data_cancer_cells[cells, i],
          breaks = quantile(data_cancer_cells[cells, i], thresholds),
          labels = c('Low', 'Medium', 'High'),
          include.lowest = TRUE
        )
      )
    }
  }

  # add to seurat object
  cancer_samples[['new_broad_type']] <- cancer_samples[['ann_coarse']]
  cancer_samples@meta.data[['new_broad_type']] <-
    as.character(cancer_samples@meta.data[['new_broad_type']])
  cancer_samples[['new_broad_type']][cancer_samples[['ann_fine']] ==
                                       'Cancer cells'] <- 'Cancer cells'

  for (i in m6a) {
    cancer_cells[[paste0(i, '_level')]] <- NA
    cancer_cells[[paste0(i, '_level')]]<-
      data_cancer_cells[[paste0(i, '_level')]]
    cancer_samples[[paste0(i, '_level')]] <- cancer_samples[['new_broad_type']]
    cancer_samples[[paste0(i, '_level')]] <-
      cancer_cells[[paste0(i, '_level')]]
  }

  return(cancer_samples)
}


# subset cancer samples ####
cancer_samples <- subset(
  NSCLC,
  platform %in% c('10x', 'Singleron') &
    origin %in% c('tumor_matastasis', 'tumor_primary')
)

# saveRDS(cancer_samples, 'outputs/NSCLC_cancer_samples.rds')
cancer_samples <- readRDS('outputs/NSCLC_cancer_samples.rds')


# run stratify_m6a_levels_in_platforms ####
cancer_samples <- stratify_m6a_levels_in_platforms(
  cancer_samples,
  platforms = c('10x', 'Singleron'),
  thresholds = c(0, 0.4, 0.6, 1)
)

saveRDS(
  cancer_samples,
  'outputs/NSCLC_cancer_samples_stratified.rds'
)
cancer_samples <-
  readRDS('outputs/NSCLC_cancer_samples_stratified.rds')


# subset primary cancer samples ####
cancer_samples_primary <- subset(cancer_samples, origin == 'tumor_primary')

saveRDS(cancer_samples_primary, 'outputs/NSCLC_cancer_primary_samples.rds')
cancer_samples_primary <-
  readRDS('outputs/NSCLC_cancer_primary_samples.rds')


# subset LUAD and LUSC samples ####
luad_samples <- subset(cancer_samples_primary, disease == 'lung adenocarcinoma')
lusc_samples <- subset(
  cancer_samples_primary,
  disease == 'squamous cell lung carcinoma'
)

saveRDS(luad_samples, 'outputs/NSCLC_luad_primary_samples.rds')
saveRDS(lusc_samples, 'outputs/NSCLC_lusc_primary_samples.rds')
luad_samples <- readRDS('outputs/NSCLC_luad_primary_samples.rds')
lusc_samples <- readRDS('outputs/NSCLC_lusc_primary_samples.rds')


# run stratify_m6a_levels_in_platforms for LUAD and LUSC ####
luad_samples <- stratify_m6a_levels_in_platforms(
  luad_samples,
  platforms = c('10x', 'Singleron'),
  thresholds = c(0, 0.4, 0.6, 1)
)

saveRDS(luad_samples, 'outputs/NSCLC_luad_primary_samples_stratified.rds')
luad_samples <- readRDS('outputs/NSCLC_luad_primary_samples_stratified.rds')

lusc_samples <- stratify_m6a_levels_in_platforms(
  lusc_samples,
  platforms = c('10x', 'Singleron'),
  thresholds = c(0, 0.4, 0.6, 1)
)

saveRDS(lusc_samples, 'outputs/NSCLC_lusc_primary_samples_stratified.rds')
lusc_samples <- readRDS('outputs/NSCLC_lusc_primary_samples_stratified.rds')
