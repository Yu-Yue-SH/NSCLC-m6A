rm(list = ls())


# load packages
library(Seurat)


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
  thresholds = c(0, 0.4, 0.6, 1),
  m6a
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


# subset smartseq2 samples ####
smart_seq2_samples <- subset(NSCLC, platform == 'Smart-seq2')
# smartseq2 has 22974 LUAD cells, 3070 LUSC cells and 8272 normal cells;
# 3636 cancer cells;
# 10051 normal adjacent cells, 10629 tumor metastasis cells and 12648
# tumor primary cells


# run stratify_m6a_levels_in_platforms ####
smart_seq2_samples <- stratify_m6a_levels_in_platforms(
  smart_seq2_samples,
  platforms = 'Smart-seq2',
  thresholds = c(0, 0.4, 0.6, 1),
  m6a = m6a
)

saveRDS(
  smart_seq2_samples,
  'outputs/NSCLC_smart_seq2_samples_stratified.rds'
)
smart_seq2_samples <-
  readRDS('outputs/NSCLC_smart_seq2_samples_stratified.rds')
