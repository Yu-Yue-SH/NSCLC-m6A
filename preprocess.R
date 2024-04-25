rm(list = ls())


# load packages
library(Seurat)


# load data
NSCLC <- readRDS('data/NSCLC_extended.rds')
NSCLC <- UpdateSeuratObject(NSCLC)


# rename features ####
counts <- GetAssayData(NSCLC, layer = 'counts')
data <- GetAssayData(NSCLC, layer = 'data')
meta_features <- NSCLC[['RNA']]@meta.features
meta_features[['ENSG_id']] <- rownames(meta_features)
rownames(meta_features) <- meta_features[['feature_name']]
rownames(counts) <- meta_features[['feature_name']]
rownames(data) <- meta_features[['feature_name']]
meta_data <- NSCLC@meta.data
re_scANVI <- NSCLC@reductions$scANVI
re_scVI <- NSCLC@reductions$scVI
re_umap <- NSCLC@reductions$umap

NSCLC <- CreateSeuratObject(
  counts = counts,
  meta.data = meta_data,
  project = 'NSCLC1m'
)

NSCLC <- SetAssayData(NSCLC, layer = 'data', data)
NSCLC[['RNA']]@meta.data <- meta_features
NSCLC@reductions$scANVI <- re_scANVI
NSCLC@reductions$scVI <- re_scVI
NSCLC@reductions$umap <- re_umap

saveRDS(NSCLC, 'outputs/NSCLC_renamed.rds')
NSCLC <- readRDS('outputs/NSCLC_renamed.rds')


# extract metadata ####
meta_data <- NSCLC@meta.data

saveRDS(meta_data, 'outputs/meta_data.rds')
meta_data <- readRDS('outputs/meta_data.rds')


# rerun umap ####
NSCLC <- RunUMAP(
  NSCLC,
  reduction = 'scANVI',
  dims = 1:10,
  reduction.name = 'new_umap'
)

saveRDS(NSCLC, 'outputs/NSCLC_new_umap.rds')
NSCLC <- readRDS('outputs/NSCLC_new_umap.rds')


# rename features in core data ####
# load data
NSCLC <- readRDS('data/NSCLC_core.rds')
NSCLC <- UpdateSeuratObject(NSCLC)


# rename features
counts <- GetAssayData(NSCLC, layer = 'counts')
data <- GetAssayData(NSCLC, layer = 'data')
meta_features <- NSCLC[['RNA']]@meta.features
meta_features[['ENSG_id']] <- rownames(meta_features)
rownames(meta_features) <- meta_features[['feature_name']]
rownames(counts) <- meta_features[['feature_name']]
rownames(data) <- meta_features[['feature_name']]
meta_data <- NSCLC@meta.data
re_scANVI <- NSCLC@reductions$scANVI
re_scVI <- NSCLC@reductions$scVI
re_umap <- NSCLC@reductions$umap

NSCLC <- CreateSeuratObject(
  counts = counts,
  meta.data = meta_data,
  project = 'NSCLC1m'
)

NSCLC <- SetAssayData(NSCLC, layer = 'data', data)
NSCLC[['RNA']]@meta.data <- meta_features
NSCLC@reductions$scANVI <- re_scANVI
NSCLC@reductions$scVI <- re_scVI
NSCLC@reductions$umap <- re_umap

saveRDS(NSCLC, 'outputs/NSCLC_core_renamed.rds')
NSCLC <- readRDS('outputs/NSCLC_core_renamed.rds')
