rm(list = ls())


# load packages
library(Seurat)
library(CellChat)
library(patchwork)


# set options
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 50 * 1024^3)


# load data
cancer_samples <-
  readRDS('outputs/NSCLC_cancer_samples_stratified_umap.rds')


# load m6a genes
writers <- c('METTL3', 'METTL5', 'METTL14', 'METTL16', 'WTAP', 'RBM15',
             'RBM15B', 'ZC3H13', 'VIRMA', 'CBLL1')
erasers <- c('FTO', 'ALKBH5')
readers <- c('YTHDF1', 'YTHDF2', 'YTHDF3', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3',
             'YTHDC1', 'YTHDC2', 'ELAVL1', 'HNRNPC', 'RBMX', 'HNRNPA2B1')
m6a <- c(writers, erasers, readers)


# function: run_cell_chat ####
run_cell_chat <- function(data_input, meta, group = 'new_broad_type', w = 8) {
  # create object
  cellchat <- createCellChat(data_input, meta, group.by = group)

  # set database
  CellChatDB <- CellChatDB.human
  CellChatDB_use <- subsetDB(CellChatDB)
  cellchat@DB <- CellChatDB_use

  # preprocess
  cellchat <- subsetData(cellchat)
  future::plan("multisession", workers = w)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # run
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  return(cellchat)
}

# run run_cell_chat ####
data_input <- GetAssayData(cancer_samples, layer = 'data')
meta <- cancer_samples@meta.data

cellchat <- run_cell_chat(data_input, meta, 'new_broad_type', 8)

saveRDS(cellchat, 'outputs/cellchat/cellchat_new_broad_type.rds')
cellchat <- readRDS('outputs/cellchat/cellchat_new_broad_type.rds')


# run for all m6a levels ####
for (i in m6a) {
  cellchat <- run_cell_chat(data_input, meta, paste0(i, '_level'), 8)
  saveRDS(cellchat, paste0('outputs/cellchat/origin/', 'cellchat_', i, '.rds'))
}


# function: cellchat_comparison ####
cellchat_comparison <- function (cellchat, hi, lo, levels_hi, levels_lo) {
  cellchat_hi <- subsetCellChat(cellchat, idents.use = hi)
  cellchat_lo <- subsetCellChat(cellchat, idents.use = lo)

  # rename idents
  cellchat_hi <- setIdent(
    cellchat_hi,
    ident.use = 'new_broad_type',
    levels = levels_hi
  )
  cellchat_lo <- setIdent(
    cellchat_lo,
    ident.use = 'new_broad_type',
    levels = levels_lo
  )

  cellchat_hi <- filterCommunication(cellchat_hi, min.cells = 10)
  cellchat_hi <- computeCommunProbPathway(cellchat_hi)
  cellchat_hi <- aggregateNet(cellchat_hi)
  cellchat_lo <- filterCommunication(cellchat_lo, min.cells = 10)
  cellchat_lo <- computeCommunProbPathway(cellchat_lo)
  cellchat_lo <- aggregateNet(cellchat_lo)

  # merge object
  object_list <- list(lo = cellchat_lo, hi = cellchat_hi)
  cellchat_diff <- mergeCellChat(
    object_list,
    add.names = names(object_list),
    cell.prefix = TRUE
  )

  return(list(object_list, cellchat_diff))
}


# run cellchat_comparison ####
hi <- c('B cell', 'High', 'cDC', 'Endothelial cell', 'Epithelial cell',
        'Macrophage/Monocyte', 'Mast cell', 'Neutrophils', 'NK cell', 'pDC',
        'Plasma cell', 'Stromal', 'T cell')
lo <- c('B cell', 'Low', 'cDC', 'Endothelial cell', 'Epithelial cell',
        'Macrophage/Monocyte', 'Mast cell', 'Neutrophils', 'NK cell', 'pDC',
        'Plasma cell', 'Stromal', 'T cell')
levels_hi <- c('B cell', 'cDC', 'Endothelial cell', 'Epithelial cell',
               'Cancer cells', 'Macrophage/Monocyte', 'Mast cell',
               'Neutrophils', 'NK cell', 'pDC', 'Plasma cell', 'Stromal',
               'T cell')
levels_lo <- c('B cell', 'cDC', 'Endothelial cell', 'Epithelial cell',
               'Cancer cells', 'Macrophage/Monocyte', 'Mast cell',
               'Neutrophils', 'NK cell', 'pDC', 'Plasma cell', 'Stromal',
               'T cell')

for (i in m6a) {
  # load data
  cellchat <- readRDS(paste0('outputs/cellchat/origin/cellchat_', i, '.rds'))

  # run cellchat_comparison
  comparison_list <- cellchat_comparison(cellchat, hi, lo, levels_hi, levels_lo)

  saveRDS(
    comparison_list[1],
    paste0('outputs/cellchat/object_list/cellchat_', i, '_object_list.rds')
  )
  saveRDS(
    comparison_list[2],
    paste0('outputs/cellchat/diff/cellchat_', i, '_diff.rds')
  )
}


# subset 10x platform only ####
cancer_samples_10x <- subset(cancer_samples, platform == '10x')

# run run_cell_chat
data_input <- GetAssayData(cancer_samples_10x, layer = 'data')
meta <- cancer_samples_10x@meta.data

cellchat <- run_cell_chat(data_input, meta, 'new_broad_type', 8)

saveRDS(cellchat, 'outputs/cellchat/10x/cellchat_new_broad_type.rds')
cellchat <- readRDS('outputs/cellchat/10x/cellchat_new_broad_type.rds')

# run for all m6a levels
for (i in m6a) {
  cellchat <- run_cell_chat(data_input, meta, paste0(i, '_level'), 8)
  saveRDS(
    cellchat,
    paste0('outputs/cellchat/10x/origin/', 'cellchat_', i, '.rds')
  )
}

# run cellchat_comparison
for (i in m6a) {
  # load data
  cellchat <-
    readRDS(paste0('outputs/cellchat/10x/origin/cellchat_', i, '.rds'))

  # run cellchat_comparison
  comparison_list <- cellchat_comparison(cellchat, hi, lo, levels_hi, levels_lo)

  saveRDS(
    comparison_list[1],
    paste0('outputs/cellchat/10x/object_list/cellchat_', i, '_object_list.rds')
  )
  saveRDS(
    comparison_list[2],
    paste0('outputs/cellchat/10x/diff/cellchat_', i, '_diff.rds')
  )
}


# run cellchat for LUAD and LUSC ####
# LUAD
luad_samples <- readRDS('outputs/NSCLC_luad_primary_samples_stratified.rds')

# run run_cell_chat
data_input <- GetAssayData(luad_samples, layer = 'data')
meta <- luad_samples@meta.data

cellchat <- run_cell_chat(data_input, meta, 'new_broad_type', 8)

saveRDS(cellchat, 'outputs/cellchat/luad/cellchat_new_broad_type.rds')
cellchat <- readRDS('outputs/cellchat/luad/cellchat_new_broad_type.rds')

# run for all m6a levels
for (i in m6a) {
  cellchat <- run_cell_chat(data_input, meta, paste0(i, '_level'), 8)

  saveRDS(
    cellchat,
    paste0('outputs/cellchat/luad/origin/', 'cellchat_', i, '.rds')
  )
}

# run cellchat_comparison
for (i in m6a) {
  # load data
  cellchat <-
    readRDS(paste0('outputs/cellchat/luad/origin/cellchat_', i, '.rds'))

  # run cellchat_comparison
  comparison_list <- cellchat_comparison(cellchat, hi, lo, levels_hi, levels_lo)

  saveRDS(
    comparison_list[[1]],
    paste0('outputs/cellchat/luad/object_list/cellchat_', i, '_object_list.rds')
  )
  saveRDS(
    comparison_list[[2]],
    paste0('outputs/cellchat/luad/diff/cellchat_', i, '_diff.rds')
  )
}

# LUSC
lusc_samples <- readRDS('outputs/NSCLC_lusc_primary_samples_stratified.rds')

# run run_cell_chat
data_input <- GetAssayData(lusc_samples, layer = 'data')
meta <- lusc_samples@meta.data

cellchat <- run_cell_chat(data_input, meta, 'new_broad_type', 8)

saveRDS(cellchat, 'outputs/cellchat/lusc/cellchat_new_broad_type.rds')
cellchat <- readRDS('outputs/cellchat/lusc/cellchat_new_broad_type.rds')

# run for all m6a levels
for (i in m6a) {
  cellchat <- run_cell_chat(data_input, meta, paste0(i, '_level'), 8)

  saveRDS(
    cellchat,
    paste0('outputs/cellchat/lusc/origin/', 'cellchat_', i, '.rds')
  )
}

# run cellchat_comparison
for (i in m6a) {
  # load data
  cellchat <-
    readRDS(paste0('outputs/cellchat/lusc/origin/cellchat_', i, '.rds'))

  # run cellchat_comparison
  comparison_list <- cellchat_comparison(cellchat, hi, lo, levels_hi, levels_lo)

  saveRDS(
    comparison_list[[1]],
    paste0('outputs/cellchat/lusc/object_list/cellchat_', i, '_object_list.rds')
  )
  saveRDS(
    comparison_list[[2]],
    paste0('outputs/cellchat/lusc/diff/cellchat_', i, '_diff.rds')
  )
}
