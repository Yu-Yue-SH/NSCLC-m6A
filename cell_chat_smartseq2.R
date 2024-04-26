rm(list = ls())


# load packages
library(Seurat)
library(CellChat)
library(patchwork)


# set options
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 50 * 1024^3)


# load data
smart_seq2_samples <-
  readRDS('outputs/NSCLC_smart_seq2_samples_stratified_04.rds')


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


# function: cellchat_comparison ####
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

cellchat_comparison <- function (
  cellchat,
  hi = hi,
  lo = lo,
  levels_hi = levels_hi,
  levels_lo = levels_hi
) {
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


# run run_cell_chat ####
data_input <- GetAssayData(smart_seq2_samples, layer = 'data')
meta <- smart_seq2_samples@meta.data

cellchat <- run_cell_chat(data_input, meta, 'new_broad_type', 8)

saveRDS(cellchat, 'outputs/cellchat/smartseq2/cellchat_new_broad_type.rds')
cellchat <-
  readRDS('outputs/cellchat/smartseq2/cellchat_new_broad_type.rds')

# run for all m6a levels
for (i in m6a) {
  cellchat <- run_cell_chat(data_input, meta, paste0(i, '_level'), 8)
  saveRDS(
    cellchat,
    paste0(
      'outputs/cellchat/smartseq2/04/origin/',
      'cellchat_',
      i,
      '.rds'
    )
  )
}


# run cellchat_comparison ####
for (i in m6a) {
  # load data
  cellchat <-
    readRDS(paste0('outputs/cellchat/smartseq2/04/origin/cellchat_', i, '.rds'))

  # run cellchat_comparison
  comparison_list <- cellchat_comparison(cellchat)

  saveRDS(
    comparison_list[[1]],
    paste0(
      'outputs/cellchat/smartseq2/04/object_list/cellchat_',
      i,
      '_object_list.rds'
    )
  )
  saveRDS(
    comparison_list[[2]],
    paste0('outputs/cellchat/smartseq2/04/diff/cellchat_', i, '_diff.rds')
  )
}


# change m6a thresholds ####
smart_seq2_samples <-
  readRDS('outputs/NSCLC_smart_seq2_samples_stratified_03.rds')

data_input <- GetAssayData(smart_seq2_samples, layer = 'data')
meta <- smart_seq2_samples@meta.data

for (i in m6a) {
  cellchat <- run_cell_chat(data_input, meta, paste0(i, '_level'), 8)
  saveRDS(
    cellchat,
    paste0(
      'outputs/cellchat/smartseq2/03/origin/',
      'cellchat_',
      i,
      '.rds'
    )
  )
}

for (i in m6a) {
  # load data
  cellchat <-
    readRDS(paste0('outputs/cellchat/smartseq2/03/origin/cellchat_', i, '.rds'))

  # run cellchat_comparison
  comparison_list <- cellchat_comparison(cellchat)

  saveRDS(
    comparison_list[[1]],
    paste0(
      'outputs/cellchat/smartseq2/03/object_list/cellchat_',
      i,
      '_object_list.rds'
    )
  )
  saveRDS(
    comparison_list[[2]],
    paste0('outputs/cellchat/smartseq2/03/diff/cellchat_', i, '_diff.rds')
  )
}
