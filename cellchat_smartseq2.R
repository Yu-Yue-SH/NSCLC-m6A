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
run_cellchat <- function(data_input, meta, group = 'new_broad_type', w = 8) {
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
cellchat_comparison <- function (cellchat_hi, cellchat_lo) {
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

cellchat <- run_cellchat(data_input, meta, 'new_broad_type', 8)

saveRDS(cellchat, 'outputs/cellchat/smartseq2/cellchat_new_broad_type.rds')
cellchat <-
  readRDS('outputs/cellchat/smartseq2/cellchat_new_broad_type.rds')

# run for all m6a levels
for (i in m6a) {
  # first high then low
  for (j in c('High', 'Low')) {
    seurat_obj <- subset(
      smart_seq2_samples,
      cells = which(
        smart_seq2_samples$new_broad_type != 'Cancer cells' |
          smart_seq2_samples[[paste0(i, '_level')]] == j
      )
    )
    data_input <- GetAssayData(seurat_obj, layer = 'data')
    meta <- seurat_obj@meta.data
    cellchat <- run_cellchat(data_input, meta, 'new_broad_type', 8)
    saveRDS(
      cellchat,
      paste0('outputs/cellchat/smartseq2/04/raw/cellchat_', i, '_', j, '.rds')
    )
  }
}


# run cellchat_comparison ####
for (i in m6a) {
  # load data
  cellchat_hi <- readRDS(
    paste0(
      'outputs/cellchat/smartseq2/04/raw/cellchat_',
      i,
      '_',
      'High',
      '.rds'
    )
  )
  cellchat_lo <- readRDS(
    paste0(
      'outputs/cellchat/smartseq2/04/raw/cellchat_',
      i,
      '_',
      'Low',
      '.rds'
    )
  )

  # run cellchat_comparison
  comparison_list <- cellchat_comparison(cellchat_hi, cellchat_lo)

  saveRDS(
    comparison_list[[1]],
    paste0(
      'outputs/cellchat/smartseq2/04/list/cellchat_',
      i,
      '_object_list.rds'
    )
  )
  saveRDS(
    comparison_list[[2]],
    paste0('outputs/cellchat/smartseq2/04/diff/cellchat_', i, '_diff.rds')
  )
}
