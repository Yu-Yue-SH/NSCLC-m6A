# NSCLC-m6a-project

code for NSCLC-m6a-project

## preprocess.R

1. rename feature names
2. extract meta data
3. rerun umap

## m6a_stratify.R

1. create function: function: stratify_m6a_levels_in_platforms()
1. run stratify_m6a_levels_in_platforms() for cancer samples
1. run stratify_m6a_levels_in_platforms()  for LUAD and LUSC samples

## m6a_stratify_smartseq2.R

1. run stratify_m6a_levels_in_platforms() for Smart-seq2 samples

## cell_chat.R

1. create function: run_cell_chat()
1. run run_cell_chat() for cancer samples and all m6a levels
1. create function: cellchat_comparison()
1. run cellchat_comparison() for all m6a levels
1. subset 10x and rerun all functions
1. run cellchat for LUAD and LUSC

## cell_chat_smartseq2.R

