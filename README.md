# NSCLC-m6a-project

code for NSCLC-m6a-project

## preprocess.R

1. rename feature names
2. extract meta data
3. rerun umap

## m6a_stratify.R

1. subset cancer samples and cancer cells in 10x and Singleron platforms, and rerun umap
1. extract m6a-gene expression data and meta data
1. stratify m6a levels in each platform into 'Not detected', 'Low', 'Medium' and 'High'
1. visualize stratification (dimplot, barplot and vlnplot)

## cell_chat.R

1. 