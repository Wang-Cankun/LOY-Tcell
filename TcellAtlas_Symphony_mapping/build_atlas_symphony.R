
# Convert Seurat object to h5ad
library(here)
library(qs)
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(symphony)
library(singlecellmethods)

# remotes::install_github("immunogenomics/singlecellmethods")

setwd('/Users/wang.13246/Documents/Project/loy/model_train')
set.seed(42)


# Atlas
sample_name <- "GC_sample.rds"
combined <- read_rds(sample_name)
DefaultAssay(combined) <- "RNA"
combined <- UpdateSeuratObject(combined) # GC_sample.rds is an old version Seurat object, need to run UpdateSeuratObject


# Run symphony
ref_exp_full <- GetAssayData(combined, slot = "data")
ref_metadata = combined@meta.data

reference = symphony::buildReference(
  ref_exp_full,                   # reference expression (genes by cells)
  ref_metadata,              # reference metadata (cells x attributes)
  vars = c('sample'),         # variable(s) to integrate over
  K = 100,                   # number of Harmony soft clusters
  verbose = TRUE,            # display verbose output
  do_umap = TRUE,            # run UMAP and save UMAP model to file
  do_normalize = FALSE,      # perform log(CP10k) normalization on reference expression
  vargenes_method = 'vst',   # variable gene selection method: 'vst' or 'mvp'
  vargenes_groups = 'sample', # metadata column specifying groups for variable gene selection within each group
  topn = 2000,               # number of variable genes (per group)
  theta = 2,                 # Harmony parameter(s) for diversity term
  d = 20,                    # number of dimensions for PCA
  save_uwot_path = './reference_Tcell_atlas', # file path to save uwot UMAP model
  additional_genes = NULL    # vector of any additional genes to force include
)

