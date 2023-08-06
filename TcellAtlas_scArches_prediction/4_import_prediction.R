
# Convert Seurat object to h5ad
library(here)
library(qs)
library(tidyverse)
library(Seurat)
library(SeuratDisk)

setwd('/bmbl_data/cankun_notebook/loss_y/')
set.seed(42)


# Query data
sample_name <- "Seurat.merge.T.excluded.qs"
combined <- qs::qread(sample_name)
DefaultAssay(combined) <- "RNA"

pred_embed <- read.csv(paste0("TcellAtlas_query_emb_full.csv"), row.names = 1)

# Check if the cell names are matched
colnames(combined)[1:50]
rownames(pred_embed)[1:50]


#rownames(pred_embed) <- colnames(combined)
combined <- AddMetaData(combined, pred_embed)

