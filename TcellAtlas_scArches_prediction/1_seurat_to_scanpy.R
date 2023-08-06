
# Convert Seurat object to h5ad
library(here)
library(qs)
library(tidyverse)
library(Seurat)
library(SeuratDisk)

setwd('/bmbl_data/cankun_notebook/loss_y/')
set.seed(42)




# Atlas
sample_name <- "GC_sample.rds"
combined <- read_rds(sample_name)
DefaultAssay(combined) <- "RNA"

combined2 <- UpdateSeuratObject(combined) # GC_sample.rds is an old version Seurat object, need to run UpdateSeuratObject
combined2@assays$RNA@scale.data <- matrix() # scaled data is not needed
#combined2@assays$RNA@data <- combined2@assays$RNA@counts
SaveH5Seurat(combined2,
             filename = paste0(sample_name, ".h5Seurat"),
             overwrite = T)
Convert(paste0(sample_name, ".h5Seurat"),
        dest = "h5ad",
        overwrite = T)





# Query data
sample_name <- "Seurat.merge.T.excluded.qs"
combined <- qs::qread(sample_name)
DefaultAssay(combined) <- "RNA"

combined2 <- UpdateSeuratObject(combined)
combined2@assays$RNA@scale.data <- matrix() # scaled data is not needed
combined2@assays$RNA@data <- combined2@assays$RNA@counts
SaveH5Seurat(combined2,
             filename = paste0(sample_name, ".h5Seurat"),
             overwrite = T)
Convert(paste0(sample_name, ".h5Seurat"),
        dest = "h5ad",
        overwrite = T)


expr <- GetAssayData(combined, slot = "counts")
write.table(expr, "expr.txt", quote = F, row.names = T, col.names = T)

