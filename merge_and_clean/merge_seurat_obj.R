library(Seurat)


merged_seurat <- Reduce(function(x,y){merge(x,y)}, list(readRDS("/scratch/negretn/epi_full_celltype_labels.rds"),
                                       readRDS("/scratch/negretn/endo_full_celltype_labels.rds"),
                                       readRDS("/scratch/negretn/meso_full_celltype_labels.rds")))


merged_seurat$cell_subtype <- Idents(merged_seurat)
merged_seurat$bulk_celltype <- as.factor(merged_seurat$bulk_celltype)
saveRDS(merged_seurat, "/scratch/negretn/merged_alldata.rds")