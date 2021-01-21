library(Seurat)

wd <- "/Users/negretn/postdoc/code/devo_scseq/"

setwd(wd)
infile <- paste0(wd, "data/meso_full_celltype_labels.rds")

data <- readRDS(infile)


scaled_data <- GetAssayData(data, slot = "data", assay = "SCT")

Matrix::writeMM(scaled_data, paste0(wd, "velocity/meso/data/exp_data.mtx"))

write.csv(Idents(data), paste0(wd, "velocity/meso/data/idents.csv"))
write.csv(Embeddings(data, reduction = "umap"), paste0(wd, "velocity/meso/data/umap.csv"))
write.csv(Embeddings(data, reduction = "pca"), paste0(wd, "velocity/meso/data/pca.csv"))
write.csv(data$timepoint, paste0(wd, "velocity/meso/data/timepoint.csv"))
write.csv(rownames(data), paste0(wd, "velocity/meso/data/genes.csv"))
write.csv(colnames(data), paste0(wd, "velocity/meso/data/cells.csv"))

write.csv(data$bulk_celltype, paste0(wd, "velocity/meso/data/bulk_celltype.csv"))
write.csv(Idents(data), paste0(wd, "velocity/meso/data/cell_subtype.csv"))



system(paste0("gzip ", wd, "velocity/meso/data/exp_data.mtx"))