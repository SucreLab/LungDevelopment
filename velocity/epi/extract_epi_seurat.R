library(Seurat)

wd <- "/Users/negretn/postdoc/code/devo_scseq/"

setwd(wd)
infile <- paste0(wd, "data/epi_full_celltype_labels.rds")

data <- readRDS(infile)


scaled_data <- GetAssayData(data, slot = "data", assay = "SCT")

Matrix::writeMM(scaled_data, paste0(wd, "velocity/epi/data/exp_data.mtx"))

write.csv(Idents(data), paste0(wd, "velocity/epi/data/idents.csv"))
write.csv(Embeddings(data, reduction = "umap"), paste0(wd, "velocity/epi/data/umap.csv"))
write.csv(Embeddings(data, reduction = "pca"), paste0(wd, "velocity/epi/data/pca.csv"))
write.csv(data$timepoint, paste0(wd, "velocity/epi/data/timepoint.csv"))
write.csv(rownames(data), paste0(wd, "velocity/epi/data/genes.csv"))
write.csv(colnames(data), paste0(wd, "velocity/epi/data/cells.csv"))

write.csv(data$bulk_celltype, paste0(wd, "velocity/epi/data/bulk_celltype.csv"))
write.csv(Idents(data), paste0(wd, "velocity/epi/data/cell_subtype.csv"))



system(paste0("gzip ", wd, "velocity/epi/data/exp_data.mtx"))