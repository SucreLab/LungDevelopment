library(Seurat)

wd <- "/scratch/neretn"

setwd(wd)
infile <- paste0(wd, "merged_alldata_noGm42418_sct_p7b_integrated_retransform.rds")

data <- readRDS(infile)
levels(Idents(data))[levels(Idents(data))=="Matrix fibroblast"] <- "Adventitial fibroblast"

scaled_data <- GetAssayData(data, slot = "data", assay = "SCT")

Matrix::writeMM(scaled_data, paste0(wd, "exp_data.mtx"))
write.csv(Idents(data), paste0(wd, "idents.csv"))
write.csv(Embeddings(data, reduction = "umap"), paste0(wd, "umap.csv"))
write.csv(data$timepoint, paste0(wd, "timepoint.csv"))
write.csv(rownames(data), paste0(wd, "genes.csv"))

write.csv(data$bulk_celltype, paste0(wd, "bulk_celltype.csv"))
write.csv(data$cell_subtype, paste0(wd, "cell_subtype.csv"))
