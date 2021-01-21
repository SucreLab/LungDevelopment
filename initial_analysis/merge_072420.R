library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(SeuratDisk)

#read in post-soupX counts
# These are all available in GSE and GSE.
e12.data <- Read10X(data.dir = "~/Dropbox/Development Post SoupX/E12/filtered_feature_bc_matrix/")
e15.data <- Read10X(data.dir = "~/Dropbox/Development Post SoupX/E15_2/filtered_feature_bc_matrix/")
e18.data <- Read10X(data.dir = "~/Dropbox/Development Post SoupX/E18/filtered_feature_bc_matrix/")
p0a.data <- Read10X(data.dir = "~/Dropbox/Development Post SoupX/P0/filtered_feature_bc_matrix/")
p0b.data <- Read10X(data.dir = "~/Dropbox/Development Post SoupX/P0_2/filtered_feature_bc_matrix/")
p3.data <- Read10X(data.dir = "~/Dropbox/Development Post SoupX/P3/filtered_feature_bc_matrix/")
p5a.data <- Read10X(data.dir = "~/Dropbox/Development Post SoupX/P5_1/filtered_feature_bc_matrix/")
p5b.data <- Read10X(data.dir = "~/Dropbox/Development Post SoupX/P5_3/filtered_feature_bc_matrix/")
p7a.data <- Read10X(data.dir = "~/Dropbox/Development Post SoupX/P7/filtered_feature_bc_matrix/")
p7b.data <- Read10X(data.dir = "~/Dropbox/Development Post SoupX/P7_3/filtered_feature_bc_matrix/")
p7c.data <- Read10X(data.dir = "~/Dropbox/Development Post SoupX/P7_enriched/filtered_feature_bc_matrix/")
p14a.data <- Read10X(data.dir = "~/Dropbox/Development Post SoupX/P14/filtered_feature_bc_matrix/")
p14b.data <- Read10X(data.dir = "~/Dropbox/Development Post SoupX/P14_2/filtered_feature_bc_matrix/")

e12.object <- CreateSeuratObject(counts = e12.data, project = "E12")
e15.object <- CreateSeuratObject(counts = e15.data, project = "E15")
e18.object <- CreateSeuratObject(counts = e18.data, project = "E18")
p0a.object <- CreateSeuratObject(counts = p0a.data, project = "P0a")
p0b.object <- CreateSeuratObject(counts = p0b.data, project = "P0b")
p3.object <- CreateSeuratObject(counts = p3.data, project = "P3")
p5a.object <- CreateSeuratObject(counts = p5a.data, project = "P5a")
p5b.object <- CreateSeuratObject(counts = p5b.data, project = "P5b")
p7a.object <- CreateSeuratObject(counts = p7a.data, project = "P7a")
p7b.object <- CreateSeuratObject(counts = p7b.data, project = "P7b")
p7c.object <- CreateSeuratObject(counts = p7c.data, project = "P7c")
p14a.object <- CreateSeuratObject(counts = p14a.data, project = "P14a")
p14b.object <- CreateSeuratObject(counts = p14b.data, project = "P14b")

#merge prenatal object
prenatal <- merge(e12.object, c(e15.object, e18.object))
prenatal <- PercentageFeatureSet(prenatal, pattern = "^mt-", col.name = "percent.mt")
prenatal <- subset(prenatal, subset = percent.mt < 10 & percent.mt >0.5 & nFeature_RNA >700)

#annotate timepoints in metadata
Idents(prenatal) <- "orig.ident"
Idents(prenatal, cells = WhichCells(prenatal, idents = c("E18"))) <- "E18"
Idents(prenatal, cells = WhichCells(prenatal, idents = c("E15"))) <- "E15"
Idents(prenatal, cells = WhichCells(prenatal, idents = c("E12"))) <- "E12"

prenatal$timepoint <- Idents(prenatal)

#normalize, scale and cluster
prenatal <- SCTransform(prenatal, variable.features.n = 2000, conserve.memory = TRUE, batch_var = "orig.ident", vars.to.regress = c("percent.mt"))
prenatal <- RunPCA(prenatal, verbose = F)
prenatal <- RunUMAP(prenatal, dims = 1:35)
prenatal <- FindNeighbors(prenatal, dims = 1:35)
prenatal <- FindClusters(prenatal, resolution = 0.8)
prenatal <- FindClusters(prenatal, resolution = 0.6)
prenatal <- FindClusters(prenatal, resolution = 1.0)

DimPlot(prenatal, group.by = "timepoint")
DimPlot(prenatal, label = T)

marker_genes =c("Hbb-bs", "Epcam", "Foxj1", "Scgb1a1", "Scgb3a2", "Abca3", "Hopx", "Col1a1", "Dcn", "Lum", "Acta2", "Wnt2", "Wnt5a", "Lgr6", "Pdgfra", "Pdgfrb", "Cspg4", "Wt1", "Pecam1", "Ccl21a", "Vwf", "Nrg1", "Plvap", "Car4", "Mki67", "Tnnt2", "Mpz")
DotPlot(prenatal, features = marker_genes) + theme(axis.text.x = element_text(angle = 45, hjust=1))

#merge postnatal
postnatal <- merge(e18.object, c(p0a.object, p0b.object, p3.object, p5a.object, p5b.object, p7a.object, p7b.object, p7c.object, p14a.object, p14b.object))
postnatal <- PercentageFeatureSet(postnatal, pattern = "^mt-", col.name = "percent.mt")
postnatal <- subset(postnatal, subset = percent.mt < 10 & percent.mt >0.5 & nFeature_RNA >700)

#annotate timepoints in metadata
Idents(postnatal) <- "orig.ident"
Idents(postnatal, cells = WhichCells(postnatal, idents = c("P14a", "P14b"))) <- "P14"
Idents(postnatal, cells = WhichCells(postnatal, idents = c("P7a", "P7b", "P7c"))) <- "P7"
Idents(postnatal, cells = WhichCells(postnatal, idents = c("P5a", "P5b"))) <- "P5"
Idents(postnatal, cells = WhichCells(postnatal, idents = c("P3"))) <- "P3"
Idents(postnatal, cells = WhichCells(postnatal, idents = c("P0a", "P0b"))) <- "P0"
Idents(postnatal, cells = WhichCells(postnatal, idents = c("E18"))) <- "E18"
postnatal$timepoint <- Idents(postnatal)

#normalize, scale and cluster
postnatal <- SCTransform(postnatal, variable.features.n = 2000, conserve.memory = TRUE, batch_var = "orig.ident", vars.to.regress = c("percent.mt"))
postnatal <- RunPCA(postnatal, verbose = F)
postnatal <- RunUMAP(postnatal, dims = 1:35)
postnatal <- FindNeighbors(postnatal, dims = 1:35)
postnatal <- FindClusters(postnatal, resolution = 0.8)
postnatal <- FindClusters(postnatal, resolution = 0.6)
postnatal <- FindClusters(postnatal, resolution = 1.0)

DimPlot(postnatal, group.by = "timepoint")
DimPlot(postnatal, label = T)

marker_genes =c("Hbb-bs", "Epcam", "Foxj1", "Scgb1a1", "Scgb3a2", "Abca3", "Hopx", "Col1a1", "Dcn", "Lum", "Acta2", "Wnt2", "Wnt5a", "Lgr6", "Pdgfra", "Pdgfrb", "Cspg4", "Wt1", "Pecam1", "Ccl21a", "Vwf", "Nrg1", "Plvap", "Car4", "Mki67", "Tnnt2", "Myl7")
DotPlot(postnatal, features = marker_genes) + theme(axis.text.x = element_text(angle = 45, hjust=1))

#save files
saveRDS(prenatal, file = "~/Dropbox/Development Post SoupX/prenatal_072420.rds")
saveRDS(postnatal, file = "~/Dropbox/Development Post SoupX/postnatal_072420.rds")

#remove RBC, doublets from prenatal
prenatal <- subset(prenatal, idents = c(1:25,27,28,31,32,34,35))
prenatal <- RunPCA(prenatal, verbose = F)
prenatal <- RunUMAP(prenatal, dims = 1:35)
prenatal <- FindNeighbors(prenatal, dims = 1:35)
prenatal <- FindClusters(prenatal, resolution = 0.8)

DimPlot(prenatal, group.by = "timepoint")
DimPlot(prenatal, label = T)
DotPlot(prenatal, features = marker_genes) + theme(axis.text.x = element_text(angle = 45, hjust=1))

#prenatal epi subset
prenatal_epi <- subset(prenatal, idents = c(17,19,25))
prenatal_epi <- SCTransform(prenatal_epi, variable.features.n = 2000, conserve.memory = TRUE, batch_var = "orig.ident", vars.to.regress = c("percent.mt"))
prenatal_epi <- RunPCA(prenatal_epi, verbose = F)
prenatal_epi <- RunUMAP(prenatal_epi, dims = 1:10)
prenatal_epi <- FindNeighbors(prenatal_epi, dims = 1:10)
prenatal_epi <- FindClusters(prenatal_epi, resolution = 0.5)

DimPlot(prenatal_epi, group.by = "timepoint")
DimPlot(prenatal_epi, label = T)

prenatal_epi_markers <- FindAllMarkers(prenatal_epi, test.use = "negbinom")

#remove doublets
prenatal_epi <- subset(prenatal_epi, idents = c(0:2, 4:7))
prenatal_epi <- SCTransform(prenatal_epi, variable.features.n = 2000, conserve.memory = TRUE, batch_var = "orig.ident", vars.to.regress = c("percent.mt"))
prenatal_epi <- RunPCA(prenatal_epi, verbose = F)
prenatal_epi <- RunUMAP(prenatal_epi, dims = 1:15)
prenatal_epi <- FindNeighbors(prenatal_epi, dims = 1:15)
prenatal_epi <- FindClusters(prenatal_epi, resolution = 0.5)

DimPlot(prenatal_epi, group.by = "timepoint")
DimPlot(prenatal_epi, label = T)
DotPlot(prenatal_epi, features = marker_genes) + theme(axis.text.x = element_text(angle = 45, hjust=1))

prenatal_epi_markers <- FindAllMarkers(prenatal_epi, test.use = "negbinom")
saveRDS(prenatal_epi, file = "~/Dropbox/Development Post SoupX/prenatal_epi_072620.rds")

#convert to h5ad
library(SeuratDisk)
#Set working directory (probably folder where you are keeping the data for the project)
setwd("~/Dropbox/Development Post SoupX/")

#Set random seed to get reproducible results
set.seed(33)

#convert to h5ad files
prenatal_epi <- UpdateSeuratObject(prenatal_epi)
SaveH5Seurat(prenatal_epi, filename = "prenatal_epi.h5Seurat", overwrite = T)
Convert("prenatal_epi.h5Seurat", dest = "h5ad")



#prenatal mesenchymal subset
prenatal_mesenchyme <- subset(prenatal, idents = c(0:4,6:8,10:16,18,20,21,23,26))
prenatal_mesenchyme <- SCTransform(prenatal_mesenchyme, variable.features.n = 2000, conserve.memory = TRUE, batch_var = "orig.ident", vars.to.regress = c("percent.mt"))
prenatal_mesenchyme <- RunPCA(prenatal_mesenchyme, verbose = F)
prenatal_mesenchyme <- RunUMAP(prenatal_mesenchyme, dims = 1:25)
prenatal_mesenchyme <- FindNeighbors(prenatal_mesenchyme, dims = 1:25)
prenatal_mesenchyme <- FindClusters(prenatal_mesenchyme, resolution = 0.4)

DimPlot(prenatal_mesenchyme, group.by = "timepoint")
DimPlot(prenatal_mesenchyme, label = T)
DotPlot(prenatal_mesenchyme, features = marker_genes) + theme(axis.text.x = element_text(angle = 45, hjust=1))

saveRDS(prenatal_mesenchyme, file = "~/Dropbox/Development Post SoupX/prenatal_mesenchyme_072620.rds")

#convert to h5ad
library(SeuratDisk)
#Set working directory (probably folder where you are keeping the data for the project)
setwd("~/Dropbox/Development Post SoupX/")

#Set random seed to get reproducible results
set.seed(33)

#convert to h5ad files
prenatal_mesenchyme <- UpdateSeuratObject(prenatal_mesenchyme)
SaveH5Seurat(prenatal_mesenchyme, filename = "prenatal_mes.h5Seurat", overwrite = T)
Convert("prenatal_mes.h5Seurat", dest = "h5ad")

#prenatal endothelial subset
prenatal_endo <- subset(prenatal, idents = c(5,9,22,24))
prenatal_endo <- SCTransform(prenatal_endo, variable.features.n = 2000, conserve.memory = TRUE, batch_var = "orig.ident", vars.to.regress = c("percent.mt"))
prenatal_endo <- RunPCA(prenatal_endo, verbose = F)
prenatal_endo <- RunUMAP(prenatal_endo, dims = 1:12)
prenatal_endo <- FindNeighbors(prenatal_endo, dims = 1:12)
prenatal_endo <- FindClusters(prenatal_endo, resolution = 0.4)

DimPlot(prenatal_endo, group.by = "timepoint")
DimPlot(prenatal_endo, label = T)

DotPlot(prenatal_endo, features = marker_genes) + theme(axis.text.x = element_text(angle = 45, hjust=1))

#remove doublets
prenatal_endo <- subset(prenatal_endo, idents = c(1:7,9,10))
prenatal_endo <- SCTransform(prenatal_endo, variable.features.n = 2000, conserve.memory = TRUE, batch_var = "orig.ident", vars.to.regress = c("percent.mt"))
prenatal_endo <- RunPCA(prenatal_endo, verbose = F)
prenatal_endo <- RunUMAP(prenatal_endo, dims = 1:15)
prenatal_endo <- FindNeighbors(prenatal_endo, dims = 1:15)
prenatal_endo <- FindClusters(prenatal_endo, resolution = 0.4)

DimPlot(prenatal_endo, group.by = "timepoint")
DimPlot(prenatal_endo, label = T)

DotPlot(prenatal_endo, features = marker_genes) + theme(axis.text.x = element_text(angle = 45, hjust=1))

saveRDS(prenatal_endo, file = "~/Dropbox/Development Post SoupX/prenatal_endo_072620.rds")


#convert to h5ad files
prenatal_endo <- UpdateSeuratObject(prenatal_endo)
SaveH5Seurat(prenatal_endo, filename = "prenatal_endo.h5Seurat", overwrite = T)
Convert("prenatal_endo.h5Seurat", dest = "h5ad")



##############Post-natal object##################

#postnatal epi object
postnatal_epi <- subset(postnatal, idents = c(10,14,18,24,29,31,34,35,37,39,40))
postnatal_epi <- SCTransform(postnatal_epi, variable.features.n = 2000, conserve.memory = TRUE, batch_var = "orig.ident", vars.to.regress = c("percent.mt"))
postnatal_epi <- RunPCA(postnatal_epi, verbose = F)
postnatal_epi <- RunUMAP(postnatal_epi, dims = 1:15)
postnatal_epi <- FindNeighbors(postnatal_epi, dims = 1:15)
postnatal_epi <- FindClusters(postnatal_epi, resolution = 0.4)

DimPlot(postnatal_epi, group.by = "timepoint")
DimPlot(postnatal_epi, label = T)

DotPlot(postnatal_epi, features = marker_genes) + theme(axis.text.x = element_text(angle = 45, hjust=1))

#remove doublets
postnatal_epi <- subset(postnatal_epi, idents = c(0:6,8,9,11,13:15))
postnatal_epi <- SCTransform(postnatal_epi, variable.features.n = 2000, conserve.memory = TRUE, batch_var = "orig.ident", vars.to.regress = c("percent.mt"))
postnatal_epi <- RunPCA(postnatal_epi, verbose = F)
postnatal_epi <- RunUMAP(postnatal_epi, dims = 1:15)
postnatal_epi <- FindNeighbors(postnatal_epi, dims = 1:15)
postnatal_epi <- FindClusters(postnatal_epi, resolution = 0.4)

DimPlot(postnatal_epi, group.by = "timepoint")
DimPlot(postnatal_epi, label = T)

DotPlot(postnatal_epi, features = marker_genes) + theme(axis.text.x = element_text(angle = 45, hjust=1))

saveRDS(postnatal_epi, file = "~/Dropbox/Development Post SoupX/postnatal_epi_072720.rds")

#convert to h5ad files
#Set working directory (probably folder where you are keeping the data for the project)
setwd("~/Dropbox/Development Post SoupX/")

#Set random seed to get reproducible results
set.seed(33)

postnatal_epi <- UpdateSeuratObject(postnatal_epi)
SaveH5Seurat(postnatal_epi, filename = "postnatal_epi_072720.h5Seurat", overwrite = T)
Convert("postnatal_epi_072720.h5Seurat", dest = "h5ad")


#postnatal mesenchyme object
postnatal_mesenchyme <- subset(postnatal, idents = c(1:6,8,9,11,13,15,16,19,20,23,25,26,28,32,33,38,41,42))
postnatal_mesenchyme <- SCTransform(postnatal_mesenchyme, variable.features.n = 2000, conserve.memory = TRUE, batch_var = "orig.ident", vars.to.regress = c("percent.mt"))
postnatal_mesenchyme <- RunPCA(postnatal_mesenchyme, verbose = F)
postnatal_mesenchyme <- RunUMAP(postnatal_mesenchyme, dims = 1:25)
postnatal_mesenchyme <- FindNeighbors(postnatal_mesenchyme, dims = 1:25)
postnatal_mesenchyme <- FindClusters(postnatal_mesenchyme, resolution = 0.4)

DimPlot(postnatal_mesenchyme, group.by = "timepoint")
DimPlot(postnatal_mesenchyme, label = T)

DotPlot(postnatal_mesenchyme, features = marker_genes) + theme(axis.text.x = element_text(angle = 45, hjust=1))


#remove doublets
postnatal_mesenchyme <- subset(postnatal_mesenchyme, idents = c(0:5,9:18))
postnatal_mesenchyme <- SCTransform(postnatal_mesenchyme, variable.features.n = 2000, conserve.memory = TRUE, batch_var = "orig.ident", vars.to.regress = c("percent.mt"))
postnatal_mesenchyme <- RunPCA(postnatal_mesenchyme, verbose = F)
postnatal_mesenchyme <- RunUMAP(postnatal_mesenchyme, dims = 1:25)
postnatal_mesenchyme <- FindNeighbors(postnatal_mesenchyme, dims = 1:25)
postnatal_mesenchyme <- FindClusters(postnatal_mesenchyme, resolution = 0.4)

DimPlot(postnatal_mesenchyme, group.by = "timepoint")
DimPlot(postnatal_mesenchyme, label = T)

DotPlot(postnatal_mesenchyme, features = marker_genes) + theme(axis.text.x = element_text(angle = 45, hjust=1))

#remove doublets
postnatal_mesenchyme <- subset(postnatal_mesenchyme, idents = c(0:16))
postnatal_mesenchyme <- RunPCA(postnatal_mesenchyme, verbose = F)
postnatal_mesenchyme <- RunUMAP(postnatal_mesenchyme, dims = 1:25)
postnatal_mesenchyme <- FindNeighbors(postnatal_mesenchyme, dims = 1:25)
postnatal_mesenchyme <- FindClusters(postnatal_mesenchyme, resolution = 0.4)

DimPlot(postnatal_mesenchyme, group.by = "timepoint")
DimPlot(postnatal_mesenchyme, label = T)

DotPlot(postnatal_mesenchyme, features = marker_genes) + theme(axis.text.x = element_text(angle = 45, hjust=1))


saveRDS(postnatal_mesenchyme, file = "~/Dropbox/Development Post SoupX/postnatal_mesenchyme_072720.rds")

#convert to h5ad files
#Set working directory (probably folder where you are keeping the data for the project)
setwd("~/Dropbox/Development Post SoupX/")

#Set random seed to get reproducible results
set.seed(33)

#convert to h5ad
postnatal_mesenchyme <- UpdateSeuratObject(postnatal_mesenchyme)
SaveH5Seurat(postnatal_mesenchyme, filename = "postnatal_mesenchyme_072720.h5Seurat", overwrite = T)
Convert("postnatal_mesenchyme_072720.h5Seurat", dest = "h5ad")

#postnatal endothelial object
postnatal_endo<- subset(postnatal, idents = c(0,7,9,12,17,21,22,29,36))
postnatal_endo <- SCTransform(postnatal_endo, variable.features.n = 2000, conserve.memory = TRUE, batch_var = "orig.ident", vars.to.regress = c("percent.mt"))
postnatal_endo <- RunPCA(postnatal_endo, verbose = F)
postnatal_endo <- RunUMAP(postnatal_endo, dims = 1:15)
postnatal_endo <- FindNeighbors(postnatal_endo, dims = 1:15)
postnatal_endo <- FindClusters(postnatal_endo, resolution = 0.4)

DimPlot(postnatal_endo, group.by = "timepoint")
DimPlot(postnatal_endo, label = T)

DotPlot(postnatal_endo, features = marker_genes) + theme(axis.text.x = element_text(angle = 45, hjust=1))

#remove doublets
postnatal_endo<- subset(postnatal_endo, idents = c(0:8,11))
postnatal_endo <- SCTransform(postnatal_endo, variable.features.n = 2000, conserve.memory = TRUE, batch_var = "orig.ident", vars.to.regress = c("percent.mt"))
postnatal_endo <- RunPCA(postnatal_endo, verbose = F)
postnatal_endo <- RunUMAP(postnatal_endo, dims = 1:15)
postnatal_endo <- FindNeighbors(postnatal_endo, dims = 1:15)
postnatal_endo <- FindClusters(postnatal_endo, resolution = 0.4)

DimPlot(postnatal_endo, group.by = "timepoint")
DimPlot(postnatal_endo, label = T)

DotPlot(postnatal_endo, features = marker_genes) + theme(axis.text.x = element_text(angle = 45, hjust=1))

#remove doublets
postnatal_endo<- subset(postnatal_endo, idents = c(0:8,10))
postnatal_endo <- SCTransform(postnatal_endo, variable.features.n = 2000, conserve.memory = TRUE, batch_var = "orig.ident", vars.to.regress = c("percent.mt"))
postnatal_endo <- RunPCA(postnatal_endo, verbose = F)
postnatal_endo <- RunUMAP(postnatal_endo, dims = 1:15)
postnatal_endo <- FindNeighbors(postnatal_endo, dims = 1:15)
postnatal_endo <- FindClusters(postnatal_endo, resolution = 0.4)

DimPlot(postnatal_endo, group.by = "timepoint")
DimPlot(postnatal_endo, label = T)

DotPlot(postnatal_endo, features = marker_genes) + theme(axis.text.x = element_text(angle = 45, hjust=1))


saveRDS(postnatal_endo, file = "~/Dropbox/Development Post SoupX/postnatal_endo_072720.rds")

#convert to h5ad files
#Set working directory (probably folder where you are keeping the data for the project)
setwd("~/Dropbox/Development Post SoupX/")

#Set random seed to get reproducible results
set.seed(33)

#convert to h5ad
postnatal_endo <- UpdateSeuratObject(postnatal_endo)
SaveH5Seurat(postnatal_endo, filename = "postnatal_endo_072720.h5Seurat", overwrite = T)
Convert("postnatal_endo_072720.h5Seurat", dest = "h5ad")









library(SeuratDisk)
#Set working directory (probably folder where you are keeping the data for the project)
setwd("~/Dropbox/Development Post SoupX/")

#Set random seed to get reproducible results
set.seed(33)

#convert to h5ad files

prenatal <- UpdateSeuratObject(prenatal)
SaveH5Seurat(prenatal, filename = "prenatal.h5Seurat", overwrite = T)
Convert("prenatal.h5Seurat", dest = "h5ad")

postnatal <- UpdateSeuratObject(postnatal)
SaveH5Seurat(postnatal, filename = "postnatal.h5Seurat", overwrite = T)
Convert("postnatal.h5Seurat", dest = "h5ad")






#re-scale with all_genes
prenatal_epi <- SCTransform(prenatal_epi, variable.features.n = 30000, conserve.memory = TRUE, batch_var = "orig.ident", vars.to.regress = c("percent.mt"))
prenatal_epi <- RunPCA(prenatal_epi, verbose = F)
prenatal_epi <- RunUMAP(prenatal_epi, dims = 1:15)
prenatal_epi <- FindNeighbors(prenatal_epi, dims = 1:15)
prenatal_epi <- FindClusters(prenatal_epi, resolution = 0.5)

DimPlot(prenatal_epi, group.by = "timepoint")
DimPlot(prenatal_epi, label = T)
DotPlot(prenatal_epi, features = marker_genes) + theme(axis.text.x = element_text(angle = 45, hjust=1))

prenatal_epi_markers <- FindAllMarkers(prenatal_epi, test.use = "negbinom")
saveRDS(prenatal_epi, file = "~/Dropbox/Development Post SoupX/prenatal_epi_081120.rds")

#convert to h5ad
library(SeuratDisk)
#Set working directory (probably folder where you are keeping the data for the project)
setwd("~/Dropbox/Development Post SoupX/")

#Set random seed to get reproducible results
set.seed(33)

#convert to h5ad files
prenatal_epi <- UpdateSeuratObject(prenatal_epi)
SaveH5Seurat(prenatal_epi, filename = "prenatal_epi_081120.h5Seurat", overwrite = T)
Convert("prenatal_epi_081120.h5Seurat", dest = "h5ad")


#prenatal_mesenchyme
postnatal_mesenchyme <- SCTransform(postnatal_mesenchyme, variable.features.n = 30000, conserve.memory = TRUE, batch_var = "orig.ident", vars.to.regress = c("percent.mt"))
postnatal_mesenchyme <- RunPCA(postnatal_mesenchyme, verbose = F)
postnatal_mesenchyme <- RunUMAP(postnatal_mesenchyme, dims = 1:25)
postnatal_mesenchyme <- FindNeighbors(postnatal_mesenchyme, dims = 1:25)
postnatal_mesenchyme <- FindClusters(postnatal_mesenchyme, resolution = 0.4)

DimPlot(postnatal_mesenchyme, group.by = "timepoint")
DimPlot(postnatal_mesenchyme, label = T)

DotPlot(postnatal_mesenchyme, features = marker_genes) + theme(axis.text.x = element_text(angle = 45, hjust=1))


saveRDS(postnatal_mesenchyme, file = "~/Dropbox/Development Post SoupX/postnatal_mesenchyme_081120.rds")

#convert to h5ad files
#Set working directory (probably folder where you are keeping the data for the project)
setwd("~/Dropbox/Development Post SoupX/")

#Set random seed to get reproducible results
set.seed(33)

#convert to h5ad
postnatal_mesenchyme <- UpdateSeuratObject(postnatal_mesenchyme)
SaveH5Seurat(postnatal_mesenchyme, filename = "postnatal_mesenchyme_081120.h5Seurat", overwrite = T)
Convert("postnatal_mesenchyme_081120.h5Seurat", dest = "h5ad")










