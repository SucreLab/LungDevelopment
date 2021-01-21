library(Seurat)
library(future)
library(sctransform)
set.seed(42)
## Run this in an R-session outside of RStudio

## This takes ~110 gb memory


args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) < 3) {
  stop("Needs to args. <input rds prefix> <int number cells expressing a gene> <n_workers>\n", call.=FALSE)
}

file_prefix <- paste0(args[1])

if (!file.exists(paste0(file_prefix, ".rds"))){
  stop("Input RDS does not exist")
}

n_drop <- args[2]
n_workers <- as.integer(args[3])

plan("multiprocess", workers = n_workers)
options(future.globals.maxSize=12*1024*1024^2) # First num in GB

cluster_pca_umap <- function(obj, dims_umap = 1:15, dims_neighbors = 1:15, k_param = 10, cluster_res = 0.3){
  obj <- RunPCA(obj, verbose = F)
  obj <- RunUMAP(obj, dims = dims_umap, verbose = F)
  obj <- FindNeighbors(obj, dims = dims_neighbors, k.param = k_param)
  obj <- FindClusters(obj, resolution = cluster_res)
  return(obj)
}


dropoutGene <- function(obj, gene){
  CreateSeuratObject(
    (GetAssayData(obj, slot = "counts", assay = "RNA")[rownames(GetAssayData(obj, slot = "counts", assay = "RNA")) != gene,]),
    project = "CreateSeuratObject",
    assay = "RNA",
    meta.data = obj@meta.data,
    min.cells = 1
  )
}
print(n_drop)

merged_alldata <- dropoutGene(readRDS(paste0(file_prefix, ".rds")), "Gm42418")

print(nrow(merged_alldata))

merged_alldata <- PercentageFeatureSet(merged_alldata, pattern = "^mt-", col.name = "percent.mt")
merged_alldata <- subset(merged_alldata, subset = percent.mt < 10 & percent.mt >0.5 & nFeature_RNA >700)

print(nrow(merged_alldata))

saveRDS(merged_alldata, paste0(file_prefix, "_subset.rds"), compress=FALSE)
