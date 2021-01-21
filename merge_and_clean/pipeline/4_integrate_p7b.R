library(Seurat)
library(future)
library(sctransform)
set.seed(42)

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

integrate <- readRDS(paste0(file_prefix, "_integrate_P7b_sct.rds"))


fn <- paste0(file_prefix, "_integrate_P7b_sct.rds")
if (file.exists(fn)) {
  file.remove(fn)
}


integrate_features <- SelectIntegrationFeatures(object.list = integrate, nfeatures = 2000)


integrate <- PrepSCTIntegration(object.list = integrate, anchor.features = integrate_features,
                                verbose = FALSE)


integrate_anchors <- FindIntegrationAnchors(object.list = integrate, normalization.method = "SCT",
                                            anchor.features = integrate_features, verbose = TRUE)


plan("multiprocess", workers = n_workers)
integrated <- IntegrateData(anchorset = integrate_anchors, normalization.method = "SCT",
                            verbose = FALSE)
## Calc UMAP clusters
plan("sequential")

rm(integrate)
gc()
# Make the tiempoints an ordered factor, so they are displayed consistantly
integrated$timepoint <- ordered(as.factor(integrated$timepoint), unique(integrated$timepoint))




cluster_pca_umap <- function(obj, dims_umap = 1:15, dims_neighbors = 1:15, k_param = 10, cluster_res = 0.3){
  obj <- RunPCA(obj, verbose = F)
  obj <- RunUMAP(obj, dims = dims_umap, verbose = F)
  obj <- FindNeighbors(obj, dims = dims_neighbors, k.param = k_param)
  obj <- FindClusters(obj, resolution = cluster_res)
  return(obj)
}
integrated <- cluster_pca_umap(integrated, k_param = 25, dims_neighbors = 1:15, cluster_res = 0.15) # Note, this breaks if future is set to plan: multiprocess

saveRDS(integrated, file = paste0(file_prefix, "_noGm42418_p7b_integrated_umap.rds"), compress=FALSE)
