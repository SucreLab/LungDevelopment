# Title     : TODO
# Objective : TODO
# Created by: Nick
# Created on: 9/29/20


#' Run the common Seurat clustring, PCA, and UMAP
#'
#' @param dims_umap dimensions to run UMAP
#' @param dims_neighbors How many dims for neighbor finding
#' @param k_param k for finding neighbors
#' @param cluster_res Resolution for clustering
#' @return Seurat object
cluster_pca_umap <- function(obj, dims_umap = 1:15, dims_neighbors = 1:15, k_param = 10, cluster_res = 0.3, umap_neighbors = 30){
  if ("integrated" %in% names(obj)){
    DefaultAssay(obj) <- "integrated"
  }
  obj <- RunPCA(obj, verbose = F)
  obj <- RunUMAP(obj, dims = dims_umap, verbose = F, n.neighbors = umap_neighbors)
  obj <- FindNeighbors(obj, dims = dims_neighbors, k.param = k_param)
  obj <- FindClusters(obj, resolution = cluster_res)
  if ("integrated" %in% names(obj)){
    DefaultAssay(obj) <- "SCT"
  }
  return(obj)
}

#' This can use a ton of memory, but it will save a lot of time
#' Important: The output is the same order as levels(Idents(obj))
#'
#' This needs to export the Seurat object to all workers - this takes a lot of RAM.
#'
#' @param obj Seurat object
#' @param n_cor number of CPU cores
#' @return List of data tables, one for each numeric cluster
parallelFindAllMarkers <- function(obj){

  all_markers <- future_lapply(levels(Idents(obj)), function(x){ # Five expression patterns
    FindMarkers(obj, ident.1 = x, ident.2 = NULL, test.use = "MAST")
  })

  return(value(all_markers))
}