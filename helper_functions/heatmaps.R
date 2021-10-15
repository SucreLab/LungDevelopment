# Title     : TODO
# Objective : TODO
# Created by: Nick
# Created on: 9/30/20

renameFactorIdent <- function(obj, start, end){
  levels(obj)[levels(obj) == start] <- end
  return(obj)
}


fortify.Seurat <- function(x){
  xy <- as.data.frame(Embeddings(x, reduction = "umap"))
  colnames(xy) <- c("x", "y")
  xy$x <- as.numeric(xy$x)
  xy$y <- as.numeric(xy$y)

  return(cbind(xy, as.data.frame(x@meta.data)))
}