# Title     : TODO
# Objective : TODO
# Created by: Nick
# Created on: 9/30/20

renameFactorIdent <- function(obj, start, end){
  levels(obj)[levels(obj) == start] <- end
  return(obj)
}


#' Make a wide DF that has some summary data about the clusters, useful for heatmaps
#'
#' @param obj seurat object to run on, needs to have a $timepoint variable
#' @param marker_genes List of marker genes to include
#' @param include_zeros Include zero expression clusters
#' @param sort_genes Use hclust to sort the order of the genes
#' @param sort_clusters Use hclust to sort the order of the clusters
#' @return Returns Wide dataframe
#' @examples
#'
#'heatmap_df <- make_heatmap_df(meso, marker_genes)
#'
#'
#'heatmap_df <- heatmap_df %>% group_by(gene,cluster) %>% summarise(expression_mean=mean(expression))
#'ggplot(heatmap_df, aes(x = gene, y = cluster, fill = expression_mean)) +
#'        geom_tile(color = "white", size = 0.1) +
#'        scale_fill_distiller(palette = "Blues", direction = 1, trans = "sqrt") +
#'        ggtitle("Expression") +
#'        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#'        coord_equal() +
#'
make_heatmap_df <- function(obj, marker_genes, include_zeros = FALSE, sort_genes = FALSE, sort_clusters = TRUE) {
  marker_genes <- marker_genes
  obj <- obj
  timepoint_factor <- factor(obj$timepoint)
  timepoint_ordered <- ordered(timepoint_factor, levels = unique(obj$timepoint))
  cell_ident_df <- data.frame(timepoint = timepoint_ordered,
                              cluster = Idents(obj))
  count_table <- table(cell_ident_df$timepoint, cell_ident_df$cluster)
  prop_celltypes <- round(prop.table(count_table, 1) * 100, 1)


  counts_markers <- t(FetchData(object = obj, slot = "data", vars = marker_genes))
  rownames(counts_markers) <- gsub("^rna_","",rownames(counts_markers))

  heatmap_df <- data.frame()
  # Todo: Replace this with an apply
  for (cluster in levels(Idents(obj))) {
    for (time in unique(obj$timepoint)) {
      if (as.data.frame.matrix(count_table)[time, cluster] > 3) {

        cell_expression <- melt(t(counts_markers[, Idents(obj) == cluster & obj$timepoint == time]))
        colnames(cell_expression) <- c("cell", "gene", "expression")
        n_rep <- nrow(cell_expression)
        cell_expression[,'cluster'] <- rep(cluster, n_rep)
        cell_expression[,'time'] <- rep(time, n_rep)
        cell_expression[,'prop'] <- rep(as.data.frame.matrix(prop_celltypes)[time, cluster] + 1, n_rep)



        heatmap_df <- rbind(heatmap_df, cell_expression)
      } else if (include_zeros){
        gene_rep <- length(marker_genes)
        zero_df <- data.frame(cell = rep("NA", gene_rep), gene = marker_genes, expression = rep(0, gene_rep),
                              cluster = rep(cluster, gene_rep), time = rep(time, gene_rep), prop = rep(0, gene_rep))
        heatmap_df <- rbind(heatmap_df, zero_df)
      }
    }
  }


  # Use hclust to sort the genes by complete linkage
  if (sort_genes) {
    expression_summary_wide <- acast(heatmap_df, cluster+time~gene, value.var = "expression")
    hc <- hclust(dist(t(expression_summary_wide)))
    heatmap_df$gene <- ordered(as.factor(heatmap_df$gene), colnames(expression_summary_wide)[hc$order])
  } else {
    #Order by defined order..
    heatmap_df$gene <- ordered(as.factor(heatmap_df$gene), marker_genes)
  }

  if (sort_clusters) {
    # Use hclust to sort the clusters by complete linkage
    expression_summary_wide_cluster <- acast(heatmap_df, gene + time ~ cluster, value.var = "expression", fun.aggregate = mean)
    expression_summary_wide_cluster[is.na(expression_summary_wide_cluster)] <- 0
    hc_cluster <- hclust(dist(t(expression_summary_wide_cluster)))
    heatmap_df$cluster <- ordered(as.factor(heatmap_df$cluster), colnames(expression_summary_wide_cluster)[hc_cluster$order])
  } else {
    heatmap_df$cluster <- ordered(as.factor(heatmap_df$cluster), unique(heatmap_df$cluster))
  }


  heatmap_df$time <- ordered(as.factor(heatmap_df$time), levels = rev(unique(obj$timepoint)))
  return(heatmap_df)
}



fortify.Seurat <- function(x){
  xy <- as.data.frame(Embeddings(x, reduction = "umap"))
  colnames(xy) <- c("x", "y")
  xy$x <- as.numeric(xy$x)
  xy$y <- as.numeric(xy$y)

  return(cbind(xy, as.data.frame(x@meta.data)))
}