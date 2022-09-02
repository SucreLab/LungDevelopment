set.seed(42) # For reproducability
N_WORKERS <- 2
n_print <- 1:20
options(future.globals.maxSize=20*1024*1024^2)

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

#' @param x A sparse matrix from the Matrix package.
#' @param file A filename that ends in ".gz".
#' From: https://slowkow.com/notes/sparse-matrix/#writemmgz
writeMMgz <- function(x, file) {
  mtype <- "real"
  if (is(x, "ngCMatrix")) {
    mtype <- "integer"
  }
  writeLines(
    c(
      sprintf("%%%%MatrixMarket matrix coordinate %s general", mtype),
      sprintf("%s %s %s", x@Dim[1], x@Dim[2], length(x@x))
    ),
    gzfile(file)
  )
  data.table::fwrite(
    x = summary(x),
    file = file,
    append = TRUE,
    sep = " ",
    row.names = FALSE,
    col.names = FALSE
  )
}

pbzip2_location <- Sys.which("pbzip2")
if(pbzip2_location != ""){
  print("pbzip2 available")
}else{
  errorCondition("Please install pbzip2 before running")
}


jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return (intersection/union)
}




saveTiff <- function(path, image, width = 5, height = 5, dpi = 300, units = "in"){
  if (!file.exists(path)){
    dir.create(dirname(path), showWarnings = FALSE)
  }

  if (Sys.info()["sysname"]=='Darwin'){
    # lzw doesn't work on mac with quartz
    tmp_path <- suppressWarnings(normalizePath(paste0(path, "_tmp.tiff")))
    out_path <- suppressWarnings(normalizePath(path))

    tiff(tmp_path, width = width, height = height, units = units, res = dpi, compression = "lzw", type = "quartz")
    print(image)
    dev.off()
    # requires imagemagick
    Sys.sleep(0.5)
    system(paste0("convert ", tmp_path, " -density ", dpi,  " -resize ", width * dpi, "x", height * dpi, "\\> -compress lzw ", out_path), ignore.stdout = TRUE)

    if (file.exists(out_path)) {
      #Delete file if it exists
      file.remove(tmp_path)
    }

  } else {
    tiff(path, width = width, height = height, units = units, res = dpi, compression = "lzw")
    print(image)
    dev.off()
  }

}

# From https://gist.github.com/retrography/359e0cc56d2cf1acd161b5645bc801a8
# The functions below use parallelized versions of gzip, xz, and bzip2 to
# improve compression/decompression performance of RDS serialization in R.
# Each function searches for the appropriate program (based on the required
# compression format) and if found, offloads the compression handling to the
# external program and therefore leaves R free to do the data import/export.
# The two main functions (saveRDS and readRDS) mask R's native read and write
# functions. The functions have been only tested on macOS, but they must work
# on any Linux/Unix.
#
# Requires the following packages: pxz, pbzip2, and pigz.
#
# Run the following line at the command prompt before using the functions.
#
#     brew install pigz pbzip2 pigz
#

#! Todo: Use saveRDS's original structure from base (create new gzfile, etc. functions instead of creating specialized saveRDS functions)
#! Todo: Rename loadRDS and writeRDS to something appropriate

#! Note: Tested on Ubuntu and it all works, but pxz doesn't compile on macOS for now while pbzip2 hangs up on macOS. No idea about Windows for now.

library(parallel)

cmdAvail <- function(cmd) as.logical(nchar(Sys.which(cmd)))

writeRDS <- function(object, con) {
  tryCatch({
    base::saveRDS(
      object,
      file = con
    )
  }, warning = function(w) {
    print(paste("WARNING: ", w))
  }, error = function(e) {
    print(paste("ERROR: ", e))
  }, finally = {
    close(con)
  })
}

loadRDS <- function(con) {
  tryCatch({
    base::readRDS(
      file = con
    )
  }, warning = function(w) {
    print(paste("WARNING: ", w))
  }, error = function(e) {
    print(paste("ERROR: ", e))
  }, finally = {
    close(con)
  })
}

saveRDS.xz <-
  function(object,
           file,
           threads = parallel::detectCores(),
           compression_level = 6) {
    if (cmdAvail("pxz")) {
      writeRDS(
        object,
        pipe(
          paste0(
            "pxz -c -k -T",
            threads,
            " -",
            compression_level,
            " > ",
            file
          ),
          "wb"
        )
      )
    } else {
      base::saveRDS(
        object,
        file = file,
        compress = "xz"
      )
    }
  }

readRDS.xz <-
  function(file,
           threads = parallel::detectCores()) {
    if (cmdAvail("pxz")) {
      object <-
        loadRDS(
          pipe(
            paste0(
              "pxz -d -k -c -T",
              threads,
              " ",
              file
            )
          )
        )
    } else {
      object <-
        base::readRDS(
          file
        )
    }
    return(object)
  }

saveRDS.gz <-
  function(object,
           file,
           threads = parallel::detectCores(),
           compression_level = 6) {
    if (cmdAvail("pigz")) {
      writeRDS(
        object,
        pipe(
          paste0(
            "pigz -c -k -p",
            threads,
            " -",
            compression_level,
            " > ",
            file
          ),
          "wb"
        )
      )
    } else {
      base::saveRDS(
        object,
        file = file,
        compress = "gzip"
      )
    }
  }

readRDS.gz <-
  function(file,
           threads = parallel::detectCores()) {
    if (cmdAvail("pigz")) {
      object <-
        loadRDS(
          pipe(
            paste0(
              "pigz -d -k -c -p",
              threads,
              " ",
              file
            )
          )
        )
    } else {
      object <-
        base::readRDS(
          file
        )
    }
    return(object)
  }

saveRDS.bz2 <-
  function(object,
           file,
           threads = parallel::detectCores(),
           compression_level = 9) {
    if (cmdAvail("pbzip2")) {
      writeRDS(
        object,
        pipe(
          paste0(
            "pbzip2 -c -k -p",
            threads,
            " -",
            compression_level,
            " > ",
            file
          ),
          "wb"
        )
      )
    } else {
      base::saveRDS(
        object,
        file = file,
        compress = "bzip2"
      )
    }
  }

readRDS.bz2 <-
  function(file,
           threads = parallel::detectCores()) {
    if (cmdAvail("pbzip2")) {
      object <-
        loadRDS(
          pipe(
            paste0(
              "pbzip2 -d -k -c -p",
              threads,
              " ",
              file
            )
          )
        )
    } else {
      object <-
        base::readRDS(
          file
        )
    }
    return(object)
  }

readRDS <-
  function(file,
           threads = parallel::detectCores()) {
    if (!file.exists(file)) {
      stop(
        paste0(
          file,
          " does not exist!"
        )
      )
    }
    fileDetails <-
      system2(
        "file",
        args = file,
        stdout = TRUE
      )
    selector <-
      sapply(
        c("gzip", "XZ", "bzip2"),
        function (x) {grepl(x, fileDetails)}
      )
    format <-
      names(selector)[selector]
    if (length(format) == 0) format <- "none"
    if (format == "gzip") {
      object <- readRDS.gz(file, threads = threads)
    } else if (format == "XZ") {
      object <- readRDS.xz(file, threads = threads)
    } else if (format == "bzip2") {
      object <- readRDS.bz2(file, threads = threads)
    } else {
      object <- force(base::readRDS(file))
    }
    return(object)
  }

saveRDS <-
  function(object,
           file = "",
           compress = TRUE) {
    if (compress %in% c(TRUE, "gz", "gzip")) {
      saveRDS.gz(object, file)
    } else if (compress %in% c("bzip", "bzip2", "bz", "bz2")) {
      saveRDS.bz2(object, file)
    } else if (compress %in% c("xz", "7zip", "7z")) {
      saveRDS.xz(object, file)
    } else if (compress == FALSE) {
      base::saveRDS(object, file)
    } else {
      stop(paste0(compress, " is not a recognized compression method!"))
    }
  }