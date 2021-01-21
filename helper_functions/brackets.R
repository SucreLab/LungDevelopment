# Title     : TODO
# Objective : TODO
# Created by: Nick
# Created on: 9/30/20

#' Make a grob from grid.brackets
#'
bracketsGrob <- function(...){
  l <- list(...)
  e <- new.env()
  e$l <- l
  grid:::recordGrob(  {
    do.call(grid.brackets, l)
  }, e)
}



addBrackets <- function(annotations){
  total_len <- sum(sapply(annotations, function(x){as.integer(x[2])}))
  scale_factor <- 1 / total_len
  start_x <- (scale_factor / 2) + 0.01
  gap_x <- scale_factor - (scale_factor * start_x)
  out_list <- list()

  for (a in annotations){
    width_x <- ((as.integer(a[2]) - 1) * scale_factor)
    out_list <- append(out_list,
                       annotation_custom(bracketsGrob(start_x, # x1
                                                      1, #y1
                                                      start_x + width_x, #x2
                                                      1, #y2
                                                      h=0.05, lwd=1, col="black", type = 4)))

    start_x <- start_x + width_x + gap_x

  }

  return(out_list)
}
addText <- function(annotations, n_clusters, fontsize = 10){
  total_len <- sum(sapply(annotations, function(x){as.integer(x[2])}))
  scale_factor <- 1 / total_len
  start_x <- (scale_factor / 2) + 0.01
  gap_x <- scale_factor - (scale_factor * start_x)
  out_list <- list()

  for (a in annotations){
    width_x <- ((as.integer(a[2]) - 1) * scale_factor)
    out_list <- append(out_list,
                       annotation_custom(textGrob(a[1], (start_x * 2 + width_x) /2 , # x1
                                                  1.1,
                                                  just = "centre",
                                                  gp = gpar(fontsize = fontsize))))

    start_x <- start_x + width_x + gap_x

  }
  return(out_list)
}