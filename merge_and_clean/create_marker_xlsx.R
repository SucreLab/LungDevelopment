setwd("~/postdoc/code/devo_scseq")

source("./helper_functions/libraries.R")

opts_knit$set(root.dir = getwd())

source("./helper_functions/trajectory.R")
source("./helper_functions/cluster.R")
source("./helper_functions/colors.R")
source("./helper_functions/brackets.R")
source("./helper_functions/heatmaps.R")
set.seed(42) # For reproducability

N_WORKERS <- 3
options(future.globals.maxSize=2048*1024^2)
plan("multiprocess", workers = N_WORKERS)


epi_marker_genes <- readRDS("./data/epi_markers_clusters_relabeled.rds")
endo_marker_genes <- readRDS("./data/endo_markers_clusters_relabeled.rds")
meso_marker_genes <- readRDS("./data/meso_markers_clusters_relabeled.rds")

epi_cluster_names <- levels(Idents(readRDS("./data/epi_full_celltype_labels.rds")))
endo_cluster_names <- levels(Idents(readRDS("./data/endo_full_celltype_labels.rds")))
meso_cluster_names <- levels(Idents(readRDS("./data/meso_full_celltype_labels.rds")))


epi_cluster_names_prefix <- paste0("Epi - ", epi_cluster_names)
endo_cluster_names_prefix <- paste0("Endo - ", endo_cluster_names)
meso_cluster_names_prefix <- paste0("Meso - ", meso_cluster_names)

names <- c(epi_cluster_names_prefix, endo_cluster_names_prefix, meso_cluster_names_prefix)
lists <- c(epi_marker_genes, endo_marker_genes, meso_marker_genes)

wb_markers <- createWorkbook()
for (i in c(1:length(names))){
  addWorksheet(wb_markers, names[i])
  writeData(wb_markers, names[i], lists[[i]], rowNames = TRUE)
}
saveWorkbook(wb_markers, file = "./figures/all/marker_genes.xlsx", overwrite = TRUE)

