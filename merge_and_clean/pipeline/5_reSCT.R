## Retransform values for better plots..
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

integrated <- SCTransform(readRDS(paste0(file_prefix, "_noGm42418_p7b_integrated_umap.rds")), variable.features.n = 2000,
                          conserve.memory = TRUE, batch_var = "orig.ident",
                          vars.to.regress = c("percent.mt"), min_cells = n_drop,
                          method = 'glmGamPoi')

fn <- paste0(file_prefix, "_noGm42418_p7b_integrated_umap.rds")
if (file.exists(fn)) {
  file.remove(fn)
}

DefaultAssay(integrated) <- "SCT"

saveRDS(integrated, file = paste0(file_prefix, "_noGm42418_sct_p7b_integrated_retransform.rds"))