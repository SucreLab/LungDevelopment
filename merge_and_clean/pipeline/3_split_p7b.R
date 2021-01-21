### Step3
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


integrate <- readRDS(paste0(file_prefix, "_integrate_P7b.rds"))

fn <- paste0(file_prefix, "_integrate_P7b.rds")
if (file.exists(fn)) {
  file.remove(fn)
}

integrate[[1]] <- suppressWarnings(SCTransform(integrate[[1]], variable.features.n = 2000, conserve.memory = TRUE, batch_var = "orig.ident", vars.to.regress = c("percent.mt"), method = 'glmGamPoi'))

integrate[[2]] <- suppressWarnings(SCTransform(integrate[[2]], variable.features.n = 2000, conserve.memory = TRUE, vars.to.regress = c("percent.mt"), method = 'glmGamPoi'))
saveRDS(integrate, paste0(file_prefix, "_integrate_P7b_sct.rds"), compress=FALSE)