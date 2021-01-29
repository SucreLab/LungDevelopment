library(reticulate)
use_python("~/anaconda3/envs/scseq/bin/python")

library(sctransform)
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(patchwork)
library(ggsignif)
library(viridis)
library(RColorBrewer)
library(reshape2)
library(RColorBrewer)
library(pBrackets)
library(future)
library(future.apply)
library(MAST)
library(grid)
library(boot)
library(openxlsx)
library(magick)

library(knitr) # for kable
options(knitr.table.format = "html")
library(kableExtra) # for pretty tables kable_styling()

# For SCTrnsform
library(devtools)
library(glmGamPoi)

library(ggnewscale)
library(ggbeeswarm)
library(ggtext)

library(anndata)
