library(Seurat)

epi <- merge(subset(readRDS("./data/input/prenatal_epi_072620.rds"), subset = timepoint != "E18"),
             readRDS(file = "./data/input/postnatal_epi_072720.rds"))
saveRDS(epi, file = "./data/not_cleaned/epi_full.rds")

endo <- merge(subset(readRDS("../data/input/prenatal_endo_072620.rds"), subset = timepoint != "E18"),
              readRDS("./data/input/postnatal_endo_072720.rds"))
saveRDS(epi, file = "./data/not_cleaned/endo_full.rds")

meso <- merge(subset(readRDS("../data/input/prenatal_mesenchyme_072620.rds"), subset = timepoint != "E18"),
              readRDS(file = "../data/input/postnatal_mesenchyme_072720.rds"))
saveRDS(epi, file = "./data/not_cleaned/meso_full.rds")