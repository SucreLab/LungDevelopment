library(SoupX)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(DropletUtils)
set.seed(33)


args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
  stop("Needs to args. <directory with counts>\n", call.=FALSE)
}

counts_dir <- paste0(args[1])


#Read matrix files into R
epi_epi.raw.matrix <- Read10X(data.dir = counts_dir, gene.column = 2, unique.features = TRUE)
epi_epi.filtered.matrix <- Read10X(data.dir = counts_dir, gene.column = 2, unique.features = TRUE)

#Data come out as a list of two matrices (i.e. 'Gene Expression' and 'Antibody Capture')
#Subset the list above to get the gene expression matrix 
epi_epi.raw.matrix.GE <- epi_epi.raw.matrix$`Gene Expression`
epi_epi.filtered.matrix.GE <- epi_epi.filtered.matrix$`Gene Expression`

#Manually create SoupChannel
epi_epi.sc <- SoupChannel(tod = epi_epi.raw.matrix.GE, toc = epi_epi.filtered.matrix.GE, metaData = NULL, soupRange = c(0, 10), keepDroplets = FALSE)

# tSNE coordinates for epi_epi data
srat_epi_epi <-  CreateSeuratObject(epi_epi.sc$toc)
srat_epi_epi <-  NormalizeData(srat_epi_epi)
srat_epi_epi <-  FindVariableFeatures(srat_epi_epi) 
srat_epi_epi <-  ScaleData(srat_epi_epi) 
srat_epi_epi <-  RunPCA(srat_epi_epi,pcs.compute=30) 
srat_epi_epi <-  RunUMAP(srat_epi_epi,dims=seq(30)) 
srat_epi_epi <-  FindNeighbors(srat_epi_epi,dims=seq(30)) 
srat_epi_epi <-  FindClusters(srat_epi_epi,resolution=1) 
epi_epi_DR <-  as.data.frame(srat_epi_epi@reductions$umap@cell.embeddings)
colnames(epi_epi_DR) = c('RD1','RD2')
epi_epi_DR$Cluster = factor(srat_epi_epi@meta.data[rownames(epi_epi_DR),'RNA_snn_res.1'])


# epi_epi Plotting of marker maps to check for RNA contamination. See if the gene of interest is expessing different clusters
epi_epi_DR$sftpc = epi_epi.sc$toc["Sftpc", rownames(epi_epi_DR)]
epi_epi_DR$scgb1a1 = epi_epi.sc$toc["Scgb1a1", rownames(epi_epi_DR)]
epi_epi_DR$ptprc = epi_epi.sc$toc["Ptprc", rownames(epi_epi_DR)]
epi_epi_DR$pdgfra = epi_epi.sc$toc["Pdgfra", rownames(epi_epi_DR)]
epi_epi_DR$pecam1 = epi_epi.sc$toc["Pecam1", rownames(epi_epi_DR)]
epi_epi_DR$dcn = epi_epi.sc$toc["Dcn", rownames(epi_epi_DR)]
ggplot(epi_epi_DR, aes(RD1, RD2)) + geom_point(aes(colour = sftpc > 0))
ggplot(epi_epi_DR, aes(RD1, RD2)) + geom_point(aes(colour = scgb1a1 > 0))
ggplot(epi_epi_DR, aes(RD1, RD2)) + geom_point(aes(colour = ptprc > 0))
ggplot(epi_epi_DR, aes(RD1, RD2)) + geom_point(aes(colour = pdgfra > 0))
ggplot(epi_epi_DR, aes(RD1, RD2)) + geom_point(aes(colour = pecam1 > 0))
ggplot(epi_epi_DR, aes(RD1, RD2)) + geom_point(aes(colour = dcn > 0))
plotMarkerMap(epi_epi.sc, "Scgb1a1", epi_epi_DR)
plotMarkerMap(epi_epi.sc, "Sftpc", epi_epi_DR)
plotMarkerMap(epi_epi.sc, "Ptprc", epi_epi_DR)
plotMarkerMap(epi_epi.sc, "Pdgfra", epi_epi_DR)
plotMarkerMap(epi_epi.sc, "Pecam1", epi_epi_DR)
plotMarkerMap(epi_epi.sc, "Dcn", epi_epi_DR)

# epi_epi, which genes are expressed most highly in the background? 
head(epi_epi.sc$soupProfile[order(epi_epi.sc$soupProfile$est, decreasing = TRUE), ], n = 20)
#plotMarkerDistribution can be used to visualise how this gene's expression is distributed across cells
plotMarkerDistribution(epi_epi.sc)

#Specify background RNA genes
background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1")
#Specify background RNA genes for epithelial specific cell sort
useToEst_epi_epi = estimateNonExpressingCells(epi_epi.sc, nonExpressedGeneList = list(SLRP = background_RNA_genes), clusters = setNames(epi_epi_DR$Cluster, rownames(epi_epi_DR)))
plotMarkerMap(epi_epi.sc, geneSet = background_RNA_genes, DR = epi_epi_DR, useToEst = useToEst_epi_epi)

## calculate contamination and adjust the counts. 
epi_epi.sc = setClusters(epi_epi.sc, epi_epi_DR$Cluster)
epi_epi.sc = calculateContaminationFraction(epi_epi.sc, list(SLRP = background_RNA_genes), useToEst = useToEst_epi_epi)
head(epi_epi.sc$metaData)

## adjust the counts based on contamination fraction
out_epi_epi = adjustCounts(epi_epi.sc)

## rechecking some plots and genes to check the genes are mostly specific to  clusters. 
cntSoggy = rowSums(epi_epi.sc$toc > 0)
cntStrained = rowSums(out_epi_epi > 0)
mostZeroed_epi_epi = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed_epi_epi
plotChangeMap(epi_epi.sc, out_epi_epi, "Cd74", epi_epi_DR)
plotChangeMap(epi_epi.sc, out_epi_epi, "Scgb3a2", epi_epi_DR)
plotChangeMap(epi_epi.sc, out_epi_epi, "Scgb1a1", epi_epi_DR)
plotChangeMap(epi_epi.sc, out_epi_epi, "Sftpc", epi_epi_DR)
plotChangeMap(epi_epi.sc, out_epi_epi, "Dcn", epi_epi_DR)

### write the file to a new count matrix. saves similar to original one. We use these for starting over in seurat. 
DropletUtils:::write10xCounts(paste0(counts_dir, "/post_soupx/"), out_epi_epi)
                        
                        