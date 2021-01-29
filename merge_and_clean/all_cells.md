scSEQ analysis of all lung cells
================
Nick Negretti
9/30/20

## Load libraries and helper functions

``` r
setwd("~/postdoc/code/devo_scseq_github")
source("./helper_functions/globals.R")
source("./helper_functions/libraries.R")

opts_knit$set(root.dir = getwd())

source("./helper_functions/trajectory.R")
source("./helper_functions/cluster.R")
source("./helper_functions/colors.R")
source("./helper_functions/brackets.R")
source("./helper_functions/heatmaps.R")

plan("multiprocess", workers = N_WORKERS)
```

## UMAP and clustering

``` r
all_cells <- readRDS("./data/merged_alldata_noGm42418_sct_p7b_integrated_retransform.rds")
all_cells$cell_subtype[all_cells$cell_subtype == "Matrix fibroblast"] <- "Adventitial fibroblast"
all_cells$cell_subtype[all_cells$cell_subtype == "Prenatal Wnt2+"] <- "Proliferating Wnt2+"
all_cells$cell_subtype[all_cells$cell_subtype == "Proliferating myofibroblast"] <- "Proliferating myofibroblast"

print(ncol(all_cells))
```

    ## [1] 92238

``` r
timepoint_factor <- factor(all_cells$timepoint)
timepoint_ordered <- ordered(timepoint_factor, levels = unique(all_cells$timepoint))

cell_ident_df <- data.frame(timepoint = timepoint_ordered,
                            cluster = Idents(all_cells))

table(cell_ident_df$timepoint)
```

    ## 
    ##   E12   E15   E18    P0    P3    P5    P7   P14 
    ##  3891 11064  9642 10187 10673 18798 21233  6750

``` r
Command(all_cells)
```

    ## [1] "FindIntegrationAnchors"       "PairwiseIntegrateReference"  
    ## [3] "RunPCA.integrated"            "RunUMAP.integrated.pca"      
    ## [5] "FindNeighbors.integrated.pca" "FindClusters"                
    ## [7] "SCTransform.RNA"

``` r
plan("sequential")

# Make the tiempoints an ordered factor, so they are displayed consistantly
all_cells$timepoint <- ordered(as.factor(all_cells$timepoint), unique(all_cells$timepoint))


#all_cells <- cluster_pca_umap(all_cells, k_param = 25, dims_neighbors = 1:15, cluster_res = 0.15) # Note, this breaks if future is set to plan: multiprocess

p_cluster <- DimPlot(all_cells) + umap_theme() +
        scale_colour_manual(name = "Cluster", values = color_scanpy_viridis20) +
        theme(aspect.ratio=1)
p_cluster
```

![](./merge_and_clean/all_cells_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
p_time <- DimPlot(all_cells, group.by = "timepoint") + umap_theme() +
        scale_colour_manual(name = "Timepoint", values = color_category_20)+
        theme(aspect.ratio=1)

# Function for going from Seurat objects to UMAP projections easily.
fortify.Seurat <- function(x){
  xy <- as.data.frame(Embeddings(x, reduction = "umap"))
  colnames(xy) <- c("x", "y")
  xy$x <- as.numeric(xy$x)
  xy$y <- as.numeric(xy$y)

  return(cbind(xy, as.data.frame(x@meta.data)))
}


all_fortify <- fortify.Seurat(all_cells)
epi_fortify <- fortify.Seurat(subset(all_cells, bulk_celltype == "Epithelium"))
endo_fortify <- fortify.Seurat(subset(all_cells, bulk_celltype == "Endothelium"))
meso_fortify <- fortify.Seurat(subset(all_cells, bulk_celltype == "Mesenchyme"))


epi_fortify$cell_subtype <- ordered(as.factor(epi_fortify$cell_subtype), c("Primordial", "Type I", "Type II", "Cilliated", "Secretory", "Transitional", "Neuroendocrine"))
endo_fortify$cell_subtype <- ordered(endo_fortify$cell_subtype, c("Arterial maEC", "Venous maEC", "miEC", "Car4+", "Proliferating", "Lymphatic"))
meso_fortify$cell_subtype <- ordered(as.factor(meso_fortify$cell_subtype), c("Proliferating Wnt2+", "Wnt2+", "Proliferating Myofibroblast", "Myofibroblast", "Adventitial fibroblast", "Pericyte", "Mesothelium", "Smooth muscle", "Neuron"))


epi_fortify$cell_subtype <- renameFactorIdent(epi_fortify$cell_subtype, "Type I", "AT1")
epi_fortify$cell_subtype <- renameFactorIdent(epi_fortify$cell_subtype, "Type II", "AT2")
meso_fortify$cell_subtype <- renameFactorIdent(meso_fortify$cell_subtype, "Proliferating Wnt2+", "Proliferating *Wnt2*+ FB")
meso_fortify$cell_subtype <- renameFactorIdent(meso_fortify$cell_subtype, "Wnt2+", "*Wnt2*+ FB")
meso_fortify$cell_subtype <- renameFactorIdent(meso_fortify$cell_subtype, "Neuron", "Other")
endo_fortify$cell_subtype <- renameFactorIdent(endo_fortify$cell_subtype, "Car4+", "*Car4*+ EC")
endo_fortify$cell_subtype <- renameFactorIdent(endo_fortify$cell_subtype, "Proliferating", "Proliferating miEC")

all_umap <- ggplot() + umap_theme() +
        geom_point(data = epi_fortify, aes(x = x, y = y, color = cell_subtype), size = .1) +
        scale_colour_manual(name = "Epithelium", values = color_category_20, guide = guide_legend(override.aes = list(size=3), order = 1)) +
        theme(aspect.ratio=1,
              legend.text=element_text(size=14),
              legend.title=element_text(size=14)) +
        new_scale_color() +
          geom_point(data = endo_fortify, aes(x = x, y = y, color = cell_subtype), size = .1) +
          scale_colour_manual(name = "Endothelium", values = color_category_20, guide = guide_legend(override.aes = list(size=3), order = 2)) +
          theme(legend.text=element_markdown(size=14),
                legend.title=element_text(size=14)) +
          guides(colour = guide_legend(override.aes = list(size=3))) +
        new_scale_color() +
          geom_point(data = meso_fortify, aes(x = x, y = y, color = cell_subtype), size = .1) +
          scale_colour_manual(name = "Mesenchyme", values = color_category_20, guide = guide_legend(override.aes = list(size=3), order = 3)) +
          theme(legend.text=element_markdown(size=14),
                legend.title=element_text(size=14)) +
          guides(colour = guide_legend(override.aes = list(size=3))) +
          labs(x = "UMAP_1", y = "UMAP_2") #+
        #new_scale_color() +
        #  stat_ellipse(data = all_fortify, aes(x = x, y = y, color = bulk_celltype), type = "norm", level = 0.99) +
        #  scale_color_manual(name = "Cell type", values = color_scanpy_default)

all_umap
```

    ## Warning: Removed 3622 rows containing missing values (geom_point).

![](./merge_and_clean/all_cells_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
saveTiff("./figures/all/UMAP_time.tiff", p_time,
         width = 11, height = 7)
```

    ## [1] TRUE

``` r
saveTiff("./figures/all/all_umap_nolines.tiff", all_umap,
         width = 11, height = 7)
```

    ## Warning: Removed 3622 rows containing missing values (geom_point).

    ## [1] TRUE

## Small plots of just individual timepoints

``` r
timepoint_str <- "P3"

for (timepoint_str in unique(as.character(all_fortify$timepoint))){
  print(timepoint_str)
all_umap <- ggplot() + umap_theme() +
        geom_point(data = all_fortify, aes(x = x, y = y, color = cell_subtype), size = .1, color = "gray") +
        scale_colour_manual(name = "Epithelium", values = color_category_20, guide = guide_legend(override.aes = list(size=3), order = 1)) +
        theme(aspect.ratio=1,
              legend.position = "none") +
        new_scale_color() +
        geom_point(data = epi_fortify[epi_fortify$timepoint == timepoint_str, ], aes(x = x, y = y, color = "#d62728"), size = .15) +
        scale_colour_manual(name = "Epithelium", values = "#d62728", guide = guide_legend(override.aes = list(size=3), order = 1)) +
        new_scale_color() +
        geom_point(data = endo_fortify[endo_fortify$timepoint == timepoint_str, ], aes(x = x, y = y, color = "#d62728"), size = .15) +
        scale_colour_manual(name = "Endothelium", values = "#d62728", guide = guide_legend(override.aes = list(size=3), order = 2)) +
        guides(colour = guide_legend(override.aes = list(size=3))) +
        new_scale_color() +
        geom_point(data = meso_fortify[meso_fortify$timepoint == timepoint_str, ], aes(x = x, y = y, color = "#d62728"), size = .15) +
        scale_colour_manual(name = "Mesenchyme", values = "#d62728", guide = guide_legend(override.aes = list(size=3), order = 3)) +
        guides(colour = guide_legend(override.aes = list(size=3))) +
        labs(x = "UMAP_1", y = "UMAP_2") +
        ggtitle(timepoint_str) +
        theme(title = element_text(size = 20))
#all_umap

  saveTiff(paste0("./figures/all/time_umap/",timepoint_str,"_umap.tiff"), all_umap,
           width = 12, height = 12, unit = "cm")
}
```

    ## [1] "E12"

    ## [1] "E15"

    ## [1] "E18"

    ## [1] "P0"

    ## [1] "P3"

    ## [1] "P5"

    ## [1] "P7"

    ## [1] "P14"

## Identify marker genes in each cluster

``` r
plan("multiprocess", workers = N_WORKERS)
filename <- "./data/all_retransform_markers.rds"
if (!file.exists(filename)) {
  markers <- parallelFindAllMarkers(all_cells)
  saveRDS(markers, filename)
} else {
  markers <- readRDS(filename)
}
```

### Cluster 0

``` r
markers[[0 + 1]][n_print,]
```

    ##        p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Tgfbi      0  3.472482 0.886 0.069         0
    ## Tagln      0  2.788986 0.889 0.128         0
    ## Acta2      0  2.388155 0.907 0.322         0
    ## Agt        0  2.127518 0.775 0.020         0
    ## Aspn       0  2.027304 0.341 0.126         0
    ## Pdlim3     0  2.017052 0.817 0.042         0
    ## Actg2      0  1.962395 0.672 0.013         0
    ## Mustn1     0  1.902027 0.668 0.036         0
    ## Htra1      0  1.888923 0.756 0.080         0
    ## Ckb        0  1.665731 0.880 0.224         0
    ## Myh11      0  1.629541 0.733 0.029         0
    ## Myl9       0  1.601325 0.857 0.257         0
    ## Spon2      0  1.599700 0.754 0.129         0
    ## Eln        0  1.597175 0.793 0.445         0
    ## Egfem1     0  1.557632 0.711 0.025         0
    ## Loxl2      0  1.421000 0.783 0.227         0
    ## Hhip       0  1.401656 0.430 0.126         0
    ## Tnc        0  1.388665 0.709 0.091         0
    ## Enpp2      0  1.378978 0.446 0.106         0
    ## Des        0  1.354576 0.675 0.092         0

### Cluster 1

``` r
markers[[1 + 1]][n_print,]
```

    ##         p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Plvap       0  2.393468 0.924 0.037         0
    ## Calcrl      0  2.043811 0.867 0.063         0
    ## Ramp2       0  2.000225 0.988 0.199         0
    ## Hpgd        0  1.941152 0.698 0.071         0
    ## Gpihbp1     0  1.937330 0.807 0.043         0
    ## Egfl7       0  1.913772 0.980 0.079         0
    ## Clec14a     0  1.855245 0.873 0.058         0
    ## Cdh5        0  1.841741 0.950 0.068         0
    ## Cd93        0  1.836656 0.886 0.040         0
    ## Ptprb       0  1.784695 0.795 0.026         0
    ## Cav1        0  1.767221 0.953 0.255         0
    ## Lyve1       0  1.758559 0.722 0.049         0
    ## Tspan7      0  1.746028 0.882 0.161         0
    ## Pecam1      0  1.745703 0.934 0.066         0
    ## Tm4sf1      0  1.736069 0.791 0.136         0
    ## Ly6e        0  1.663498 0.914 0.176         0
    ## Tmem100     0  1.636846 0.760 0.128         0
    ## Cldn5       0  1.607527 0.926 0.067         0
    ## Gja4        0  1.605778 0.604 0.047         0
    ## Esam        0  1.560258 0.905 0.133         0

### Cluster 2

``` r
markers[[2 + 1]][n_print,]
```

    ##         p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Mfap4       0  2.746768 0.999 0.440         0
    ## Gyg         0  2.302386 0.980 0.200         0
    ## Tcf21       0  2.263327 0.959 0.190         0
    ## Fhl1        0  2.092471 0.987 0.456         0
    ## Fn1         0  2.063467 0.947 0.243         0
    ## Macf1       0  1.919245 0.929 0.305         0
    ## Adh1        0  1.852411 0.958 0.170         0
    ## Cebpb       0  1.811627 0.932 0.278         0
    ## Limch1      0  1.752607 0.932 0.135         0
    ## Npnt        0  1.638729 0.872 0.145         0
    ## Enpep       0  1.631089 0.811 0.094         0
    ## G0s2        0  1.625537 0.798 0.121         0
    ## Lbh         0  1.429065 0.896 0.146         0
    ## Col13a1     0  1.345592 0.881 0.108         0
    ## Meox2       0  1.333883 0.870 0.109         0
    ## Cyr61       0  1.331456 0.822 0.432         0
    ## Olfml3      0  1.282885 0.950 0.342         0
    ## Clec3b      0  1.276261 0.829 0.107         0
    ## Plxdc2      0  1.230623 0.902 0.250         0
    ## Cdh11       0  1.212182 0.923 0.280         0

### Cluster 3

``` r
markers[[3 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Dcn          0  2.787013 0.660 0.072         0
    ## Col1a1       0  2.449949 0.970 0.493         0
    ## Col1a2       0  2.229123 0.975 0.541         0
    ## Col3a1       0  2.110465 0.940 0.513         0
    ## Lum          0  2.054010 0.471 0.026         0
    ## Mfap5        0  1.941410 0.681 0.093         0
    ## Gsn          0  1.859978 0.896 0.334         0
    ## Dlk1         0  1.768924 0.540 0.165         0
    ## Serpinf1     0  1.764062 0.680 0.061         0
    ## Ptn          0  1.682809 0.555 0.086         0
    ## Clec3b       0  1.580207 0.665 0.163         0
    ## Pcolce       0  1.490010 0.882 0.335         0
    ## Meg3         0  1.437541 0.837 0.317         0
    ## Col14a1      0  1.430573 0.725 0.149         0
    ## Mgp          0  1.336125 0.890 0.527         0
    ## Itm2a        0  1.319069 0.753 0.208         0
    ## Gpc3         0  1.312669 0.920 0.330         0
    ## Ogn          0  1.306638 0.781 0.193         0
    ## Igfbp4       0  1.299719 0.959 0.645         0
    ## Sparc        0  1.266589 0.997 0.889         0

### Cluster 4

``` r
markers[[4 + 1]][n_print,]
```

    ##           p_val avg_logFC pct.1 pct.2 p_val_adj
    ## H2afz         0 1.2148378 0.898 0.614         0
    ## Mif           0 1.1906708 0.875 0.446         0
    ## Hist1h2ap     0 1.0652510 0.421 0.087         0
    ## Npm1          0 1.0502341 0.963 0.769         0
    ## Hmgb2         0 0.9821023 0.674 0.313         0
    ## Pclaf         0 0.9712766 0.515 0.080         0
    ## Ran           0 0.9694105 0.893 0.532         0
    ## Tuba1b        0 0.9115440 0.948 0.730         0
    ## Hist1h2ae     0 0.9088379 0.405 0.069         0
    ## Hmgn2         0 0.8797202 0.938 0.626         0
    ## Gapdh         0 0.8595688 0.965 0.836         0
    ## Hnrnpa1       0 0.8512648 0.937 0.624         0
    ## Hspd1         0 0.8340794 0.799 0.419         0
    ## Eno1          0 0.8270343 0.892 0.556         0
    ## Anp32b        0 0.8167900 0.860 0.461         0
    ## Hmgb1         0 0.8149364 0.979 0.904         0
    ## Rps2          0 0.8024562 0.992 0.941         0
    ## Rps26         0 0.8021146 0.979 0.940         0
    ## Rpsa          0 0.7972935 0.992 0.950         0
    ## Hist1h1b      0 0.7881522 0.337 0.051         0

### Cluster 5

``` r
markers[[5 + 1]][n_print,]
```

    ##         p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Sftpc       0  6.559549 0.999 0.014         0
    ## Sftpa1      0  4.542499 0.968 0.004         0
    ## Lyz2        0  4.441398 0.809 0.005         0
    ## Sftpb       0  4.377879 0.986 0.005         0
    ## Cxcl15      0  4.163704 0.974 0.006         0
    ## Lyz1        0  3.683135 0.379 0.001         0
    ## Slc34a2     0  3.166963 0.967 0.001         0
    ## Napsa       0  3.162398 0.988 0.005         0
    ## Chil1       0  3.086301 0.903 0.004         0
    ## Cd74        0  3.030069 0.909 0.021         0
    ## Lcn2        0  2.683242 0.773 0.003         0
    ## Wfdc2       0  2.671769 0.974 0.040         0
    ## Sftpd       0  2.645028 0.938 0.007         0
    ## Npc2        0  2.548953 0.990 0.501         0
    ## Cbr2        0  2.337539 0.868 0.020         0
    ## Dram1       0  2.274046 0.972 0.020         0
    ## Lpcat1      0  2.246523 0.964 0.106         0
    ## Fabp5       0  2.086988 0.872 0.120         0
    ## S100g       0  2.062303 0.865 0.003         0
    ## Muc1        0  2.049581 0.971 0.011         0

### Cluster 6

``` r
markers[[6 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Cox4i2       0  3.045889 0.977 0.030         0
    ## Postn        0  2.951741 0.989 0.119         0
    ## Gucy1a1      0  2.893545 0.989 0.173         0
    ## Gucy1b1      0  2.593526 0.979 0.133         0
    ## Ndufa4l2     0  2.254671 0.952 0.074         0
    ## Higd1b       0  2.183155 0.886 0.009         0
    ## Mfge8        0  2.171739 0.991 0.556         0
    ## Pde5a        0  2.166465 0.956 0.240         0
    ## Pdgfrb       0  2.117362 0.955 0.093         0
    ## Itga1        0  2.069755 0.955 0.173         0
    ## Itm2a        0  2.040222 0.958 0.235         0
    ## Emid1        0  1.936840 0.922 0.142         0
    ## Pdzd2        0  1.848383 0.837 0.051         0
    ## Pcdh18       0  1.831570 0.910 0.065         0
    ## Ebf1         0  1.775160 0.923 0.064         0
    ## Lhfp         0  1.763600 0.975 0.464         0
    ## Mcam         0  1.712261 0.936 0.186         0
    ## Trpc6        0  1.696690 0.848 0.012         0
    ## Lmcd1        0  1.611507 0.839 0.145         0
    ## Fermt2       0  1.576541 0.969 0.560         0

### Cluster 7

``` r
markers[[7 + 1]][n_print,]
```

    ##         p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Car4        0  3.242374 0.983 0.035         0
    ## Ly6a        0  2.368690 0.967 0.144         0
    ## Cyp4b1      0  2.309966 0.957 0.099         0
    ## Cldn5       0  2.250818 0.987 0.201         0
    ## Kdr         0  2.230425 0.970 0.150         0
    ## Kitl        0  2.121693 0.959 0.245         0
    ## Icam2       0  2.092494 0.977 0.199         0
    ## Igfbp7      0  2.073530 0.995 0.481         0
    ## Pmp22       0  1.971059 0.991 0.651         0
    ## Scn7a       0  1.908134 0.964 0.179         0
    ## Tspan13     0  1.906367 0.964 0.200         0
    ## Ly6c1       0  1.857931 0.911 0.148         0
    ## Emp2        0  1.852273 0.944 0.371         0
    ## Ednrb       0  1.829835 0.911 0.138         0
    ## Ramp2       0  1.823048 0.997 0.324         0
    ## Acvrl1      0  1.811737 0.979 0.246         0
    ## Apln        0  1.758652 0.827 0.039         0
    ## Egfl7       0  1.699209 0.994 0.221         0
    ## Ecscr       0  1.694582 0.946 0.202         0
    ## Nrp1        0  1.672349 0.960 0.380         0

### Cluster 8

``` r
markers[[8 + 1]][n_print,]
```

    ##           p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Hist1h2ap     0  2.078320 0.767 0.105         0
    ## Hmgb2         0  1.736637 0.978 0.334         0
    ## Top2a         0  1.523614 0.815 0.083         0
    ## Ube2c         0  1.511580 0.666 0.072         0
    ## Stmn1         0  1.414663 0.936 0.299         0
    ## Egfl7         0  1.385533 0.984 0.228         0
    ## Tuba1b        0  1.385145 0.989 0.747         0
    ## Rrm2          0  1.384600 0.729 0.086         0
    ## H2afz         0  1.359487 0.987 0.635         0
    ## H2afx         0  1.332676 0.807 0.175         0
    ## Plvap         0  1.322847 0.957 0.183         0
    ## Arl6ip1       0  1.294890 0.822 0.467         0
    ## Lmnb1         0  1.284930 0.910 0.234         0
    ## Birc5         0  1.282021 0.805 0.094         0
    ## Hist1h2ae     0  1.200004 0.630 0.090         0
    ## Cdk1          0  1.199272 0.758 0.090         0
    ## Tubb5         0  1.193092 0.986 0.785         0
    ## Lyve1         0  1.180325 0.841 0.157         0
    ## Cenpf         0  1.171917 0.636 0.066         0
    ## Prc1          0  1.160593 0.666 0.058         0

### Cluster 9

``` r
markers[[9 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Upk3b        0  3.190390 0.994 0.002         0
    ## Msln         0  3.065888 0.704 0.016         0
    ## C3           0  3.036491 0.758 0.016         0
    ## Igfbp6       0  2.624369 0.928 0.052         0
    ## Igfbp5       0  2.523861 0.934 0.109         0
    ## Sfrp2        0  2.253994 0.796 0.089         0
    ## Clu          0  2.226662 0.906 0.091         0
    ## Aldh1a2      0  2.168606 0.954 0.022         0
    ## Hspb1        0  2.075570 0.960 0.228         0
    ## Crip1        0  2.067946 0.914 0.451         0
    ## Upk1b        0  2.007552 0.939 0.058         0
    ## Rarres2      0  1.993372 0.884 0.374         0
    ## Slurp1       0  1.968390 0.292 0.000         0
    ## Slc9a3r1     0  1.877850 0.948 0.060         0
    ## Krt19        0  1.858813 0.914 0.047         0
    ## Aebp1        0  1.728846 0.923 0.188         0
    ## C2           0  1.704690 0.791 0.038         0
    ## S100a6       0  1.702469 0.840 0.216         0
    ## Gpc3         0  1.689129 0.983 0.383         0
    ## Gbp2         0  1.652648 0.580 0.058         0

### Cluster 10

``` r
markers[[10 + 1]][n_print,]
```

    ##        p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Ager       0  3.871486 0.997 0.087         0
    ## Aqp5       0  3.062096 0.977 0.047         0
    ## Cldn18     0  3.052421 0.993 0.093         0
    ## S100a6     0  2.720783 0.984 0.215         0
    ## Krt7       0  2.675420 0.984 0.060         0
    ## Gprc5a     0  2.372778 0.960 0.064         0
    ## Emp2       0  2.308915 0.995 0.378         0
    ## Clic3      0  2.305374 0.962 0.069         0
    ## Hopx       0  2.289320 0.954 0.057         0
    ## Ndnf       0  2.254551 0.919 0.044         0
    ## Bcam       0  2.244909 0.979 0.321         0
    ## Akap5      0  2.073886 0.938 0.017         0
    ## Col4a4     0  2.060922 0.931 0.020         0
    ## Rtkn2      0  2.021920 0.912 0.010         0
    ## Icam1      0  1.998542 0.956 0.155         0
    ## Igfbp2     0  1.916123 0.657 0.050         0
    ## Col4a3     0  1.887803 0.883 0.012         0
    ## Vegfa      0  1.848410 0.906 0.119         0
    ## Scnn1g     0  1.802897 0.857 0.006         0
    ## Anxa3      0  1.760991 0.973 0.280         0

### Cluster 11

``` r
markers[[11 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Scgb1a1      0  5.852857 0.327 0.010         0
    ## Scgb3a2      0  3.959930 0.261 0.004         0
    ## Cyp2f2       0  2.981545 0.662 0.004         0
    ## Cbr2         0  2.946341 0.735 0.083         0
    ## Sec14l3      0  2.901414 0.719 0.032         0
    ## Hp           0  2.843043 0.640 0.033         0
    ## AU040972     0  2.678771 0.445 0.001         0
    ## Wfdc2        0  2.263642 0.936 0.107         0
    ## Dynlrb2      0  2.183624 0.521 0.002         0
    ## Mt1          0  2.098952 0.570 0.118         0
    ## Cd24a        0  2.059444 0.860 0.147         0
    ## Tmem212      0  1.972326 0.496 0.001         0
    ## Tppp3        0  1.945488 0.520 0.108         0
    ## Cyp2s1       0  1.931094 0.528 0.010         0
    ## Foxj1        0  1.918776 0.517 0.001         0
    ## Tubb4b       0  1.885993 0.801 0.454         0
    ## Cldn3        0  1.826430 0.931 0.100         0
    ## Krt8         0  1.787449 0.916 0.071         0
    ## Chchd10      0  1.775294 0.773 0.096         0
    ## Anxa1        0  1.712751 0.730 0.127         0

``` r
# UMAP:
# All cells grayed out except E15 / another except P3.


ggplot_data <- as.data.frame(Embeddings(all_cells, reduction = "umap"))
ggplot_data$timepoint <- all_cells$timepoint
ggplot_data$cell_subtype <- all_cells$cell_subtype
ggplot_data$bulk_celltype <- all_cells$bulk_celltype


ggdata_epi <- subset(ggplot_data, bulk_celltype == "Epithelium")
ggdata_endo <- subset(ggplot_data, bulk_celltype == "Endothelium")
ggdata_meso <- subset(ggplot_data, bulk_celltype == "Mesenchyme")
ggdata_meso$cell_subtype[ggdata_meso$cell_subtype == "Matrix fibroblast"] <- "Adventitial fibroblast"


ggdata_epi$cell_subtype <- ordered(as.factor(ggdata_epi$cell_subtype), c("Primordial", "Transitional", "AT1", "AT2", "Cilliated", "Secretory", "Neuroendocrine"))
ggdata_endo$cell_subtype <- ordered(as.factor(ggdata_endo$cell_subtype), c("Arterial maEC", "Venous maEC", "miEC", "Car4+", "Proliferating", "Lymphatic"))
# Note - the next few lines combines the prenatal Wnt2+ and the Proliferating myofibroblasts
ggmeso_subtype_factor <- as.factor(ggdata_meso$cell_subtype)
levels(ggmeso_subtype_factor) <- levels(ggmeso_subtype_factor)[c(1:5,3,9,8:9)]
ggdata_meso$cell_subtype <- ordered(ggmeso_subtype_factor, c("Wnt2+", "Myofibroblast", "Adventitial fibroblast", "Pericyte", "Mesothelium", "Smooth muscle", "Neuron"))


ggdata_meso$cell_subtype <- renameFactorIdent(ggdata_meso$cell_subtype, "Wnt2+", "*Wnt2*+ FB")

ggdata_epi <- ggdata_epi[ggdata_epi$cell_subtype %in% c("Primordial", "Transitional"),]
epi_colors <- c("#1f77b4", "#aec7e8", "#2ca02c", "#98df8a")
ggdata_endo <- ggdata_endo[ggdata_endo$cell_subtype %in% c("miEC", "Car4+"),]
endo_colors <- c("#ff7f0e", "#ffbb78")
ggdata_meso <- ggdata_meso[ggdata_meso$cell_subtype %in% c("Myofibroblast", "*Wnt2*+ FB"),]
meso_colors <- c("#d62728", "#ff9896")



e15_highlight <- ggplot() +
        geom_point(data = ggplot_data[!(ggplot_data$timepoint %in% c("E15")),],
                   aes(x = UMAP_1, y = UMAP_2), size = 0.1, alpha = 0.2, color = "gray") +
        umap_theme() +
        scale_colour_manual(name = "Timepoint", values = color_scanpy_default) +
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3))) +
        new_scale_color() +

        geom_point(data = ggdata_epi[(ggdata_epi$timepoint %in% c("E15")),], aes(x = UMAP_1, y = UMAP_2, color = cell_subtype), size = .1) +
        scale_colour_manual(name = "Epithelium", values = epi_colors, guide = guide_legend(override.aes = list(size=3), order = 1)) +
        theme(aspect.ratio=1,
              legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        new_scale_color() +

        geom_point(data = ggdata_endo[(ggdata_endo$timepoint %in% c("E15")),], aes(x = UMAP_1, y = UMAP_2, color = cell_subtype), size = .1) +
        scale_colour_manual(name = "Endothelium", values = endo_colors, guide = guide_legend(override.aes = list(size=3), order = 2)) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3))) +
        new_scale_color() +

        geom_point(data = ggdata_meso[(ggdata_endo$timepoint %in% c("E15")),], aes(x = UMAP_1, y = UMAP_2, color = cell_subtype), size = .1) +
        scale_colour_manual(name = "Mesenchyme", values = meso_colors, guide = guide_legend(override.aes = list(size=3), order = 3)) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3))) +
        labs(x = "UMAP_1", y = "UMAP_2") +

        ggtitle("E15 cells")

e15_highlight
```

![](./merge_and_clean/all_cells_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
saveTiff("./figures/all/e15_highlight_UMAP_revised.tiff", e15_highlight,
         width = 9, height = 7)
```

    ## [1] TRUE

``` r
ggdata_epi <- subset(ggplot_data, bulk_celltype == "Epithelium")
ggdata_endo <- subset(ggplot_data, bulk_celltype == "Endothelium")
ggdata_meso <- subset(ggplot_data, bulk_celltype == "Mesenchyme")
ggdata_meso$cell_subtype[ggdata_meso$cell_subtype == "Matrix fibroblast"] <- "Adventitial fibroblast"


ggdata_epi$cell_subtype <- ordered(as.factor(ggdata_epi$cell_subtype), c("Primordial", "Transitional", "Type I", "Type II", "Cilliated", "Secretory", "Neuroendocrine"))
ggdata_endo$cell_subtype <- ordered(as.factor(ggdata_endo$cell_subtype), c("Arterial maEC", "Venous maEC", "miEC", "Car4+", "Proliferating", "Lymphatic"))
# Note - the next few lines combines the prenatal Wnt2+ and the Proliferating myofibroblasts
ggmeso_subtype_factor <- as.factor(ggdata_meso$cell_subtype)
levels(ggmeso_subtype_factor) <- levels(ggmeso_subtype_factor)[c(1:5,3,9,8:9)]
ggdata_meso$cell_subtype <- ordered(ggmeso_subtype_factor, c("Wnt2+", "Myofibroblast", "Adventitial fibroblast", "Pericyte", "Mesothelium", "Smooth muscle", "Neuron"))



ggdata_epi <- ggdata_epi[ggdata_epi$cell_subtype %in% c("Primordial", "Transitional", "Type I", "Type II"),]
epi_colors <- c("#1f77b4", "#aec7e8", "#2ca02c", "#98df8a")
ggdata_endo <- ggdata_endo[ggdata_endo$cell_subtype %in% c("miEC", "Car4+"),]
endo_colors <- c("#ff7f0e", "#ffbb78")
ggdata_meso <- ggdata_meso[ggdata_meso$cell_subtype %in% c("Myofibroblast", "*Wnt2*+ FB"),]
meso_colors <- c("#d62728", "#ff9896")

ggdata_endo$cell_subtype <- renameFactorIdent(ggdata_endo$cell_subtype, "Car4+", "*Car4*+")
ggdata_epi$cell_subtype <- renameFactorIdent(ggdata_epi$cell_subtype, "Type I", "AT1")
ggdata_epi$cell_subtype <- renameFactorIdent(ggdata_epi$cell_subtype, "Type II", "AT2")

p3_highlight <- ggplot() +
        geom_point(data = ggplot_data[!(ggplot_data$timepoint %in% c("P3")),],
                   aes(x = UMAP_1, y = UMAP_2), size = 0.1, alpha = 0.2, color = "gray") +
        umap_theme() +
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3))) +
        new_scale_color() +

        geom_point(data = ggdata_epi[(ggdata_epi$timepoint %in% c("P3")),], aes(x = UMAP_1, y = UMAP_2, color = cell_subtype), size = .1) +
        scale_colour_manual(name = "Epithelium", values = epi_colors, guide = guide_legend(override.aes = list(size=3), order = 1)) +
        theme(aspect.ratio=1,
              legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        new_scale_color() +

        geom_point(data = ggdata_endo[(ggdata_endo$timepoint %in% c("P3")),], aes(x = UMAP_1, y = UMAP_2, color = cell_subtype), size = .1) +
        scale_colour_manual(name = "Endothelium", values = endo_colors, guide = guide_legend(override.aes = list(size=3), order = 2)) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3))) +
        new_scale_color() +

        geom_point(data = ggdata_meso[(ggdata_endo$timepoint %in% c("P3")),], aes(x = UMAP_1, y = UMAP_2, color = cell_subtype), size = .1) +
        scale_colour_manual(name = "Mesenchyme", values = meso_colors, guide = guide_legend(override.aes = list(size=3), order = 3)) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3))) +
        labs(x = "UMAP_1", y = "UMAP_2") +

        ggtitle("P3 cells")

p3_highlight
```

![](./merge_and_clean/all_cells_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
saveTiff("./figures/all/p3_highlight_UMAP_revised.tiff", p3_highlight,
         width = 9, height = 7)
```

    ## [1] TRUE

# Create a UMAP by time

``` r
all_umap_bytime <- ggplot() + umap_theme() +
        geom_point(data = all_cells, aes(x = x, y = y, color = timepoint), size = .1) +
        scale_colour_manual(name = "Timepoint", values = color_category_20, guide = guide_legend(override.aes = list(size=3), order = 1)) +
        theme(aspect.ratio=1,
              legend.text=element_text(size=14),
              legend.title=element_text(size=14))+
        xlab("UMAP_1") +
        ylab("UMAP_2")

#new_scale_color() +
#  stat_ellipse(data = all_fortify, aes(x = x, y = y, color = bulk_celltype), type = "norm", level = 0.99) +
#  scale_color_manual(name = "Cell type", values = color_scanpy_default)
all_umap_bytime
```

![](./merge_and_clean/all_cells_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
saveTiff("./figures/all/p3_highlight_UMAP_revised.tiff", all_umap_bytime,
         width = 9, height = 7)
```

    ## [1] TRUE

``` r
all_cells$leiden <- Idents(all_cells)
Idents(all_cells) <- all_cells$cell_subtype
Idents(all_cells) <- renameFactorIdent(Idents(all_cells), "Neuron", "Other")
Idents(all_cells) <- renameFactorIdent(Idents(all_cells), "Wnt2+", "*Wnt2*+ FB")
Idents(all_cells) <- renameFactorIdent(Idents(all_cells), "Proliferating Wnt2+", "Proliferating *Wnt2*+ FB")
Idents(all_cells) <- renameFactorIdent(Idents(all_cells), "Prenatal Myofibroblast", "Proliferating myofibroblast")
Idents(all_cells) <- renameFactorIdent(Idents(all_cells), "Type I", "AT1")
Idents(all_cells) <- renameFactorIdent(Idents(all_cells), "Type II", "AT2")
Idents(all_cells) <- renameFactorIdent(Idents(all_cells), "Car4+", "*Car4*+")
Idents(all_cells) <- renameFactorIdent(Idents(all_cells), "Proliferating", "Proliferating miEC")
```

# Create a giant heatmap of every cluster..

``` r
marker_genes <- c("Epcam", "Foxj1", "Scgb1a1", "Scgb3a2", "Abca3", "Hopx",
                  "Pecam1", "Ccl21a", "Vwf", "Nrg1", "Plvap", "Car4",
                  "Col1a1", "Dcn", "Lum", "Acta2", "Wnt2", "Wnt5a", "Lgr6", "Pdgfra", "Pdgfrb", "Cspg4", "Wt1",
                  "Mki67", "Tnnt2", "Mpz")


heatmap_df <- make_heatmap_df(all_cells, marker_genes, sort_clusters = FALSE)

heatmap_df <- heatmap_df %>% group_by(gene,cluster) %>% summarise(expression_mean=mean(expression))
```

    ## `summarise()` regrouping output by 'gene' (override with `.groups` argument)

``` r
# For brackets and text (see ./helper_functions/brackets.R)
n_clusters <- max(as.numeric(Idents(all_cells)))
text_vec <- list(c("Epithelial", 6),
                 c("Endothelial", 6),
                 c("Mesenchymal", 11),
                 c("Other", 3)
)

all_heatmap <- ggplot(heatmap_df, aes(x = gene, y = cluster, fill = expression_mean)) +
        geom_tile(color = "white", size = 0.1) +
        scale_fill_distiller(palette = "Blues", direction = 1, trans = "sqrt", name = "Expression") +
        coord_fixed(ratio = 1, xlim = NULL, ylim = c(1,n_clusters), expand = TRUE, clip = "off") +
        theme(plot.margin=unit(c(1.5,1,1,1),"cm")) +
        addBrackets(text_vec) +
        addText(text_vec, n_clusters) +
        theme(panel.background = element_rect(fill = "white", colour = "black", size = 0),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title.x  = element_blank(),
              axis.text.x = element_text(face = "italic", size = 10),
              axis.text.y = element_markdown(size = 10)
              ) +
        labs(y = "Cluster") +
        theme(axis.text.x = element_text(angle = 45, hjust=1))

all_heatmap
```

![](./merge_and_clean/all_cells_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
saveTiff("./figures/all/heatmap_allgenes.tiff", all_heatmap,
         width = 9, height = 7)
```

    ## [1] TRUE

# Elastin Heatmaps\!\!

``` r
Idents(all_cells) <- all_cells$cell_subtype
marker_genes <- c('Eln', 'Fbn1', 'Fbn2', 'Fbln1', 'Fbln2', 'Fbln3', 'Fbln4', 'Fbln5',
                  'Emilin1', 'Magp1', 'Magp2', 'Ltbp1', 'Ltbp2', 'Ltbp3', 'Ltbp4', 'Lox',
                  'Loxl1', 'Loxl2', 'Loxl3', 'Loxl4')

# Absent: Fbln3, Fbln4, Magp1, Magp2
# todo: Make these zeros

heatmap_df <- make_heatmap_df(subset(all_cells, idents = c("Prenatal Myofibroblast", "Myofibroblast")), marker_genes, sort_clusters = FALSE)
```

    ## Warning in FetchData(object = obj, slot = "data", vars = marker_genes): The
    ## following requested variables were not found: Fbln3, Fbln4, Magp1, Magp2

``` r
heatmap_df <- heatmap_df %>% group_by(gene,time) %>% summarise(expression_mean=mean(expression))
```

    ## `summarise()` regrouping output by 'gene' (override with `.groups` argument)

``` r
addBackZeros <- function(df, markers){
  out <- do.call(rbind, lapply(markers[!(markers %in% unique(df$gene))], function(x){
    t(do.call(cbind, list(sapply(unique(levels(df$time)), function(y){
      c(x, y, 0)
    }))))
  }))
  colnames(out) <- colnames(df)
  out <- as.data.frame(out)
  out$expression_mean <- as.double(out$expression_mean)
  out$time <- ordered(out$time, levels(df$time))
  rbind(df, out)
}

n_clusters <- max(as.numeric(Idents(subset(all_cells, idents = c("Prenatal Myofibroblast", "Myofibroblast")))))
#grid.brackets()
# Look in to grid.brackets to customize the brackets. Can prob make a function for this.
myo_matrix_plot <- ggplot(addBackZeros(heatmap_df, marker_genes), aes(x = gene, y = time, fill = expression_mean)) +
        geom_tile(color = "white", size = 0.1) +
        scale_fill_distiller(palette = "Blues", direction = 1, trans = "sqrt", name = "Expression", limits = c(0,5)) +
        #ggtitle("Expression") +
        coord_fixed(ratio = 0.6) +
        theme(panel.background = element_rect(fill = "white"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title.x  = element_blank(),
              legend.text=element_text(size=14),
              legend.title=element_text(size=14),
              axis.text.x  = element_text(size=14, angle = 45, vjust = 1, hjust=1, face = "italic"),
              axis.text.y = element_text(size=14),
              axis.title.y = element_text(size=14),
              legend.key.size = unit(0.8, "cm"),
              plot.title = element_text(size=14)
        ) +
        labs(y = "Timepoint") +
        ggtitle("Expression patterns in Myofibroblasts")


myo_matrix_plot
```

![](./merge_and_clean/all_cells_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
saveTiff("./figures/meso/myo_matrix.tiff", myo_matrix_plot,
         width = 12, height = 4)
```

    ## [1] TRUE

``` r
#marker_genes <- c('Eln', 'Fbn1', 'Fbn2', 'Fbln1', 'Fbln2', 'Fbln3', 'Fbln4', 'Fbln5',
#                  'Emilin1', 'Magp1', 'Magp2', 'Ltbp1', 'Ltbp2', 'Ltbp3', 'Ltbp4', 'Lox',
#                  'Loxl1', 'Loxl2', 'Loxl3', 'Loxl4')

marker_genes <- c('Eln', 'Fbln5',
                  'Col4a1', 'Col4a2', 'Col4a3', 'Col4a4', 'Col4a5', 'Col4a6', 'Hspg2', 'Lama1', 'Lama3', 'Lama5', 'Lamb1', 'Lamb2', 'Lamb3', 'Lamc1', 'Lamc2')

# Absent: Fbln3, Fbln4, Magp1, Magp2
# todo: Make these zeros

heatmap_df <- make_heatmap_df(subset(subset(all_cells, idents = "Type I"), timepoint %in% c("E12", "E15"), invert = TRUE), marker_genes, sort_clusters = FALSE)

heatmap_df <- heatmap_df %>% group_by(gene,time) %>% summarise(expression_mean=mean(expression))
```

    ## `summarise()` regrouping output by 'gene' (override with `.groups` argument)

``` r
addBackZeros <- function(df, markers){
  out <- do.call(rbind, lapply(markers[!(markers %in% unique(df$gene))], function(x){
    t(do.call(cbind, list(sapply(unique(levels(df$time)), function(y){
      c(x, y, 0)
    }))))
  }))
  colnames(out) <- colnames(df)
  out <- as.data.frame(out)
  out$expression_mean <- as.double(out$expression_mean)
  out$time <- ordered(out$time, levels(df$time))
  rbind(df, out)
}

n_clusters <- max(as.numeric(Idents(subset(all_cells, idents = "Type I"))))
#grid.brackets()
# Look in to grid.brackets to customize the brackets. Can prob make a function for this.
#head(heatmap_df)
at1_matrix_plot <- ggplot(heatmap_df, aes(x = gene, y = time, fill = expression_mean)) +
        geom_tile(color = "white", size = 0.05) +
        scale_fill_distiller(palette = "Blues", direction = 1, trans = "sqrt", name = "Expression", limits = c(0,4)) +
        #ggtitle("Expression") +
        coord_fixed(ratio = 0.6) +
        theme(panel.background = element_rect(fill = "white"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title.x  = element_blank(),
              legend.text=element_text(size=14),
              legend.title=element_text(size=14),
              axis.text.x  = element_text(size=14, angle = 45, vjust = 1, hjust=1, face = "italic"),
              axis.text.y = element_text(size=14),
              axis.title.y = element_text(size=14),
              legend.key.size = unit(0.8, "cm"),
              plot.title = element_text(size=14)
        ) +
        labs(y = "Timepoint") +
        ggtitle("Expression patterns in AT1 cells")


at1_matrix_plot
```

![](./merge_and_clean/all_cells_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
saveTiff("./figures/epi/AT1_matrix.tiff", at1_matrix_plot,
         width = 12, height = 4)
```

    ## [1] TRUE

``` r
marker_genes <- c('Eln', 'Fbln5',
                  'Col4a1', 'Col4a2', 'Col4a3', 'Col4a4', 'Col4a5', 'Col4a6', 'Hspg2', 'Lama1', 'Lama3', 'Lama5', 'Lamb1', 'Lamb2', 'Lamb3', 'Lamc1', 'Lamc2',
                  'Vegfa', 'Pdgfa', 'Ctgf', 'F3', 'Tgfb1', 'Tgfb2')

# Absent: Fbln3, Fbln4, Magp1, Magp2
# todo: Make these zeros

heatmap_df <- make_heatmap_df(all_cells, marker_genes, sort_clusters = FALSE)

heatmap_df <- heatmap_df %>% group_by(gene,cluster) %>% summarise(expression_mean=mean(expression),
                                                                  frac_pos=(sum(expression > 0) * 100 / length(expression)))
```

    ## `summarise()` regrouping output by 'gene' (override with `.groups` argument)

``` r
n_clusters <- max(as.numeric(Idents(all_cells)))

# For brackets and text (see ./helper_functions/brackets.R)

all_mat_dot_heatmap <- ggplot(heatmap_df, aes(x = gene, y = cluster, color = expression_mean, size = frac_pos)) +
        geom_point(shape = 16) +
        scale_size(range = c(1, 7), name = "% expressing") +
        scale_color_distiller(palette = "Blues", direction = 1, trans = "sqrt", name = "Expression", limits = c(0,5)) +
        coord_fixed(ratio = 1, xlim = NULL, ylim = c(1,n_clusters), expand = TRUE, clip = "off") +
        theme(plot.margin=unit(c(1.5,1,1,1),"cm")) +
        theme(panel.background = element_rect(fill = "white", colour = "black", size = 0),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title.x  = element_blank(),
              axis.text.x = element_text(face = "italic", size = 10),
              axis.text.y = element_markdown(size = 10)
        ) +
        labs(y = "Cluster") +
        theme(axis.text.x = element_text(angle = 45, hjust=1))

all_mat_dot_heatmap
```

![](./merge_and_clean/all_cells_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

``` r
saveTiff("./figures/all/dotmap_matrix_all.tiff", all_mat_dot_heatmap,
         width = 9, height = 7)
```

    ## [1] TRUE

``` r
timepoint_factor <- factor(subset(all_cells, bulk_celltype == "Mesenchyme")$timepoint)
timepoint_ordered <- ordered(timepoint_factor, levels = unique(all_cells$timepoint))

cell_ident_df <- data.frame(timepoint = timepoint_ordered,
                            cluster = subset(all_cells, bulk_celltype == "Mesenchyme")$cell_subtype)

prop_celltypes <- round(prop.table(table(cell_ident_df$timepoint, cell_ident_df$cluster), 1) * 100, 1)


meso_melt_prop_celltypes <- melt(prop_celltypes)

#unique(meso_melt_prop_celltypes$Var2)

levels(meso_melt_prop_celltypes$Var2)[levels(meso_melt_prop_celltypes$Var2) == "Wnt2+"] <- "*Wnt2*+ FB"
levels(meso_melt_prop_celltypes$Var2)[levels(meso_melt_prop_celltypes$Var2) == "Proliferating Wnt2+"] <- "Proliferating *Wnt2*+ FB"
levels(meso_melt_prop_celltypes$Var2)[levels(meso_melt_prop_celltypes$Var2) == "Neuron"] <- "Other"
unique(meso_melt_prop_celltypes$Var2)
```

    ## [1] Adventitial fibroblast   Mesothelium              Myofibroblast           
    ## [4] Other                    Pericyte                 Prenatal Myofibroblast  
    ## [7] Proliferating *Wnt2*+ FB Smooth muscle            *Wnt2*+ FB              
    ## 9 Levels: Adventitial fibroblast Mesothelium Myofibroblast Other ... *Wnt2*+ FB

``` r
meso_melt_prop_celltypes$Var2 <- ordered(as.factor(meso_melt_prop_celltypes$Var2), c("Proliferating *Wnt2*+ FB", "*Wnt2*+ FB", "Proliferating Myofibroblast", "Myofibroblast", "Adventitial fibroblast", "Pericyte", "Mesothelium", "Smooth muscle", "Other"))

print(prop_celltypes)
```

    ##      
    ##       Adventitial fibroblast Mesothelium Myofibroblast Neuron Pericyte
    ##   E12                   21.7         3.0           0.6    1.0      0.4
    ##   E15                   11.3         3.0           4.6    0.2      4.8
    ##   E18                   11.6         0.7           3.0    0.1      2.3
    ##   P0                    18.9         3.5           3.7    0.9      5.7
    ##   P3                    28.4         5.5          36.8    1.1      7.0
    ##   P5                    12.8         2.6          36.8    0.6      9.2
    ##   P7                    25.3         3.9          40.9    1.1     11.9
    ##   P14                   21.4         7.8          37.8    1.6     17.7
    ##      
    ##       Prenatal Myofibroblast Proliferating Wnt2+ Smooth muscle Wnt2+
    ##   E12                   10.9                51.6           0.2  10.7
    ##   E15                    8.3                41.2           0.7  26.0
    ##   E18                   25.4                 2.7           0.1  54.0
    ##   P0                     3.9                 4.3           1.7  57.3
    ##   P3                     1.7                 2.6           5.0  11.9
    ##   P5                     0.5                 1.4           6.0  30.0
    ##   P7                     1.1                 1.2           2.1  12.6
    ##   P14                    0.5                 0.8           1.1  11.4

``` r
#colnames(prop_celltypes)[1:2] <- c("Prenatal *Wnt2*+", "*Wnt2*+")
meso_mountain <- ggplot(meso_melt_prop_celltypes, aes(y=value, x=Var1)) +
        geom_area(aes(color = as.factor(Var2), group=as.factor(Var2), fill = as.factor(Var2)), alpha = 1) +
        xlab("Day") +
        ylab("% of cells") +
        labs(fill = "Cluster", color = "Cluster", group = "Cluster") +
        scale_colour_manual(aesthetics = c("color", "fill"), values= color_category_20) +
        theme(legend.key = element_blank(),
              legend.text = element_markdown(size=14),
              legend.title = element_text(size=14),
              axis.text.x  = element_text(size=14),
              axis.text.y = element_text(size=14),
              axis.title.x = element_text(size=14),
              axis.title.y = element_text(size=14),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "white"),
              aspect.ratio = 0.9,
              legend.justification = "left"
        ) +
        scale_x_discrete(expand = c(.01, .01)) +
        scale_y_continuous(expand = c(.01,.01)) +
        coord_fixed(0.05)


meso_mountain
```

![](./merge_and_clean/all_cells_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

``` r
saveTiff("./figures/all/meso_mountain.tiff", meso_mountain,
         width = 9, height = 4)
```

    ## [1] TRUE

``` r
timepoint_factor <- factor(subset(all_cells, bulk_celltype == "Epithelium")$timepoint)
timepoint_ordered <- ordered(timepoint_factor, levels = unique(all_cells$timepoint))

cell_ident_df <- data.frame(timepoint = timepoint_ordered,
                            cluster = subset(all_cells, bulk_celltype == "Epithelium")$cell_subtype)

prop_celltypes <- round(prop.table(table(cell_ident_df$timepoint, cell_ident_df$cluster), 1) * 100, 1)


epi_melt_prop_celltypes <- melt(prop_celltypes)

unique(epi_melt_prop_celltypes$Var2)
```

    ## [1] Cilliated      Neuroendocrine Primordial     Secretory      Transitional  
    ## [6] Type I         Type II       
    ## 7 Levels: Cilliated Neuroendocrine Primordial Secretory ... Type II

``` r
levels(epi_melt_prop_celltypes$Var2)[levels(epi_melt_prop_celltypes$Var2) == "Type I"] <- "AT1"
levels(epi_melt_prop_celltypes$Var2)[levels(epi_melt_prop_celltypes$Var2) == "Type II"] <- "AT2"
unique(epi_melt_prop_celltypes$Var2)
```

    ## [1] Cilliated      Neuroendocrine Primordial     Secretory      Transitional  
    ## [6] AT1            AT2           
    ## 7 Levels: Cilliated Neuroendocrine Primordial Secretory Transitional ... AT2

``` r
epi_melt_prop_celltypes$Var2 <- ordered(as.factor(epi_melt_prop_celltypes$Var2), c("Primordial", "Transitional", "AT1", "AT2", "Cilliated", "Secretory", "Neuroendocrine"))

print(prop_celltypes)
```

    ##      
    ##       Cilliated Neuroendocrine Primordial Secretory Transitional Type I Type II
    ##   E12       0.0            0.0       94.1       0.0          0.0    0.7     5.2
    ##   E15       6.9            0.0       66.3       1.0          5.0    5.0    15.8
    ##   E18       2.1            2.1       17.6      11.1          3.9   40.0    23.2
    ##   P0        4.3            0.4        0.2       1.3          0.1   35.8    58.1
    ##   P3        3.1            0.3        0.2       1.7          0.1   27.2    67.3
    ##   P5        4.8            0.2        0.8       1.1          0.4   15.1    77.6
    ##   P7       10.8            0.6        0.5       4.9          0.4    7.9    74.7
    ##   P14      11.9            3.0        0.2      15.5         12.8    6.4    50.2

``` r
#colnames(prop_celltypes)[1:2] <- c("Prenatal *Wnt2*+", "*Wnt2*+")
epi_mountain <- ggplot(epi_melt_prop_celltypes, aes(y=value, x=Var1)) +
        geom_area(aes(color = as.factor(Var2), group=as.factor(Var2), fill = as.factor(Var2)), alpha = 1) +
        xlab("Day") +
        ylab("% of cells") +
        labs(fill = "Cluster", color = "Cluster", group = "Cluster") +
        scale_colour_manual(aesthetics = c("color", "fill"), values= color_category_20) +
        theme(legend.key = element_blank(),
              legend.text = element_markdown(size=14),
              legend.title = element_text(size=14),
              axis.text.x  = element_text(size=14),
              axis.text.y = element_text(size=14),
              axis.title.x = element_text(size=14),
              axis.title.y = element_text(size=14),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "white"),
              aspect.ratio = 0.9,
              legend.justification = "left"
        ) +
        scale_x_discrete(expand = c(.01, .01)) +
        scale_y_continuous(expand = c(.01,.01)) +
        coord_fixed(0.05)


epi_mountain
```

![](./merge_and_clean/all_cells_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
saveTiff("./figures/all/epi_mountain.tiff", epi_mountain,
         width = 9, height = 4)
```

    ## [1] TRUE

``` r
timepoint_factor <- factor(subset(all_cells, bulk_celltype == "Endothelium")$timepoint)
timepoint_ordered <- ordered(timepoint_factor, levels = unique(all_cells$timepoint))

cell_ident_df <- data.frame(timepoint = timepoint_ordered,
                            cluster = subset(all_cells, bulk_celltype == "Endothelium")$cell_subtype)

prop_celltypes <- round(prop.table(table(cell_ident_df$timepoint, cell_ident_df$cluster), 1) * 100, 1)


endo_melt_prop_celltypes <- melt(prop_celltypes)

unique(endo_melt_prop_celltypes$Var2)
```

    ## [1] Arterial maEC Car4+         Lymphatic     Proliferating Venous maEC  
    ## [6] miEC         
    ## Levels: Arterial maEC Car4+ Lymphatic Proliferating Venous maEC miEC

``` r
levels(endo_melt_prop_celltypes$Var2)[levels(endo_melt_prop_celltypes$Var2) == "Car4+"] <- "*Car4*+ EC"
levels(endo_melt_prop_celltypes$Var2)[levels(endo_melt_prop_celltypes$Var2) == "Proliferating"] <- "Proliferating miEC"
unique(endo_melt_prop_celltypes$Var2)
```

    ## [1] Arterial maEC      *Car4*+ EC         Lymphatic          Proliferating miEC
    ## [5] Venous maEC        miEC              
    ## 6 Levels: Arterial maEC *Car4*+ EC Lymphatic ... miEC

``` r
endo_melt_prop_celltypes$Var2 <- ordered(as.factor(endo_melt_prop_celltypes$Var2), c("Arterial maEC", "Venous maEC", "miEC", "*Car4*+ EC", "Proliferating miEC", "Lymphatic"))

print(prop_celltypes)
```

    ##      
    ##       Arterial maEC Car4+ Lymphatic Proliferating Venous maEC miEC
    ##   E12           4.8   0.0       3.3           0.5        11.0 80.4
    ##   E15           2.8   0.0       1.6           2.3         5.5 87.7
    ##   E18           6.3   5.8       2.9          12.5         5.9 66.6
    ##   P0            8.4  17.3       2.7           7.8         5.1 58.6
    ##   P3           10.6  15.9       1.8          11.5         9.4 50.8
    ##   P5            5.9  12.5       1.5          10.8         7.8 61.4
    ##   P7            7.2  11.8       3.5           7.6         7.6 62.4
    ##   P14           3.0  17.9       1.2           3.0         4.9 70.1

``` r
#colnames(prop_celltypes)[1:2] <- c("Prenatal *Wnt2*+", "*Wnt2*+")
endo_mountain <- ggplot(endo_melt_prop_celltypes, aes(y=value, x=Var1)) +
        geom_area(aes(color = as.factor(Var2), group=as.factor(Var2), fill = as.factor(Var2)), alpha = 1) +
        xlab("Day") +
        ylab("% of cells") +
        labs(fill = "Cluster", color = "Cluster", group = "Cluster") +
        scale_colour_manual(aesthetics = c("color", "fill"), values= color_category_20) +
        theme(legend.key = element_blank(),
              legend.text = element_markdown(size=14),
              legend.title = element_text(size=14),
              axis.text.x  = element_text(size=14),
              axis.text.y = element_text(size=14),
              axis.title.x = element_text(size=14),
              axis.title.y = element_text(size=14),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "white"),
              aspect.ratio = 0.9
        ) +
        scale_x_discrete(expand = c(.01, .01)) +
        scale_y_continuous(expand = c(.01,.01)) +
        coord_fixed(0.05)


endo_mountain
```

![](./merge_and_clean/all_cells_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

``` r
saveTiff("./figures/all/endo_mountain.tiff", endo_mountain,
         width = 9, height = 4)
```

    ## [1] TRUE

``` r
saveTiff("./figures/all/all_mountains.tiff",
         epi_mountain + endo_mountain + meso_mountain + plot_layout(ncol = 1, nrow = 3) & theme(legend.justification = "left", plot.margin = unit(c(5.5,5.5,35,5.5), "pt")),
         width = 8, height = 12)
```

    ## [1] TRUE
