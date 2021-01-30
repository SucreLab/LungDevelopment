scSEQ analysis of the developing mesenchyme
================
Nick Negretti
11/17/20

# Analysis of the lung mesenchyme

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

## Load data from merge and clean pipeline

``` r
meso <- readRDS("./data/meso_full_noGm42418_sct_p7b_integrated_retransform.rds")
```

## UMAP and clustering

``` r
plan("sequential")

# Make the tiempoints an ordered factor, so they are displayed consistantly
meso$timepoint <- ordered(as.factor(meso$timepoint), unique(meso$timepoint))

# cluster and UMAP from the previous processing seems fine
# k_param = 25, dims_neighbors = 1:15, cluster_res = 0.15  - is what is after the integration pipeline

p_cluster <- DimPlot(meso) + umap_theme() +
        scale_colour_manual(name = "Cluster", values = color_category_20) +
        theme(aspect.ratio=1) +
        theme(legend.text=element_text(size=14),
              legend.title=element_text(size=14))
p_time <- DimPlot(meso, group.by = "timepoint") + umap_theme() +
        scale_colour_manual(name = "Timepoint", values = color_scanpy_default)+
        theme(aspect.ratio=1) +
        theme(legend.text=element_text(size=14),
              legend.title=element_text(size=14))

p_cluster + p_time + plot_annotation("Mesenchyme")
```

![](/Users/negretn/postdoc/code/devo_scseq_github/mesenchyme/all_meso_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
saveTiff("./figures/meso/supp_umap.tiff",
         p_cluster + p_time + plot_annotation("Mesenchyme"),
         width = 12, height = 4.5)
```

    ## [1] TRUE

## Identify marker genes in each cluster

``` r
plan("multiprocess", workers = N_WORKERS)
filename <- "./data/meso_markers_clusters.rds"
if (!file.exists(filename)) {
  meso_markers <- parallelFindAllMarkers(meso)
  saveRDS(meso_markers, filename)
} else {
  meso_markers <- readRDS(filename)
}
```

### Cluster 0

``` r
meso_markers[[0 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Mfap4        0  2.098592 0.987 0.668         0
    ## Macf1        0  2.017740 0.909 0.318         0
    ## Gyg          0  1.988729 0.966 0.347         0
    ## Tcf21        0  1.864860 0.912 0.382         0
    ## Limch1       0  1.724045 0.857 0.148         0
    ## Fn1          0  1.692297 0.920 0.465         0
    ## Enpep        0  1.667648 0.772 0.097         0
    ## Fhl1         0  1.628276 0.951 0.724         0
    ## Npnt         0  1.626439 0.836 0.182         0
    ## Cebpb        0  1.550074 0.909 0.444         0
    ## Adh1         0  1.539196 0.941 0.335         0
    ## G0s2         0  1.486004 0.756 0.210         0
    ## Col13a1      0  1.303174 0.859 0.176         0
    ## Lbh          0  1.280005 0.845 0.242         0
    ## Slc27a6      0  1.239157 0.692 0.062         0
    ## Meox2        0  1.222511 0.828 0.183         0
    ## Nebl         0  1.219599 0.789 0.126         0
    ## Cpm          0  1.164033 0.479 0.036         0
    ## Slit2        0  1.161618 0.680 0.085         0
    ## Selenbp1     0  1.145685 0.855 0.364         0

### Cluster 1

``` r
meso_markers[[1 + 1]][n_print,]
```

    ##         p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Tgfbi       0  2.511384 0.981 0.242         0
    ## Aspn        0  2.107025 0.420 0.230         0
    ## Pdlim3      0  2.060527 0.918 0.115         0
    ## Mustn1      0  1.918069 0.746 0.085         0
    ## Agt         0  1.878755 0.915 0.091         0
    ## Htra1       0  1.833499 0.838 0.204         0
    ## Tagln       0  1.805600 0.919 0.346         0
    ## Egfem1      0  1.624733 0.842 0.091         0
    ## Spon2       0  1.570425 0.883 0.293         0
    ## Prss35      0  1.495087 0.706 0.164         0
    ## Tnc         0  1.487866 0.853 0.201         0
    ## P2ry14      0  1.479401 0.834 0.120         0
    ## Ckb         0  1.435215 0.941 0.396         0
    ## Loxl2       0  1.396974 0.855 0.422         0
    ## Net1        0  1.357518 0.610 0.152         0
    ## Actg2       0  1.341935 0.694 0.095         0
    ## Filip1l     0  1.242676 0.737 0.140         0
    ## Des         0  1.205209 0.758 0.212         0
    ## Wnt5a       0  1.204440 0.753 0.214         0
    ## Acta2       0  1.167611 0.937 0.611         0

### Cluster 2

``` r
meso_markers[[2 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Dcn          0  2.547963 0.645 0.126         0
    ## Lum          0  2.054821 0.461 0.046         0
    ## Col1a1       0  1.972118 0.963 0.773         0
    ## Serpinf1     0  1.862776 0.675 0.056         0
    ## Col1a2       0  1.766497 0.966 0.843         0
    ## Mfap5        0  1.745047 0.673 0.155         0
    ## Ptn          0  1.642166 0.554 0.119         0
    ## Dlk1         0  1.628111 0.551 0.220         0
    ## Col3a1       0  1.609137 0.941 0.799         0
    ## Gsn          0  1.497312 0.907 0.550         0
    ## Igfbp7       0  1.462536 0.860 0.454         0
    ## Meg3         0  1.300578 0.848 0.425         0
    ## Igfbp4       0  1.285860 0.970 0.741         0
    ## Col14a1      0  1.276534 0.707 0.246         0
    ## S100a6       0  1.268439 0.741 0.197         0
    ## C1qtnf3      0  1.240650 0.277 0.013         0
    ## Cpz          0  1.220920 0.501 0.046         0
    ## Cthrc1       0  1.215085 0.378 0.025         0
    ## Clec3b       0  1.207051 0.673 0.308         0
    ## Pcolce       0  1.182981 0.896 0.614         0

### Cluster 3

``` r
meso_markers[[3 + 1]][n_print,]
```

    ##           p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Hist1h2ap     0 1.4568609 0.496 0.076         0
    ## H2afz         0 1.4370052 0.965 0.657         0
    ## Hmgb2         0 1.2528511 0.781 0.329         0
    ## Pclaf         0 1.1600875 0.611 0.084         0
    ## Mif           0 1.1509764 0.934 0.593         0
    ## Ran           0 1.1052982 0.965 0.619         0
    ## Ube2c         0 1.0575719 0.417 0.055         0
    ## Tuba1b        0 1.0480118 0.981 0.818         0
    ## H2afx         0 1.0350882 0.605 0.153         0
    ## Hist1h2ae     0 1.0059444 0.453 0.071         0
    ## Hmgn2         0 0.9857195 0.993 0.746         0
    ## Stmn1         0 0.9496317 0.778 0.205         0
    ## Gapdh         0 0.9410165 0.981 0.840         0
    ## Birc5         0 0.9225034 0.567 0.071         0
    ## Npm1          0 0.9172386 0.984 0.891         0
    ## Selenoh       0 0.9029108 0.835 0.322         0
    ## Anp32b        0 0.8861377 0.942 0.544         0
    ## Hmgb1         0 0.8784957 0.998 0.941         0
    ## Hist1h1b      0 0.8704823 0.387 0.051         0
    ## Cenpa         0 0.8635937 0.453 0.071         0

### Cluster 4

``` r
meso_markers[[4 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Cox4i2       0  3.162804 0.977 0.039         0
    ## Gucy1a1      0  2.861935 0.989 0.287         0
    ## Postn        0  2.643348 0.990 0.211         0
    ## Gucy1b1      0  2.619673 0.979 0.214         0
    ## Itga1        0  2.349548 0.958 0.138         0
    ## Higd1b       0  2.320877 0.889 0.006         0
    ## Ndufa4l2     0  2.273249 0.953 0.114         0
    ## Pdgfrb       0  2.162783 0.958 0.160         0
    ## Mfge8        0  2.162532 0.997 0.687         0
    ## Mcam         0  2.151101 0.941 0.051         0
    ## Emid1        0  2.137377 0.922 0.111         0
    ## Pde5a        0  2.033240 0.961 0.398         0
    ## Pdzd2        0  2.022086 0.843 0.025         0
    ## Col4a1       0  2.013812 0.986 0.515         0
    ## Col4a2       0  1.952886 0.960 0.370         0
    ## Pcdh18       0  1.904550 0.912 0.109         0
    ## Ebf1         0  1.845813 0.926 0.088         0
    ## Trpc6        0  1.823760 0.853 0.017         0
    ## Itm2a        0  1.798879 0.963 0.381         0
    ## Fermt2       0  1.646052 0.981 0.637         0

### Cluster 5

``` r
meso_markers[[5 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Ednrb        0  2.052480 0.901 0.124         0
    ## Enpp2        0  1.946599 0.617 0.208         0
    ## Crh          0  1.846816 0.427 0.003         0
    ## Nnat         0  1.721900 0.854 0.306         0
    ## Hhip         0  1.611922 0.816 0.265         0
    ## H19          0  1.536411 0.787 0.223         0
    ## Fam105a      0  1.512982 0.627 0.125         0
    ## Igf1         0  1.439690 0.894 0.342         0
    ## Crispld2     0  1.396503 0.699 0.110         0
    ## Actg2        0  1.266296 0.622 0.203         0
    ## Pdgfra       0  1.256950 0.723 0.302         0
    ## Rgs2         0  1.251188 0.662 0.190         0
    ## Myh11        0  1.214048 0.791 0.237         0
    ## Txnip        0  1.191003 0.634 0.329         0
    ## Tsc22d3      0  1.154017 0.654 0.323         0
    ## Penk         0  1.131542 0.436 0.140         0
    ## Aldh2        0  1.079622 0.874 0.825         0
    ## Gm13889      0  1.075241 0.669 0.207         0
    ## Mylk         0  1.059601 0.925 0.657         0
    ## Mettl7a1     0  1.015284 0.542 0.157         0

### Cluster 6

``` r
meso_markers[[6 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Upk3b        0  3.286119 0.992 0.003         0
    ## Msln         0  3.286119 0.706 0.001         0
    ## C3           0  3.146365 0.758 0.008         0
    ## Igfbp6       0  2.631859 0.933 0.072         0
    ## Bcam         0  2.522937 0.867 0.066         0
    ## Igfbp5       0  2.430135 0.935 0.157         0
    ## Clu          0  2.419644 0.913 0.046         0
    ## Aldh1a2      0  2.229059 0.955 0.037         0
    ## Krt19        0  2.172003 0.916 0.001         0
    ## Crip1        0  2.140877 0.920 0.397         0
    ## Sfrp2        0  2.139940 0.803 0.146         0
    ## Ezr          0  2.064742 0.924 0.050         0
    ## Slurp1       0  2.063462 0.297 0.000         0
    ## Cav1         0  2.057731 0.949 0.197         0
    ## Slc9a3r1     0  2.033189 0.948 0.022         0
    ## Upk1b        0  2.001888 0.940 0.094         0
    ## Krt7         0  1.987846 0.803 0.002         0
    ## Hspb1        0  1.961283 0.962 0.299         0
    ## Gbp2         0  1.900108 0.598 0.030         0
    ## S100a6       0  1.892695 0.848 0.276         0

### Cluster 7

``` r
meso_markers[[7 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Acta2        0  1.858180 0.999 0.676         0
    ## Eln          0  1.826326 0.996 0.695         0
    ## Tagln        0  1.686080 0.992 0.461         0
    ## Map3k7cl     0  1.615958 0.802 0.013         0
    ## Actc1        0  1.586931 0.478 0.021         0
    ## Crip1        0  1.582157 0.935 0.401         0
    ## Cxcl12       0  1.563795 0.707 0.205         0
    ## Myl9         0  1.524350 0.993 0.618         0
    ## Cav1         0  1.297555 0.833 0.207         0
    ## Gpc6         0  1.273014 0.715 0.116         0
    ## Wisp2        0  1.272134 0.594 0.013         0
    ## Tm4sf1       0  1.270081 0.674 0.111         0
    ## Itga4        0  1.216852 0.617 0.064         0
    ## Csrp2        0  1.214317 0.756 0.433         0
    ## Myh11        0  1.210371 0.808 0.258         0
    ## Parm1        0  1.203789 0.615 0.093         0
    ## Tpm2         0  1.133883 0.936 0.601         0
    ## Gadd45g      0  1.113346 0.587 0.153         0
    ## Aoc3         0  1.109874 0.718 0.057         0
    ## Igfbp7       0  1.086200 0.970 0.517         0

### Cluster 8

``` r
meso_markers[[8 + 1]][n_print,]
```

    ##                 p_val avg_logFC pct.1 pct.2     p_val_adj
    ## Mpz      0.000000e+00 2.8172901 0.454 0.001  0.000000e+00
    ## Fabp3    0.000000e+00 2.2536888 0.290 0.043  0.000000e+00
    ## Dbi      0.000000e+00 2.2208684 0.925 0.597  0.000000e+00
    ## Plp1     0.000000e+00 2.1718686 0.650 0.007  0.000000e+00
    ## Cnp      0.000000e+00 1.6614826 0.601 0.053  0.000000e+00
    ## Sox10    0.000000e+00 1.6116184 0.618 0.002  0.000000e+00
    ## Mal      0.000000e+00 1.5885925 0.486 0.001  0.000000e+00
    ## Mbp      0.000000e+00 1.5641929 0.461 0.012  0.000000e+00
    ## Cryab    0.000000e+00 1.5382939 0.857 0.448  0.000000e+00
    ## Arpc1b   0.000000e+00 1.4591853 0.744 0.603  0.000000e+00
    ## Egfl8    0.000000e+00 1.4242010 0.464 0.003  0.000000e+00
    ## Plekhb1  0.000000e+00 1.1542162 0.531 0.003  0.000000e+00
    ## Sema3b   0.000000e+00 1.1317556 0.493 0.023  0.000000e+00
    ## Ttyh1    0.000000e+00 1.1200421 0.553 0.007  0.000000e+00
    ## Cmtm5    0.000000e+00 1.1123070 0.507 0.004  0.000000e+00
    ## Fbxo7    0.000000e+00 0.9079919 0.517 0.085  0.000000e+00
    ## Metrn    0.000000e+00 0.8812798 0.507 0.097  0.000000e+00
    ## Gjc3     0.000000e+00 0.8245382 0.420 0.001  0.000000e+00
    ## L1cam    0.000000e+00 0.7898473 0.425 0.001  0.000000e+00
    ## Ldhb    3.550919e-318 1.0847139 0.710 0.168 7.739583e-314

## Use known marker genes to determine the identity of the clusters

``` r
marker_genes <- c("Wnt2", "Macf1", # Wnt2 fibroblasts 0
                  "Tgfbi", "Wnt5a", # Myofibroblasts 1
                  "Dcn", "Col1a1", # Adventitial fibroblast 2
                  "Top2a", "Mki67", # Prenatal Wnt2 3
                  "Cspg4", "Cox4i2", # Pericyte 4
                  "Ednrb", "Pdgfra", # Prenatal myo 5
                  "Wt1", "Upk3b", # Mesothelium 6
                  "Eln", "Acta2", # Smooth muscle 7
                  "Mpz", "Mal", # Neuron 8
                  "Tnnt2", "Actc1" # Cardiomyocyte
)


heatmap_df <- make_heatmap_df(meso, marker_genes, sort_clusters = FALSE)

heatmap_df <- heatmap_df %>% group_by(gene,cluster) %>% summarise(expression_mean=mean(expression))
```

    ## `summarise()` regrouping output by 'gene' (override with `.groups` argument)

``` r
# For brackets and text (see ./helper_functions/brackets.R)
n_clusters <- max(as.numeric(Idents(meso)))
text_vec <- list(c("Wnt2+ FB", 2),
                 c("Myofibroblast", 2),
                 c("Adventitial FB", 2),
                 c("Prenatal\nWnt2+ FB", 2),
                 c("Pericyte", 2),
                 c("Prenatal\nMyofibroblast", 2),
                 c("Mesothelium", 2),
                 c("Smooth\nmuscle", 2),
                 c("Neuron", 2),
                 c("Cardiomyocyte", 2)
)

meso_heatmap <- ggplot(heatmap_df, aes(x = gene, y = cluster, fill = expression_mean)) +
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
              axis.text.x = element_text(size = 12, face = "italic"),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 14),
              legend.title = element_text(size = 14),
              legend.text = element_text(size = 12)
        ) +
        labs(y = "Cluster") +
        theme(axis.text.x = element_text(angle = 45, hjust=1))

meso_heatmap
```

![](/Users/negretn/postdoc/code/devo_scseq_github/mesenchyme/all_meso_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
saveTiff("./figures/meso/supp_cluster_markers.tiff",
         meso_heatmap,
         width = 12, height = 5.5)
```

    ## [1] TRUE

## Fraction of the cells over time

``` r
timepoint_factor <- factor(meso$timepoint)
timepoint_ordered <- ordered(timepoint_factor, levels = unique(meso$timepoint))

cell_ident_df <- data.frame(timepoint = timepoint_ordered,
                            cluster = Idents(meso))

prop_celltypes <- round(prop.table(table(cell_ident_df$timepoint, cell_ident_df$cluster), 1) * 100, 1)


ggplot(melt(prop_celltypes), aes(y=value, x=Var1)) +
        geom_area(aes(color = as.factor(Var2), group=as.factor(Var2), fill = as.factor(Var2))) +
        xlab("Day") +
        ylab("% of cells") +
        labs(fill = "Cluster", color = "Cluster", group = "Cluster") +
        scale_colour_manual(aesthetics = c("color", "fill"), values=color_category_20)
```

![](/Users/negretn/postdoc/code/devo_scseq_github/mesenchyme/all_meso_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

## Apply labels

``` r
meso_relabel <- meso
# Label the cell-types
text_vec <- list(c("Wnt2+ FB", 2),
                 c("Myofibroblast", 2),
                 c("Adventitial", 2),
                 c("Prolif.\nWnt2+ FB", 2),
                 c("Pericyte", 2),
                 c("Prolif.\nMyofibroblast", 2),
                 c("Mesothelium", 2),
                 c("Smooth\nmuscle", 2),
                 c("Neuron", 2))
levels(Idents(meso_relabel)) <- c("Wnt2+ FB", #0
                                  "Myofibroblast", #1
                                  "Adventitial fibroblast", #2
                                  "Proliferating Wnt2+ FB", #3
                                  "Pericyte", #4
                                  "Proliferating myofibroblast", #5
                                  "Mesothelium", #6
                                  "Smooth muscle", #7
                                  "Neuron" #8

)
Idents(meso_relabel) <- ordered(Idents(meso_relabel), c("Proliferating Wnt2+ FB", "Wnt2+ FB", "Proliferating myofibroblast", "Myofibroblast", "Adventitial fibroblast", "Pericyte", "Mesothelium", "Smooth muscle", "Neuron"))

Ident_md <- Idents(meso_relabel)
levels(Ident_md)[1:2] <- c("Prenatal *Wnt2*+ FB", "*Wnt2*+ FB")

meso_relabel$Ident_md <- Ident_md

meso_relabel$bulk_celltype <- "Mesenchyme"

# Save this for use by other stages of the analysis
filename <- "./data/meso_full_celltype_labels.rds"
if (!file.exists(filename)) {
  saveRDS(meso_relabel, filename)
}
```

# Relabel clusters based on marker gene expression

``` r
p_rename_cluster <- DimPlot(meso_relabel, group.by = "Ident_md") + umap_theme() +
        scale_colour_manual(name = "Cluster", values = color_category_20) +
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14),
              plot.title = element_text(size=14)) +
        ggtitle("Mesenchyme")


p_time <- DimPlot(meso_relabel, group.by = "timepoint") + umap_theme() +
        scale_colour_manual(name = "Timepoint", values = color_category_20) +
        theme(aspect.ratio=1) +
        theme(legend.text=element_text(size=14),
              legend.title=element_text(size=14))

p_time + p_rename_cluster + plot_annotation("Mesenchyme")
```

![](/Users/negretn/postdoc/code/devo_scseq_github/mesenchyme/all_meso_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
saveTiff("./figures/meso/umap.tiff",
         p_rename_cluster,
         width = 6, height = 4.5)
```

    ## [1] TRUE

## Re-confirm the markers in the relabeled clusters

``` r
plan("multiprocess", workers = N_WORKERS)
filename <- "./data/meso_markers_clusters_relabeled.rds"
if (!file.exists(filename)) {
  meso_relabeled_markers <- parallelFindAllMarkers(meso_relabel)
  saveRDS(meso_relabeled_markers, filename)
} else {
  meso_relabeled_markers <- readRDS(filename)
}
```

### Prenatal Wnt2+

``` r
meso_relabeled_markers[[0 + 1]][n_print,]
```

    ##           p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Hist1h2ap     0 1.4568609 0.496 0.076         0
    ## H2afz         0 1.4370052 0.965 0.657         0
    ## Hmgb2         0 1.2528511 0.781 0.329         0
    ## Pclaf         0 1.1600875 0.611 0.084         0
    ## Mif           0 1.1509764 0.934 0.593         0
    ## Ran           0 1.1052982 0.965 0.619         0
    ## Ube2c         0 1.0575719 0.417 0.055         0
    ## Tuba1b        0 1.0480118 0.981 0.818         0
    ## H2afx         0 1.0350882 0.605 0.153         0
    ## Hist1h2ae     0 1.0059444 0.453 0.071         0
    ## Hmgn2         0 0.9857195 0.993 0.746         0
    ## Stmn1         0 0.9496317 0.778 0.205         0
    ## Gapdh         0 0.9410165 0.981 0.840         0
    ## Birc5         0 0.9225034 0.567 0.071         0
    ## Npm1          0 0.9172386 0.984 0.891         0
    ## Selenoh       0 0.9029108 0.835 0.322         0
    ## Anp32b        0 0.8861377 0.942 0.544         0
    ## Hmgb1         0 0.8784957 0.998 0.941         0
    ## Hist1h1b      0 0.8704823 0.387 0.051         0
    ## Cenpa         0 0.8635937 0.453 0.071         0

### Wnt2+

``` r
meso_relabeled_markers[[1 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Mfap4        0  2.098592 0.987 0.668         0
    ## Macf1        0  2.017740 0.909 0.318         0
    ## Gyg          0  1.988729 0.966 0.347         0
    ## Tcf21        0  1.864860 0.912 0.382         0
    ## Limch1       0  1.724045 0.857 0.148         0
    ## Fn1          0  1.692297 0.920 0.465         0
    ## Enpep        0  1.667648 0.772 0.097         0
    ## Fhl1         0  1.628276 0.951 0.724         0
    ## Npnt         0  1.626439 0.836 0.182         0
    ## Cebpb        0  1.550074 0.909 0.444         0
    ## Adh1         0  1.539196 0.941 0.335         0
    ## G0s2         0  1.486004 0.756 0.210         0
    ## Col13a1      0  1.303174 0.859 0.176         0
    ## Lbh          0  1.280005 0.845 0.242         0
    ## Slc27a6      0  1.239157 0.692 0.062         0
    ## Meox2        0  1.222511 0.828 0.183         0
    ## Nebl         0  1.219599 0.789 0.126         0
    ## Cpm          0  1.164033 0.479 0.036         0
    ## Slit2        0  1.161618 0.680 0.085         0
    ## Selenbp1     0  1.145685 0.855 0.364         0

### Prenatal Myofibroblast

``` r
meso_relabeled_markers[[2 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Ednrb        0  2.052480 0.901 0.124         0
    ## Enpp2        0  1.946599 0.617 0.208         0
    ## Crh          0  1.846816 0.427 0.003         0
    ## Nnat         0  1.721900 0.854 0.306         0
    ## Hhip         0  1.611922 0.816 0.265         0
    ## H19          0  1.536411 0.787 0.223         0
    ## Fam105a      0  1.512982 0.627 0.125         0
    ## Igf1         0  1.439690 0.894 0.342         0
    ## Crispld2     0  1.396503 0.699 0.110         0
    ## Actg2        0  1.266296 0.622 0.203         0
    ## Pdgfra       0  1.256950 0.723 0.302         0
    ## Rgs2         0  1.251188 0.662 0.190         0
    ## Myh11        0  1.214048 0.791 0.237         0
    ## Txnip        0  1.191003 0.634 0.329         0
    ## Tsc22d3      0  1.154017 0.654 0.323         0
    ## Penk         0  1.131542 0.436 0.140         0
    ## Aldh2        0  1.079622 0.874 0.825         0
    ## Gm13889      0  1.075241 0.669 0.207         0
    ## Mylk         0  1.059601 0.925 0.657         0
    ## Mettl7a1     0  1.015284 0.542 0.157         0

### Myofibroblast

``` r
meso_relabeled_markers[[3 + 1]][n_print,]
```

    ##         p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Tgfbi       0  2.511384 0.981 0.242         0
    ## Aspn        0  2.107025 0.420 0.230         0
    ## Pdlim3      0  2.060527 0.918 0.115         0
    ## Mustn1      0  1.918069 0.746 0.085         0
    ## Agt         0  1.878755 0.915 0.091         0
    ## Htra1       0  1.833499 0.838 0.204         0
    ## Tagln       0  1.805600 0.919 0.346         0
    ## Egfem1      0  1.624733 0.842 0.091         0
    ## Spon2       0  1.570425 0.883 0.293         0
    ## Prss35      0  1.495087 0.706 0.164         0
    ## Tnc         0  1.487866 0.853 0.201         0
    ## P2ry14      0  1.479401 0.834 0.120         0
    ## Ckb         0  1.435215 0.941 0.396         0
    ## Loxl2       0  1.396974 0.855 0.422         0
    ## Net1        0  1.357518 0.610 0.152         0
    ## Actg2       0  1.341935 0.694 0.095         0
    ## Filip1l     0  1.242676 0.737 0.140         0
    ## Des         0  1.205209 0.758 0.212         0
    ## Wnt5a       0  1.204440 0.753 0.214         0
    ## Acta2       0  1.167611 0.937 0.611         0

### Adventitial fibroblast

``` r
meso_relabeled_markers[[4 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Dcn          0  2.547963 0.645 0.126         0
    ## Lum          0  2.054821 0.461 0.046         0
    ## Col1a1       0  1.972118 0.963 0.773         0
    ## Serpinf1     0  1.862776 0.675 0.056         0
    ## Col1a2       0  1.766497 0.966 0.843         0
    ## Mfap5        0  1.745047 0.673 0.155         0
    ## Ptn          0  1.642166 0.554 0.119         0
    ## Dlk1         0  1.628111 0.551 0.220         0
    ## Col3a1       0  1.609137 0.941 0.799         0
    ## Gsn          0  1.497312 0.907 0.550         0
    ## Igfbp7       0  1.462536 0.860 0.454         0
    ## Meg3         0  1.300578 0.848 0.425         0
    ## Igfbp4       0  1.285860 0.970 0.741         0
    ## Col14a1      0  1.276534 0.707 0.246         0
    ## S100a6       0  1.268439 0.741 0.197         0
    ## C1qtnf3      0  1.240650 0.277 0.013         0
    ## Cpz          0  1.220920 0.501 0.046         0
    ## Cthrc1       0  1.215085 0.378 0.025         0
    ## Clec3b       0  1.207051 0.673 0.308         0
    ## Pcolce       0  1.182981 0.896 0.614         0

### Pericyte

``` r
meso_relabeled_markers[[5 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Cox4i2       0  3.162804 0.977 0.039         0
    ## Gucy1a1      0  2.861935 0.989 0.287         0
    ## Postn        0  2.643348 0.990 0.211         0
    ## Gucy1b1      0  2.619673 0.979 0.214         0
    ## Itga1        0  2.349548 0.958 0.138         0
    ## Higd1b       0  2.320877 0.889 0.006         0
    ## Ndufa4l2     0  2.273249 0.953 0.114         0
    ## Pdgfrb       0  2.162783 0.958 0.160         0
    ## Mfge8        0  2.162532 0.997 0.687         0
    ## Mcam         0  2.151101 0.941 0.051         0
    ## Emid1        0  2.137377 0.922 0.111         0
    ## Pde5a        0  2.033240 0.961 0.398         0
    ## Pdzd2        0  2.022086 0.843 0.025         0
    ## Col4a1       0  2.013812 0.986 0.515         0
    ## Col4a2       0  1.952886 0.960 0.370         0
    ## Pcdh18       0  1.904550 0.912 0.109         0
    ## Ebf1         0  1.845813 0.926 0.088         0
    ## Trpc6        0  1.823760 0.853 0.017         0
    ## Itm2a        0  1.798879 0.963 0.381         0
    ## Fermt2       0  1.646052 0.981 0.637         0

### Mesothelium

``` r
meso_relabeled_markers[[6 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Upk3b        0  3.286119 0.992 0.003         0
    ## Msln         0  3.286119 0.706 0.001         0
    ## C3           0  3.146365 0.758 0.008         0
    ## Igfbp6       0  2.631859 0.933 0.072         0
    ## Bcam         0  2.522937 0.867 0.066         0
    ## Igfbp5       0  2.430135 0.935 0.157         0
    ## Clu          0  2.419644 0.913 0.046         0
    ## Aldh1a2      0  2.229059 0.955 0.037         0
    ## Krt19        0  2.172003 0.916 0.001         0
    ## Crip1        0  2.140877 0.920 0.397         0
    ## Sfrp2        0  2.139940 0.803 0.146         0
    ## Ezr          0  2.064742 0.924 0.050         0
    ## Slurp1       0  2.063462 0.297 0.000         0
    ## Cav1         0  2.057731 0.949 0.197         0
    ## Slc9a3r1     0  2.033189 0.948 0.022         0
    ## Upk1b        0  2.001888 0.940 0.094         0
    ## Krt7         0  1.987846 0.803 0.002         0
    ## Hspb1        0  1.961283 0.962 0.299         0
    ## Gbp2         0  1.900108 0.598 0.030         0
    ## S100a6       0  1.892695 0.848 0.276         0

### Neuron

``` r
meso_relabeled_markers[[7 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Acta2        0  1.858180 0.999 0.676         0
    ## Eln          0  1.826326 0.996 0.695         0
    ## Tagln        0  1.686080 0.992 0.461         0
    ## Map3k7cl     0  1.615958 0.802 0.013         0
    ## Actc1        0  1.586931 0.478 0.021         0
    ## Crip1        0  1.582157 0.935 0.401         0
    ## Cxcl12       0  1.563795 0.707 0.205         0
    ## Myl9         0  1.524350 0.993 0.618         0
    ## Cav1         0  1.297555 0.833 0.207         0
    ## Gpc6         0  1.273014 0.715 0.116         0
    ## Wisp2        0  1.272134 0.594 0.013         0
    ## Tm4sf1       0  1.270081 0.674 0.111         0
    ## Itga4        0  1.216852 0.617 0.064         0
    ## Csrp2        0  1.214317 0.756 0.433         0
    ## Myh11        0  1.210371 0.808 0.258         0
    ## Parm1        0  1.203789 0.615 0.093         0
    ## Tpm2         0  1.133883 0.936 0.601         0
    ## Gadd45g      0  1.113346 0.587 0.153         0
    ## Aoc3         0  1.109874 0.718 0.057         0
    ## Igfbp7       0  1.086200 0.970 0.517         0

``` r
timepoint_factor <- factor(meso_relabel$timepoint)
timepoint_ordered <- ordered(timepoint_factor, levels = unique(meso_relabel$timepoint))

cell_ident_df <- data.frame(timepoint = timepoint_ordered,
                            cluster = Idents(meso_relabel))

prop_celltypes <- round(prop.table(table(cell_ident_df$timepoint, cell_ident_df$cluster), 1) * 100, 1)

print(prop_celltypes)
```

    ##      
    ##       Proliferating Wnt2+ FB Wnt2+ FB Proliferating myofibroblast Myofibroblast
    ##   E12                   51.6     10.7                        10.9           0.6
    ##   E15                   41.2     26.0                         8.3           4.6
    ##   E18                    2.8     54.0                        25.4           3.0
    ##   P0                     4.3     57.3                         3.9           3.7
    ##   P3                     2.6     11.9                         1.7          36.8
    ##   P5                     1.4     30.0                         0.5          36.8
    ##   P7                     1.2     12.6                         1.1          40.8
    ##   P14                    0.8     11.4                         0.5          37.8
    ##      
    ##       Adventitial fibroblast Pericyte Mesothelium Smooth muscle Neuron
    ##   E12                   21.7      0.4         3.0           0.2    1.0
    ##   E15                   11.3      4.8         3.0           0.7    0.2
    ##   E18                   11.6      2.3         0.7           0.1    0.1
    ##   P0                    18.9      5.7         3.5           1.7    0.9
    ##   P3                    28.4      7.0         5.5           5.0    1.1
    ##   P5                    12.8      9.2         2.6           6.0    0.6
    ##   P7                    25.2     11.9         3.9           2.1    1.1
    ##   P14                   21.4     17.7         7.8           1.1    1.6

``` r
colnames(prop_celltypes)[1:2] <- c("Prenatal *Wnt2*+ FB", "*Wnt2*+ FB")
meso_mountain <- ggplot(melt(prop_celltypes), aes(y=value, x=Var1)) +
        geom_area(aes(color = as.factor(Var2), group=as.factor(Var2), fill = as.factor(Var2))) +
        xlab("Day") +
        ylab("% of cells") +
        labs(fill = "Cluster", color = "Cluster", group = "Cluster") +
        scale_colour_manual(aesthetics = c("color", "fill"), values=color_category_20) +
        theme(legend.key = element_blank(),
              legend.text = element_markdown(size=14),
              legend.title = element_text(size=14),
              axis.text.x  = element_text(size=14),
              axis.text.y = element_text(size=14),
              axis.title.x = element_text(size=14),
              axis.title.y = element_text(size=14),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "white")
        ) +
        scale_x_discrete(expand = c(.01, .01)) +
        scale_y_continuous(expand = c(.01,.01)) +
        coord_fixed(0.05)


meso_mountain
```

![](/Users/negretn/postdoc/code/devo_scseq_github/mesenchyme/all_meso_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

``` r
saveTiff("./figures/meso/labeled_mountain.tiff",
         meso_mountain,
         width = 9, height = 4)
```

    ## [1] TRUE

## Wnt2 and Wnt5a expression

``` r
wnt2_meso_plot <- FeaturePlot(meso, "Wnt2") +
        featureplot_theme() +
        theme(aspect.ratio=1) +
        scale_color_viridis(name = "Expression", direction = -1)
```

    ## Scale for 'colour' is already present. Adding another scale for 'colour',
    ## which will replace the existing scale.

``` r
wnt5a_meso_plot <- FeaturePlot(meso, "Wnt5a") +
        featureplot_theme() +
        theme(aspect.ratio=1) +
        scale_color_viridis(name = "Expression", direction = -1)
```

    ## Scale for 'colour' is already present. Adding another scale for 'colour',
    ## which will replace the existing scale.

``` r
wnt2_meso_plot + wnt5a_meso_plot + plot_layout(ncol = 2)
```

![](/Users/negretn/postdoc/code/devo_scseq_github/mesenchyme/all_meso_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
saveTiff("./figures/meso/meso_wnt_plot.tiff",
         wnt2_meso_plot + wnt5a_meso_plot + plot_layout(ncol = 2),
         width = 12, height = 4.5)
```

    ## [1] TRUE

# Higher res clustering

## UMAP and clustering with increased granularity

``` r
plan("sequential")

# cluster and UMAP from the previous processing seems fine
# k_param = 25, dims_neighbors = 1:15, cluster_res = 0.15  - is what is after the integration pipeline
meso_highres <- cluster_pca_umap(meso, k_param = 25, dims_umap = 1:15, dims_neighbors = 1:25, cluster_res = 0.2) # Note, this breaks if future is set to plan: multiprocess
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

    ## Computing nearest neighbor graph

    ## Computing SNN

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 58562
    ## Number of edges: 2082165
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9515
    ## Number of communities: 13
    ## Elapsed time: 15 seconds

``` r
p_cluster <- DimPlot(meso_highres) + umap_theme() +
        scale_colour_manual(name = "Cluster", values = color_category_20) +
        theme(aspect.ratio=1)
p_time <- DimPlot(meso_highres, group.by = "timepoint") + umap_theme() +
        scale_colour_manual(name = "Timepoint", values = color_scanpy_default)+
        theme(aspect.ratio=1)

p_cluster + p_time + plot_annotation("Mesenchyme")
```

![](/Users/negretn/postdoc/code/devo_scseq_github/mesenchyme/all_meso_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

``` r
#pdf("./figures/meso/supp_umap_highres.pdf", width = 12, height = 4.5)
#p_cluster + p_time + plot_annotation("Mesenchyme")
#dev.off()
```

``` r
ggplot_data <- as.data.frame(Embeddings(meso_highres, reduction = "umap"))
ggplot_data$named_clusters <- Idents(meso_relabel)
ggplot_data$timepoint <- meso_relabel$timepoint
ggplot_data$subclusters <- Idents(meso_highres)

myo_highlight <- ggplot() +
        geom_point(data = ggplot_data[!(ggplot_data$named_clusters %in% c("Proliferating myofibroblast", "Myofibroblast")),],
                   aes(x = UMAP_1, y = UMAP_2), size = 0.1, alpha = 0.2, color = "gray") +
        geom_point(data = ggplot_data[(ggplot_data$named_clusters %in% c("Proliferating myofibroblast", "Myofibroblast")),],
                   aes(x = UMAP_1, y = UMAP_2, color = timepoint), size = 0.1) +
        umap_theme() +
        scale_colour_manual(name = "Timepoint", values = color_scanpy_default) +
        theme(aspect.ratio=1) +
        theme(legend.text=element_text(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3))) +
        ggtitle("Myofibroblasts by time")

myo_highlight
```

![](/Users/negretn/postdoc/code/devo_scseq_github/mesenchyme/all_meso_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

``` r
myo_supp_highlight <- ggplot() +
        geom_point(data = ggplot_data[!(ggplot_data$subclusters %in% c(0,6,8,10)),],
                   aes(x = UMAP_1, y = UMAP_2), size = 0.1, alpha = 0.2, color = "gray") +
        geom_point(data = ggplot_data[(ggplot_data$subclusters %in% c(0,6,8,10)),],
                   aes(x = UMAP_1, y = UMAP_2, color = subclusters), size = 0.1) +
        umap_theme() +
        scale_colour_manual(name = "Cluster", values = color_scanpy_default, labels = c("E12/E15 Myofibroblast",
                                                                                        "E18 Myofibroblast",
                                                                                        "Replicating myofibroblast",
                                                                                        "Mature myofibroblast"),
        breaks = c(10,6,8,0)) +
        theme(aspect.ratio=1) +
        theme(legend.text=element_text(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3))) +
        ggtitle("Myofibroblast subclusters")

myo_supp_highlight
```

![](/Users/negretn/postdoc/code/devo_scseq_github/mesenchyme/all_meso_files/figure-gfm/unnamed-chunk-30-2.png)<!-- -->

``` r
saveTiff("./figures/meso/umap_myo_timepoint.tiff",
         myo_highlight,
         width = 5, height = 4.5)
```

    ## [1] TRUE

``` r
saveTiff("./figures/meso/umap_myo_supp_cluster.tiff",
         myo_supp_highlight,
         width = 6, height = 4.5)
```

    ## [1] TRUE

## Identify marker genes in each cluster

``` r
plan("multiprocess", workers = N_WORKERS)
filename <- "./data/meso_markers_clusters_highres.rds"
if (!file.exists(filename)) {
  meso_markers_highres <- parallelFindAllMarkers(meso_highres)
  saveRDS(meso_markers_highres, filename)
} else {
  meso_markers_highres <- readRDS(filename)
}
```

### Cluster 0

``` r
meso_markers_highres[[0 + 1]][n_print,]
```

    ##         p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Tgfbi       0  2.376759 0.995 0.265         0
    ## Aspn        0  2.265244 0.482 0.222         0
    ## Mustn1      0  2.010895 0.799 0.096         0
    ## Pdlim3      0  1.970410 0.935 0.140         0
    ## Tagln       0  1.848136 0.946 0.359         0
    ## Agt         0  1.785046 0.923 0.118         0
    ## Spon2       0  1.586474 0.899 0.310         0
    ## Egfem1      0  1.493564 0.839 0.119         0
    ## Net1        0  1.449415 0.660 0.156         0
    ## Prss35      0  1.438266 0.708 0.183         0
    ## Ckb         0  1.424501 0.955 0.412         0
    ## Htra1       0  1.424291 0.833 0.228         0
    ## P2ry14      0  1.421968 0.838 0.145         0
    ## Actg2       0  1.406082 0.721 0.110         0
    ## Loxl2       0  1.398569 0.867 0.434         0
    ## Tnc         0  1.365862 0.855 0.224         0
    ## Filip1l     0  1.274808 0.762 0.156         0
    ## Des         0  1.252093 0.785 0.225         0
    ## Acta2       0  1.213925 0.946 0.620         0
    ## Wnt5a       0  1.208393 0.766 0.230         0

### Cluster 1

``` r
meso_markers_highres[[1 + 1]][n_print,]
```

    ##               p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Tcf21             0 1.8021338 0.972 0.429         0
    ## Gyg               0 1.7391961 0.989 0.412         0
    ## Mfap4             0 1.5380255 0.998 0.702         0
    ## Cebpb             0 1.4974228 0.957 0.486         0
    ## G0s2              0 1.4546530 0.821 0.258         0
    ## Adh1              0 1.3486055 0.969 0.399         0
    ## Fhl1              0 1.3040038 0.980 0.744         0
    ## Limch1            0 1.1457287 0.901 0.219         0
    ## Selenbp1          0 1.0916150 0.915 0.407         0
    ## Gpx3              0 0.9502986 0.764 0.361         0
    ## Lbh               0 0.9480626 0.878 0.303         0
    ## 6030408B16Rik     0 0.9422001 0.851 0.236         0
    ## Col13a1           0 0.9414232 0.862 0.254         0
    ## Meox2             0 0.9396783 0.862 0.250         0
    ## Tagln2            0 0.9339594 0.987 0.692         0
    ## Macf1             0 0.9229008 0.900 0.387         0
    ## Ppp1r14a          0 0.9221137 0.932 0.494         0
    ## Enpep             0 0.9170086 0.754 0.179         0
    ## Cebpd             0 0.9015292 0.766 0.472         0
    ## Ifitm1            0 0.8751800 0.693 0.303         0

### Cluster 2

``` r
meso_markers_highres[[2 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Dcn          0  2.590543 0.668 0.127         0
    ## Lum          0  2.096589 0.479 0.047         0
    ## Col1a1       0  2.014426 0.984 0.770         0
    ## Serpinf1     0  1.906197 0.706 0.057         0
    ## Col1a2       0  1.811145 0.985 0.840         0
    ## Mfap5        0  1.795427 0.707 0.154         0
    ## Col3a1       0  1.659136 0.961 0.796         0
    ## Ptn          0  1.654842 0.568 0.121         0
    ## Dlk1         0  1.598978 0.544 0.225         0
    ## Gsn          0  1.536382 0.931 0.549         0
    ## Igfbp7       0  1.503902 0.892 0.452         0
    ## Igfbp4       0  1.324224 0.981 0.742         0
    ## Col14a1      0  1.319812 0.737 0.245         0
    ## S100a6       0  1.306620 0.766 0.198         0
    ## Meg3         0  1.274346 0.853 0.428         0
    ## Cpz          0  1.258053 0.527 0.046         0
    ## C1qtnf3      0  1.257884 0.281 0.015         0
    ## Clec3b       0  1.257150 0.705 0.305         0
    ## Cthrc1       0  1.254659 0.393 0.026         0
    ## Pcolce       0  1.227176 0.913 0.613         0

### Cluster 3

``` r
meso_markers_highres[[3 + 1]][n_print,]
```

    ##           p_val avg_logFC pct.1 pct.2 p_val_adj
    ## H2afz         0 1.4134423 0.956 0.657         0
    ## Hist1h2ap     0 1.3788171 0.477 0.077         0
    ## Hmgb2         0 1.2195780 0.765 0.329         0
    ## Mif           0 1.1473591 0.934 0.591         0
    ## Pclaf         0 1.1211947 0.590 0.085         0
    ## Ran           0 1.0875979 0.960 0.618         0
    ## Ube2c         0 1.0376455 0.403 0.055         0
    ## Tuba1b        0 1.0304509 0.980 0.817         0
    ## H2afx         0 1.0107236 0.585 0.154         0
    ## Hist1h2ae     0 0.9795241 0.435 0.072         0
    ## Hmgn2         0 0.9781727 0.993 0.744         0
    ## Gapdh         0 0.9382266 0.980 0.840         0
    ## Npm1          0 0.9256433 0.985 0.891         0
    ## Stmn1         0 0.9163134 0.763 0.205         0
    ## Birc5         0 0.9019783 0.550 0.071         0
    ## Selenoh       0 0.8860963 0.827 0.321         0
    ## Anp32b        0 0.8735221 0.936 0.543         0
    ## Hmgb1         0 0.8724610 0.997 0.941         0
    ## Cenpa         0 0.8501519 0.438 0.072         0
    ## Cks2          0 0.8360244 0.521 0.084         0

### Cluster 4

``` r
meso_markers_highres[[4 + 1]][n_print,]
```

    ##         p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Macf1       0 1.7700599 0.918 0.434         0
    ## Npnt        0 1.6726925 0.836 0.312         0
    ## Fn1         0 1.5807968 0.915 0.556         0
    ## Enpep       0 1.4523326 0.785 0.230         0
    ## Cdh11       0 1.2728185 0.973 0.610         0
    ## Limch1      0 1.2505241 0.751 0.300         0
    ## Slit2       0 1.2386763 0.706 0.200         0
    ## Nebl        0 1.1776991 0.791 0.257         0
    ## Cpm         0 1.1749249 0.527 0.119         0
    ## Slc27a6     0 1.1720122 0.687 0.188         0
    ## Mfap4       0 1.1471252 0.959 0.734         0
    ## Colec12     0 1.0882855 0.905 0.494         0
    ## Plxdc2      0 1.0615415 0.894 0.550         0
    ## Itga8       0 1.0232868 0.726 0.248         0
    ## Frem1       0 1.0175317 0.691 0.194         0
    ## Prex2       0 1.0079951 0.723 0.282         0
    ## Col13a1     0 0.9962600 0.844 0.313         0
    ## Slc38a5     0 0.9652047 0.800 0.307         0
    ## Lbh         0 0.9466547 0.771 0.369         0
    ## Dock4       0 0.9391997 0.606 0.163         0

### Cluster 5

``` r
meso_markers_highres[[5 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Cox4i2       0  3.163604 0.978 0.039         0
    ## Gucy1a1      0  2.866603 0.990 0.287         0
    ## Postn        0  2.650208 0.993 0.211         0
    ## Gucy1b1      0  2.624873 0.980 0.214         0
    ## Itga1        0  2.353818 0.960 0.138         0
    ## Higd1b       0  2.328165 0.896 0.007         0
    ## Ndufa4l2     0  2.277715 0.954 0.115         0
    ## Mfge8        0  2.169474 0.997 0.688         0
    ## Pdgfrb       0  2.168905 0.960 0.161         0
    ## Mcam         0  2.155575 0.942 0.052         0
    ## Emid1        0  2.145206 0.927 0.111         0
    ## Pde5a        0  2.041829 0.965 0.398         0
    ## Pdzd2        0  2.030477 0.849 0.026         0
    ## Col4a1       0  2.021871 0.986 0.515         0
    ## Col4a2       0  1.960193 0.962 0.370         0
    ## Pcdh18       0  1.912913 0.916 0.110         0
    ## Ebf1         0  1.848619 0.927 0.089         0
    ## Trpc6        0  1.832859 0.861 0.017         0
    ## Itm2a        0  1.807121 0.967 0.382         0
    ## Fermt2       0  1.651868 0.982 0.637         0

### Cluster 6

``` r
meso_markers_highres[[6 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Enpp2        0  2.448455 0.972 0.206         0
    ## Crh          0  2.072185 0.508 0.011         0
    ## Hhip         0  2.017895 0.956 0.274         0
    ## Ednrb        0  1.983063 0.911 0.145         0
    ## Fam105a      0  1.939520 0.935 0.127         0
    ## H19          0  1.744705 0.826 0.236         0
    ## Crispld2     0  1.711006 0.821 0.121         0
    ## Pdgfra       0  1.637922 0.915 0.306         0
    ## Txnip        0  1.598657 0.861 0.328         0
    ## Rgs2         0  1.594711 0.837 0.196         0
    ## Penk         0  1.577194 0.717 0.137         0
    ## Tsc22d3      0  1.556421 0.850 0.324         0
    ## Aldh2        0  1.525440 0.974 0.823         0
    ## Nnat         0  1.493042 0.856 0.321         0
    ## Eln          0  1.438281 0.983 0.692         0
    ## Agt          0  1.428151 0.979 0.249         0
    ## Gm13889      0  1.416430 0.856 0.213         0
    ## Igf1         0  1.401142 0.925 0.355         0
    ## Mettl7a1     0  1.398585 0.805 0.158         0
    ## Glul         0  1.380768 0.882 0.407         0

### Cluster 7

``` r
meso_markers_highres[[7 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Upk3b        0  3.286810 0.993 0.003         0
    ## Msln         0  3.285169 0.706 0.001         0
    ## C3           0  3.147639 0.760 0.008         0
    ## Igfbp6       0  2.633647 0.935 0.072         0
    ## Bcam         0  2.520723 0.866 0.066         0
    ## Igfbp5       0  2.431049 0.937 0.157         0
    ## Clu          0  2.419801 0.913 0.046         0
    ## Aldh1a2      0  2.229059 0.956 0.037         0
    ## Krt19        0  2.171626 0.916 0.001         0
    ## Sfrp2        0  2.141841 0.803 0.146         0
    ## Crip1        0  2.140859 0.921 0.397         0
    ## Ezr          0  2.064209 0.924 0.050         0
    ## Slurp1       0  2.063462 0.297 0.000         0
    ## Cav1         0  2.056195 0.950 0.197         0
    ## Slc9a3r1     0  2.034493 0.950 0.022         0
    ## Upk1b        0  2.003140 0.942 0.094         0
    ## Krt7         0  1.986889 0.803 0.002         0
    ## Hspb1        0  1.962249 0.963 0.299         0
    ## Gbp2         0  1.900275 0.598 0.030         0
    ## S100a6       0  1.890141 0.849 0.276         0

### Cluster 8

``` r
meso_markers_highres[[8 + 1]][n_print,]
```

    ##           p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Hist1h2ap     0  1.871968 0.601 0.111         0
    ## Hmgb2         0  1.697567 0.966 0.365         0
    ## Htra1         0  1.648883 0.879 0.329         0
    ## Ube2c         0  1.503592 0.697 0.079         0
    ## Stmn1         0  1.497443 0.945 0.253         0
    ## Top2a         0  1.492381 0.724 0.092         0
    ## Cenpa         0  1.451805 0.707 0.099         0
    ## Pclaf         0  1.413736 0.866 0.125         0
    ## H2afz         0  1.411883 0.993 0.684         0
    ## Cks2          0  1.334873 0.831 0.117         0
    ## Ccna2         0  1.333444 0.829 0.098         0
    ## Rrm2          0  1.321967 0.672 0.098         0
    ## Tubb4b        0  1.285789 0.874 0.440         0
    ## Lmnb1         0  1.281101 0.948 0.270         0
    ## Mki67         0  1.271615 0.814 0.088         0
    ## Prc1          0  1.257905 0.674 0.063         0
    ## Cenpf         0  1.239652 0.688 0.075         0
    ## Ccnb2         0  1.233960 0.786 0.100         0
    ## Birc5         0  1.231995 0.850 0.107         0
    ## Cdk1          0  1.200289 0.752 0.102         0

### Cluster 9

``` r
meso_markers_highres[[9 + 1]][n_print,]
```

    ##          p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Acta2        0  1.865169 1.000 0.676         0
    ## Eln          0  1.819380 0.985 0.695         0
    ## Tagln        0  1.704592 0.994 0.461         0
    ## Map3k7cl     0  1.626769 0.818 0.013         0
    ## Crip1        0  1.602831 0.946 0.401         0
    ## Actc1        0  1.571231 0.478 0.021         0
    ## Cxcl12       0  1.567188 0.706 0.205         0
    ## Myl9         0  1.538171 0.995 0.618         0
    ## Cav1         0  1.292805 0.823 0.207         0
    ## Tm4sf1       0  1.282246 0.672 0.111         0
    ## Wisp2        0  1.277377 0.596 0.013         0
    ## Gpc6         0  1.276395 0.718 0.116         0
    ## Itga4        0  1.229091 0.633 0.064         0
    ## Parm1        0  1.222428 0.635 0.092         0
    ## Csrp2        0  1.221446 0.761 0.433         0
    ## Myh11        0  1.216695 0.825 0.257         0
    ## Tpm2         0  1.150314 0.941 0.601         0
    ## Aoc3         0  1.119091 0.730 0.056         0
    ## Gadd45g      0  1.111235 0.583 0.154         0
    ## Heyl         0  1.090956 0.723 0.123         0

### Cluster 10

``` r
meso_markers_highres[[10 + 1]][n_print,]
```

    ##         p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Actg2       0 2.0769038 0.918 0.213         0
    ## Cnn1        0 1.8415319 0.800 0.103         0
    ## Myh11       0 1.8324136 0.888 0.257         0
    ## Nnat        0 1.7407203 0.844 0.328         0
    ## Ednrb       0 1.7070704 0.884 0.156         0
    ## Acta2       0 1.6763454 0.998 0.676         0
    ## Mylk        0 1.6441212 0.963 0.667         0
    ## Myl9        0 1.6047302 0.979 0.619         0
    ## Tpm2        0 1.4011013 0.931 0.602         0
    ## Igf1        0 1.3051025 0.842 0.365         0
    ## Tmem158     0 1.2458742 0.630 0.256         0
    ## Fam213a     0 1.2126161 0.776 0.212         0
    ## Dstn        0 1.1835183 0.953 0.743         0
    ## Ckb         0 1.1729968 0.911 0.509         0
    ## Crh         0 1.0566464 0.355 0.022         0
    ## Tagln       0 1.0557718 0.958 0.463         0
    ## Tpm1        0 1.0403484 0.992 0.846         0
    ## Myl6        0 0.9782969 0.994 0.947         0
    ## Fbxl22      0 0.9631445 0.498 0.034         0
    ## Svil        0 0.9468287 0.903 0.489         0

### Cluster 11

``` r
meso_markers_highres[[11 + 1]][n_print,]
```

    ##                  p_val avg_logFC pct.1 pct.2     p_val_adj
    ## Myl7      0.000000e+00 4.0432928 0.376 0.012  0.000000e+00
    ## Myl4      0.000000e+00 3.3805649 0.426 0.006  0.000000e+00
    ## Tnnt2     0.000000e+00 2.9100809 0.442 0.015  0.000000e+00
    ## Tnnc1     0.000000e+00 2.4675903 0.383 0.002  0.000000e+00
    ## Fabp3     0.000000e+00 2.4293149 0.391 0.042  0.000000e+00
    ## Slc25a4   0.000000e+00 1.4566140 0.967 0.860  0.000000e+00
    ## Tnnt1     0.000000e+00 1.3199855 0.635 0.073  0.000000e+00
    ## Eno3      0.000000e+00 1.2098463 0.462 0.039  0.000000e+00
    ## Cdh15     0.000000e+00 1.0701165 0.543 0.002  0.000000e+00
    ## Fitm1     0.000000e+00 0.6461060 0.495 0.003  0.000000e+00
    ## Actc1    2.519735e-322 3.3193472 0.647 0.028 5.492014e-318
    ## Ppp1r14b 5.088876e-322 0.8618547 0.599 0.207 1.109171e-317
    ## Csrp3    6.620480e-322 1.9875667 0.386 0.001 1.443000e-317
    ## Tnni3    4.901131e-321 2.4546008 0.378 0.002 1.068251e-316
    ## Ank1     3.099076e-319 0.7319070 0.536 0.009 6.754746e-315
    ## Pgam2    7.251387e-314 1.6274782 0.383 0.002 1.580512e-309
    ## Sln      2.993310e-306 2.0006361 0.371 0.001 6.524218e-302
    ## Cycs     1.262775e-305 1.2594365 0.794 0.531 2.752344e-301
    ## Cox6a2   6.936893e-304 1.7092175 0.371 0.001 1.511965e-299
    ## Pdlim4   2.684886e-300 1.1929857 0.657 0.110 5.851977e-296

### Cluster 12

``` r
meso_markers_highres[[12 + 1]][n_print,]
```

    ##         p_val avg_logFC pct.1 pct.2 p_val_adj
    ## Mpz         0  3.045782 0.614 0.001         0
    ## Dbi         0  2.456146 0.995 0.597         0
    ## Plp1        0  2.429933 0.945 0.006         0
    ## Sostdc1     0  2.082210 0.496 0.017         0
    ## Sox10       0  1.914786 0.921 0.000         0
    ## Cd9         0  1.883186 0.819 0.106         0
    ## Cnp         0  1.879504 0.838 0.052         0
    ## Mal         0  1.821225 0.668 0.000         0
    ## Matn2       0  1.815978 0.797 0.170         0
    ## Mbp         0  1.751278 0.627 0.012         0
    ## Arpc1b      0  1.692785 0.929 0.602         0
    ## Egfl8       0  1.621841 0.638 0.002         0
    ## Gatm        0  1.546819 0.523 0.009         0
    ## Fabp7       0  1.516997 0.649 0.011         0
    ## Ptprz1      0  1.509560 0.704 0.016         0
    ## Gfra3       0  1.409245 0.496 0.001         0
    ## Plekhb1     0  1.348575 0.721 0.002         0
    ## Sema3b      0  1.348103 0.710 0.022         0
    ## Ncam1       0  1.341968 0.811 0.146         0
    ## Ttyh1       0  1.335385 0.784 0.006         0

# Save markers for detailed meso clusters

``` r
meso_highres_myo_clusters <- c(0,6,8,10)
meso_highres_myo_names <- c("Mature myofibroblasts", "E18 myofibroblasts", "Proliferating myofibroblasts", "E12-E15 myofibroblasts")

wb_myo_markers <- createWorkbook()
for (idx in seq_along(meso_highres_myo_clusters)){
  i <- meso_highres_myo_clusters[idx] + 1
  addWorksheet(wb_myo_markers, meso_highres_myo_names[idx])
  writeData(wb_myo_markers, meso_highres_myo_names[idx], meso_markers_highres[[i]], rowNames = TRUE)
}
saveWorkbook(wb_myo_markers, file = "./figures/meso/myo_highres_marker_genes.xlsx", overwrite = TRUE)
```

``` r
wnt5a_all_myo <- VlnPlot(subset(meso_relabel, ident = c("Myofibroblast", "Proliferating myofibroblast")), "Wnt5a", group.by = "timepoint") +ggtitle("All myofibroblasts")
wnt5a_mature_myo <-  VlnPlot(subset(meso_highres, ident = 0), "Wnt5a", group.by = "timepoint") + ggtitle("Mature myofibroblasts")
wnt5a_proliferating_myo <- VlnPlot(subset(meso_highres, ident = 8), "Wnt5a", group.by = "timepoint") + ggtitle("Proliferating myofibroblasts")
wnt5a_e18_myo <- VlnPlot(subset(meso_highres, ident = 6), "Wnt5a", group.by = "timepoint") + ggtitle("'E18' myofibroblasts")
```
