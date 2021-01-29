Plottings outputs of the cellrank analysis
================
Nick Negretti
12/07/20

# Plot cellrank outputs of the lung endothelium

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
cellrank_out <- anndata::read_h5ad("./velocity/endo/only_metadata.h5ad")
#cellrank_out <- anndata::read_h5ad("./velocity/endo/after_cellrank.h5ad")
```

``` r
cellrank_out
```

    ## AnnData object with n_obs × n_vars = 22267 × 18783
    ##     obs: 'timepoint', 'bulk_cellype', 'cell_subtype', 'sample_batch', 'Clusters', '_X', '_Y', 'initial_size_spliced', 'initial_size_unspliced', 'initial_size', 'n_counts', 'velocity_self_transition', 'clusters_gradients', 'terminal_states', 'terminal_states_probs', 'to_terminal_states_dp', 'initial_states', 'initial_states_probs', 'from_initial_states_dp', 'velocity_pseudotime', 'latent_time', 'dpt_pseudotime', 'S_score', 'G2M_score', 'phase', 'velocity_length', 'velocity_confidence', 'velocity_confidence_transition', 'root_cells', 'end_points'
    ##     var: 'fit_alpha', 'fit_beta', 'fit_gamma', 'fit_t_', 'fit_scaling', 'fit_alignment_scaling', 'fit_std_u', 'fit_std_s', 'fit_likelihood', 'fit_u0', 'fit_s0', 'fit_pval_steady', 'fit_steady_u', 'fit_steady_s', 'velocity_genes', 'fit_r2', 'fit_diff_kinetics', 'fit_pval_kinetics', 'to miEC_1 corr', 'to Lymphatic_1 corr', 'to Lymphatic_2 corr', 'to Arterial maEC corr', 'to Venous maEC corr', 'to miEC_2 corr', 'to Proliferating\nmiEC corr', 'to Car4+ corr', 'to miEC_3 corr', 'to miEC_1 qval', 'to Lymphatic_1 qval', 'to Lymphatic_2 qval', 'to Arterial maEC qval', 'to Venous maEC qval', 'to miEC_2 qval', 'to Proliferating\nmiEC qval', 'to Car4+ qval', 'to miEC_3 qval'
    ##     uns: 'from_initial_states_names', 'initial_states_names', 'to_terminal_states_names'
    ##     obsm: 'X_diffmap', 'X_pca', 'X_umap', 'from_initial_states', 'macrostates_bwd', 'macrostates_fwd', 'to_terminal_states', 'velocity_umap'
    ##     varm: 'PCs', 'fit_pvals_kinetics', 'loss'

``` r
endo_ms <- Matrix::readMM("./velocity/endo/endo_expression_Ms.mtx")
```

``` r
gene <- "Car4"
include_celltypes <- c("Proliferating miEC", "miEC", "*Car4+* EC")
exp_over_time <- data.frame(expression = endo_ms[,colnames(cellrank_out) == gene],
                            Cluster = cellrank_out$obs$cell_subtype,
                            latent_time = cellrank_out$obs$latent_time
)

unique(exp_over_time$Cluster)
```

    ## [1] miEC                Venous maEC         Lymphatic          
    ## [4] Arterial maEC       Proliferating\nmiEC Car4+              
    ## 6 Levels: Arterial maEC Car4+ Lymphatic Proliferating\nmiEC ... miEC

``` r
levels(exp_over_time$Cluster)[levels(exp_over_time$Cluster) == "Car4+"] <- "*Car4+* EC"
levels(exp_over_time$Cluster)[levels(exp_over_time$Cluster) == "Proliferating\nmiEC"] <- "Proliferating miEC"
exp_over_time$Cluster <- ordered(as.factor(exp_over_time$Cluster), c("Arterial maEC", "Venous maEC", "miEC", "*Car4+* EC", "Proliferating miEC", "Lymphatic"))

exp_over_time <- exp_over_time[exp_over_time$Cluster %in% include_celltypes,]

car4_exp <- ggplot(exp_over_time, aes(y = expression, x = latent_time, color = Cluster)) +
        geom_point(alpha = 0.5) +
        scale_color_manual(aesthetics = c("color", "fill"), values= color_category_20[3:5], guide = guide_legend(override.aes = list(size=3), order = 1)) +
        theme(legend.key = element_blank(),
              legend.text = element_markdown(size=14),
              legend.title = element_markdown(size=16),
              #axis.text.x  = element_text(size=14),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size=14),
              axis.title.x = element_text(size=14, vjust = -1),
              axis.title.y = element_text(size=14),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "white", color = "black"),
              aspect.ratio = 0.8,
              plot.title = element_markdown(size = 16)
        ) +
        scale_x_continuous(breaks = c(0,1)) +
        ylab("Expression") +
        xlab("Latent time") +
        ggtitle(paste0("*", gene, "*"))

car4_exp
```

![](./velocity/endo/endo_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
gene <- "Kdr"
include_celltypes <- c("Proliferating miEC", "miEC", "*Car4+* EC")
exp_over_time <- data.frame(expression = endo_ms[,colnames(cellrank_out) == gene],
                            Cluster = cellrank_out$obs$cell_subtype,
                            latent_time = cellrank_out$obs$latent_time
)

unique(exp_over_time$Cluster)
```

    ## [1] miEC                Venous maEC         Lymphatic          
    ## [4] Arterial maEC       Proliferating\nmiEC Car4+              
    ## 6 Levels: Arterial maEC Car4+ Lymphatic Proliferating\nmiEC ... miEC

``` r
levels(exp_over_time$Cluster)[levels(exp_over_time$Cluster) == "Car4+"] <- "*Car4+* EC"
levels(exp_over_time$Cluster)[levels(exp_over_time$Cluster) == "Proliferating\nmiEC"] <- "Proliferating miEC"
exp_over_time$Cluster <- ordered(as.factor(exp_over_time$Cluster), c("Arterial maEC", "Venous maEC", "miEC", "*Car4+* EC", "Proliferating miEC", "Lymphatic"))

exp_over_time <- exp_over_time[exp_over_time$Cluster %in% include_celltypes,]

kdr_exp <- ggplot(exp_over_time, aes(y = expression, x = latent_time, color = Cluster)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(aesthetics = c("color", "fill"), values= color_category_20[3:5], guide = guide_legend(override.aes = list(size=3), order = 1)) +
  theme(legend.key = element_blank(),
        legend.text = element_markdown(size=14),
        legend.title = element_markdown(size=16),
        #axis.text.x  = element_text(size=14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=14, vjust = -1),
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        aspect.ratio = 0.8,
        plot.title = element_markdown(size = 16)
  ) +
  scale_x_continuous(breaks = c(0,1)) +
  ylab("Expression") +
  xlab("Latent time") +
  ggtitle(paste0("*", gene, "*"))

kdr_exp
```

![](./velocity/endo/endo_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
saveTiff("./figures/endo/latent_time/car4_kdr.tiff",
         car4_exp + kdr_exp + plot_layout(ncol = 2, nrow = 1, guides = "collect"),
         width = 15, height = 4)
```

    ## Warning in normalizePath(paste0(path, "_tmp.tiff")): path[1]="./figures/endo/
    ## latent_time/car4_kdr.tiff_tmp.tiff": No such file or directory

    ## [1] TRUE

``` r
ggdata_endo <- data.frame(x = cellrank_out$obsm$X_umap[,1],
                     y = cellrank_out$obsm$X_umap[,2],
                          latent_time = cellrank_out$obs$latent_time,
                     group = as.character(cellrank_out$obs$cell_subtype),
                     timepoint = as.character(cellrank_out$obs$timepoint))
head(ggdata_endo)
```

    ##           x          y latent_time       group timepoint
    ## 1 -4.074359 -1.9203233  0.02768251        miEC       E12
    ## 2 -5.113158 -1.5107509  0.05834358        miEC       E12
    ## 3 -3.916332 -1.8675787  0.04629836        miEC       E12
    ## 4  2.025012 -3.0183446  0.04420196        miEC       E12
    ## 5  5.131426 -0.5862549  0.33998461        miEC       E12
    ## 6 -1.381284 -5.4132347  0.06785460 Venous maEC       E12

``` r
colnames(ggdata_endo) <- c("x", "y", "latent_time", "group", "timepoint")
mean(ggdata_endo$latent_time)
```

    ## [1] 0.5328144

``` r
end_latent_time <- ggplot(ggdata_endo %>% arrange(latent_time), aes(x = x, y = y, color = latent_time)) +
  geom_point(size = 0.5) +
  umap_theme() +
  theme(aspect.ratio = 1,
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
  plot.title = element_text(size = 18)) +
  scale_color_viridis(name = "Latent time", direction = -1) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle("Latent time")

end_latent_time
```

![](./velocity/endo/endo_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
saveTiff("./figures/endo/latent_time_UMAP.tiff",
         end_latent_time,
         width = 6, height = 4.5)
```

    ## Warning in normalizePath(paste0(path, "_tmp.tiff")): path[1]="./figures/endo/
    ## latent_time_UMAP.tiff_tmp.tiff": No such file or directory

    ## [1] TRUE

``` r
to_term <- cellrank_out$obsm$to_terminal_states

colnames(to_term) <- cellrank_out$uns$to_terminal_states_names

# Combine similar terminal types
to_term_unique <- to_term[,c(4,5,7,8)]

head(to_term_unique)
```

    ##      Arterial maEC Venous maEC Proliferating\nmiEC       Car4+
    ## [1,]   0.011960215 0.019704247        0.0017609942 0.025674085
    ## [2,]   0.019250436 0.030857740        0.0021600244 0.023113562
    ## [3,]   0.013471921 0.021899814        0.0019379260 0.023770796
    ## [4,]   0.009447754 0.016621258        0.0021754607 0.027077385
    ## [5,]   0.002235276 0.006700978        0.0000954821 0.009095503
    ## [6,]   0.025645083 0.039427649        0.0025135529 0.019334186

``` r
to_term <- cbind(to_term_unique, rowSums(to_term[,c(1, 6, 9)]), rowSums(to_term[,c(2, 3)]))


colnames(to_term) <- c(colnames(to_term_unique), "miEC", "Lymphatic")
```

``` r
to_term <- as.data.frame(to_term)
to_term$timepoint <- ordered(as.factor(cellrank_out$obs$timepoint), c("E12", "E15", "E18", "P0", "P3", "P5", "P7", "P14"))

head(to_term)
```

    ##   Arterial maEC Venous maEC Proliferating\nmiEC       Car4+      miEC
    ## 1   0.011960215 0.019704247        0.0017609942 0.025674085 0.9087324
    ## 2   0.019250436 0.030857740        0.0021600244 0.023113562 0.8840505
    ## 3   0.013471921 0.021899814        0.0019379260 0.023770796 0.9043489
    ## 4   0.009447754 0.016621258        0.0021754607 0.027077385 0.9231763
    ## 5   0.002235276 0.006700978        0.0000954821 0.009095503 0.9773286
    ## 6   0.025645083 0.039427649        0.0025135529 0.019334186 0.8681315
    ##     Lymphatic timepoint
    ## 1 0.032168065       E12
    ## 2 0.040567785       E12
    ## 3 0.034570598       E12
    ## 4 0.021501796       E12
    ## 5 0.004544193       E12
    ## 6 0.044948067       E12

``` r
head(to_term[cellrank_out$obs$cell_subtype == "Car4+",])
```

    ##     Arterial maEC Venous maEC Proliferating\nmiEC      Car4+      miEC
    ## 775  0.0006197730 0.001012945        2.951419e-06 0.06086449 0.9373938
    ## 796  0.0002748540 0.001237382        3.028182e-06 0.36251190 0.6358924
    ## 814  0.0019920317 0.007077044        8.928033e-04 0.05172400 0.9358478
    ## 836  0.0002224017 0.001672800        2.110708e-05 0.53364125 0.4642028
    ## 842  0.0003801710 0.001286513        2.313533e-06 0.41194386 0.5862970
    ## 862  0.0004194676 0.001255191        1.359633e-05 0.20635573 0.7917574
    ##        Lymphatic timepoint
    ## 775 1.060400e-04       E18
    ## 796 8.045579e-05       E18
    ## 814 2.466332e-03       E18
    ## 836 2.396656e-04       E18
    ## 842 9.013576e-05       E18
    ## 862 1.986499e-04       E18

``` r
ggplot(melt(to_term[cellrank_out$obs$cell_subtype == "Car4+",]), aes(x = value, color = variable, fill = variable)) +
        geom_density(alpha = 0.3) +
        scale_x_continuous(limits = c(0.2,1)) +
        facet_grid(vars(timepoint), scales = "free")
```

    ## Using timepoint as id variables

    ## Warning: Removed 12016 rows containing non-finite values (stat_density).

![](./velocity/endo/endo_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
car4_plot <- ggplot(melt(to_term[cellrank_out$obs$cell_subtype == "Car4+",]), aes(x = value, color = variable, fill = variable)) +
  geom_density(alpha = 0.3) +
  scale_x_continuous(limits = c(0.2,1)) +
        ggtitle("Fate probability of 'Car4+' cells")
```

    ## Using timepoint as id variables

``` r
miEC_plot <- ggplot(melt(to_term[cellrank_out$obs$cell_subtype == "miEC",]), aes(x = value, color = variable, fill = variable)) +
        geom_density(alpha = 0.3) +
        scale_x_continuous(limits = c(0.2,1)) +
        ggtitle("Fate probability of 'miEC' cells")
```

    ## Using timepoint as id variables

``` r
car4_plot + miEC_plot + patchwork::plot_layout(ncol = 1)
```

    ## Warning: Removed 12016 rows containing non-finite values (stat_density).

    ## Warning: Removed 69713 rows containing non-finite values (stat_density).

![](./velocity/endo/endo_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

``` r
summary_stats <- melt(to_term[cellrank_out$obs$cell_subtype == "miEC",]) %>% group_by(variable) %>% summarize(mean = mean(value) * 100,
                                                                                                              median = median(value) * 100,
                                                                                                              sd = sd(value) * 100,
                                                                                                              sem = (sd(value) / sqrt(length(value))) * 100,
                                                                                                              q1 = quantile(value, 1/4) * 100,
                                                                                                              q2 = quantile(value, 3/4) * 100,
                                                                                                              #ci1 = get_lower_ci(value) * 100,
                                                                                                              #ci2 = get_upper_ci(value) * 100
)
```

    ## Using timepoint as id variables

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
miEC_barplot <- ggplot(summary_stats, aes(y = mean, x = variable)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  ggtitle("Fate of 'miEC' cells") +
  prism_theme() +
  theme(aspect.ratio = 1.6) +
  ylab("Transition probability")

miEC_barplot_noself <- ggplot(summary_stats[summary_stats$variable != "miEC",], aes(y = mean * 100, x = variable)) +
        geom_col() +
        geom_errorbar(aes(ymin = mean * 100 - sem * 100, ymax = mean * 100 + sem * 100), width = 0.2) +
        ggtitle("Fate of miEC cells") +
        prism_theme() +
        theme(aspect.ratio = 1.6) +
        ylab("Transition probability")

miEC_barplot + miEC_barplot_noself
```

![](./velocity/endo/endo_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
to_term_long <- melt(to_term[cellrank_out$obs$cell_subtype == "miEC",])
```

    ## Using timepoint as id variables

``` r
to_term_long$value <- to_term_long$value * 100
to_term_long_above_zero <- to_term_long[to_term_long$value > 1,]


get_early_pct <- function(data, summary_stats, cell_group){
  timepoints <- data$obs$timepoint[data$obs$cell_subtype == cell_group]
  early <- sum(timepoints %in% c("E12", "E15", "E18")) * 100 / length(timepoints)
  return(early * summary_stats[summary_stats$variable == cell_group, "median"] / 100) # scale vals to max bar
}

for (type in unique(summary_stats$variable)){
  summary_stats[summary_stats$variable == type, "early"] <- get_early_pct(cellrank_out, summary_stats, type)
}


fate_miec_bar <- ggplot() +
        geom_col(data = summary_stats, aes(y = median, x = variable), color = "black", fill = "white") +
        #geom_col(data = summary_stats, aes(x = variable, y = early), fill = "#CDCDCD", width = 0.85) + # color early here, overlay on the other line
        geom_errorbar(data = summary_stats, aes(y = median, x = variable, ymin = q1, ymax = q2), width = 0.2) +
        #geom_boxplot(data = to_term_long, aes(y = value, x = variable), outlier.shape = NA, lwd=0.25) +
        #ggbeeswarm::geom_quasirandom(data = to_term_long, aes(y = value, x = variable), size = 0.8, shape = 1, alpha = 0.6) +
        ggtitle("Fate of miEC cells") +
        prism_theme() +
        theme(aspect.ratio = 1.6,
              plot.title = element_markdown(size = 14),
              axis.text.x = element_markdown(size = 14)
        ) +
        theme(aspect.ratio = 1.6) +
        ylab("Transition probability (%)")
fate_miec_bar
```

![](./velocity/endo/endo_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
saveTiff("./figures/endo/fate_miEC.tiff",
         fate_miec_bar,
         width = 3, height = 5)
```

    ## Warning in normalizePath(paste0(path, "_tmp.tiff")): path[1]="./figures/endo/
    ## fate_miEC.tiff_tmp.tiff": No such file or directory

    ## Warning in normalizePath(path): path[1]="./figures/endo/fate_miEC.tiff": No such
    ## file or directory

    ## [1] TRUE

``` r
to_term_orig <- cellrank_out$obsm$to_terminal_states

colnames(to_term_orig) <- cellrank_out$uns$to_terminal_states_names

ggdata_endo_to_Car4 <- data.frame(x = cellrank_out$obsm$X_umap[,1],
                             y = cellrank_out$obsm$X_umap[,2],
                             to_Car4 = to_term_orig[,"Car4+"],
                             cell_id = rownames(cellrank_out$obs))


to_Car4_time_plot <- ggplot(ggdata_endo_to_Car4 %>% arrange(to_Car4), aes(x = x, y = y, color = to_Car4)) + geom_point(alpha = 0.6, shape = 16) + umap_theme() + theme(aspect.ratio = 1) + scale_color_viridis(limits = c(0,max(ggdata_endo_to_Car4$to_Car4)))
to_Car4_time_plot
```

![](./velocity/endo/endo_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
saveRDS(ggdata_endo_to_Car4, "./velocity/endo/data/endo_car4_prob.rds")
```

``` r
from_initial <- cellrank_out$obsm$from_initial_states

colnames(from_initial) <- cellrank_out$uns$from_initial_states_names

head(from_initial)
```

    ##         Lymphatic Proliferating\nmiEC Arterial maEC      miEC        Car4+
    ## [1,] 0.0007350025          0.05611097  0.0001468950 0.9426707 0.0003364242
    ## [2,] 0.0010248540          0.07983325  0.0002091305 0.9184531 0.0004797040
    ## [3,] 0.0005110572          0.03740273  0.0001010959 0.9617545 0.0002305707
    ## [4,] 0.0007974362          0.06492642  0.0001716621 0.9337107 0.0003937840
    ## [5,] 0.0013954276          0.11156992  0.0003771833 0.8858651 0.0007923923
    ## [6,] 0.0009970946          0.07079114  0.0001906157 0.9275893 0.0004318048

``` r
# Combine similar terminal types
#from_initial_no_T2 <- from_initial[,c(1,2,3,5)]

#from_initial <- cbind(from_initial_no_T2, rowSums(from_initial[,c(4,6,7)])) # Merge the Type II fates
#colnames(from_initial) <- c(colnames(from_initial_no_T2), "Type II")
from_initial <- as.data.frame(from_initial)
from_initial$timepoint <- ordered(as.factor(cellrank_out$obs$timepoint), c("E12", "E15", "E18", "P0", "P3", "P5", "P7", "P14"))

head(from_initial[cellrank_out$obs$cell_subtype == "Car4+",])
```

    ##       Lymphatic Proliferating\nmiEC Arterial maEC      miEC        Car4+
    ## 775 0.001772509           0.1550197  0.0011538667 0.8402722 0.0017816695
    ## 796 0.001860410           0.1987458  0.0012026782 0.7863850 0.0118061465
    ## 814 0.001236359           0.1094682  0.0003037521 0.8882972 0.0006944475
    ## 836 0.001937752           0.2001877  0.0012585656 0.7860964 0.0105195862
    ## 842 0.001856539           0.1951749  0.0011821461 0.7926166 0.0091697451
    ## 862 0.001756838           0.1655073  0.0010292002 0.8296676 0.0020389930
    ##     timepoint
    ## 775       E18
    ## 796       E18
    ## 814       E18
    ## 836       E18
    ## 842       E18
    ## 862       E18

``` r
to_car4_plot <- ggplot(melt(from_initial[cellrank_out$obs$cell_subtype == "Car4+",]), aes(x = value, color = variable, fill = variable)) +
  geom_density(alpha = 0.3) +
  scale_x_continuous(limits = c(0.05,1)) +
  ggtitle("Fate probability towards 'Car4+' cells")
```

    ## Using timepoint as id variables

``` r
to_car4_plot
```

    ## Warning: Removed 7266 rows containing non-finite values (stat_density).

    ## Warning: Groups with fewer than two data points have been dropped.

    ## Warning in max(ids, na.rm = TRUE): no non-missing arguments to max; returning -
    ## Inf

![](./velocity/endo/endo_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
colnames(from_initial)[2] <- "Proliferating miEC"
colnames(from_initial)[5] <- "*Car4*+"
summary_stats <- melt(from_initial[cellrank_out$obs$cell_subtype == "Car4+",]) %>% group_by(variable) %>% summarize(mean = mean(value) * 100,
                                                                                                                    median = median(value) * 100,
                                                                                                                    sd = sd(value) * 100,
                                                                                                                    sem = (sd(value) / sqrt(length(value))) * 100,
                                                                                                                    q1 = quantile(value, 1/4) * 100,
                                                                                                                    q2 = quantile(value, 3/4) * 100,
)
```

    ## Using timepoint as id variables

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
transitional_barplot <- ggplot(summary_stats, aes(y = mean, x = variable)) +
        geom_col() +
        geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
        ggtitle("Towards 'Car4+' cells") +
        prism_theme() +
        theme(aspect.ratio = 1.6) +
        ylab("Transition probability")

transitional_barplot
```

![](./velocity/endo/endo_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
from_initial_long <- melt(from_initial[cellrank_out$obs$cell_subtype == "Car4+",])
```

    ## Using timepoint as id variables

``` r
from_initial_long$value <- from_initial_long$value * 100
from_initial_long_above_zero <- from_initial_long[from_initial_long$value > 1,]

get_early_pct <- function(data, summary_stats, cell_group){
  timepoints <- data$obs$timepoint[data$obs$cell_subtype == cell_group]
  early <- sum(timepoints %in% c("E12", "E15", "E18")) * 100 / length(timepoints)
  return(early * summary_stats[summary_stats$variable == cell_group, "median"] / 100) # scale vals to max bar
}

for (type in unique(summary_stats$variable)){
  summary_stats[summary_stats$variable == type, "early"] <- get_early_pct(cellrank_out, summary_stats, type)
}

source_car4_bar <- ggplot() +
        geom_col(data = summary_stats, aes(y = median, x = variable), color = "black", fill = "white") +

        #geom_col(data = summary_stats, aes(x = variable, y = early), fill = "#CDCDCD", width = 0.84) + # color early here, overlay on the other line
        geom_errorbar(data = summary_stats, aes(y = median, x = variable, ymin = q1, ymax = q2), width = 0.6) +
        #geom_boxplot(data = from_initial_long, aes(y = value, x = variable), outlier.shape = NA, lwd=0.25) +
        #ggbeeswarm::geom_quasirandom(data = from_initial_long, aes(y = value, x = variable), size = 0.8, shape = 1, alpha = 0.6) +
        ggtitle("Source of *Car4*+ cells") +
        prism_theme() +
        theme(aspect.ratio = 1.6,
              plot.title = element_markdown(size = 14),
              axis.text.x = element_markdown(size = 14)
        ) +
        ylab("Transition probability (%)")

source_car4_bar
```

![](./velocity/endo/endo_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

``` r
saveTiff("./figures/endo/to_car4_probability.tiff",
         source_car4_bar,
         width = 3, height = 5)
```

    ## Warning in normalizePath(paste0(path, "_tmp.tiff")): path[1]="./figures/endo/
    ## to_car4_probability.tiff_tmp.tiff": No such file or directory

    ## Warning in normalizePath(path): path[1]="./figures/endo/
    ## to_car4_probability.tiff": No such file or directory

    ## [1] TRUE

``` r
saveTiff("./figures/endo/fate_miEC_source_car4.tiff",
         fate_miec_bar + source_car4_bar + plot_layout(ncol = 2, nrow = 1) & theme(plot.margin = margin(5.5,105.5,5.5,5.5)),
         width = 15, height = 4)
```

    ## Warning in normalizePath(paste0(path, "_tmp.tiff")): path[1]="./figures/endo/
    ## fate_miEC_source_car4.tiff_tmp.tiff": No such file or directory

    ## Warning in normalizePath(path): path[1]="./figures/endo/
    ## fate_miEC_source_car4.tiff": No such file or directory

    ## [1] TRUE
