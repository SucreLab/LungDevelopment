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
cellrank_out <- anndata::read_h5ad("./velocity/meso/with_latent_time.h5ad")
ggdata_meso <- data.frame(x = cellrank_out$obsm$X_umap[,1],
                          y = cellrank_out$obsm$X_umap[,2],
                          velocity_length = cellrank_out$obs$velocity_length,
                          group = as.character(cellrank_out$obs$cell_subtype),
                          timepoint = as.character(cellrank_out$obs$timepoint))
```

``` r
cellrank_out
```

    ## AnnData object with n_obs × n_vars = 58560 × 21781
    ##     obs: 'timepoint', 'bulk_cellype', 'cell_subtype', 'sample_batch', 'Clusters', '_X', '_Y', 'initial_size_spliced', 'initial_size_unspliced', 'initial_size', 'n_counts', 'velocity_self_transition', 'clusters_gradients', 'terminal_states', 'terminal_states_probs', 'to_terminal_states_dp', 'initial_states', 'initial_states_probs', 'from_initial_states_dp', 'velocity_length', 'velocity_confidence', 'velocity_confidence_transition', 'velocity_pseudotime', 'latent_time'
    ##     var: 'fit_alpha', 'fit_beta', 'fit_gamma', 'fit_t_', 'fit_scaling', 'fit_alignment_scaling', 'fit_std_u', 'fit_std_s', 'fit_likelihood', 'fit_u0', 'fit_s0', 'fit_pval_steady', 'fit_steady_u', 'fit_steady_s', 'velocity_genes', 'fit_r2', 'fit_diff_kinetics', 'fit_pval_kinetics'
    ##     uns: 'T_bwd_params', 'T_fwd_params', 'cell_subtype_colors', 'clusters_gradients_colors', 'eig_bwd', 'eig_fwd', 'from_initial_states_colors', 'from_initial_states_names', 'initial_states_colors', 'initial_states_names', 'neighbors', 'pca', 'recover_dynamics', 'terminal_states_colors', 'terminal_states_names', 'timepoint_colors', 'to_terminal_states_colors', 'to_terminal_states_names', 'velocity_graph', 'velocity_graph_neg'
    ##     obsm: 'X_pca', 'X_umap', 'from_initial_states', 'macrostates_bwd', 'macrostates_fwd', 'to_terminal_states', 'velocity_umap'
    ##     varm: 'PCs', 'fit_pvals_kinetics', 'loss'
    ##     layers: 'fit_t', 'velocity'
    ##     obsp: 'T_bwd', 'T_fwd', 'connectivities', 'distances'

``` r
meso_ms <- Matrix::readMM("./velocity/meso/meso_expression_Ms.mtx.gz")
```

``` r
gene <- "Plin2"
include_celltypes <- c("Proliferating *Wnt2*+ FB", "*Wnt2*+ FB", "Proliferating myofibroblast", "Myofibroblast")
exp_over_time <- data.frame(expression = meso_ms[,colnames(cellrank_out) == gene],
                            Cluster = cellrank_out$obs$cell_subtype,
                            latent_time = cellrank_out$obs$latent_time
)

unique(exp_over_time$Cluster)
```

    ## [1] Adventitial fibroblast Prenatal Wnt2+         Pericyte              
    ## [4] Prenatal Myofibroblast Wnt2+                  Mesothelium           
    ## [7] Smooth muscle          Myofibroblast          Neuron                
    ## 9 Levels: Adventitial fibroblast Mesothelium Myofibroblast Neuron ... Wnt2+

``` r
levels(exp_over_time$Cluster)[levels(exp_over_time$Cluster) == "Prenatal Wnt2+"] <- "Proliferating *Wnt2*+ FB"
levels(exp_over_time$Cluster)[levels(exp_over_time$Cluster) == "Wnt2+"] <- "*Wnt2*+ FB"
levels(exp_over_time$Cluster)[levels(exp_over_time$Cluster) == "Prenatal Myofibroblast"] <- "Proliferating myofibroblast"

exp_over_time$Cluster <- ordered(as.factor(exp_over_time$Cluster), c("Proliferating *Wnt2*+ FB", "*Wnt2*+ FB", "Proliferating myofibroblast", "Myofibroblast", "Adventitial fibroblast", "Pericyte", "Mesothelium", "Smooth muscle", "Neuron"))

exp_over_time <- exp_over_time[exp_over_time$Cluster %in% include_celltypes,]

plin2_exp <- ggplot(exp_over_time, aes(y = expression, x = latent_time, color = Cluster)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(aesthetics = c("color", "fill"), values= color_category_20[1:4], guide = guide_legend(override.aes = list(size=3), order = 1)) +
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

plin2_exp
```

![](./velocity/meso/meso_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
gene <- "Tgfbi"
include_celltypes <- c("Proliferating *Wnt2*+ FB", "*Wnt2*+ FB", "Proliferating myofibroblast", "Myofibroblast")
exp_over_time <- data.frame(expression = meso_ms[,colnames(cellrank_out) == gene],
                            Cluster = cellrank_out$obs$cell_subtype,
                            latent_time = cellrank_out$obs$latent_time
)

unique(exp_over_time$Cluster)
```

    ## [1] Adventitial fibroblast Prenatal Wnt2+         Pericyte              
    ## [4] Prenatal Myofibroblast Wnt2+                  Mesothelium           
    ## [7] Smooth muscle          Myofibroblast          Neuron                
    ## 9 Levels: Adventitial fibroblast Mesothelium Myofibroblast Neuron ... Wnt2+

``` r
levels(exp_over_time$Cluster)[levels(exp_over_time$Cluster) == "Prenatal Wnt2+"] <- "Proliferating *Wnt2*+ FB"
levels(exp_over_time$Cluster)[levels(exp_over_time$Cluster) == "Wnt2+"] <- "*Wnt2*+ FB"
levels(exp_over_time$Cluster)[levels(exp_over_time$Cluster) == "Prenatal Myofibroblast"] <- "Proliferating myofibroblast"

exp_over_time$Cluster <- ordered(as.factor(exp_over_time$Cluster), c("Proliferating *Wnt2*+ FB", "*Wnt2*+ FB", "Proliferating myofibroblast", "Myofibroblast", "Adventitial fibroblast", "Pericyte", "Mesothelium", "Smooth muscle", "Neuron"))

exp_over_time <- exp_over_time[exp_over_time$Cluster %in% include_celltypes,]

tgfbi_exp <- ggplot(exp_over_time, aes(y = expression, x = latent_time, color = Cluster)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(aesthetics = c("color", "fill"), values= color_category_20[1:4], guide = guide_legend(override.aes = list(size=3), order = 1)) +
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

tgfbi_exp
```

![](./velocity/meso/meso_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
saveTiff("./figures/meso/plin2_tgfbi.tiff",
         plin2_exp + tgfbi_exp + plot_layout(ncol = 2, nrow = 1, guides = "collect"),
         width = 15, height = 4)
```

    ## Warning in normalizePath(paste0(path, "_tmp.tiff")): path[1]="./figures/meso/
    ## plin2_tgfbi.tiff_tmp.tiff": No such file or directory

    ## Warning in normalizePath(path): path[1]="./figures/meso/plin2_tgfbi.tiff": No
    ## such file or directory

    ## [1] TRUE

``` r
ggdata_meso <- data.frame(x = cellrank_out$obsm$X_umap[,1],
                     y = cellrank_out$obsm$X_umap[,2],
                     velocity_length = cellrank_out$obs$velocity_length,
                     group = as.character(cellrank_out$obs$cell_subtype),
                     timepoint = as.character(cellrank_out$obs$timepoint))
head(ggdata_meso)
```

    ##            x          y velocity_length                  group timepoint
    ## 1  0.1800986  1.6966943           14.79 Adventitial fibroblast       E12
    ## 2 -0.6809665  0.1726872           14.79         Prenatal Wnt2+       E12
    ## 3  2.2270545  2.1155751           14.02 Adventitial fibroblast       E12
    ## 4 -9.6195186 -0.6785916           15.86               Pericyte       E12
    ## 5 -0.2681494  7.6894823           17.03 Adventitial fibroblast       E12
    ## 6 -3.2464713 -8.1006193           14.35 Prenatal Myofibroblast       E12

``` r
colnames(ggdata_meso) <- c("x", "y", "velocity_length", "group", "timepoint")
mean(ggdata_meso$velocity_length)
```

    ## [1] 17.18124

``` r
mean(ggdata_meso$velocity_length[cellrank_out$obs$cell_subtype == "Prenatal Wnt2+"])
```

    ## [1] 15.69845

``` r
ggplot(ggdata_meso %>% arrange(velocity_length), aes(x = x, y = y, color = velocity_length)) + geom_point() + umap_theme() + theme(aspect.ratio = 1) + scale_color_viridis()
```

![](./velocity/meso/meso_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
to_term <- cellrank_out$obsm$to_terminal_states

colnames(to_term) <- cellrank_out$uns$to_terminal_states_names
head(to_term)
```

    ##       Pericyte_1 Myofibroblast_1 Prenatal Myofibroblast_1
    ## [1,] 0.011270299     0.050363731               0.09519438
    ## [2,] 0.012406293     0.039895762               0.14178947
    ## [3,] 0.006151310     0.039782381               0.14835289
    ## [4,] 0.017688361     0.049250249               0.09892252
    ## [5,] 0.005146711     0.039248168               0.07971109
    ## [6,] 0.000900429     0.006396649               0.88388480
    ##      Prenatal Myofibroblast_2   Pericyte_2    Wnt2+_1    Wnt2+_2     Wnt2+_3
    ## [1,]              0.011841697 7.005173e-04 0.53742788 0.09163262 0.037682306
    ## [2,]              0.008290472 8.133982e-04 0.51268626 0.08579478 0.052109291
    ## [3,]              0.008668414 3.987316e-04 0.51850791 0.08155978 0.046911452
    ## [4,]              0.007532461 1.157461e-03 0.53278246 0.09633190 0.041214608
    ## [5,]              0.006935428 3.341290e-04 0.54172863 0.08685698 0.044982369
    ## [6,]              0.007060686 5.885997e-05 0.06581638 0.01120488 0.005238058
    ##      Myofibroblast_2 Adventitial fibroblast
    ## [1,]     0.012023512             0.15186305
    ## [2,]     0.009232252             0.13698201
    ## [3,]     0.008763507             0.14090362
    ## [4,]     0.012641185             0.14247879
    ## [5,]     0.007442288             0.18761421
    ## [6,]     0.002381798             0.01705747

``` r
# De-duplicate. Other 'lineages' may be of interset later, but aren't a good broad overview.
to_term_dedup <- to_term[,c(10)]

to_term <- cbind(to_term_dedup,
                 rowSums(to_term[,c(1, 5)]), # Pericyte
                 rowSums(to_term[,c(2, 9)]), # Myofibroblast
                 rowSums(to_term[,c(3, 4)]), # Prenatal myofibroblast
                 rowSums(to_term[,c(6, 7, 8)]) # Wnt2+
)
colnames(to_term) <- c("Adventitial fibroblast",
                       "Pericyte",
                       "Myofibroblast",
                       "Proliferating Myofibroblast",
                       "*Wnt2*+ FB"
)

to_term <- as.data.frame(to_term)
to_term$timepoint <- ordered(as.factor(cellrank_out$obs$timepoint), c("E12", "E15", "E18", "P0", "P3", "P5", "P7", "P14"))
```

``` r
summary_stats <- melt(to_term[cellrank_out$obs$cell_subtype == "Prenatal Wnt2+",]) %>% group_by(variable) %>% summarize(mean = mean(value) * 100,
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
aspect_ratio <- 1.4 # other charts are 1.6


to_term_long <- melt(to_term[cellrank_out$obs$cell_subtype == "Prenatal Wnt2+",])
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


transitional_jitter_and_bar <- ggplot() +
        geom_col(data = summary_stats, aes(y = median, x = variable), color = "black", fill = "white") +
        #geom_col(data = summary_stats, aes(x = variable, y = early), fill = "#CDCDCD", width = 0.85) + # color early here, overlay on the other line
        geom_errorbar(data = summary_stats, aes(y = median, x = variable, ymin = q1, ymax = q2), width = 0.2) +
        #geom_boxplot(data = to_term_long, aes(y = value, x = variable), outlier.shape = NA, lwd=0.25) +
        #ggbeeswarm::geom_quasirandom(data = to_term_long, aes(y = value, x = variable), size = 0.8, shape = 1, alpha = 0.6) +
        ggtitle("Fate of<br>proliferating *Wnt2*+ FB") +
        prism_theme() +
        theme(aspect.ratio = aspect_ratio,
              plot.title = element_markdown(size = 14),
              axis.text.x = element_markdown(size = 12)
        ) +
        ylab("Transition probability (%)")
transitional_jitter_and_bar
```

![](./velocity/meso/meso_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
saveTiff("./figures/meso/fate_prenatal_wnt2.tiff",
         transitional_jitter_and_bar,
         width = 4, height = 5)
```

    ## Warning in normalizePath(paste0(path, "_tmp.tiff")): path[1]="./figures/meso/
    ## fate_prenatal_wnt2.tiff_tmp.tiff": No such file or directory

    ## Warning in normalizePath(path): path[1]="./figures/meso/
    ## fate_prenatal_wnt2.tiff": No such file or directory

    ## [1] TRUE

``` r
summary_stats <- melt(to_term[cellrank_out$obs$cell_subtype == "Prenatal Myofibroblast",]) %>% group_by(variable) %>% summarize(mean = mean(value) * 100,
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
aspect_ratio <- 1.4 # other charts are 1.6

to_term_long <- melt(to_term[cellrank_out$obs$cell_subtype == "Prenatal Myofibroblast",])
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


prenatal_myo_jitter_and_bar <- ggplot() +
        geom_col(data = summary_stats, aes(y = median, x = variable), color = "black", fill = "white") +
        #geom_col(data = summary_stats, aes(x = variable, y = early), fill = "#CDCDCD", width = 0.85) + # color early here, overlay on the other line
        geom_errorbar(data = summary_stats, aes(y = median, x = variable, ymin = q1, ymax = q2), width = 0.2) +
        #geom_boxplot(data = to_term_long, aes(y = value, x = variable), outlier.shape = NA, lwd=0.25) +
        #ggbeeswarm::geom_quasirandom(data = to_term_long, aes(y = value, x = variable), size = 0.8, shape = 1, alpha = 0.6) +
        ggtitle("Fate of<br>proliferating myofibroblasts") +
        prism_theme() +
  theme(aspect.ratio = aspect_ratio,
        plot.title = element_markdown(size = 14),
        axis.text.x = element_markdown(size = 12)
  ) +
        ylab("Transition probability (%)")
prenatal_myo_jitter_and_bar
```

![](./velocity/meso/meso_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
saveTiff("./figures/meso/fate_prenatal_myo.tiff",
         transitional_jitter_and_bar,
         width = 4, height = 5)
```

    ## Warning in normalizePath(paste0(path, "_tmp.tiff")): path[1]="./figures/meso/
    ## fate_prenatal_myo.tiff_tmp.tiff": No such file or directory

    ## Warning in normalizePath(path): path[1]="./figures/meso/fate_prenatal_myo.tiff":
    ## No such file or directory

    ## [1] TRUE

``` r
from_initial <- cellrank_out$obsm$from_initial_states #absorbtion from X
colnames(from_initial) <- cellrank_out$uns$from_initial_states_names
head(from_initial)
```

    ##      Pericyte_1 Prenatal Myofibroblast_1 Smooth muscle Prenatal Myofibroblast_2
    ## [1,] 0.02158014              0.015121016   0.015637004              0.009181764
    ## [2,] 0.01503855              0.010311835   0.010543063              0.006078948
    ## [3,] 0.02115275              0.014782641   0.015246841              0.008866057
    ## [4,] 0.01338573              0.009062634   0.009282653              0.005325300
    ## [5,] 0.02224369              0.015627428   0.017390107              0.009543524
    ## [6,] 0.02071325              0.015509465   0.015286890              0.009048846
    ##      Mesothelium Pericyte_2       Wnt2+ Prenatal Wnt2+ Myofibroblast
    ## [1,]  0.11120129 0.02611937 0.016387391      0.7733464   0.011425650
    ## [2,]  0.05457693 0.01964872 0.011049239      0.8647177   0.008034985
    ## [3,]  0.10443003 0.02567313 0.015967312      0.7827365   0.011144708
    ## [4,]  0.04498640 0.01785182 0.009636186      0.8833810   0.007088271
    ## [5,]  0.14097661 0.02762748 0.018278121      0.7362658   0.012047196
    ## [6,]  0.07794668 0.02576398 0.015129980      0.8091149   0.011485997

``` r
# De-duplicate. Other 'lineages' may be of interset later, but aren't a good broad overview.
from_initial_dedup <- from_initial[,c(3, 5, 7, 8, 9)]

from_initial <- cbind(from_initial_dedup,
                 rowSums(from_initial[,c(2, 4)]), # Prenatal Myofibroblast
                 rowSums(from_initial[,c(1, 6)]) # Pericyte

)
colnames(from_initial) <- c("Smooth muscle",
                            "Mesothelium",
                            "*Wnt2*+ FB",
                            "Proliferating *Wnt2*+ FB",
                            "Myofibroblast",
                            "Proliferating myofibroblast",
                            "Pericyte"
)
```

``` r
from_initial <- as.data.frame(from_initial)
from_initial$timepoint <- ordered(as.factor(cellrank_out$obs$timepoint), c("E12", "E15", "E18", "P0", "P3", "P5", "P7", "P14"))

head(from_initial[cellrank_out$obs$cell_subtype == "Proliferating Wnt2+",])
```

    ## [1] Smooth muscle               Mesothelium                
    ## [3] *Wnt2*+ FB                  Proliferating *Wnt2*+ FB   
    ## [5] Myofibroblast               Proliferating myofibroblast
    ## [7] Pericyte                    timepoint                  
    ## <0 rows> (or 0-length row.names)

``` r
to_wnt2_plot <- ggplot(melt(from_initial[cellrank_out$obs$cell_subtype == "Prenatal Wnt2+",]), aes(x = value, color = variable, fill = variable)) +
        geom_density(alpha = 0.3) +
        scale_x_continuous(limits = c(0.05,1)) +
        ggtitle("Fate probability towards 'Proliferating Wnt2+' cells")
```

    ## Using timepoint as id variables

``` r
to_wnt2_plot
```

    ## Warning: Removed 31217 rows containing non-finite values (stat_density).

![](./velocity/meso/meso_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
summary_stats <- melt(from_initial[cellrank_out$obs$cell_subtype == "Wnt2+",]) %>% group_by(variable) %>% summarize(mean = mean(value) * 100,
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
to_wnt2_long <- melt(from_initial[cellrank_out$obs$cell_subtype == "Wnt2+",])
```

    ## Using timepoint as id variables

``` r
to_wnt2_long$value <- to_wnt2_long$value * 100
to_wnt2_long <- to_wnt2_long[to_wnt2_long$value > 1,]


transitional_barplot <- ggplot() +
  geom_col(data = summary_stats, aes(y = median, x = variable), color = "black", fill = "white") +

  #geom_col(data = summary_stats, aes(x = variable, y = early), fill = "#CDCDCD", width = 0.84) + # color early here, overlay on the other line
  geom_errorbar(data = summary_stats, aes(y = median, x = variable, ymin = q1, ymax = q2), width = 0.6) +
  #geom_boxplot(data = from_initial_long, aes(y = value, x = variable), outlier.shape = NA, lwd=0.25) +
  #ggbeeswarm::geom_quasirandom(data = from_initial_long, aes(y = value, x = variable), size = 0.8, shape = 1, alpha = 0.6) +
  ggtitle("Source of *Wnt2*+ FB") +
  prism_theme() +
  theme(aspect.ratio = 1.6,
        plot.title = element_markdown(size = 14),
        axis.text.x = element_markdown(size = 14)
  ) +
  ylab("Transition probability (%)")

transitional_barplot
```

![](./velocity/meso/meso_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
saveTiff("./figures/meso/to_wnt2_probability.tiff",
         transitional_barplot,
         width = 3, height = 5)
```

    ## Warning in normalizePath(paste0(path, "_tmp.tiff")): path[1]="./figures/meso/
    ## to_wnt2_probability.tiff_tmp.tiff": No such file or directory

    ## Warning in normalizePath(path): path[1]="./figures/meso/
    ## to_wnt2_probability.tiff": No such file or directory

    ## [1] TRUE

``` r
summary_stats <- melt(from_initial[cellrank_out$obs$cell_subtype == "Myofibroblast",]) %>% group_by(variable) %>% summarize(mean = mean(value) * 100,
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
to_wnt2_long <- melt(from_initial[cellrank_out$obs$cell_subtype == "Myofibroblast",])
```

    ## Using timepoint as id variables

``` r
to_wnt2_long$value <- to_wnt2_long$value * 100
to_wnt2_long <- to_wnt2_long[to_wnt2_long$value > 1,]


transitional_barplot <- ggplot() +
  geom_col(data = summary_stats, aes(y = median, x = variable), color = "black", fill = "white") +

  #geom_col(data = summary_stats, aes(x = variable, y = early), fill = "#CDCDCD", width = 0.84) + # color early here, overlay on the other line
  geom_errorbar(data = summary_stats, aes(y = median, x = variable, ymin = q1, ymax = q2), width = 0.6) +
  #geom_boxplot(data = from_initial_long, aes(y = value, x = variable), outlier.shape = NA, lwd=0.25) +
  #ggbeeswarm::geom_quasirandom(data = from_initial_long, aes(y = value, x = variable), size = 0.8, shape = 1, alpha = 0.6) +
  ggtitle("Source of myofibroblast") +
  prism_theme() +
  theme(aspect.ratio = 1.6,
        plot.title = element_markdown(size = 14),
        axis.text.x = element_markdown(size = 14)
  ) +
  ylab("Transition probability (%)")

transitional_barplot
```

![](./velocity/meso/meso_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
saveTiff("./figures/meso/to_myo_probability.tiff",
         transitional_barplot,
         width = 3, height = 5)
```

    ## Warning in normalizePath(paste0(path, "_tmp.tiff")): path[1]="./figures/meso/
    ## to_myo_probability.tiff_tmp.tiff": No such file or directory

    ## Warning in normalizePath(path): path[1]="./figures/meso/
    ## to_myo_probability.tiff": No such file or directory

    ## [1] TRUE
