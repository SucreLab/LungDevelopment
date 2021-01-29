Plottings outputs of the cellrank analysis
================
Nick Negretti
12/07/20

# Plot cellrank outputs of the lung epithelium

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
cellrank_out <- anndata::read_h5ad("./velocity/epi/only_metadata.h5ad")
cellrank_out
```

    ## AnnData object with n_obs × n_vars = 11472 × 19455
    ##     obs: 'timepoint', 'bulk_cellype', 'cell_subtype', 'sample_batch', 'Clusters', '_X', '_Y', 'initial_size_spliced', 'initial_size_unspliced', 'initial_size', 'n_counts', 'velocity_self_transition', 'velocity_length', 'velocity_confidence', 'velocity_confidence_transition', 'clusters_gradients', 'terminal_states', 'terminal_states_probs', 'to_terminal_states_dp', 'initial_states', 'initial_states_probs', 'from_initial_states_dp', 'velocity_pseudotime', 'latent_time', 'dpt_pseudotime', 'S_score', 'G2M_score', 'phase', 'root_cells', 'end_points'
    ##     var: 'fit_alpha', 'fit_beta', 'fit_gamma', 'fit_t_', 'fit_scaling', 'fit_alignment_scaling', 'fit_std_u', 'fit_std_s', 'fit_likelihood', 'fit_u0', 'fit_s0', 'fit_pval_steady', 'fit_steady_u', 'fit_steady_s', 'velocity_genes', 'fit_r2', 'fit_diff_kinetics', 'fit_pval_kinetics', 'spearmans_score', 'velocity_score'
    ##     uns: 'from_initial_states_names', 'initial_states_names', 'to_terminal_states_names'
    ##     obsm: 'X_diffmap', 'X_pca', 'X_umap', 'from_initial_states', 'macrostates_bwd', 'macrostates_fwd', 'to_terminal_states', 'velocity_umap'
    ##     varm: 'PCs', 'fit_pvals_kinetics', 'loss'

``` r
cellrank_out$uns$initial_states_names
```

    ## [1] "Secretory"  "Primordial" "Cilliated"  "Type II_1"  "Type I"    
    ## [6] "Type II_2"  "Type II_3"

``` r
ggdata_epi <- data.frame(x = cellrank_out$obsm$X_umap[,1],
                     y = cellrank_out$obsm$X_umap[,2],
                     velocity_length = cellrank_out$obs$velocity_length,
                     group = as.character(cellrank_out$obs$cell_subtype),
                     timepoint = as.character(cellrank_out$obs$timepoint))
head(ggdata_epi)
```

    ##            x          y velocity_length      group timepoint
    ## 1 -0.9757824 -1.9436623           15.58 Primordial       E12
    ## 2 -2.7454511 -0.4228817           10.24 Primordial       E12
    ## 3 -1.4709311 -2.4200785           14.36 Primordial       E12
    ## 4 -2.4639739  0.6939415           12.16 Primordial       E12
    ## 5 -2.8190245  0.2094124            9.35 Primordial       E12
    ## 6 -2.2454604 -1.1322407           13.84 Primordial       E12

``` r
colnames(ggdata_epi) <- c("x", "y", "velocity_length", "group", "timepoint")
mean(ggdata_epi$velocity_length)
```

    ## [1] 14.63572

``` r
ggplot(ggdata_epi %>% arrange(velocity_length), aes(x = x, y = y, color = velocity_length)) + geom_point() + umap_theme() + theme(aspect.ratio = 1) + scale_color_viridis()
```

![](./velocity/epi/epi_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
epi_ms <- Matrix::readMM("./velocity/epi/epi_expression_Ms.mtx.gz")
```

``` r
gene <- "Hopx"
include_celltypes <- c("AT1", "AT2", "Primordial", "Transitional")
exp_over_time <- data.frame(expression = epi_ms[,colnames(cellrank_out) == gene],
                            Cluster = cellrank_out$obs$cell_subtype,
                            latent_time = cellrank_out$obs$latent_time
)

unique(exp_over_time$Cluster)
```

    ## [1] Primordial     Type I         Type II        Cilliated      Transitional  
    ## [6] Secretory      Neuroendocrine
    ## 7 Levels: Cilliated Neuroendocrine Primordial Secretory ... Type II

``` r
levels(exp_over_time$Cluster)[levels(exp_over_time$Cluster) == "Type I"] <- "AT1"
levels(exp_over_time$Cluster)[levels(exp_over_time$Cluster) == "Type II"] <- "AT2"



exp_over_time$Cluster <- ordered(as.factor(exp_over_time$Cluster), c("Primordial", "Transitional", "AT1", "AT2", "Cilliated", "Secretory", "Neuroendocrine"))

exp_over_time <- exp_over_time[exp_over_time$Cluster %in% include_celltypes,]

hopx_exp <- ggplot(exp_over_time, aes(y = expression, x = latent_time, color = Cluster)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(aesthetics = c("color", "fill"), values= color_category_20, guide = guide_legend(override.aes = list(size=3), order = 1)) +
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

hopx_exp
```

![](./velocity/epi/epi_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
gene <- "Sftpc"
include_celltypes <- c("AT1", "AT2", "Primordial", "Transitional")
exp_over_time <- data.frame(expression = epi_ms[,colnames(cellrank_out) == gene],
                            Cluster = cellrank_out$obs$cell_subtype,
                            latent_time = cellrank_out$obs$latent_time
)

unique(exp_over_time$Cluster)
```

    ## [1] Primordial     Type I         Type II        Cilliated      Transitional  
    ## [6] Secretory      Neuroendocrine
    ## 7 Levels: Cilliated Neuroendocrine Primordial Secretory ... Type II

``` r
levels(exp_over_time$Cluster)[levels(exp_over_time$Cluster) == "Type I"] <- "AT1"
levels(exp_over_time$Cluster)[levels(exp_over_time$Cluster) == "Type II"] <- "AT2"



exp_over_time$Cluster <- ordered(as.factor(exp_over_time$Cluster), c("Primordial", "Transitional", "AT1", "AT2", "Cilliated", "Secretory", "Neuroendocrine"))

exp_over_time <- exp_over_time[exp_over_time$Cluster %in% include_celltypes,]

sftpc_exp <- ggplot(exp_over_time, aes(y = expression, x = latent_time, color = Cluster)) +
        geom_point(alpha = 0.5) +
        scale_color_manual(aesthetics = c("color", "fill"), values= color_category_20, guide = guide_legend(override.aes = list(size=3), order = 1)) +
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

sftpc_exp
```

![](./velocity/epi/epi_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
gene <- "Cdkn1a"
include_celltypes <- c("AT1", "AT2", "Primordial", "Transitional")
exp_over_time <- data.frame(expression = epi_ms[,colnames(cellrank_out) == gene],
                            Cluster = cellrank_out$obs$cell_subtype,
                            latent_time = cellrank_out$obs$latent_time
)

unique(exp_over_time$Cluster)
```

    ## [1] Primordial     Type I         Type II        Cilliated      Transitional  
    ## [6] Secretory      Neuroendocrine
    ## 7 Levels: Cilliated Neuroendocrine Primordial Secretory ... Type II

``` r
levels(exp_over_time$Cluster)[levels(exp_over_time$Cluster) == "Type I"] <- "AT1"
levels(exp_over_time$Cluster)[levels(exp_over_time$Cluster) == "Type II"] <- "AT2"



exp_over_time$Cluster <- ordered(as.factor(exp_over_time$Cluster), c("Primordial", "Transitional", "AT1", "AT2", "Cilliated", "Secretory", "Neuroendocrine"))

exp_over_time <- exp_over_time[exp_over_time$Cluster %in% include_celltypes,]

cdkn1a_exp <- ggplot(exp_over_time, aes(y = expression, x = latent_time, color = Cluster)) +
        geom_point(alpha = 0.5) +
        scale_color_manual(aesthetics = c("color", "fill"), values= color_category_20, guide = guide_legend(override.aes = list(size=3), order = 1)) +
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

cdkn1a_exp
```

![](./velocity/epi/epi_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
saveTiff("./figures/epi/latent_time/sftpc_hopx.tiff",
         hopx_exp + sftpc_exp +  plot_layout(ncol = 2, nrow = 1, guides = "collect"),
         width = 13, height = 4.5)
```

    ## [1] TRUE

``` r
saveTiff("./figures/epi/latent_time/cdkn1a.tiff",
         cdkn1a_exp,
         width = 6, height = 4.5)
```

    ## [1] TRUE

``` r
gene <- "Fbln5"
include_celltypes <- c("AT1", "AT2", "Primordial", "Transitional")
exp_over_time <- data.frame(expression = epi_ms[,colnames(cellrank_out) == gene],
                            Cluster = cellrank_out$obs$cell_subtype,
                            latent_time = cellrank_out$obs$latent_time
)

unique(exp_over_time$Cluster)
```

    ## [1] Primordial     Type I         Type II        Cilliated      Transitional  
    ## [6] Secretory      Neuroendocrine
    ## 7 Levels: Cilliated Neuroendocrine Primordial Secretory ... Type II

``` r
levels(exp_over_time$Cluster)[levels(exp_over_time$Cluster) == "Type I"] <- "AT1"
levels(exp_over_time$Cluster)[levels(exp_over_time$Cluster) == "Type II"] <- "AT2"



exp_over_time$Cluster <- ordered(as.factor(exp_over_time$Cluster), c("Primordial", "Transitional", "AT1", "AT2", "Cilliated", "Secretory", "Neuroendocrine"))

exp_over_time <- exp_over_time[exp_over_time$Cluster %in% include_celltypes,]

fbln5_exp <- ggplot(exp_over_time, aes(y = expression, x = latent_time, color = Cluster)) +
        geom_point(alpha = 0.5) +
        scale_color_manual(values= color_category_20, guide = guide_legend(override.aes = list(size=3), order = 1)) +
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

fbln5_exp
```

![](./velocity/epi/epi_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
saveTiff("./figures/epi/latent_time/fbln5.tiff",
         fbln5_exp,
         width = 6, height = 4.5)
```

    ## [1] TRUE

``` r
to_term_nomerge <- cellrank_out$obsm$to_terminal_states

colnames(to_term_nomerge) <- cellrank_out$uns$to_terminal_states_names


# Combine similar terminal types
to_term <- cbind(to_term_nomerge[,c(1,2)], rowSums(to_term_nomerge[,c(3,4)]), rowSums(to_term_nomerge[,c(5,6,8)]), rowSums(to_term_nomerge[,c(7,9)]))
colnames(to_term) <- c("Neuroendocrine", "Secretory", "Cilliated", "Type II", "Type I")
to_term <- as.data.frame(to_term)
to_term$timepoint <- ordered(as.factor(cellrank_out$obs$timepoint), c("E12", "E15", "E18", "P0", "P3", "P5", "P7", "P14"))
```

``` r
colnames(to_term)[4:5] <- c("AT2", "AT1")

summary_stats <- melt(to_term[cellrank_out$obs$cell_subtype == "Transitional",]) %>% group_by(variable) %>% summarize(mean = mean(value) * 100,
                                                                                                                      median = median(value) * 100,
                                                                                                                      sd = sd(value) * 100,
                                                                                                                      sem = (sd(value) / sqrt(length(value))) * 100,
                                                                                                                      q1 = quantile(value, 1/4) * 100,
                                                                                                                      q2 = quantile(value, 3/4) * 100,
                                                                                                                      ci1 = get_lower_ci(value) * 100,
                                                                                                                      ci2 = get_upper_ci(value) * 100)
```

    ## Using timepoint as id variables

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
to_term_long <- melt(to_term[cellrank_out$obs$cell_subtype == "Transitional",])
```

    ## Using timepoint as id variables

``` r
to_term_long$value <- to_term_long$value * 100
to_term_long_above_zero <- to_term_long[to_term_long$value > 1,]

transitional_jitter <- ggplot(to_term_long, aes(y = value, x = variable)) +
        ggbeeswarm::geom_quasirandom(size = 0.8, shape = 1, alpha = 0.6, bandwidth = 1) +
        ggtitle("Fate of 'transitional' cells") +
        prism_theme() +
        theme(aspect.ratio = 1.6) +
        ylab("Transition probability (%)")

transitional_violin <- ggplot(to_term_long_above_zero, aes(y = value, x = variable)) +
        geom_violin() +
        ggtitle("Fate of 'transitional' cells") +
        prism_theme() +
        theme(aspect.ratio = 1.6) +
        ylab("Transition probability (%)")

transitional_barplot <- ggplot(summary_stats, aes(y = median, x = variable)) +
  geom_col() +
  geom_errorbar(aes(ymin = q1, ymax = q2), width = 0.2) +
  ggtitle("Fate of 'transitional' cells") +
  prism_theme() +
  theme(aspect.ratio = 1.6) +
  ylab("Transition probability (%)")

get_early_pct <- function(data, summary_stats, cell_group){
  timepoints <- data$obs$timepoint[data$obs$cell_subtype == cell_group]
  early <- sum(timepoints %in% c("E12", "E15", "E18")) * 100 / length(timepoints)
  return(early * summary_stats[summary_stats$variable == cell_group, "median"] / 100) # scale vals to max bar
}

for (type in unique(summary_stats$variable)){
  summary_stats[summary_stats$variable == type, "early"] <- get_early_pct(cellrank_out, summary_stats, type)
}

summary_stats$variable <- ordered(summary_stats$variable, levels(summary_stats$variable)[order(as.character(levels(summary_stats$variable)))])
transitional_jitter_and_bar <- ggplot() +
        geom_col(data = summary_stats, aes(y = median, x = variable), color = "black", fill = "white") +
        geom_col(data = summary_stats, aes(x = variable, y = early), fill = "#CDCDCD", width = 0.85) + # color early here, overlay on the other line
        geom_errorbar(data = summary_stats, aes(y = median, x = variable, ymin = q1, ymax = q2), width = 0.6) +
        #geom_boxplot(data = to_term_long, aes(y = value, x = variable), outlier.shape = NA, lwd=0.25) +
        #ggbeeswarm::geom_quasirandom(data = to_term_long, aes(y = value, x = variable), size = 0.8, shape = 1, alpha = 0.6) +
        ggtitle("Fate of transitional cells") +
        prism_theme() +
        theme(aspect.ratio = 1.6) +
        ylab("Transition probability (%)")

transitional_jitter_and_bar
```

    ## Warning: Removed 2 rows containing missing values (position_stack).

![](./velocity/epi/epi_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
transitional_violin + transitional_jitter + transitional_barplot + transitional_jitter_and_bar + plot_layout(nrow = 1)
```

    ## Warning: Removed 2 rows containing missing values (position_stack).

![](./velocity/epi/epi_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

``` r
saveTiff("./figures/epi/from_transitional_probability.tiff",
         transitional_jitter_and_bar,
         width = 3.5, height = 5.3)
```

    ## Warning: Removed 2 rows containing missing values (position_stack).

    ## [1] TRUE

``` r
summary_stats <- melt(to_term[cellrank_out$obs$cell_subtype == "Type II",]) %>% group_by(variable) %>% summarize(mean = mean(value) * 100,
                                                                                                                      median = median(value) * 100,
                                                                                                                      sd = sd(value) * 100,
                                                                                                                      sem = (sd(value) / sqrt(length(value))) * 100,
                                                                                                                      q1 = quantile(value, 1/4) * 100,
                                                                                                                      q2 = quantile(value, 3/4) * 100
                                                                                                                      #ci1 = get_lower_ci(value) * 100,
                                                                                                                      #ci2 = get_upper_ci(value) * 100
                                                                                                                      )
```

    ## Using timepoint as id variables

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
get_early_pct <- function(data, summary_stats, cell_group){
  timepoints <- data$obs$timepoint[data$obs$cell_subtype == cell_group]
  early <- sum(timepoints %in% c("E12", "E15", "E18")) * 100 / length(timepoints)
  return(early * summary_stats[summary_stats$variable == cell_group, "median"] / 100) # scale vals to max bar
}

for (type in unique(summary_stats$variable)){
  summary_stats[summary_stats$variable == type, "early"] <- get_early_pct(cellrank_out, summary_stats, type)
}


to_term_long <- melt(to_term[cellrank_out$obs$cell_subtype == "Type II",])
```

    ## Using timepoint as id variables

``` r
to_term_long$value <- to_term_long$value * 100
to_term_long_above_zero <- to_term_long[to_term_long$value > 1,]

at2_1_transitional_jitter_and_bar <- ggplot() +
        geom_col(data = summary_stats, aes(y = median, x = variable), color = "black", fill = "white") +
        geom_col(data = summary_stats, aes(x = variable, y = early), fill = "#CDCDCD", width = 0.85) + # color early here, overlay on the other line
        geom_errorbar(data = summary_stats, aes(y = median, x = variable, ymin = q1, ymax = q2), width = 0.6) +
        #geom_boxplot(data = to_term_long, aes(y = value, x = variable), outlier.shape = NA, lwd=0.25) +
        #ggbeeswarm::geom_quasirandom(data = to_term_long, aes(y = value, x = variable), size = 0.8, shape = 1, alpha = 0.6) +
        ggtitle("Fate of AT2 cells") +
        prism_theme() +
        theme(aspect.ratio = 1.6,
              plot.title = element_markdown(size = 14),
              axis.text.x = element_markdown(size = 14)) +
        ylab("Transition probability (%)")

at2_1_transitional_jitter_and_bar
```

    ## Warning: Removed 2 rows containing missing values (position_stack).

![](./velocity/epi/epi_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
saveTiff("./figures/epi/fate_AT2.tiff",
         at2_1_transitional_jitter_and_bar,
         width = 3.5, height = 5.3)
```

    ## Warning: Removed 2 rows containing missing values (position_stack).

    ## [1] TRUE

``` r
from_initial <- cellrank_out$obsm$from_initial_states

colnames(from_initial) <- cellrank_out$uns$from_initial_states_names

head(from_initial)
```

    ##        Secretory Primordial    Cilliated    Type II_1       Type I    Type II_2
    ## [1,] 0.005713296  0.9839824 0.0004660468 0.0021542966 0.0012070385 0.0013901565
    ## [2,] 0.003746096  0.9912577 0.0003159508 0.0011424542 0.0001945432 0.0007141155
    ## [3,] 0.006911929  0.9838334 0.0008493751 0.0018939761 0.0008305552 0.0012079991
    ## [4,] 0.002101344  0.9949269 0.0001584361 0.0007014545 0.0001187391 0.0004170327
    ## [5,] 0.000000000  1.0000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000
    ## [6,] 0.008964042  0.9789582 0.0006712044 0.0027411185 0.0006612547 0.0017024564
    ##        Type II_3
    ## [1,] 0.005086777
    ## [2,] 0.002629140
    ## [3,] 0.004472747
    ## [4,] 0.001576047
    ## [5,] 0.000000000
    ## [6,] 0.006301773

``` r
# Combine similar terminal types

from_initial <- cbind(from_initial[,c(1:3, 5)], rowSums(from_initial[,c(4,6,7)])) # Merge the Type II fates
colnames(from_initial) <- c(colnames(from_initial)[c(1:4)], "Type II")
from_initial <- as.data.frame(from_initial)
from_initial$timepoint <- ordered(as.factor(cellrank_out$obs$timepoint), c("E12", "E15", "E18", "P0", "P3", "P5", "P7", "P14"))

head(from_initial[cellrank_out$obs$cell_subtype == "Transitional",])
```

    ##      Secretory Primordial   Cilliated      Type I    Type II timepoint
    ## 291 0.02989208  0.9342984 0.001150781 0.002160138 0.03249855       E15
    ## 302 0.06808330  0.7469612 0.002236005 0.009183589 0.17353590       E15
    ## 342 0.02264489  0.9464021 0.001109534 0.001273976 0.02856949       E15
    ## 343 0.06257360  0.8002389 0.002051083 0.006784564 0.12835187       E15
    ## 363 0.02608466  0.9405292 0.001103450 0.001360134 0.03092254       E15
    ## 389 0.07617952  0.7120763 0.002145189 0.012110885 0.19748813       E18

``` r
colnames(from_initial)[4:5] <- c("AT1", "AT2")

transitional_plot <- ggplot(melt(from_initial[cellrank_out$obs$cell_subtype == "Transitional",]), aes(x = value, color = variable, fill = variable)) +
  geom_density(alpha = 0.3) +
  scale_x_continuous(limits = c(0.05,1)) +
  ggtitle("Fate probability towards 'transitional' cells")
```

    ## Using timepoint as id variables

``` r
transitional_plot
```

    ## Warning: Removed 249 rows containing non-finite values (stat_density).

    ## Warning: Groups with fewer than two data points have been dropped.

    ## Warning in max(ids, na.rm = TRUE): no non-missing arguments to max; returning -
    ## Inf

![](./velocity/epi/epi_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
summary_stats <- melt(from_initial[cellrank_out$obs$cell_subtype == "Type II",]) %>% group_by(variable) %>% summarize(mean = mean(value) * 100,
                                                                                                                           median = median(value) * 100,
                                                                                                                           sd = sd(value) * 100,
                                                                                                                           sem = (sd(value) / sqrt(length(value))) * 100,
                                                                                                                           q1 = quantile(value, 1/4) * 100,
                                                                                                                           q2 = quantile(value, 3/4) * 100)
```

    ## Using timepoint as id variables

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
from_initial_long <- melt(from_initial[cellrank_out$obs$cell_subtype == "Type II",])
```

    ## Using timepoint as id variables

``` r
from_initial_long$value <- from_initial_long$value * 100
#from_initial_long_above_zero <- from_initial_long[to_term_long$value > 1,]


summary_stats$early <- 0

get_early_pct <- function(data, summary_stats, cell_group){
  timepoints <- data$obs$timepoint[data$obs$cell_subtype == cell_group]
  early <- sum(timepoints %in% c("E12", "E15", "E18")) * 100 / length(timepoints)
  return(early * summary_stats[summary_stats$variable == cell_group, "median"] / 100) # scale vals to max bar
}

for (type in unique(summary_stats$variable)){
  summary_stats[summary_stats$variable == type, "early"] <- get_early_pct(cellrank_out, summary_stats, type)
}

AT2_jitter_and_bar <- ggplot() +
        geom_col(data = summary_stats, aes(y = median, x = variable), color = "black", fill = "white") +
        geom_col(data = summary_stats, aes(x = variable, y = early), fill = "#CDCDCD", width = 0.85) + # color early here, overlay on the other line
        geom_errorbar(data = summary_stats, aes(y = median, x = variable, ymin = q1, ymax = q2), width = 0.6) +
        #geom_boxplot(data = from_initial_long, aes(y = value, x = variable), outlier.shape = NA, lwd=0.25) +
        #ggbeeswarm::geom_quasirandom(data = from_initial_long, aes(y = value, x = variable), size = 0.8, shape = 1, alpha = 0.6) +
        ggtitle("Source of AT2 cells") +
        prism_theme() +
        theme(aspect.ratio = 1.6,
              plot.title = element_markdown(size = 16),
              axis.text.x = element_markdown(size = 14)) +
        ylab("Transition probability (%)")

AT2_jitter_and_bar
```

    ## Warning: Removed 2 rows containing missing values (position_stack).

![](./velocity/epi/epi_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
saveTiff("./figures/epi/to_AT2.tiff",
         AT2_jitter_and_bar,
         width = 3.2, height = 5.2)
```

    ## Warning: Removed 2 rows containing missing values (position_stack).

    ## [1] TRUE

``` r
summary_stats <- melt(from_initial[cellrank_out$obs$cell_subtype == "Type I",]) %>% group_by(variable) %>% summarize(mean = mean(value) * 100,
                                                                                                                      median = median(value) * 100,
                                                                                                                      sd = sd(value) * 100,
                                                                                                                      sem = (sd(value) / sqrt(length(value))) * 100,
                                                                                                                      q1 = quantile(value, 1/4) * 100,
                                                                                                                      q2 = quantile(value, 3/4) * 100)
```

    ## Using timepoint as id variables

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
from_initial_long <- melt(from_initial[cellrank_out$obs$cell_subtype == "Type I",])
```

    ## Using timepoint as id variables

``` r
from_initial_long$value <- from_initial_long$value * 100
#from_initial_long_above_zero <- from_initial_long[to_term_long$value > 1,]


summary_stats$early <- 0

get_early_pct <- function(data, summary_stats, cell_group){
  timepoints <- data$obs$timepoint[data$obs$cell_subtype == cell_group]
  early <- sum(timepoints %in% c("E12", "E15", "E18")) * 100 / length(timepoints)
  return(early * summary_stats[summary_stats$variable == cell_group, "median"] / 100) # scale vals to max bar
}

for (type in unique(summary_stats$variable)){
  summary_stats[summary_stats$variable == type, "early"] <- get_early_pct(cellrank_out, summary_stats, type)
}

AT1_jitter_and_bar <- ggplot() +
        geom_col(data = summary_stats, aes(y = median, x = variable), color = "black", fill = "white") +
        geom_col(data = summary_stats, aes(x = variable, y = early), fill = "#CDCDCD", width = 0.85) + # color early here, overlay on the other line
        geom_errorbar(data = summary_stats, aes(y = median, x = variable, ymin = q1, ymax = q2), width = 0.6) +
        #geom_boxplot(data = from_initial_long, aes(y = value, x = variable), outlier.shape = NA, lwd=0.25) +
        #ggbeeswarm::geom_quasirandom(data = from_initial_long, aes(y = value, x = variable), size = 0.8, shape = 1, alpha = 0.6) +
        ggtitle("Source of AT1 cells") +
        prism_theme() +
        theme(aspect.ratio = 1.6,
              plot.title = element_markdown(size = 16),
              axis.text.x = element_markdown(size = 14)) +
        ylab("Transition probability (%)")

AT1_jitter_and_bar
```

    ## Warning: Removed 2 rows containing missing values (position_stack).

![](./velocity/epi/epi_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
saveTiff("./figures/epi/to_AT1.tiff",
         AT1_jitter_and_bar,
         width = 3.2, height = 5.2)
```

    ## Warning: Removed 2 rows containing missing values (position_stack).

    ## [1] TRUE

``` r
summary_stats <- melt(from_initial[cellrank_out$obs$cell_subtype == "Transitional",]) %>% group_by(variable) %>% summarize(mean = mean(value) * 100,
                                                                                                                     median = median(value) * 100,
                                                                                                                     sd = sd(value) * 100,
                                                                                                                     sem = (sd(value) / sqrt(length(value))) * 100,
                                                                                                                     q1 = quantile(value, 1/4) * 100,
                                                                                                                     q2 = quantile(value, 3/4) * 100)
```

    ## Using timepoint as id variables

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
from_initial_long <- melt(from_initial[cellrank_out$obs$cell_subtype == "Transitional",])
```

    ## Using timepoint as id variables

``` r
from_initial_long$value <- from_initial_long$value * 100
#from_initial_long_above_zero <- from_initial_long[to_term_long$value > 1,]


summary_stats$early <- 0

get_early_pct <- function(data, summary_stats, cell_group){
  timepoints <- data$obs$timepoint[data$obs$cell_subtype == cell_group]
  early <- sum(timepoints %in% c("E12", "E15", "E18")) * 100 / length(timepoints)
  return(early * summary_stats[summary_stats$variable == cell_group, "median"] / 100) # scale vals to max bar
}

for (type in unique(summary_stats$variable)){
  summary_stats[summary_stats$variable == type, "early"] <- get_early_pct(cellrank_out, summary_stats, type)
}

summary_stats$variable <- ordered(summary_stats$variable, levels(summary_stats$variable)[order(as.character(levels(summary_stats$variable)))])

to_transitional_jitter_and_bar <- ggplot() +
        geom_col(data = summary_stats, aes(y = median, x = variable), color = "black", fill = "white") +
        geom_col(data = summary_stats, aes(x = variable, y = early), fill = "#CDCDCD", width = 0.85) + # color early here, overlay on the other line
        geom_errorbar(data = summary_stats, aes(y = median, x = variable, ymin = q1, ymax = q2), width = 0.6) +
        #geom_boxplot(data = from_initial_long, aes(y = value, x = variable), outlier.shape = NA, lwd=0.25) +
        #ggbeeswarm::geom_quasirandom(data = from_initial_long, aes(y = value, x = variable), size = 0.8, shape = 1, alpha = 0.6) +
        ggtitle("Source of transitional cells") +
        prism_theme() +
        theme(aspect.ratio = 1.6) +
        ylab("Transition probability (%)")

to_transitional_jitter_and_bar
```

    ## Warning: Removed 2 rows containing missing values (position_stack).

![](./velocity/epi/epi_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
saveTiff("./figures/epi/to_transitional_probability.tiff",
         to_transitional_jitter_and_bar,
         width = 3.2, height = 5.2)
```

    ## Warning: Removed 2 rows containing missing values (position_stack).

    ## [1] TRUE

``` r
to_transitional_jitter_and_bar + transitional_jitter_and_bar + plot_layout(ncol = 2, nrow = 1) & theme(plot.margin = margin(5.5,105.5,5.5,5.5))
```

    ## Warning: Removed 2 rows containing missing values (position_stack).
    
    ## Warning: Removed 2 rows containing missing values (position_stack).

![](./velocity/epi/epi_plot_cellrank_outs_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
saveTiff("./figures/epi/main_transitions.pdf",
         to_transitional_jitter_and_bar + transitional_jitter_and_bar + plot_layout(ncol = 2, nrow = 1) & theme(plot.margin = margin(5.5,100.5,5.5,5.5)),
         height = 5, width = 10,)
```

    ## Warning: Removed 2 rows containing missing values (position_stack).
    
    ## Warning: Removed 2 rows containing missing values (position_stack).

    ## [1] TRUE
