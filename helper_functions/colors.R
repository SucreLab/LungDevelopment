# Title     : A few vectors with improved color palletes, and some themes for pretty plots

color_scanpy_default <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
color_scanpy_viridis20 <- c("#023fa5", "#7d87b9", "#bec1d4", "#d6bcc0", "#bb7784", "#8e063b", "#4a6fe3",
                             "#8595e1", "#b5bbe3", "#e6afb9", "#e07b91", "#d33f6a", "#11c638", "#8dd593",
                             "#c6dec7", "#ead3c6", "#f0b98d", "#ef9708", "#0fcfc0", "#9cded6", "#d5eae7",
                             "#f3e1eb", "#f6c4e1", "#f79cd4",
                             "#7f7f7f", "#c7c7c7", "#1CE6FF", "#336600")
color_category_20 <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728", "#ff9896", "#9467bd",
                       "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5")


umap_theme <- function(){
  theme_grey() %+replace%
    theme(panel.background = element_rect(fill = "white", colour = "black", size = 1),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key=element_blank())
}


featureplot_theme <- function(){
  theme_grey() %+replace%
    theme(panel.background = element_rect(fill = "white", colour = "black", size = 1),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}


prism_theme <- function(){
  theme_grey() %+replace%
    theme(panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(),
          panel.grid.minor = element_blank(),
          legend.key=element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 0.95),
          axis.text.y = element_text(size = 14),
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(0.25, "cm"),
          axis.title.y = element_text(size = 16, vjust = 4, angle = 90),
          axis.title.x = element_blank(),
          plot.margin = margin(t = 5, r = 5, b = 5, l = 20),
          legend.text = element_text(size = 14),
          legend.key.size = unit(2,"line"),
          legend.title = element_blank())
}