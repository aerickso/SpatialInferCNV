#' Plotting the spatial distribution of genes with an inferred copy number alteration from an underlying matrix
#'
#' Plot_PGA_Visualization_Matrix()

Plot_PGA_Visualization_Matrix <- function(SectionName, InputMatrix, MaxValInput) {
  ggplot(InputMatrix, aes(x = x, y = y)) +
    geom_raster(aes(fill=PercentageGenomeAltered)) +
    scale_fill_gradient(limits = c(0, MaxValInput), low="blue", high="yellow", na.value = "grey50") +
    labs(x="X-coord", y="Y-coord") +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11)) +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          line = element_blank(),
          title = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          legend.position = "none",
          plot.margin=grid::unit(c(0,0,0,0), "mm"))

  ggsave(paste0(SectionName, "_", "PGA_SpatialVisualization_", Sys.Date(), ".png"),
         plot = last_plot() + labs(x=NULL, y=NULL),
         device = NULL,
         path = NULL,
         scale = 1,
         width = NA,
         height = NA,
         dpi = 300,
         limitsize = TRUE)
}
