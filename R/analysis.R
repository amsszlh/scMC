#' Stacked bar plot showing the proportion of cells across certrain cell groups
#'
#' @param object Seurat object
#' @param x Name of one metadata column to show on the x-axis
#' @param fill Name of one metadata column to compare the proportion of cells
#' @param facet Name of one metadata column defining faceting groups
#' @param colors.use defining the color of stacked bar plot; either a char vector defining a color for each cell group or a palette name from brewer.pal
#' @param n.colors Number of colors when setting colors.use to be a palette name from brewer.pal
#' @param n.row Number of rows in facet_grid()
#' @param title.name Name of the main title
#' @param legend.title Name of legend
#' @param xlabel Name of x label
#' @param ylabel Name of y label
#' @param width bar width
#' @param show.legend Whether show the legend
#' @param x.lab.rot Whether rorate the xtick labels
#' @param text.size font size
#' @param flip Whether flip the cartesian coordinates so that horizontal becomes vertical
#' @return ggplot2 object
#' @export
#'
#' @import ggplot2
#' @importFrom plyr ddply as.quoted
computeProportion <- function(object, x = "cellType", fill = "condition", facet = NULL, colors.use = NULL, n.colors = 8, n.row = 1, title.name = NULL, legend.title = NULL,
                              xlabel = NULL, ylabel = "Cellular composition (%)", width = 0.8,
                              show.legend = TRUE, x.lab.rot = TRUE, text.size = 10, flip = FALSE) {

  df <- plyr::ddply(object@meta.data, plyr::as.quoted(c(x,fill)), nrow)

  gg <- ggplot(df, aes_string(x = x, y = "V1", fill = fill)) +
    geom_bar(stat="identity", position="fill", width = width) +
    scale_y_continuous(name = ylabel, labels = c(0,25,50,75,100))
  if (!is.null(facet)) {
    gg <- gg + facet_wrap(facet, nrow = n.row)
  }

  gg <- gg + theme_classic() + ylab(ylabel) + xlab(xlabel) +
    labs(title = title.name) +  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(text = element_text(size = text.size), axis.text = element_text(colour="black"))

  if (!is.null(colors.use)) {
    if (length(colors.use) == 1) {
      colors.use <- RColorBrewer::brewer.pal(n.colors, colors.use)[1:length(unique(df[, fill]))]
    }
    gg <- gg + scale_fill_manual(values = alpha(colors.use, alpha = 1), drop = FALSE)
    #   gg <- gg + scale_color_manual(values = alpha(color.use, alpha = 1), drop = FALSE) + guides(colour = FALSE)
  }
  if (is.null(legend.title)) {
    gg <- gg + theme(legend.title = element_blank())
  } else {
    gg <- gg + guides(fill=guide_legend(legend.title))
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    gg <- gg + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=text.size))
  }
  if (flip) {
    gg <- gg + coord_flip()
  }
  gg
  return(gg)
}
