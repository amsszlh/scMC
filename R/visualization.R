#' ggplot theme in scMC
#'
#' @export
#'
#' @examples
#' @importFrom ggplot2 theme_classic element_rect theme element_blank element_line element_text
scMC_theme_opts <- function() {
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
    theme_classic() +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(color = "black")) +
    theme(axis.line.y = element_line(color = "black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(legend.key = element_blank()) + theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
}


#' Generate ggplot2 colors
#'
#' @param n number of colors to generate
#' @importFrom grDevices hcl
#' @export
#'
ggPalette <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Generate colors from a customed color palette
#'
#' @param n number of colors
#'
#' @return A color palette for plotting
#' @importFrom grDevices colorRampPalette
#'
#' @export
#'
scPalette <- function(n) {
  colorSpace <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
                  '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69','#8DD3C7','#999999')
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  } else {
    colors <- grDevices::colorRampPalette(colorSpace)(n)
  }
  return(colors)
}


#' Stacked Violin plot
#'
#' @param object Seurat object
#' @param features Features to plot (gene expression, metrics)
#' @param colors.use defining the color for each cell group
#' @param colors.ggplot whether use ggplot color scheme; default: colors.ggplot = FALSE
#' @param split.by Name of a metadata column to split plot by;
#' @param idents Which classes to include in the plot (default is all)
#' @param show.text.y whther show y-axis text
#' @param line.size line width in the violin plot
#' @param pt.size size of the dots
#' @param plot.margin adjust the white space between each plot
#' @param angle.x angle for x-axis text rotation
#' @param vjust.x adjust x axis text
#' @param hjust.x adjust x axis text
#' @param ... Extra parameters passed to VlnPlot from Seurat package
#' @return ggplot2 object
#' @export
#' @import ggplot2
#' @importFrom  patchwork wrap_plots
#' @importFrom purrr map map2 map_dbl
#' @importFrom Seurat VlnPlot
#'
# #' This is modified from Ming Tang's post at https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
#'
StackedVlnPlot<- function(object, features, idents = NULL, split.by = NULL,
                          colors.use = NULL, colors.ggplot = FALSE,
                          angle.x = 90, vjust.x = NULL, hjust.x = NULL, show.text.y = TRUE, line.size = NULL,
                          pt.size = 0,
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          ...) {
  options(warn=-1)
  if (is.null(colors.use)) {
    numCluster <- length(levels(Idents(object)))
    if (colors.ggplot) {
      colors.use <- NULL
    } else {
      colors.use <- scPalette(numCluster)
    }
  }
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle=c(0, 45, 90)
    hjust=c(0, 0.5, 1)
    vjust=c(0, 0.5, 0.5)
    vjust.x = vjust[angle == angle.x]
    hjust.x = hjust[angle == angle.x]
  }

  plot_list<- purrr::map(features, function(x) modify_vlnplot(object = object, features = x, idents = idents, split.by = split.by, cols = colors.use, pt.size = pt.size,
                                                              show.text.y = show.text.y, line.size = line.size, ...))

  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line()) +
    theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x, vjust = vjust.x)) +
    theme(axis.text.x = element_text(size = 10))

  # change the y-axis tick to only max value
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x +
                            scale_y_continuous(breaks = c(y)) +
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

#' modified vlnplot
#' @param object Seurat object
#' @param features Features to plot (gene expression, metrics)
#' @param split.by Name of a metadata column to split plot by;
#' @param idents Which classes to include in the plot (default is all)
#' @param cols defining the color for each cell group
#' @param show.text.y whther show y-axis text
#' @param line.size line width in the violin plot
#' @param pt.size size of the dots
#' @param plot.margin adjust the white space between each plot
#' @param ... pass any arguments to VlnPlot in Seurat
#' @import ggplot2
#' @importFrom Seurat VlnPlot
#'
modify_vlnplot<- function(object,
                          features,
                          idents = NULL,
                          split.by = NULL,
                          cols = NULL,
                          show.text.y = TRUE,
                          line.size = NULL,
                          pt.size = 0,
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          ...) {
  options(warn=-1)
  p<- VlnPlot(object, features = features, cols = cols, pt.size = pt.size, idents = idents, split.by = split.by,  ... )  +
    xlab("") + ylab(features) + ggtitle("")
  p <- p + theme(text = element_text(size = 10)) + theme(axis.line = element_line(size=line.size)) +
    theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 8), axis.line.x = element_line(colour = 'black', size=line.size),axis.line.y = element_line(colour = 'black', size= line.size))
  # theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
  p <- p + theme(legend.position = "none",
                 plot.title= element_blank(),
                 axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_text(size = rel(1), angle = 0),
                 axis.text.y = element_text(size = rel(1)),
                 plot.margin = plot.margin ) +
    theme(axis.text.y = element_text(size = 8))
  p <- p + theme(element_line(size=line.size))

  if (!show.text.y) {
    p <- p + theme(axis.ticks.y=element_blank(), axis.text.y=element_blank())
  }
  return(p)
}

#' extract the max value of the y axis
#' @param p ggplot object
#' @importFrom  ggplot2 ggplot_build
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}



