##' Plot data with archetypes in 2D and 3D
##' @rdname plot_arc
##' @name plot_arc
##' @author Vitalii Kleshchevnikov
##' @description \code{plot_arc()} plot data with polytope representing the Pareto front, where vertices are archetypes (dots connected with lines). When archetype data is "r_pch_fit" all archetype locations from each subsample are shown with lines connecting the average location (type "average"); or lines connecting archetypes in each of the experiments (colored differently, type "all").
##' @param arc_data list of matrices storing the position of archetypes,  dim(dimensions, archetypes), class "pch_fit", "r_pch_fit". Each element of a list represents an independent run of the polytope fitting algorithm
##' @param data matrix of data in which archetypes/polytope were found, dim(variables/dimentions, examples)
##' @param which_dimensions indices or character vector specifying dimension names
##' @param type used when arc_data is "r_pch_fit", one of "average", "all"
##' @param average_func used when arc_data is "r_pch_fit", function telling how to find average position of vertices
##' @param geom plotting function to plot data in 2D, useful options are ggplot2::geom_point (scatterplot) and ggplot2::geom_bin2d (density)
##' @param colors character vector giving color palette for different archetype fits and the data (both 3D and 2D plot)
##' @param arch_size size of archetype point
##' @param line_size width of lines connecting archetypes
##' @return \code{plot_arc()} ggplot2 (2D) or plotly (3D) plot
##' @export plot_arc
##' @import plotly
##' @seealso \code{\link[ParetoTI]{fit_pch}}, \code{\link[ParetoTI]{arch_dist}}
##' @examples
##' library(ParetoTI)
##' library(ggplot2)
##' # Random data that fits into the triangle (2D)
##' set.seed(4355)
##' archetypes = generate_arc(arc_coord = list(c(5, 0), c(-10, 15), c(-30, -20)),
##'                           mean = 0, sd = 1, N_dim = 2)
##' data = generate_data(archetypes$XC, N_examples = 1e4, jiiter = 0.04, size = 0.9)
##' plot_arc(arch_data = archetypes, data = data,
##'     which_dimensions = 1:2) +
##'     theme_bw()
##' # Plot data as 2D density rather than points
##' plot_arc(arch_data = archetypes, data = data,
##'     which_dimensions = 1:2, geom = ggplot2::geom_bin2d)
##'
##' # Random data that fits into the triangle (3D)
##' set.seed(4355)
##' archetypes = generate_arc(arc_coord = list(c(5, 0, 4), c(-10, 15, 0), c(-30, -20, -5)),
##'                           mean = 0, sd = 1, N_dim = 3)
##' data = generate_data(archetypes$XC, N_examples = 1e4, jiiter = 0.04, size = 0.9)
##' plot_arc(arch_data = archetypes, data = data,
##'     which_dimensions = 1:3) +
##'     theme_bw()
plot_arc = function(arch_data, data, which_dimensions = as.integer(1:2),
                    type = c("average", "all")[1], average_func = mean,
                    geom = list(ggplot2::geom_point, ggplot2::geom_bin2d)[[1]],
                    colors = c("#D62728", "#1F77B4", "#2CA02C", "#17BED0", "#006400", "#FF7E0F"),
                    arch_size = NULL, line_size = NULL) {
  for_plot = ParetoTI:::.arc_data_table(arch_data, data)
  lines_for_plot = ParetoTI:::.archLines(for_plot, label = "archetypes", type, average_func)
  if(is(arch_data, "r_pch_fit") & type == "average"){
    for_plot[grepl("archetypes", lab), lab := "archetypes"]
    ly_arch_size = 2
    ly_line_size = 5
    gg_arch_size = 2
    gg_line_size = 1.5
  } else {
    ly_arch_size = 10
    ly_line_size = 5
    gg_arch_size = 5
    gg_line_size = 1.5
  }
  # set sizes to manual when provided
  if(!is.null(arch_size)){
    ly_arch_size = arch_size
    gg_arch_size = arch_size
  }
  if(!is.null(line_size)){
    ly_line_size = line_size
    gg_line_size = line_size
  }
  if(is.integer(which_dimensions)){
    x = colnames(for_plot)[which_dimensions[1]]
    y = colnames(for_plot)[which_dimensions[2]]
  } else if (is.character(which_dimensions)) {
    x = which_dimensions[1]
    y = which_dimensions[2]
  } else stop("which_dimensions is neither integer nor character vector")
  ## 2D plot ===================================================================##
  if(length(which_dimensions) == 2){
    plot = ggplot2::ggplot(for_plot, ggplot2::aes(x = get(x), y = get(y),
                                                  color = lab, group = lab)) +
      geom() +
      ggplot2::geom_path(data = lines_for_plot, inherit.aes = FALSE,
                         ggplot2::aes(x = get(x), y = get(y),
                                      color = lab, group = lab), size = gg_line_size) +
      ggplot2::geom_point(data = lines_for_plot, inherit.aes = FALSE,
                          ggplot2::aes(x = get(x), y = get(y),
                                       color = lab, group = lab), size = gg_arch_size) +
      ggplot2::xlab(colnames(for_plot)[which_dimensions[1]]) +
      ggplot2::ylab(colnames(for_plot)[which_dimensions[2]]) +
      ggplot2::scale_color_manual(values = colors)
    ## 3D plot ===================================================================##
  } else if(length(which_dimensions) == 3 & nrow(data) >= 3) {
    if(is.integer(which_dimensions)){
      z = colnames(for_plot)[which_dimensions[3]]
    } else if (is.character(which_dimensions)) {
      z = which_dimensions[3]
    } else stop("which_dimensions is neither integer nor character vector")
    x = as.formula(paste0("~", x))
    y = as.formula(paste0("~", y))
    z = as.formula(paste0("~", z))
    plot = plot_ly(for_plot, x = x, y = y, z = z,
                   color = ~ lab, colors = colors,
                   marker = list(size = 2)) %>%
      add_markers()
    plot = add_trace(p = plot, x = x, y = y, z = z, mode = 'lines',
                     data = lines_for_plot,
                     marker = list(size = ly_arch_size),
                     line = list(width = ly_line_size),
                     inherit = TRUE)
  } else stop("asked to plot < 2 or > 3 dimensions, or dataset has less dimensions than specified by which_dimensions")
  plot
}
