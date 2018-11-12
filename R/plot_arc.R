##' Plot data with archetypes in 2D and 3D
##' @rdname plot_arc
##' @name plot_arc
##' @author Vitalii Kleshchevnikov
##' @description \code{plot_arc()} plot data with polytope representing the Pareto front, where vertices are archetypes (dots connected with lines). When archetype data is "b_pch_fit" all archetype locations from each subsample are shown with lines connecting the average location (type "average"); or lines connecting archetypes in each of the experiments (colored differently, type "all").
##' @param arc_data objects of class "pch_fit", "b_pch_fit", "k_pch_fit" storing the position of archetypes, and other data from \code{\link[ParetoTI]{fit_pch}}() run. arc_data$XC is matrix of dim(dimensions, archetypes) or list where each element is XC matrix from an independent run of the polytope fitting algorithm. Set to NULL if you want to show data alone.
##' @param data matrix of data in which archetypes/polytope were found, dim(variables/dimentions, examples)
##' @param which_dimensions indices or character vector specifying dimension names. When \code{which_dimensions} exceeds the number of dimensions in arc_data these archetypes will be omitted. This can happen when fitting simplexes: lines and triangles are only 2D, so will be omitted from 3D plots.
##' @param type used when arc_data is "b_pch_fit", one of "average", "all"
##' @param average_func used when arc_data is "b_pch_fit", function telling how to find average position of vertices
##' @param geom plotting function to plot data in 2D, useful options are ggplot2::geom_point (scatterplot) and ggplot2::geom_bin2d (density)
##' @param colors character vector giving color palette for different archetype fits and the data (both 3D and 2D plot)
##' @param arch_size size of archetype points
##' @param line_size width of lines connecting archetypes
##' @param data_size size of data points in plotly. Values for ggplot are 1/2 of data_size.
##' @param arch_alpha opacity of archetype points
##' @param data_lab vector, 1L or length of data, label data points (examples) with a qualitative or quantitative label
##' @param arc_lab vector, 1L or nrow(arc_data$XC)/noc, label vertices/archetypes (points) with a qualitative or quantitative label
##' @param legend_name name to display on legend, e.g. gene name in data_lab
##' @param text_size archetype label text size
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
plot_arc = function(arch_data = NULL, data, which_dimensions = as.integer(1:2),
                    type = c("average", "all")[1], average_func = mean,
                    geom = list(ggplot2::geom_point, ggplot2::geom_bin2d)[[1]],
                    colors = c("#1F77B4", "#D62728", "#2CA02C", "#17BED0", "#006400", "#FF7E0F"),
                    arch_size = NULL, line_size = NULL,
                    data_size = 2, arch_alpha = 0.4,
                    data_lab = "data", arc_lab = "archetypes",
                    legend_name = "data",
                    text_size = NULL, nudge = c(0.05, 0.1)) {
  if(!is.null(arch_data)){
    for_plot = ParetoTI:::.arc_data_table(arch_data, data, data_lab = data_lab,
                                          which_dimensions = which_dimensions)
    lines_for_plot = ParetoTI:::.archLines(for_plot$arc_data, arc_lab = arc_lab,
                                           type, average_func)
  } else {
    for_plot = list()
    for_plot$data = as.data.table(t(data))
    for_plot$data$lab = data_lab
  }
  if(is(arch_data, "b_pch_fit") & type == "average"){
    for_plot$arc_data[grepl("archetypes", lab), lab := "archetypes"]
    for_plot$arc_data[, lab := factor(lab, levels = sort(unique(lab), decreasing = TRUE))]
    setorder(for_plot$arc_data, lab)
    ly_arch_size = 2
    ly_line_size = 5
    gg_arch_size = 2
    gg_line_size = 1.5
  } else {
    for_plot$arc_data[, lab := factor(lab, levels = sort(unique(lab), decreasing = TRUE))]
    setorder(for_plot$arc_data, lab)
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
  if(!is.null(text_size)){
    ly_text_size = text_size
    gg_text_size = text_size
  } else {
    ly_text_size = 20
    gg_text_size = 8
  }
  # get column names for corresponding dimensions
  if(is.integer(which_dimensions)){
    x = colnames(for_plot$arc_data)[1]
    y = colnames(for_plot$arc_data)[2]
  } else if (is.character(which_dimensions)) {
    x = which_dimensions[1]
    y = which_dimensions[2]
  } else stop("which_dimensions is neither integer nor character vector")

  # assign color to data
  if(is.numeric(for_plot$data$lab)){
    n_data_lab = 1
    names_data_lab = NULL
  } else {
    n_data_lab = uniqueN(for_plot$data$lab)
    names_data_lab = as.character(unique(for_plot$data$lab))
  }
  data_colors = colors[seq(1, n_data_lab)]
  names(data_colors) = names_data_lab
  arc_colors = colors[seq(n_data_lab + 1,
                          n_data_lab + uniqueN(for_plot$arc_data$lab))]
  names(arc_colors) = as.character(unique(for_plot$arc_data$lab))

  ## 2D plot ===================================================================##
  if(length(which_dimensions) == 2){
    plot = ggplot2::ggplot(for_plot$data, ggplot2::aes(x = get(x), y = get(y),
                                                       color = lab))
    if(is.numeric(for_plot$data$lab)){
      plot = plot + geom(size = data_size/2)
    } else {
      plot = plot + geom(color = data_colors, size = data_size/2)
    }
    if("lines_for_plot" %in% ls()) {
      # calculate nudge by distance for text labels
      nd_x = for_plot$arc_data[, range(get(x))]
      nd_x = abs(nd_x[1] - nd_x[2]) * nudge[1]
      nd_y = for_plot$arc_data[, range(get(y))]
      nd_y = abs(nd_y[1] - nd_y[2]) * nudge[2]
      # plot archetypes
      plot = plot +
        ggplot2::geom_point(data = for_plot$arc_data, inherit.aes = FALSE,
                            ggplot2::aes(x = get(x), y = get(y),
                                         group = lab), size = gg_arch_size,
                            color = arc_colors, alpha = arch_alpha) +
        ggplot2::geom_path(data = lines_for_plot, inherit.aes = FALSE,
                           ggplot2::aes(x = get(x), y = get(y),
                                        group = lab), size = gg_line_size,
                           color = arc_colors) +
        #ggplot2::scale_color_manual(values = colors) +
        ggplot2::geom_point(data = lines_for_plot, inherit.aes = FALSE,
                            ggplot2::aes(x = get(x), y = get(y),
                                         group = lab), size = gg_arch_size,
                            color = arc_colors) +
        ggplot2::geom_text(data = unique(lines_for_plot), inherit.aes = FALSE,
                           ggplot2::aes(x = get(x), y = get(y), label = arch_id,
                                        group = lab),
                           color = arc_colors,
                           show.legend = FALSE, size = gg_text_size,
                           nudge_x = nd_x, nudge_y = nd_y) # add separate step for plotting archetype positions
    }
    plot = plot + ggplot2::xlab(x) + ggplot2::ylab(y)
    ## 3D plot ===================================================================##
  } else if(length(which_dimensions) == 3 & nrow(data) >= 3) {
    if(is.integer(which_dimensions)){
      z = colnames(for_plot$arc_data)[3]
    } else if (is.character(which_dimensions)) {
      z = which_dimensions[3]
    } else stop("which_dimensions is neither integer nor character vector")
    x = as.formula(paste0("~", x))
    y = as.formula(paste0("~", y))
    z = as.formula(paste0("~", z))

    setorder(for_plot$data, lab)
    plot = plot_ly(for_plot$data)

    if(is.numeric(for_plot$data$lab)){
      plot = add_markers(p = plot, x = x, y = y, z = z, showlegend = TRUE, mode = "markers",
                         color = ~ lab, name = 'data',
                         marker = list(size = data_size,
                                       colorbar = list(title = legend_name)))
    } else {
      plot = add_markers(p = plot, x = x, y = y, z = z, mode = "markers",
                         colors = ~ lab, name = ~ lab,
                         marker = list(size = data_size,
                                       color = data_colors[for_plot$data$lab]))
    }

    if("lines_for_plot" %in% ls()) {
      setorder(for_plot$arc_data, lab)
      setorder(lines_for_plot, lab)
      arc_colors = arc_colors[for_plot$arc_data$lab]
      arc_line_colors = arc_colors[lines_for_plot$lab]
      arc_text_colors = arc_colors[unique(lines_for_plot)$lab]
      plot = add_markers(p = plot, x = x, y = y, z = z, showlegend = FALSE,
                         mode = "markers", colors = ~ lab, name = ~ lab,
                         marker = list(color = arc_colors, size = ly_arch_size,
                                       opacity = arch_alpha),
                         data = for_plot$arc_data,
                         inherit = TRUE) %>%
        add_trace(x = x, y = y, z = z, mode = 'lines',
                  data = lines_for_plot, colors = ~ lab,
                  showlegend = TRUE, name = ~ lab,
                  marker = list(size = ly_arch_size, color = arc_line_colors),
                  line = list(width = ly_line_size, color = arc_line_colors),
                  inherit = TRUE) %>%
        add_text(mode = "markers", showlegend = FALSE,
                 x = x, y = y, z = z, textposition = "top center",
                 colors = ~ lab, name = ~ lab,
                 data = unique(lines_for_plot), inherit = TRUE,
                 marker = list(size = 0, color = arc_colors),
                 textfont = list(color = arc_colors, size = ly_text_size),
                 text = ~ arch_id)
    }
  } else stop("asked to plot < 2 or > 3 dimensions, or dataset has less dimensions than specified by which_dimensions")
  plot
}
