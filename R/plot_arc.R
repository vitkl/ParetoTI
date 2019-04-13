##' Plot data with archetypes in 2D, 3D and a panel of 2D
##' @rdname plot_arc
##' @name plot_arc
##' @author Vitalii Kleshchevnikov
##' @description \code{plot_arc()} plot data with polytope representing the Pareto front, where vertices are archetypes (dots connected with lines). When archetype data is "b_pch_fit" all archetype locations from each subsample are shown with lines connecting the average location (type "average"); or lines connecting archetypes in each of the experiments (colored differently, type "all").
##' @param arc_data objects of class "pch_fit", "b_pch_fit", "k_pch_fit" storing the position of archetypes, and other data from \code{\link[ParetoTI]{fit_pch}}() run. arc_data$XC is matrix of dim(dimensions, archetypes) or list where each element is XC matrix from an independent run of the archetypal analysis. Set to NULL if you want to show data alone.
##' @param data matrix of data in which archetypes/polytope were found, dim(variables/dimentions, examples)
##' @param which_dimensions indices or character vector specifying dimension names. 2D plot, 3D plot or a panel for 2D plots when more than 3 dimensions. For \code{arch_to_tsne()} use 1:2 or 1:3. When \code{which_dimensions} exceeds the number of dimensions in arc_data these archetypes will be omitted. This can happen when fitting simplexes: lines and triangles are only 2D, so will be omitted from 3D plots.
##' @param type used when arc_data is "b_pch_fit", one of "average", "all"
##' @param average_func used when arc_data is "b_pch_fit", function telling how to find average position of vertices
##' @param geom plotting function to plot data in 2D, useful options are ggplot2::geom_point (scatterplot) and ggplot2::geom_bin2d (density)
##' @param colors character vector giving color palette for different archetype fits and the data (both 3D and 2D plot)
##' @param arch_size size of archetype points
##' @param line_size width of lines connecting archetypes
##' @param data_size size of data points in plotly. Values for ggplot are 1/2 of data_size.
##' @param arch_alpha opacity of archetype points
##' @param data_lab vector, 1L or length of data, label data points (examples) with a qualitative or quantitative label
##' @param arc_lab vector, 1L or nrow(arc_data$XC)/noc, label vertices/archetypes (points) with a categorical. Only used when looking at a single fit (pch_fit).
##' @param arc_names_num logical, when archetypes are named, use numbers (default, TRUE), or names (FALSE, produces cluttered plot)?
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
##'                           mean = 0, sd = 1)
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
##'                           mean = 0, sd = 1)
##' data = generate_data(archetypes$XC, N_examples = 1e4, jiiter = 0.04, size = 0.9)
##'
##' plot_arc(arch_data = archetypes, data = data,
##'     which_dimensions = 1:3)
##'
##' # Project to tSNE coordinates (from 3D to 2D)
##' arc_tsne = arch_to_tsne(archetypes, data, which_dimensions = 1:2)
##' plot_arc(arch_data = arc_tsne$arch_data, data = arc_tsne$data,
##'     which_dimensions = 1:2) +
##'     theme_bw()
##'
##' # Project to UMAP representation
##' arc_umap = arch_to_umap(archetypes, data, which_dimensions = 1:2,
##'                         method = c("naive", # implemented in R and slow
##'                                    "umap-learn")) # requires python module
##' plot_arc(arch_data = arc_umap$arch_data, data = arc_umap$data,
##'     which_dimensions = 1:2) +
##'     theme_bw()
plot_arc = function(arch_data = NULL, data, which_dimensions = as.integer(1:2),
                    type = c("average", "all")[1], average_func = mean,
                    geom = list(ggplot2::geom_point, ggplot2::geom_bin2d)[[1]],
                    colors = c("#1F77B4", "#D62728", "#2CA02C", "#17BED0", "#006400", "#FF7E0F"),
                    arch_size = NULL, line_size = NULL,
                    data_size = 4, arch_alpha = 0.4,
                    data_lab = "data", arc_lab = "archetypes", arc_names_num = TRUE,
                    legend_name = "data",
                    text_size = NULL, nudge = c(0.05, 0.1)) {


  if((uniqueN(data_lab) + uniqueN(arc_lab)) > uniqueN(colors) & !is.integer(data_lab) & ! is.numeric(data_lab)) {
    stop("uniqueN(data_lab) > colors, please add more colors")
  }
  if(length(arc_lab) > 1 & !(is(arch_data, "pch_fit") |
                             is(arch_data, "random_arc") | is.null(arch_data))) {
    stop("Archetype labels can be used only with single fit of the model
          - class(arch_data) == 'pch_fit'")
  }

  if(!is.null(arch_data)) {

    for_plot = ParetoTI:::.arc_data_table(arch_data, data,
                                          data_lab = data_lab, arc_lab = arc_lab,
                                          which_dimensions = which_dimensions)

    # convert names to numbers
    if(arc_names_num) for_plot$arc_data[, arch_id := gsub("^.+_|$", "", arch_id)]

    lines_for_plot = ParetoTI:::.archLines(for_plot$arc_data, arc_lab = arc_lab,
                                           pch_fit = is(arch_data, "pch_fit") | is(arch_data, "random_arc"),
                                           type, average_func)

  } else {

    for_plot = list()
    for_plot$data = as.data.table(Matrix::t(data))
    for_plot$data$lab = data_lab

  }
  if(is(arch_data, "b_pch_fit") & type == "average"){

    for_plot$arc_data[grepl("archetypes", lab), lab := "archetypes"]
    for_plot$arc_data[, lab := factor(lab, levels = sort(unique(lab), decreasing = TRUE))]
    setorder(for_plot$arc_data, lab)
    ly_arch_size = 3
    ly_line_size = 5
    gg_arch_size = 2
    gg_line_size = 1.5

  } else {
    if(!is.null(arch_data)) {

      for_plot$arc_data[, lab := factor(lab, levels = sort(unique(lab), decreasing = TRUE))]
      setorder(for_plot$arc_data, lab)

    }
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
    x = colnames(for_plot$data)[1]
    y = colnames(for_plot$data)[2]
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
  # assign color to data
  data_colors = colors[seq(1, n_data_lab)]
  names(data_colors) = names_data_lab
  # assign remaining colors to archetypes and lines connecting them
  arc_colors = colors[seq(n_data_lab + 1,
                          n_data_lab + uniqueN(for_plot$arc_data$lab))]
  names(arc_colors) = as.character(unique(for_plot$arc_data$lab))

  ## 2D plot ===================================================================##
  if(length(which_dimensions) == 2){

    setorder(for_plot$data, lab)

    plot = ggplot2::ggplot(for_plot$data, ggplot2::aes(x = get(x), y = get(y),
                                                       color = lab))
    if(is.numeric(for_plot$data$lab)){

      plot = plot + geom(size = data_size/2)

    } else if(identical(geom, geom_bin2d)) {

      plot = plot + geom()

    } else {

      plot = plot + geom(size = data_size/2) +
        scale_color_manual(aesthetics = "color", values = data_colors[for_plot$data$lab]) +
        guides(color = guide_legend(title="data"))

    }
    if("lines_for_plot" %in% ls()) {
      # calculate nudge by distance for text labels
      nd_x = for_plot$arc_data[, range(get(x))]
      nd_x = abs(nd_x[1] - nd_x[2]) * nudge[1]
      nd_y = for_plot$arc_data[, range(get(y))]
      nd_y = abs(nd_y[1] - nd_y[2]) * nudge[2]

      setorder(for_plot$arc_data, lab)
      setorder(lines_for_plot, lab)
      # generate archetype colors of the same length as data
      arc_colors = arc_colors[for_plot$arc_data$lab]
      arc_line_colors = arc_colors[lines_for_plot$lab]
      arc_text_colors = arc_colors[unique(lines_for_plot)$lab]

      # plot archetypes
      plot = plot +
        ggplot2::geom_point(data = for_plot$arc_data, inherit.aes = FALSE,
                            ggplot2::aes(x = get(x), y = get(y),
                                         group = lab), size = gg_arch_size,
                            color = arc_colors, alpha = arch_alpha) +
        ggplot2::geom_path(data = lines_for_plot, inherit.aes = FALSE,
                           ggplot2::aes(x = get(x), y = get(y),
                                        group = lab), size = gg_line_size,
                           color = arc_line_colors)
      #ggplot2::scale_color_manual(values = colors) +
      plot = plot + ggplot2::geom_point(data = lines_for_plot, inherit.aes = FALSE,
                                        ggplot2::aes(x = get(x), y = get(y),
                                                     group = lab, colour = lab), size = gg_arch_size,
                                        color = arc_line_colors)

      if(uniqueN(arch_data$summary$k) == 1){
        plot = plot +
          ggplot2::geom_text(data = unique(lines_for_plot), inherit.aes = FALSE,
                             ggplot2::aes(x = get(x), y = get(y), label = arch_id,
                                          group = lab),
                             color = arc_text_colors,
                             show.legend = FALSE, size = gg_text_size,
                             nudge_x = nd_x, nudge_y = nd_y)
      }
    }
    plot = plot + ggplot2::xlab(x) + ggplot2::ylab(y)

    if(is.numeric(for_plot$data$lab)){
      # if cells are colored on a gradient - add nice palette
      plot = plot + ggplot2::scale_color_viridis_c()
    }

    ## 3D plot ===================================================================##

  } else if(length(which_dimensions) == 3 & nrow(data) >= 3) {

    if(is.integer(which_dimensions)){
      z = colnames(for_plot$data)[3]
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

      data_colors_2 = data_colors[as.character(for_plot$data$lab)]
      #data_colors_2 = factor(data_colors_2, levels = data_colors)
      plot = add_markers(p = plot, x = x, y = y, z = z, mode = "markers",
                         colors = ~ lab, name = ~ lab,
                         marker = list(size = data_size,
                                       colors = data_colors_2))
    }

    if("lines_for_plot" %in% ls()) {
      setorder(for_plot$arc_data, lab)
      setorder(lines_for_plot, lab)
      # generate archetype colors of the same length as data
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
                  inherit = TRUE)
      if(uniqueN(arch_data$summary$k) == 1){
        plot = add_text(p = plot, mode = "markers", showlegend = FALSE,
                        x = x, y = y, z = z, textposition = "top center",
                        colors = ~ lab, name = ~ lab,
                        data = unique(lines_for_plot), inherit = TRUE,
                        marker = list(size = 0, color = arc_colors),
                        textfont = list(color = arc_colors, size = ly_text_size),
                        text = ~ arch_id)
      }
    }
  } else if (length(which_dimensions) > 3 & nrow(data) >= length(which_dimensions)) {

    # create a matrix of possible pairwise combinations
    combs = expand.grid(which_dimensions, which_dimensions)
    # remove the same dimension
    combs = combs[combs[, 1] != combs[, 2],]
    # remove the same pair in reverse order
    pair_id = apply(combs, 1, function(x) paste0(sort(x), collapse = ""))
    combs = split(combs, pair_id)
    combs = t(vapply(combs, function(x) as.integer(x[1,]), integer(2)))

    # remove row names that screw up subsetting
    rownames(combs) = NULL

    plot = list()
    for (i in seq_len(nrow(combs))) {
      dims = as.integer(combs[i, ])

      p_pca = plot_arc(arch_data = arch_data, data = data,
                       which_dimensions = dims,
                       type = type, average_func = average_func,
                       geom = geom, colors = colors,
                       arch_size = arch_size, line_size = line_size,
                       data_size = data_size, arch_alpha = arch_alpha,
                       data_lab = data_lab, arc_lab = arc_lab,
                       legend_name = legend_name,
                       text_size = text_size, nudge = nudge)

      plot = c(plot, list(p_pca + theme(legend.position = "none")))
    }
    plot = cowplot::plot_grid(plotlist = plot)

  } else stop("dataset has less dimensions than specified by which_dimensions")

  plot

}

##' @rdname plot_arc
##' @name arch_to_tsne
##' @description arch_to_tsne() Project archetype positions to tSNE coordinates (2D or 3D) using \code{\link[Rtsne]{Rtsne}}.
##' @param pca perform PCA? Argument to \code{\link[Rtsne]{Rtsne}}.
##' @param partial_pca perform partial PCA? Argument to \code{\link[Rtsne]{Rtsne}}.
##' @param ... additional arguments to \code{\link[Rtsne]{Rtsne}} and \code{\link[umap]{umap}}.
##' @return arch_to_tsne() list with: arch_data containing archetype positions in tSNE coordinates, and data positions in tSNE coordinates
##' @export arch_to_tsne
arch_to_tsne = function(arch_data, data, which_dimensions = 1:2,
                        pca = FALSE, partial_pca = FALSE, ...) {

  if(!(is(arch_data, "pch_fit") | is(arch_data, "random_arc"))) {
    arch_data = average_pch_fits(arch_data)
  }

  colnames(arch_data$XC) = paste0("archetype", seq_len(ncol(arch_data$XC)))
  for_tnse = t(cbind(data, arch_data$XC))

  tnse = Rtsne::Rtsne(for_tnse, which_dimensions[2], pca = pca,
                      partial_pca = partial_pca, ...)
  tnse = tnse$Y
  colnames(tnse) = paste0("TSNE", seq_len(ncol(tnse)))
  rownames(tnse) = rownames(for_tnse)

  arch_data$XC = t(tnse[colnames(arch_data$XC),])

  list(arch_data = arch_data, data = t(tnse[seq_len(ncol(data)),]))
}

##' @rdname plot_arc
##' @name arch_to_umap
##' @description arch_to_umap() Project archetype positions to UMAP coordinates using \code{\link[umap]{umap}}.
##' @param method Method for finding UMAP representation. Available methods are 'naive' (an implementation written in pure R) and 'umap-learn' (requires python package 'umap-learn'). See \code{\link[umap]{umap}} for details.
##' @param n_neighbors sensible default for \code{\link[umap]{umap}}, pass other parameters via ...
##' @param min_dist sensible default for \code{\link[umap]{umap}}
##' @param metric sensible default for \code{\link[umap]{umap}}
##' @return arch_to_umap() list with: arch_data containing archetype positions in UMAP coordinates, data positions in UMAP coordinates, and umap_config parameters used to find this representation.
##' @export arch_to_umap
arch_to_umap = function(arch_data, data, which_dimensions = 1:2,
                        method = c("naive", "umap-learn")[1],
                        n_neighbors = 30L, min_dist = 0.3,
                        metric = ifelse(method[1] == "umap-learn",
                                        "correlation", "euclidean"), ...) {

  if(!(is(arch_data, "pch_fit") | is(arch_data, "random_arc"))) {
    arch_data = average_pch_fits(arch_data)
  }

  colnames(arch_data$XC) = paste0("archetype", seq_len(ncol(arch_data$XC)))
  for_umap = t(cbind(data, arch_data$XC))

  umap_out = umap::umap(for_umap, n_components = which_dimensions[2],
                        method = method[1],
                        n_neighbors = n_neighbors, min_dist = min_dist,
                        metric = metric , ...)
  umap_config = umap_out$config
  umap_out = umap_out$layout
  colnames(umap_out) = paste0("UMAP", seq_len(ncol(umap_out)))
  rownames(umap_out) = rownames(for_umap)

  arch_data$XC = t(umap_out[colnames(arch_data$XC),])

  list(arch_data = arch_data, data = t(umap_out[seq_len(ncol(data)),]),
       umap_config = umap_config)
}

