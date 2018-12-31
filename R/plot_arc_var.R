##' Plot SSE and variance explained by polytope fit
##' @rdname plot_arc_var
##' @name plot_arc_var
##' @author Vitalii Kleshchevnikov
##' @description \code{plot_arc_var()} shows SSE and variance explained by polytope models with different number of vertices (k)
##' @param arc_data object of class "k_pch_fit", storing the position of vertices of best fit polytopes with different k, and other data from \code{\link[ParetoTI]{fit_pch}}() run. arc_data$XC is a list where each element is XC matrix of dim(dimensions, archetypes) storing positions of vertices with different k.
##' @param type which measure to plot as a function of k, one of "varexpl", "SSE", "res_varexpl", "dim". Use dim to plot variance in position in each dimension.
##' @param arch_size size of archetype point
##' @param line_size width of lines connecting archetypes
##' @param reorder reorder dimensions based on variance in position (type = "dim).
##' @return \code{plot_arc_var()} ggplot2 (2D) plot
##' @export plot_arc_var
##' @seealso \code{\link[ParetoTI]{fit_pch}}, \code{\link[ParetoTI]{plot_arc}}
##' @examples
##' library(ParetoTI)
##' library(ggplot2)
##' # Random data that fits into the triangle (2D)
##' set.seed(4355)
##' archetypes = generate_arc(arc_coord = list(c(5, 0), c(-10, 15), c(-30, -20)),
##'                           mean = 0, sd = 1, N_dim = 2)
##' data = generate_data(archetypes$XC, N_examples = 1e4, jiiter = 0.04, size = 0.9)
##' arc_data = k_fit_pch(data, ks = 1:4, check_installed = T, delta=0.1, order_by_side = F)
##' # Show polytopes and the data
##' plot_arc(arch_data = arc_data, data = data,
##'          which_dimensions = 1:2, type = "all", arch_size = 2,
##'          colors = c("#D62728", "#1F77B4", "#2CA02C", "#17BED0", "grey")) +
##'   theme_bw()
##' # Show variance explained by a polytope with each k
##' plot_arc_var(arc_data, type = c("varexpl", "SSE", "res_varexpl")[1],
##'              point_size = 2, line_size = 1.5) + theme_bw()
plot_arc_var = function(arc_data, type = c("varexpl", "SSE", "res_varexpl",
                                           "total_var", "dim")[1],
                        point_size = 2, line_size = 1.5, reorder = FALSE){
  if(!(is(arc_data, "k_pch_fit") | is(arc_data, "b_pch_fit") |
       is(arc_data, "pch_fit"))) {
    stop("arc_data should be k_pch_fit, b_pch_fit or pch_fit")
  }

  if(type != "dim"){
    type_lab = .type_lab(type)
    k_var = data.table(varexpl = arc_data$pch_fits$varexpl, SSE = arc_data$pch_fits$SSE,
                       k = sapply(arc_data$pch_fits$XC, ncol),
                       total_var = arc_data$pch_fits$total_var,
                       t_ratio = arc_data$pch_fits$t_ratio)
    setorder(k_var, k)
    k_var[k == min(k), res_varexpl := varexpl]
    for (i in seq(min(k_var$k)+1, max(k_var$k))) {
      k_var[k == i, res_varexpl := varexpl - k_var[k == i - 1, varexpl]]
    }
    ggplot2::ggplot(k_var, ggplot2::aes(x = k, y = get(type))) +
      ggplot2::geom_path(size = line_size) + ggplot2::geom_point(size = point_size) +
      ggplot2::xlab("k, number of vertices/archetypes") +
      ggplot2::ylab(type_lab)
  } else {
    if(is(arc_data, "k_pch_fit")){
      var_dim = arc_data$pch_fits$var_dim
    } else {
      var_dim = arc_data$var_dim
    }

    var_dim = melt.data.table(var_dim, id.vars = "k")
    setorder(var_dim, k, value)
    if(reorder){
      var_dim[, variable := factor(variable,
                                   levels = unique(variable[order(value)]))]
    }

    ggplot2::ggplot(var_dim, ggplot2::aes(x = variable, y = value, group = k)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~ k) +
      ggplot2::xlab("dimension name") +
      ggplot2::ylab("Variance in position across vertices")
  }

}

.type_lab = function(type, short = FALSE){
  if(short) {
    type_lab = c("Variance explained", "Sum of squared errors", "Residual \nvariance explained", "Variance in positions", "Volume of polytope /\n convex hull")
  } else {
    type_lab = c("Variance explained", "Sum of squared errors", "Variance explained on top of k-1 model", "Mean variance in position of vertices", "T-ratio of volume of polytope by volume of convex hull")
  }
  names_type_lab = c("varexpl", "SSE", "res_varexpl", "total_var", "t_ratio")
  type_lab[names_type_lab %in% type]
}
