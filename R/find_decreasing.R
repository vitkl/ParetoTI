##' Find features that are decreasing functions of distance from archetype
##' @rdname find_decreasing
##' @name find_decreasing
##' @description \code{find_decreasing()} Fits gam models to find if features are a decreasing function of distance from archetype. Both gam functions and first derivatives can be visualised using plot() method.
##' @param data_attr data.table dim(examples, dimensions) that includes distance of each example to archetype in columns given by \code{arc_col} and feature values given by \code{features}
##' @param arc_col character vector, columns that give distance to archetypes (column per archetype)
##' @param features character vector (1L), column than containg feature values
##' @param min.sp lower bound for the smoothing parameter, details: \link[mgcv]{gam}. Default value of 60 works well to stabilise curve shape near min and max distance
##' @param N_smooths number of bases used to represent the smooth term (\link[mgcv]{s}), 4 for cubic splines
##' @param n_points number of points at which to evaluate derivative
##' @param d numeric vector (1L), finite difference interval
##' @param weights how to weight points along x axis when calculating mean (integral) probability. Useful if you care that the function is decreasing near the archetype but not far away. Two defaults suggest to weight point equally or discard bottom 50 percent.
##' @param return_only_summary set to TRUE when using inside other function to fit multiple features
##' @param stop_at_10 prevents \code{find_decreasing()} from fitting too many features
##' @param one_arc_per_model If TRUE fit separate gam models for each archetype. If FALSE combine all archetypes in one model: feature ~ s(arc1) + s(arc2) + ... + s(arcN).
##' @param ... arguments passed to \link[mgcv]{gam}
##' @return \code{find_decreasing()} list (S3 object, gam_deriv) containing summary p-values for features and each archetype, function call and (optionally) a data.table with values of the first derivative
##' @export find_decreasing
##' @import data.table
find_decreasing = function(data_attr, arc_col,
                           features = c("Gpx1", "Alb", "Cyp2e1", "Apoa2")[3],
                           min.sp = c(60), N_smooths = 4,
                           n_points = 200, d = 1 / n_points,
                           weights = c(rep(1, each = n_points),
                                       rep(c(1, 0), each = n_points / 2))[1],
                           return_only_summary = FALSE, stop_at_10 = TRUE,
                           one_arc_per_model = TRUE,
                           ...){
  if(length(features) > 10 & isTRUE(stop_at_10) & !isTRUE(return_only_summary)) stop("Trying to fit more than 10 features requesting values of 1st derivative. return_only_summary = FALSE option is designed for visualising few features rather than bulk processing. Please use return_only_summary = TRUE or set stop_at_10 = FALSE")
  # start loop
  res = lapply(seq_len(length(features)), function(i){
    feature = features[i]
    if(one_arc_per_model) {
      # if fit a model for each archetype, iterate over archetypes
      derivs = lapply(arc_col, function(col){
        fit_arc_gam_1(feature = feature, col = col, N_smooths = N_smooths,
                      data_attr = data_attr, min.sp = min.sp, ..., d = d,
                      n_points = n_points, weights = weights)
      })
      # combine results
      res = list()
      res$call = match.call()
      res$derivs = rbindlist(lapply(derivs, function(x) x$derivs))
      res$gam_fit = lapply(derivs, function(x) x$gam_fit)
      res$gam_sm = rbindlist(lapply(derivs, function(x) x$gam_sm))
      class(res) = "gam_deriv"
    } else {
      res = fit_arc_gam_1(feature = feature, col = arc_col, N_smooths = N_smooths,
                          data_attr = data_attr, min.sp = min.sp, ..., d = d,
                          n_points = n_points, weights = weights)
      res$gam_fit = list(res$gam_fit)
    }
    # generate summary of the derivative
    summary = summary.gam_deriv(res)
    if(return_only_summary) return(summary) else {
      list(summary = summary, derivs = res$derivs, gam_fit = res$gam_fit)
    }
  })
  if(return_only_summary){
    return(rbindlist(res))
  } else {
    res_n = list()
    res_n$call = match.call()
    res_n$summary = rbindlist(lapply(res, function(x) x$summary))
    res_n$derivs = rbindlist(lapply(res, function(x) x$derivs))
    res_n$gam_fit = lapply(res, function(x) x$gam_fit)
    res_n$gam_fit = unlist(res_n$gam_fit, recursive = FALSE)
    res_n$gam_sm = NA
    class(res_n) = "gam_deriv"
    res_n
  }
}


##' @rdname find_decreasing
##' @name fit_arc_gam_1
##' @description \code{fit_arc_gam_1()} Finds single GAM model fit and it's first derivative for a single feature and one or several archetypes.
##' @return \code{fit_arc_gam_1()} list containing function call, 1st derivative values of GAM model (derivs), summary of GAM model (p-value and r^2, gam_sm)
##' @export fit_arc_gam_1
##' @import data.table
fit_arc_gam_1 = function(feature, col, N_smooths, data_attr, min.sp, ..., d,
                         n_points, weights){
  # generate formula specifying the model
  form = paste0(feature," ~ ", paste0("s(", col, ", bs = \"cr\", k = ", N_smooths,")", collapse = " + "))
  # fit gam model
  gam_fit = mgcv::gam(as.formula(form), data = data_attr,
                      min.sp = min.sp, ...)
  # find 1st derivative by finite differencing
  derivs = find_gam_deriv(gam_fit, N_smooths = N_smooths, d = d,
                          n_points = n_points, weights = weights,
                          return_gam = FALSE)
  derivs = list(call = derivs$call, derivs = derivs$derivs,
                gam_fit = gam_fit, gam_sm = derivs$gam_sm,
                summary = NA)
  class(derivs) = "gam_deriv"
  derivs
}
