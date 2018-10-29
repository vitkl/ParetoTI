##' Find features that are decreasing functions of distance from archetype
##' @rdname find_decreasing
##' @name find_decreasing
##' @description \code{find_decreasing()} Finds
##' @param gam_fit a fitted gam object as produced by gam().
##' @param n_points number of points at which to evaluate derivative
##' @param x numeric vector where to evaluate derivatives
##' @param N_smooths number of curves for each predictor in the model
##' @param weights how to weight points along x axis when calculating mean (integral) probability
##' @param d numeric vector (1L), finite difference interval
##' @param return_derivs return derivative and SD of derivative values? By default only summary p-values are returned
##' @param return_gam return gam model as well? By default only summary p-values are returned
##' @return \code{find_decreasing()}
##' @export find_decreasing
##' @import data.table
find_decreasing = function(){

}
##' @rdname find_decreasing
##' @name fit_arc_gam1
##' @description \code{fit_arc_gam1()} Fits gam models one by one to find if a small number of features are a decreasing function of distance from archetype. Both gam functions and first derivatives can be visualised using plot() method.
##' @param data_attr data.table dim(examples, dimensions) that includes distance of each example to archetype in columns given by \code{arc_col} and feature values given by \code{features}
##' @param arc_col character vector, columns that give distance to archetypes (column per archetype)
##' @param features character vector (1L), column than containg feature values
##' @param min.sp lower bound for the smoothing parameter, details: \link[mgcv]{gam}. Default value of 60 works well to stabilise curve shape near min and max distance
##' @param N_smooths number of bases used to represent the smooth term (\link[mgcv]{s}), 4 for cubic splines
##' @param n_points number of points at which to evaluate derivative
##' @param d numeric vector (1L), finite difference interval
##' @param weights how to weight points along x axis when calculating mean (integral) probability. Useful if you care that the function is decreasing near the archetype but not far away. Two defaults suggest to weight point equally or discard bottom 50%.
##' @param return_only_summary set to TRUE when using inside other function to fit multiple features
##' @param stop_at_10 prevents \code{fit_arc_gam1()} from fitting too many features
##' @param ... arguments passed to \link[mgcv]{gam}
##' @return \code{fit_arc_gam1()} list (S3 object, gam_deriv) containing summary p-values for features and each archetype, function call and (optionally) a data.table with values of the first derivative and gam model fit
##' @export fit_arc_gam1
##' @import data.table
fit_arc_gam1 = function(data_attr, arc_col,
                        features = c("Gpx1", "Alb", "Cyp2e1", "Apoa2")[3],
                        min.sp = c(60), N_smooths = 4,
                        n_points = 200, d = 1 / n_points,
                        weights = c(rep(1, each = n_points),
                                    rep(c(1, 0), each = n_points / 2))[1],
                        return_only_summary = FALSE, stop_at_10 = TRUE,
                        ...){
  if(length(features) > 10 & isTRUE(stop_at_10)) stop("Trying to fit more than 10 features. fit_arc_gam1() is designed for visualising few features rather than bulk processing. Please use fit_arc_gam() or set stop_at_10 = FALSE")
  # encode combinations of features and columns
  combs = as.matrix(expand.grid(features, arc_col))
  # start loop
  res = lapply(seq_len(nrow(combs)), function(i){
    # decode combinations of feature and column
    feature = combs[i, 1]
    col = combs[i, 2]
    # generate formula specifying the model
    form = paste0(feature," ~ ", paste0("s(", col, ", bs = \"cr\", k = ", N_smooths,")", collapse = " + "))
    # fit gam model
    gam_fit = mgcv::gam(as.formula(form), data = data_attr,
                        min.sp = min.sp, ...)
    # find 1st derivative by finite differencing
    derivs = find_gam_deriv(gam_fit, N_smooths = N_smooths, d = d,
                            n_points = n_points, weights = weights,
                            return_derivs = !return_only_summary, return_gam = FALSE)
    # generate and returm summary of gam fit and the derivative
    summary = derivs$summary
    sm = summary(gam_fit)
    summary$smooth_term_p_val = sm$s.pv
    summary$r_sq =  sm$r.sq

    if(return_only_summary) return(summary) else {
      list(summary = summary, derivs = derivs$derivs, gam_fit = gam_fit)
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
    class(res_n) = "gam_deriv"
    res_n
  }
}
