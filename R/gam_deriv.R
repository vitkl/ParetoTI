##' Find derivative of smooth curves from gam model fit (Generalised Additive Model, mgcv)
##' @rdname gam_deriv
##' @name gam_deriv
##' @description \code{gam_deriv()} Finds a derivative with standard errors of smooth curves from gam model fit by finite differencing. Crucially it computes the probability that derivative is negative at each point.
##' @details Plot method uses ggplot() to show mean derivative and confidence intervals (2 SD from the mean). This code is a modification of last example in \link[mgcv]{predict.gam}.
##' @param gam_fit a fitted gam object as produced by gam().
##' @param n_points number of points at which to evaluate derivative
##' @param x numeric vector where to evaluate derivatives
##' @param N_smooths number of curves for each predictor in the model
##' @param weights how to weight points along x axis when calculating mean (integral) probability
##' @param d numeric vector (1L), finite difference interval
##' @return \code{gam_deriv()} list (S3 object, gam_deriv) containing data.table-s for each predictor in gam model
##' @export gam_deriv
##' @import data.table
gam_deriv = function(gam_fit,
                     n_points = 200, x = NULL,
                     N_smooths = 9,
                     weights = list(rep(1, 200), sqrt(seq(1, 0, length.out = 200)),
                                 1 / (1 + seq(1, 0, length.out = 200)^3))[[1]],
                     d = 1 / n_points){

  ## now evaluate derivatives of smooths with associated standard
  ## errors, by finite differencing...
  cols = colnames(gam_fit$model)
  cols = cols[!cols %in% as.character(gam_fit$formula[[2]])]
  # space mesh within min and max of the data or create new x mesh
  min_val = min(gam_fit$model[, cols])
  max_val = max(gam_fit$model[, cols])
  if(is.null(x)) {
    x = seq(min_val, max_val, length = n_points)
  } else n_points = length(x)
  newd = matrix(rep(x, length(cols)), n_points, length(cols), byrow = FALSE)
  colnames(newd) = cols

  X0 = predict(gam_fit, as.data.frame(newd), type="lpmatrix")
  newd2 = newd + d ## shift the evaluation mesh
  X1 = predict(gam_fit, as.data.frame(newd2), type="lpmatrix")
  Xp = (X1-X0)/d ## maps coefficients to (fd approx.) derivatives

  N_param = length(gam_fit$coefficients)-1
  # substract 1 for intercept that is shared for all archtypes and removed in the next step
  N_smooths = N_smooths - 1
  N_covar = N_param/N_smooths
  derivs = lapply(seq(1, N_covar), function(i) {
    # compute derivatives and corresponding sd
    Xi = Xp*0
    # Xi%*%coef(b) = smooth deriv i
    Xi[, (i-1)*N_smooths+1:N_smooths+1] = Xp[, (i-1)*N_smooths+1:N_smooths+1]
    # ith smooth derivative
    df = Xi %*% coef(gam_fit)
    df_sd = rowSums(Xi %*% gam_fit$Vp * Xi)^.5
    # calculate the probability that derivative is negative
    prob_neg = pnorm(0, as.numeric(df), df_sd, lower.tail = TRUE)
    x_name = unique(gsub("s\\(|\\)\\.[[:digit:]]+", "",
                         colnames(Xp)[(i-1)*N_smooths+1:N_smooths+1]))
    data.table(deriv = as.numeric(df), deriv_sd = df_sd,
               x = as.numeric(newd[, x_name]),
               x_name = x_name,
               x_lab = paste0(x_name,",\ntotal p(der < 0): ", signif(mean(prob_neg), 2)),
               prob_neg = prob_neg, weights = weights)
  })
  derivs = list(call = match.call(),
                derivs = rbindlist(derivs))
  # add summary with mean probability for each archetype
  mean_prob = derivs$derivs[, .(p = weighted.mean(prob_neg, weights),
                                metric = "mean_prob"), by = .(x_name)]
  # add summary with product of probability for each archetype
  prod_prob = derivs$derivs[, .(p = prod(prob_neg), metric = "prod_prob"), by = .(x_name)]
  # calculate p_arc1(der < 0) and p_arc2(der > 0):
  derivs_temp = dcast.data.table(derivs$derivs, x ~ x_name, value.var = "prob_neg")
  summary = vapply(cols, function(arc, derivs_temp) {
    not_arc = cols[!cols == arc]
    derivs_temp[, prob_prod := prod(1 - .SD) * get(arc), by = x, .SDcols = not_arc]
    summary = c(derivs_temp[, weighted.mean(prob_prod, weights)],
                derivs_temp[, prod(prob_prod)])
    names(summary) = c("mean_prob_excl", "prod_prob_excl")
    summary
  }, FUN.VALUE = numeric(2), derivs_temp)
  summary = melt.data.table(as.data.table(summary, keep.rownames = "metric"),
                            id.vars = "metric", variable.name = "x_name", value.name = "p")
  # merge summaries together
  summary = rbindlist(list(mean_prob, prod_prob, summary), use.names = TRUE)
  summary$y_name = as.character(gam_fit$formula[[2]])
  derivs$summary = summary
  class(derivs) = "gam_deriv"
  derivs
}

##' @rdname gam_deriv
##' @name plot.gam_deriv
##' @export plot.gam_deriv
plot.gam_deriv = function(derivs) {
  n_der = uniqueN(derivs$derivs)
  der = derivs$derivs
  ggplot2::ggplot(der, ggplot2::aes(x, deriv)) +
    ggplot2::geom_line(ggplot2::aes(x, deriv)) +
    ggplot2::geom_line(ggplot2::aes(x, y = deriv + 2 * deriv_sd), linetype = 2) +
    ggplot2::geom_line(ggplot2::aes(x, y =  deriv - 2 * deriv_sd), linetype = 2) +
    ggplot2::ylim(min(der$deriv - 2 * der$deriv_sd),
                  max(der$deriv + 2 * der$deriv_sd)) +
    ggplot2::facet_wrap(~ x_name)
}
