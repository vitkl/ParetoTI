##' Find derivative of smooth curves from gam model fit (Generalised Additive Model, mgcv)
##' @rdname find_gam_deriv
##' @name find_gam_deriv
##' @description \code{find_gam_deriv()} Finds a derivative with standard errors of smooth curves from gam model fit by finite differencing. Crucially it computes the probability that derivative is negative at each point.
##' @details Plot method uses ggplot() to show mean derivative and confidence intervals (2 SD from the mean). This code is a modification of last example in \link[mgcv]{predict.gam}.
##' @param gam_fit a fitted gam object as produced by gam().
##' @param n_points number of points at which to evaluate derivative
##' @param x numeric vector where to evaluate derivatives
##' @param N_smooths number of bases used to represent the smooth term (\link[mgcv]{s}), 4 for cubic splines
##' @param weights how to weight points along x axis when calculating mean (integral) probability
##' @param d numeric vector (1L), finite difference interval
##' @param return_gam return gam model as well? By default only only the first derivative at n_points
##' @return \code{find_gam_deriv()} list (S3 object, find_gam_deriv) containing function call, a data.table with values of the first derivative, GAM model fit summary (gam_sm) and (optionally) gam model fit; summary entry is NA.
##' @export find_gam_deriv
##' @import data.table
find_gam_deriv = function(gam_fit,
                          n_points = 200, x = NULL,
                          N_smooths = 9,
                          weights = list(rep(1, 200), sqrt(seq(1, 0, length.out = 200)),
                                         1 / (1 + seq(1, 0, length.out = 200)^3))[[1]],
                          d = 1 / n_points, return_gam = FALSE){

  ## evaluate derivatives of smooths with associated standard
  ## errors, by finite differencing...
  # find column names for each archetype
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
  # Combine results of iterations over predictors of GAM model (distance from vertices)
  derivs = rbindlist(derivs)
  # generate summary of gam fit
  gam_sm = summary(gam_fit)

  # find function value at 0 and 1 distance and compute f(0) / f(1) ratio
  zero_one_d = as.data.frame(newd[c(1, n_points),])
  colnames(zero_one_d) = colnames(newd)
  zero_one = predict(gam_fit, zero_one_d, type="link")

  gam_sm = data.table(smooth_term_p_val = gam_sm$s.pv, r_sq =  gam_sm$r.sq,
                      y_name = as.character(gam_fit$formula[[2]]),
                      x_name = unique(derivs$x_name),
                      min_max_ratio = zero_one[1] / zero_one[2])
  # combine results into list are return
  if(!isTRUE(return_gam)) {
    derivs = list(call = match.call(),
                  derivs = derivs,
                  gam_fit = NA,
                  gam_sm = gam_sm, summary = NA)
  } else {
    derivs = list(call = match.call(),
                  derivs = derivs,
                  gam_fit = gam_fit,
                  gam_sm = gam_sm, summary = NA)
  }
  derivs$derivs$y_name = as.character(gam_fit$formula[[2]])
  class(derivs) = "gam_deriv"
  derivs
}

##' @rdname find_gam_deriv
##' @name summary.gam_deriv
##' @export summary.gam_deriv
##' @import data.table
summary.gam_deriv = function(derivs){
  summary = lapply(unique(derivs$derivs$y_name), function(yname){
    der = derivs$derivs[y_name %in% yname]

    # find average derivative values for 100%, top 50% and top 20% data points
    effect_size = der[, .(deriv100 = weighted.mean(deriv, weights),
                          deriv50 = weighted.mean(deriv,
                                                  weights * c(rep(1, .N/2),
                                                              rep(0, .N/2))),
                          deriv20 = weighted.mean(deriv,
                                                  weights * c(rep(1, .N/5),
                                                              rep(0, .N/5*4)))), by = .(x_name)]
    # add summary with mean probability for each archetype
    mean_prob = der[, .(p = weighted.mean(prob_neg, weights),
                        metric = "mean_prob"), by = .(x_name)]
    # add summary with product of probability for each archetype
    prod_prob = der[, .(p = prod(prob_neg), metric = "prod_prob"), by = .(x_name)]
    # calculate p(p_arc1(der < 0) and p_arc2(der > 0)):
    # generate a matrix of p_arc1(der < 0) and p_arc2(der > 0) ... p_arcN(der > 0)
    mean_prob_excl = copy(dcast.data.table(mean_prob, x_name ~ x_name, value.var = "p"))
    cols = unique(der$x_name)
    for (col in cols) {
      mean_prob_excl[is.na(get(col)), c(col) :=
                       mean_prob_excl[x_name == col, 1 - get(col)]]
    }
    mean_prob_excl = as.matrix(mean_prob_excl, rownames = "x_name")
    # multiply and average probabilities
    prod_prob_excl = apply(mean_prob_excl, 1, prod)
    mean_prob_excl = rowMeans(mean_prob_excl)
    prod_prob_excl = data.table(p = prod_prob_excl, x_name = names(prod_prob_excl),
                                metric = "prod_prob_excl")
    mean_prob_excl = data.table(p = mean_prob_excl, x_name = names(mean_prob_excl),
                                metric = "mean_prob_excl")
    # merge summaries together
    summary = rbindlist(list(mean_prob, prod_prob,
                             mean_prob_excl, prod_prob_excl), use.names = TRUE)
    summary$y_name = unique(der$y_name)
    merge(summary, effect_size, by = "x_name")
  })
  # rbind iterations over y_name
  summary = rbindlist(summary)
  # merge gam fit summary
  if(!isTRUE(is.na(derivs$gam_sm))) {
    summary = merge(summary, derivs$gam_sm, by = c("y_name", "x_name"), all = TRUE)
  }
  summary
}


##' @rdname find_gam_deriv
##' @name plot.gam_deriv
##' @export plot.gam_deriv
##' @import data.table
plot.gam_deriv = function(derivs, features = derivs$derivs$y_name,
                          title = "First derivative of GAM model: gene expression = function(distance from archetype)") {
  n_der = uniqueN(derivs$derivs$x_name)
  der = derivs$derivs
  der = der[y_name %in% features]
  ggplot2::ggplot(der, ggplot2::aes(x, deriv)) +
    ggplot2::geom_line(ggplot2::aes(x, deriv)) +
    ggplot2::geom_line(ggplot2::aes(x, y = deriv + 2 * deriv_sd), linetype = 2) +
    ggplot2::geom_line(ggplot2::aes(x, y =  deriv - 2 * deriv_sd), linetype = 2) +
    ggplot2::ylim(min(der$deriv - 2 * der$deriv_sd),
                  max(der$deriv + 2 * der$deriv_sd)) +
    ggplot2::facet_grid(y_name ~ x_name) +
    ggplot2::geom_hline(yintercept = 0, color = "grey80") +
    ggtitle(title)
}
