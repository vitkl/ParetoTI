##' Plot gam model fit for multiple genes and archetypes (Generalised Additive Model, mgcv)
##' @rdname plot_gam
##' @name plot_gam
##' @description \code{plot_gam()} Finds a derivative with standard errors of smooth curves from gam model fit by finite differencing. Crucially it computes the probability that derivative is negative at each point.
##' @details Plotting of multiple gam fits is inspired by https://stackoverflow.com/questions/49471300/gam-plots-with-ggplot
##' @param gam_fit a fitted gam object as produced by gam().
##' @param x numeric vector where to evaluate derivatives
##' @param N_smooths number of curves for each predictor in the model
##' @param d numeric vector (1L), finite difference interval
##' @return \code{plot_gam()} list (S3 object, plot_gam) containing data.table-s for each predictor in gam model
##' @export plot_gam
##' @import data.table
plot_gam = function(gam_fit, data, y = "y", smooth.cov = "x", groupCovs = NULL,
                    color = NULL, color_val = c("#868686FF", "#0073C2FF")){
  voxel::plotGAM(gam_fit, smooth.cov = smooth.cov, groupCovs = groupCovs) +
    geom_point(data = data, aes_string(y = y, x = smooth.cov), alpha = 0.2) + # , color = color
    geom_rug(data = data, aes_string(y = y, x = smooth.cov), alpha = 0.2) #+
    #scale_color_manual(color, values = color_val)
}

#plot_gam(gam_fit, data_gam, y = "y", smooth.cov = "x", groupCovs = NULL)

