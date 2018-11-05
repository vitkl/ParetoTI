##' Plot gam model fit for multiple genes and archetypes (Generalised Additive Model, mgcv)
##' @rdname plot_gam
##' @name plot_gam
##' @description \code{plot_gam()} plots gam model fit for multiple features (genes) and archetypes
##' @details Plotting of multiple gam fits is inspired by https://stackoverflow.com/questions/49471300/gam-plots-with-ggplot
##' @param gam_deriv object with multiple gam fits (many arhetypes, many features). Either gam_deriv or gam_fit must be provided.
##' @param gam_fit a fitted gam object as produced by gam().
##' @param data data.table dim(examples, dimensions) that includes distance of each example to archetype in columns given by \code{arc_col} and feature values given by \code{feature}. Must be provided.
##' @param feature character vector (1L), column than containg feature values
##' @param groupCovs argument for \link[voxel]{plotGAM}
##' @return \code{plot_gam()} list (S3 object, plot_gam) containing data.table-s for each predictor in gam model
##' @export plot_gam
##' @import data.table
plot_gam = function(gam_deriv = NULL, gam_fit, data, feature = NULL,
                    groupCovs = NULL, title_size = 10){
  if(!is.null(gam_deriv)){
    ns = seq_len(length(gam_deriv$gam_fit))
    plots = lapply(ns, function(n){
      gam_fit = gam_deriv$gam_fit[[n]]
      # find column names for each archetype
      cols = colnames(gam_fit$model)
      cols = cols[!cols %in% as.character(gam_fit$formula[[2]])]
      lapply(cols, function(col){
        voxel::plotGAM(gam_fit, smooth.cov = col, groupCovs = groupCovs) +
          geom_point(data = data, aes_string(y = as.character(gam_fit$formula[[2]]),
                                             x = col), alpha = 0.2) +
          geom_rug(data = data, aes_string(y = as.character(gam_fit$formula[[2]]),
                                           x = col), alpha = 0.2) +
          theme(legend.position = "none",
                title = element_text(size = title_size))
      })
    })
    plots = unlist(plots, recursive = FALSE)
    features = vapply(ns, function(n){
      as.character(gam_deriv$gam_fit[[n]]$formula[[2]])
    }, character(1))
    plots = plots[order(features)]
    features = features[order(features)]
    if(!is.null(feature)) {
      plots = plots[features %in% feature]
      features = features[features %in% feature]
    }
    cowplot::plot_grid(plotlist = plots, nrow = uniqueN(features))
  } else {
    # find column names for each archetype
    cols = colnames(gam_fit$model)
    cols = cols[!cols %in% as.character(gam_fit$formula[[2]])]
    voxel::plotGAM(gam_fit, smooth.cov = cols, groupCovs = groupCovs) +
      geom_point(data = data, aes_string(y = as.character(gam_fit$formula[[2]]),
                                         x = cols), alpha = 0.2) +
      geom_rug(data = data, aes_string(y = as.character(gam_fit$formula[[2]]),
                                       x = cols), alpha = 0.2)
  }
}

