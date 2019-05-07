##' Plot gam model fit for multiple genes and archetypes (Generalised Additive Model, mgcv)
##' @rdname plot_gam
##' @name plot_gam
##' @description \code{plot_gam()} plots gam model fit for multiple features (genes) and archetypes
##' @details Plotting of multiple gam fits is inspired by https://stackoverflow.com/questions/49471300/gam-plots-with-ggplot
##' @param gam_deriv object with multiple gam fits (many arhetypes, many features). Either gam_deriv or gam_fit must be provided.
##' @param gam_fit a fitted gam object as produced by gam().
##' @param data data.table dim(examples, dimensions) that includes distance of each example to archetype in columns given by \code{arc_col} and feature values given by \code{feature}. Must be provided.
##' @param feature character vector, which features (e.g. genes) to plot?
##' @param archetype character vector, which archetypes to plot?
##' @param groupCovs argument for \link[voxel]{plotGAM}
##' @param title_size size of individial plot titles and axis titles
##' @param point_size/rug_size/point_alpha/rug_alpha size and alpha controls for points and geom_rug.
##' @param point_position point position option, a way to add jitter but "identity" by default
##' @param title main title for all features and archetypes
##' @return \code{plot_gam()} list (S3 object, plot_gam) containing data.table-s for each predictor in gam model
##' @export plot_gam
##' @import data.table
plot_gam = function(gam_deriv = NULL, gam_fit, data,
                    feature = NULL, archetype = NULL,
                    groupCovs = NULL, title_size = 10,
                    point_size = 1, point_alpha = 0.1,
                    point_position = c("identity", "jitter"),
                    rug_size = 1, rug_alpha = 0.1,
                    title = "GAM model: gene expression = function(distance from archetype)"){
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
                                             x = col),
                     alpha = point_alpha, size = point_size,
                     position = point_position[1]) +
          geom_rug(data = data, aes_string(y = as.character(gam_fit$formula[[2]]),
                                           x = col),
                   alpha = rug_alpha, size = rug_size) +
          theme(legend.position = "none",
                title = element_text(size = title_size))
      })
    })

    plots = unlist(plots, recursive = FALSE)

    # Find and filter features & archetypes
    features = vapply(ns, function(n){
      as.character(gam_deriv$gam_fit[[n]]$formula[[2]])
    }, character(1))
    archetypes = vapply(ns, function(n){
      as.character(gam_deriv$gam_fit[[n]]$formula[[3]][[2]])
    }, character(1))
    # reorder by feature name
    order_features = order(features)
    features = features[order_features]
    archetypes = archetypes[order_features]
    plots = plots[order_features]
    if(!is.null(feature)) {
      if(!is.null(archetype)) {
        # filter archetypes if provided
        plots = plots[features %in% feature &
                        archetypes %in% archetype]
        features = features[features %in% feature &
                              archetypes %in% archetype]
      } else {
        plots = plots[features %in% feature]
        features = features[features %in% feature]
      }
    }

    # combine into one plot
    end_plot = cowplot::plot_grid(plotlist = plots, nrow = uniqueN(features))

  } else {
    # find column names for each archetype
    cols = colnames(gam_fit$model)
    cols = cols[!cols %in% as.character(gam_fit$formula[[2]])]
    end_plot = voxel::plotGAM(gam_fit, smooth.cov = cols, groupCovs = groupCovs) +
      geom_point(data = data, aes_string(y = as.character(gam_fit$formula[[2]]),
                                         x = cols),
                 alpha = point_alpha, size = point_size,
                 position = point_position[1]) +
      geom_rug(data = data, aes_string(y = as.character(gam_fit$formula[[2]]),
                                       x = cols),
               alpha = rug_alpha, size = rug_size)
  }

  if(title == "") return(end_plot) # if title empty return ggplot output
  # otherwise add the main plot title
  title = cowplot::ggdraw() +
    cowplot::draw_label(title,
                        fontface = "bold")
  cowplot::plot_grid(title, end_plot,
                     ncol = 1, rel_heights = c(0.08, 1))
}

