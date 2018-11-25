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
##' @param type one of s, m, cmq. s means single core processing using lapply. m means multi-core parallel procession using parLapply. cmq means multi-node parallel processing on a computing cluster using clustermq package.
##' @param clust_options list of options for parallel processing. The default for "m" is list(cores = parallel::detectCores()-1, cluster_type = "PSOCK"). The default for "cmq" is list(memory = 2000, template = list(), n_jobs = 10, fail_on_error = FALSE). Change these options as required.
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
                           type = c("s", "m", "cmq")[1], clust_options = list(),
                           ...){
  if(length(features) > 10 & isTRUE(stop_at_10) & !isTRUE(return_only_summary)) stop("Trying to fit more than 10 features requesting values of 1st derivative. return_only_summary = FALSE option is designed for visualising few features rather than bulk processing. Please use return_only_summary = TRUE or set stop_at_10 = FALSE")
  # start loop

  # single process -------------------------------------------------------------
  if(type == "s"){
    res = lapply(seq_len(length(features)), .find_decreasing_1,
                 features = features, arc_col = arc_col, N_smooths = N_smooths,
                 data_attr = data_attr, min.sp = min.sp, ..., d = d,
                 n_points = n_points, weights = weights,
                 return_only_summary = return_only_summary,
                 one_arc_per_model = one_arc_per_model)
  }
  # multi-process --------------------------------------------------------------
  if(type == "m"){
    # set defaults or replace them with provided options
    default = list(cores = parallel::detectCores()-1, cluster_type = "PSOCK")
    default_retain = !names(default) %in% names(clust_options)
    options = c(default[default_retain], clust_options)
    # create cluster
    cl = parallel::makeCluster(options$cores, type = options$cluster_type)
    # get library support needed to run the code
    parallel::clusterEvalQ(cl, {library(ParetoTI)})

    res = parallel::parLapply(cl, seq_len(length(features)),
                              ParetoTI:::.find_decreasing_1,
                              features = features, arc_col = arc_col,
                              N_smooths = N_smooths,
                              data_attr = data_attr, min.sp = min.sp, ..., d = d,
                              n_points = n_points, weights = weights,
                              return_only_summary = return_only_summary,
                              one_arc_per_model = one_arc_per_model)
    # stop cluster
    parallel::stopCluster(cl)
  }
  # clustermq ------------------------------------------------------------------
  if(type == "cmq"){
    # set defaults or replace them with provided options
    default = list(memory = 2000, template = list(), n_jobs = 10,
                   fail_on_error = FALSE, timeout = Inf)
    default_retain = !names(default) %in% names(clust_options)
    options = c(default[default_retain], clust_options)

    # run analysis
    suppressWarnings({ # hide "NA introduced by coersion" warning specific to cmq implementation
      suppressPackageStartupMessages({ # hide package startup warnings on each cluster
        res = clustermq::Q(fun = ParetoTI:::.find_decreasing_1,
                           i = seq_len(length(features)),
                           const = list(features = features, arc_col = arc_col,
                                        N_smooths = N_smooths,
                                        data_attr = data_attr, min.sp = min.sp,
                                        ..., d = d,
                                        n_points = n_points, weights = weights,
                                        return_only_summary = return_only_summary,
                                        one_arc_per_model = one_arc_per_model),
                           memory = options$memory, template = options$template,
                           n_jobs = options$n_jobs, rettype = "list",
                           fail_on_error = options$fail_on_error,
                           timeout = options$timeout)
      })
    })
  }

  # combine results ------------------------------------------------------------
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
  derivs = ParetoTI::find_gam_deriv(gam_fit, N_smooths = N_smooths, d = d,
                                    n_points = n_points, weights = weights,
                                    return_gam = FALSE)
  derivs = list(call = derivs$call, derivs = derivs$derivs,
                gam_fit = gam_fit, gam_sm = derivs$gam_sm,
                summary = NA)
  class(derivs) = "gam_deriv"
  derivs
}

.find_decreasing_1 = function(i, features, arc_col, N_smooths, data_attr, min.sp,
                              ..., d, n_points, weights, return_only_summary,
                              one_arc_per_model){
  feature = features[i]
  if(one_arc_per_model) {
    # if fit a model for each archetype, iterate over archetypes
    derivs = lapply(arc_col, function(col){
      ParetoTI::fit_arc_gam_1(feature = feature, col = col, N_smooths = N_smooths,
                              data_attr = data_attr, min.sp = min.sp, ..., d = d,
                              n_points = n_points, weights = weights)
    })
    # combine results
    res = list()
    res$call = match.call()
    res$derivs = data.table::rbindlist(lapply(derivs, function(x) x$derivs))
    res$gam_fit = lapply(derivs, function(x) x$gam_fit)
    res$gam_sm = data.table::rbindlist(lapply(derivs, function(x) x$gam_sm))
    class(res) = "gam_deriv"
  } else {
    res = ParetoTI::fit_arc_gam_1(feature = feature, col = arc_col, N_smooths = N_smooths,
                                  data_attr = data_attr, min.sp = min.sp, ..., d = d,
                                  n_points = n_points, weights = weights)
    res$gam_fit = list(res$gam_fit)
  }
  # generate summary of the derivative
  summary = ParetoTI::summary.gam_deriv(res)
  if(return_only_summary) return(summary) else {
    list(summary = summary, derivs = res$derivs, gam_fit = res$gam_fit)
  }
}


##' @rdname find_decreasing
##' @name get_top_decreasing
##' @description \code{get_top_decreasing()} Get top-12 genes and top-3 gene sets for each archetype.
##' @param summary_genes gam_deriv summary data.table for decreasing genes
##' @param summary_sets gam_deriv summary data.table for decreasing gene sets
##' @param cutoff_genes value of cutoff_metric (lower bound) for genes
##' @param cutoff_sets value of cutoff_metric (lower bound) for gene sets
##' @param cutoff_metric probability metric for selecting decreasing genes: mean_prob, prod_prob, mean_prob_excl or prod_prob_excl
##' @param gam_fit_pval smooth term probability in gam fit (upper bound)
##' @param order_by order decreasing feature list by measure in summary sets. By default is deriv20, the average value of derivative at 20 of point closest to vertex.
##' @param min_max_ratio what should be the ratio of gene expression at the point closest to archetype compared to point furthest from archetype? By default, at least 1.3 or closest point is 30 percent higher than furthest point.
##' @return \code{get_top_decreasing()} return character vector with one element for each archetype and print to output
##' @export get_top_decreasing
##' @import data.table
get_top_decreasing = function(summary_genes, summary_sets = NULL,
                              cutoff_genes = 0.99, cutoff_sets = 0.99,
                              cutoff_metric = "mean_prob",
                              gam_fit_pval = 0.01,
                              order_by = "deriv20", order_decreasing = FALSE,
                              min_max_ratio_cutoff = 1.3) {

  enriched = summary_genes[metric == cutoff_metric][p > cutoff_genes &
                                                      smooth_term_p_val < gam_fit_pval &
                                                      min_max_ratio > min_max_ratio_cutoff]
  # add enriched genes and sets
  enriched = enriched[order(get(order_by), decreasing = order_decreasing),
                      .(arch_name = x_name, y_name)]
  enriched[, arch_lab := paste0(arch_name, "\n\n",
                                paste0(y_name[1:4][!is.na(y_name[1:4])],
                                       collapse = ", "), "\n",
                                paste0(y_name[5:8][!is.na(y_name[5:8])],
                                       collapse = ", "), "\n",
                                paste0(y_name[9:12][!is.na(y_name[9:12])],
                                       collapse = ", ")),
           by = arch_name]
  enriched = unique(enriched[, .(arch_name, arch_lab)])
  # get and merge enriched sets
  if(!is.null(summary_sets)) {
    enriched_sets = summary_sets[metric == cutoff_metric][p > cutoff_sets &
                                                            smooth_term_p_val < gam_fit_pval &
                                                            min_max_ratio > min_max_ratio_cutoff]

    enriched_sets = enriched_sets[order(get(order_by), decreasing = order_decreasing),
                                  .(arch_name = x_name, y_name_set = y_name)]
    enriched = merge(enriched, enriched_sets, by = "arch_name",
                     all.x = T, all.y = F)
    enriched[, arch_lab := paste0(arch_lab, "\n\n",
                                  paste0(y_name_set[1][!is.na(y_name_set[1])],
                                         collapse = ", "), "\n",
                                  paste0(y_name_set[2][!is.na(y_name_set[2])],
                                         collapse = ", "), "\n",
                                  paste0(y_name_set[3][!is.na(y_name_set[3])],
                                         collapse = ", ")),
             by = arch_name]
    enriched = unique(enriched[, .(arch_name, arch_lab)])
    while(sum(grepl("\n\n\n", enriched$arch_lab))){
      enriched$arch_lab = gsub("\n\n\n", "\n\n", enriched$arch_lab)
    }
  }
  for (i in seq_len(nrow(enriched))) {
    cat(" -- ", enriched$arch_lab[i], "\n\n", sep = " ")
  }
  enriched
}
