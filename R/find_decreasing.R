##' Find features that are decreasing functions of distance from archetype
##' @rdname find_decreasing
##' @name find_decreasing
##' @description \code{find_decreasing()} Fits gam models to find features that are a decreasing function of distance from archetype. Both gam functions and first derivatives can be visualised using plot() method.
##' @param data_attr data.table dim(examples, dimensions) that includes distance of each example to archetype in columns given by \code{arc_col} and feature values given by \code{features}
##' @param arc_col character vector, columns that give distance to archetypes (column per archetype)
##' @param features character vector (1L), column than containg feature values
##' @param min.sp lower bound for the smoothing parameter, details: \link[mgcv]{gam}. Default value of 60 works well to stabilise curve shape near min and max distance
##' @param N_smooths number of bases used to represent the smooth term (\link[mgcv]{s}), 4 for cubic splines
##' @param n_points number of points at which to evaluate derivative
##' @param d numeric vector (1L), finite difference interval
##' @param weights how to weight points along x axis when calculating mean (integral) probability. Useful if you care that the function is decreasing near the archetype but not far away. Two defaults suggest to weight point equally or discard bottom 50 percent.
##' @param return_only_summary return only summary data.table containing p-values for each feature at each archetype and effect-size measures (average derivative).
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
##' @param p.adjust.method choose method for correcting p-value for multiple hypothesis testing. See p.adjust.methods and \link[stats]{p.adjust} for details.
##' @param gam_fit_pval smooth term probability in gam fit (upper bound)
##' @param invert_cutoff invert cutoff for genes and sets. If FALSE p < cutoff_genes, if TRUE p > cutoff_genes.
##' @param order_by order decreasing feature list by measure in summary sets. By default is deriv20, the average value of derivative at 20 of point closest to vertex.
##' @param min_max_diff_cutoff_g what should be the mean difference (log-ratio, when y is log-space) of gene expression at the point closest to archetype compared to point furthest from archetype? When Wilcox method was used it is difference between mean of bin closest to archetype and all other cells. By default, at least 0.3 for genes 0.1 for functions.
##' @param min_max_diff_cutoff_f see min_max_diff_cutoff_g
##' @return \code{get_top_decreasing()} print summary to output, and return list with character vector with one element for each archetype, and 2 data.table- with selection of enriched genes and functions.
##' @export get_top_decreasing
##' @import data.table
get_top_decreasing = function(summary_genes, summary_sets = NULL,
                              cutoff_genes = 0.01, cutoff_sets = 0.01,
                              cutoff_metric = "wilcoxon_p_val",
                              p.adjust.method = c("fdr", "none")[1],
                              gam_fit_pval = 0.01, invert_cutoff = FALSE,
                              order_by = "mean_diff", order_decreasing = FALSE,
                              min_max_diff_cutoff_g = 0.3, min_max_diff_cutoff_f = 0.1) {

  enriched = summary_genes[metric == cutoff_metric]
  # do fdr correction of p-value
  enriched[, p := p.adjust(p, method = p.adjust.method)]
  # look at derivatives only when the model fit itself is good
  if("gam_fit_pval" %in% colnames(enriched)) {
    enriched = enriched[smooth_term_p_val < gam_fit_pval]
  }

  if(isTRUE(invert_cutoff)) {
    enriched = enriched[p >= cutoff_genes]
  } else {
    enriched = enriched[p < cutoff_genes]
  }

  # filter by difference between min and max value (GAM) or by log-ratio
  if("min_max_diff" %in% colnames(enriched)) { # min and max value of GAM
    if(isTRUE(invert_cutoff)) {
      enriched = enriched[min_max_diff < min_max_diff_cutoff_g]
    } else {
      enriched = enriched[min_max_diff > min_max_diff_cutoff_g]
    }
  } else if("mean_diff" %in% colnames(enriched)) { # wilcox test mean log-ratio
    if(isTRUE(invert_cutoff)) {
      enriched = enriched[mean_diff < min_max_diff_cutoff_g]
    } else {
      enriched = enriched[mean_diff > min_max_diff_cutoff_g]
    }
  }
  # add enriched genes and sets
  enriched_genes = copy(enriched) # copy filtered results
  setnames(enriched_genes, c("x_name", "y_name"), c("arch_name", "genes"))
  enriched_genes = enriched_genes[, rank_by_vertex := frank(-get(order_by)),
                                  by = .(arch_name)]
  setorder(enriched_genes, arch_name, rank_by_vertex)
  # generate labels
  enriched = enriched[order(get(order_by), decreasing = order_decreasing),
                      .(arch_name = x_name, y_name)]
  enriched[, arch_lab := paste0(paste0(y_name[1:4][!is.na(y_name[1:4])],
                                       collapse = ", "), "\n",
                                paste0(y_name[5:8][!is.na(y_name[5:8])],
                                       collapse = ", "), "\n",
                                paste0(y_name[9:12][!is.na(y_name[9:12])],
                                       collapse = ", ")),
           by = arch_name]
  enriched_lab = unique(enriched[, .(arch_name, arch_lab)])
  # get and merge enriched sets
  enriched_sets = NULL # set to null when nothing provided
  if(!is.null(summary_sets)) {
    enriched_sets = summary_sets[metric == cutoff_metric]
    # do fdr correction of p-value
    enriched_sets[, p := p.adjust(p, method = p.adjust.method)]
    # look at derivatives only when the model fit itself is good
    if("gam_fit_pval"  %in% colnames(enriched_sets)) {
      enriched_sets = enriched_sets[smooth_term_p_val < gam_fit_pval]
    }

    if(isTRUE(invert_cutoff)) {
      enriched_sets = enriched_sets[p >= cutoff_genes]
    } else {
      enriched_sets = enriched_sets[p < cutoff_genes]
    }

    # filter by difference between min and max value (GAM) or by log-ratio
    if("min_max_diff" %in% colnames(enriched_sets)) { # min and max value of GAM
      if(isTRUE(invert_cutoff)) {
        enriched_sets = enriched_sets[min_max_diff < min_max_diff_cutoff_f]
      } else {
        enriched_sets = enriched_sets[min_max_diff > min_max_diff_cutoff_f]
      }
    } else if("mean_diff" %in% colnames(enriched_sets)) { # wilcox test mean log-ratio
      if(isTRUE(invert_cutoff)) {
        enriched_sets = enriched_sets[mean_diff < min_max_diff_cutoff_f]
      } else {
        enriched_sets = enriched_sets[mean_diff > min_max_diff_cutoff_f]
      }
    }

    enriched_sets_lab = enriched_sets[order(get(order_by), decreasing = order_decreasing),
                                      .(arch_name = x_name, y_name_set = y_name)]
    enriched_lab = merge(enriched_lab, enriched_sets_lab, by = "arch_name",
                         all.x = T, all.y = T)
    enriched_lab[, arch_lab := paste0(arch_lab, "\n\n",
                                      paste0(y_name_set[1][!is.na(y_name_set[1])],
                                             collapse = ", "), "\n",
                                      paste0(y_name_set[2][!is.na(y_name_set[2])],
                                             collapse = ", "), "\n",
                                      paste0(y_name_set[3][!is.na(y_name_set[3])],
                                             collapse = ", ")),
                 by = arch_name]
    enriched_lab = unique(enriched_lab[, .(arch_name, arch_lab)])
  }

  # add vertex label
  enriched_lab[, arch_lab := paste0(arch_name, "\n\n", arch_lab), by = arch_name]
  # remove excessive empty lines
  while(sum(grepl("\n\n\n", enriched_lab$arch_lab))){
    enriched_lab$arch_lab = gsub("\n\n\n", "\n\n", enriched_lab$arch_lab)
  }

  for (i in seq_len(nrow(enriched_lab))) {
    cat(" -- ", enriched_lab$arch_lab[i], "\n\n", sep = " ")
  }
  list(lab = enriched_lab, enriched = enriched,
       enriched_genes = enriched_genes, enriched_sets = enriched_sets)
}


##' @rdname find_decreasing
##' @name find_decreasing_wilcox
##' @description \code{find_decreasing_wilcox()} find features that are a decreasing function of distance from archetype by finding features with highest value (median) in bin closest to vertex (1 vs all Wilcox test).
##' @param bin_prop proportion of data to put in bin closest to archetype
##' @param method how to find_decreasing_wilcox()? Use \link[BioQC]{wmwTest} or \link[stats]{wilcox.test}. BioQC::wmwTest can be up to 1000 times faster, so it is default.
##' @return \code{find_decreasing_wilcox()} data.table containing p-values for each feature at each archetype and effect-size measures (average difference between bins). When log(counts) was used mean_diff reflects log-fold change.
##' @export find_decreasing_wilcox
##' @import data.table
find_decreasing_wilcox = function(data_attr, arc_col,
                                  features = c("Gpx1", "Alb", "Cyp2e1", "Apoa2")[3],
                                  bin_prop = 0.1, na.rm = FALSE,
                                  type = c("s", "m", "cmq")[1],
                                  clust_options = list(),
                                  method = c("BioQC", "r_stats")[1]) {

  # find which cells are in bin closest to each archetype
  arch_bin = bin_cells_by_arch(data_attr, arc_col, bin_prop, return_names = FALSE)

  if(method == "BioQC"){
    ## create gene sets using bin_prop,
    ## define gene set as 0.1 points closest to archetype
    ## use BioQC::wmwTest to do Wilcoxon tests

    # extract relevant features to matrix (cells in rows, features in columns)
    feature_mat = as.matrix(data_attr[, c(features, "sample_id"), with = FALSE],
                            rownames = "sample_id")

    # run Wilcox tests
    decreasing = BioQC::wmwTest(x = feature_mat, indexList = arch_bin,
                                col = "GeneSymbol",
                                valType = c("p.greater"), simplify = TRUE)
    decreasing = as.data.table(decreasing, keep.rownames = "x_name")
    decreasing = melt.data.table(decreasing, id.vars = "x_name",
                                 value.name = "p", variable.name = "y_name")
    # find mean and median difference between bin closest to vertex vs other bins
    decreasing[, c("median_diff", "mean_diff", "top_bin_mean") := .({
      as.numeric(median(feature_mat[arch_bin[x_name][[1]], y_name], na.rm = na.rm) -
                   median(feature_mat[-arch_bin[x_name][[1]], y_name], na.rm = na.rm))
    }, {
      as.numeric(mean(feature_mat[arch_bin[x_name][[1]], y_name], na.rm = na.rm) -
                   mean(feature_mat[-arch_bin[x_name][[1]], y_name], na.rm = na.rm))
    },
    as.numeric(mean(feature_mat[arch_bin[x_name][[1]], y_name], na.rm = na.rm))),
    by = .(y_name, x_name)]

  } else if(method == "r_stats"){

    # extract relevant features to list
    feature_list = lapply(features, function(feature) {
      m = matrix(data_attr[, get(feature)], nrow = nrow(data_attr), ncol = 1)
      rownames(m) = data_attr$sample_id
      colnames(m) = feature
      m
    })

    # create clusters with provided options --------------------------------------
    if(type == "m"){
      # set default options or replace them with provided options
      default = list(cores = parallel::detectCores()-1, cluster_type = "PSOCK")
      default_retain = !names(default) %in% names(clust_options)
      options = c(default[default_retain], clust_options)

      # create and register cluster
      cl = parallel::makePSOCKcluster(2)
      doParallel::registerDoParallel(cl)
    } else if(type == "cmq") {
      # set defaults or replace them with provided options
      default = list(memory = 2000, template = list(), n_jobs = 5,
                     fail_on_error = FALSE, timeout = Inf)
      default_retain = !names(default) %in% names(clust_options)
      options = c(default[default_retain], clust_options)

      # register cluster
      clustermq::register_dopar_cmq(n_jobs = options$n_jobs,
                                    memory = options$memory,
                                    template = options$template,
                                    fail_on_error = options$fail_on_error,
                                    timeout = options$timeout)
    }

    # run wilcox test for features -----------------------------------------------
    # define how to iterate over feature list with foreach syntax
    fr_obj = foreach::foreach(feature_mat = feature_list,
                              .combine = rbind)

    # run tests
    if(type %in% c("m", "cmq")){
      # multiple cores locally or computing cluster with clustermq
      decreasing = foreach::`%dopar%`(fr_obj,
                                      ParetoTI:::.find_decreasing_wilcox_1(feature_mat, arch_bin,
                                                                           arc_col, na.rm))
      if(type == "m") parallel::stopCluster(cl) # stop local cluster
    } else {
      # single core locally
      decreasing = foreach::`%do%`(fr_obj,
                                   ParetoTI:::.find_decreasing_wilcox_1(feature_mat, arch_bin,
                                                                        arc_col, na.rm))
    }
  }
  setorder(decreasing, x_name, p)
  decreasing$metric = "wilcoxon_p_val"
  decreasing[,.(x_name, y_name, p, median_diff, mean_diff, top_bin_mean, metric)]
}

.find_decreasing_wilcox_1 = function(feature_mat, arch_bin,
                                     arc_col, na.rm) {

  decreasing = lapply(arch_bin, function(arc, feature_mat) {
    # select 1st bin cells with indices
    y1 = feature_mat[arc, 1]
    y0 = feature_mat[-arc, 1]
    y1_mean = mean(y1, na.rm = na.rm)

    data.table(p = wilcox.test(x = y1, y = y0,
                               alternative = "greater")$p.value,
               median_diff = as.numeric(median(y1, na.rm = na.rm) -
                                          median(y0, na.rm = na.rm)),
               mean_diff = as.numeric(y1_mean -
                                        mean(y0, na.rm = na.rm)),
               top_bin_mean = as.numeric(y1_mean))
  }, feature_mat)
  decreasing = rbindlist(decreasing)
  decreasing$x_name = names(arch_bin)
  decreasing$y_name = colnames(feature_mat)
  decreasing
}

##' @rdname find_decreasing
##' @name bin_cells_by_arch
##' @description \code{bin_cells_by_arch()} find which cells are in bin closest to vertex.
##' @param return_names return list of indices of cells or names of cells?
##' @return \code{bin_cells_by_arch()} list of indices of cells or names of cells that are in bin closest to each vertex
##' @export bin_cells_by_arch
##' @import data.table
bin_cells_by_arch = function(data_attr, arc_col, bin_prop = 0.1, return_names = FALSE){

  # extract distance to archetypes into matrix,
  dist_to_arch = as.matrix(data_attr[, arc_col, with = FALSE])
  # convert distances to order
  dist_to_arch = apply(dist_to_arch, MARGIN = 2, order, decreasing = FALSE)

  # find how many points fit into bin closest to archetype
  bin1_length = round(nrow(dist_to_arch) * bin_prop, 0)
  # pick indices of cells in 1st bin
  dist_to_arch = dist_to_arch[seq_len(bin1_length), ]
  arch_bin = lapply(seq(1, ncol(dist_to_arch)), function(i) dist_to_arch[,i])
  names(arch_bin) = colnames(dist_to_arch)

  # optionally: convert indices to cells
  if(return_names) arch_bin = lapply(arch_bin, function(arch_ind) {
    data_attr$sample_id[arch_ind]
  })
  arch_bin
}
