##' Find the smallest polytope (Principal Convex Hull) that contains most of the data
##' @rdname fit_pch
##' @name fit_pch
##' @author Vitalii Kleshchevnikov
##' @description \code{fit_pch()} fits a polytope (Principal Convex Hull) to data using PCHA algorithm. All of the listed functions take input matrix dim(variables/dimentions, examples) and return archetype positions (XC) of dim(variables/dimentions, archetypes).
##' @details \code{fit_pch()} provides an R interface to python implementation of PCHA algorithm (Principal Convex Hull Analysis) by Ulf Aslak (https://github.com/ulfaslak/py_pcha) which was originally developed for Archetypal Analysis by Mørup et. al.
##' @param data numeric matrix in which to find archetypes, dim(variables/dimentions, examples)
##' @param noc integer, number of archetypes to find
##' @param I vector, entries of data to use for dictionary in C (optional)
##' @param U vector, entries of data to model in S (optional)
##' @param delta parameter that inflates original polytope(simplex) fit such that it may contain more points of the dataset
##' @param verbose if TRUE display messages
##' @param conv_crit The convergence criteria (default: 10^-6 relative change in SSE). Reduce to 1e-4 or 1e-3 for reduced computation time by more approximate results. 1e-4 gives very similar results to 1e-5 or 1e-6 on simulated example and STARmap 1020 gene data.
##' @param maxiter maximum number of iterations (default: 500 iterations)
##' @param check_installed if TRUE, check if python module py_pcha is found. Useful to set to FALSE for running analysis or within other functions
##' @param order_by integer, dimensions to be used for ordering vertices/archetypes. Vertices are ordered by angle (cosine) between c(1, 1) vector and a vector pointing to that vertex. Additional step finds when vertex vector is to the left (counter-clockwise) of the c(1, 1) vector. When bootstraping vertices can be aligned to reference and ordered (order_type == "align") by these dimensions.
##' @param order_type order archetypes by: cosine distance from c(1,1, ..., 1) vector ("cosine"), dot product that measures position to each side of the c(1,1) vector ("side"), align positions to reference when bootstraping(fit using all data, "align"). See \link[ParetoTI]{align_arc} When order_type is "align" vertices are ordered by "cosine" first.
##' @param volume_ratio find volume of the convex hull of the data and the t-ratio ("t_ratio") or ratio of variance of vertex positions to variance of data ("variance_ratio") or do nothing ("none")? Caution! Geometric figure should be at least simplex: qhull algorhirm will fail to find convex hull of flat 2D shapes in 3D, 3D shapes in 4D and so on. So, for calculating this dimentions are selected based order of rows in data. Makes sense for principal components but not for original data. Caution 2! Computation time and memory use increse very quickly with dimensions. Do not use for more than 7-8 dimentions.
##' @param converge_else_fail throw an error and stop execution if PCHA did not converge in \code{maxiter} steps.
##' @return \code{fit_pch()}: object of class pch_fit (list) containing the following elements:
##' XC - numeric matrix, dim(I, noc)/dim(dimensions, archetypes) feature matrix (i.e. XC=data[,I]*C forming the archetypes);
##' S - numeric matrix, dim(noc, length(U)) matrix, S>=0 |S_j|_1=1;
##' C - numeric matrix, dim(noc, length(U)) matrix, S>=0 |S_j|_1=1;
##' SSE - numeric vector (1L), Sum of Squared Errors;
##' varexpl - numeric vector (1L), Percent variation explained by the model.
##' hull_vol - numeric vector (1L), Volume of convex hull of the data.
##' arc_vol - numeric vector (1L), Volume of polytope enclosed by vertices.
##' t_ratio - numeric vector (1L), Ratio of \code{arc_vol} to \code{hull_vol}
##' var_vert - numeric vector (noc L), Variance in position of each vertex.
##' var_dim - numeric vector (noc L), Variance in positions of vertices in each dimension.
##' total_var - numeric vector (1L), Mean variance in position of all vertices.
##' summary - data.table with varexpl, t_ratio, total_var and noc for ease of combining multiple shape fits.
##' call - function call.
##' @export fit_pch
##' @export py_PCHA
##' @seealso \code{\link[parallel]{parLapply}}, \code{\link[base]{lapply}}, \code{\link[clustermq]{Q}}
##' @examples
##' set.seed(4355)
##' archetypes = generate_arc(arc_coord = list(c(5, 0), c(-10, 15), c(-30, -20)),
##'                           mean = 0, sd = 1, N_dim = 2)
##' data = generate_data(archetypes, N_examples = 1e4, jiiter = 0.04, size = 0.9)
##' dim(data)
##' # Fit a polytope with 3 vertices to data matrix
##' arc = fit_pch(data, noc=as.integer(3), delta=0)
##' # Fit the same polytope 3 times without resampling to test convergence of the algorithm.
##' arc_rob = fit_pch_bootstrap(data, n = 3, sample_prop = NULL,
##'                          noc=as.integer(3), delta=0)
##' # Fit the 10 polytopes to resampled datasets each time looking at 70% of examples.
##' arc_data = fit_pch_bootstrap(data, n = 10, sample_prop = 0.7,
##'                          noc=as.integer(3), delta=0)
##'
##' # Use local parallel processing to fit the 10 polytopes to resampled datasets each time looking at 70% of examples.
##' arc_data = fit_pch_bootstrap(data, n = 10, sample_prop = 0.7,
##'                          noc=as.integer(3), delta=0, type = "m")
##'
##' # Fit polytopes with 2-4 vertices
##' arc_ks = k_fit_pch(data, ks = 2:4, check_installed = T, delta=0)
##'
##' # Evaluate how much vertices vary in randomised data, (variable shuffled
##' # without replacement)
##' rand3 = randomise_fit_pch1(i = 1, data, true_fit = NULL,
##'     replace = FALSE, bootstrap_N = 200, seed = 2543,
##'     return_data = T, return_arc = T, sample_prop = 0.65,
##'     order_type = "align", noc = as.integer(3),
##'     delta = 0.1, bootstrap_type = "cmq")
##' p = plot_arc(arc_data, data, which_dimensions = 1:2, line_size = 1)
##' # plot shapes of observed data and randomised data side-by-side using cowplot
##' library(cowplot)
##' # create compound plot
##' p_all = plot_grid(plotlist = list(p + theme(legend.position = "none"),
##'                                   plot_arc(rand3$arc_data, rand3$data,
##'                                     which_dimensions = 1:2, line_size = 1) +
##'                                     theme(legend.position = "none")))
##' legend = get_legend(p)
##' # add legend to plot
##' plot_grid(p_all, legend, rel_widths = c(3.2, .6))
fit_pch = function(data, noc = as.integer(3), I = NULL, U = NULL,
                   delta = 0, verbose = FALSE, conv_crit = 1e-6,
                   maxiter = 500, check_installed = T,
                   order_by = seq(1, nrow(data)),
                   order_type = c("cosine", "side", "align")[3],
                   volume_ratio = c("t_ratio", "variance_ratio", "none")[1],
                   converge_else_fail = TRUE,
                   var_in_dims = FALSE, normalise_var = TRUE) {
  if(check_installed) .py_pcha_installed()
  # run PCHA
  res = tryCatch({
    res = py_PCHA$PCHA(X = data, noc = as.integer(noc), I = I, U = U,
                       delta = delta, verbose = verbose,
                       conv_crit = conv_crit, maxiter = maxiter)
    names(res) = c("XC", "S", "C", "SSE", "varexpl")
    if(!is.null(rownames(data))) rownames(res$XC) = rownames(data)
    class(res) = "pch_fit"
    res
  }, error = function(err) {
    if(isTRUE(converge_else_fail)) stop(paste0("fit_pch(noc = ",noc,") error: ", err))
    return(NULL)
  })
  if(is.null(res)) return(NULL)

  # step sorting archetypes to improve reproducibility
  # (archetypes are inherently exchangeable)
  if(noc > 1) {
    XC2 = res$XC[order_by,]
    arch_order = .find_vertex_order(XC2, noc, order_type = order_type)
    res$XC = res$XC[, arch_order]
    res$S = res$S[arch_order, ]
    res$C = res$C[arch_order, ]
  } else {
    # when only one vertex make sure data is still in the matrix form
    res$XC = matrix(res$XC, length(res$XC), 1)
    res$S = matrix(res$S, 1, length(res$S))
    res$C = matrix(res$C, 1, 1)
  }

  # Calculate ratio of volumes to measure goodness of fit ----------------------
  # 1 t_ratio  -----------------------------------------------------------------
  # when calculating convex hull and volume of the polytope
  # adjust number of dimensions to noc
  data_dim = seq(1, noc-1)
  if(volume_ratio == "t_ratio" & nrow(data) >= length(data_dim) & noc > 2){
    # calculate volume or area of the polytope only when number of vertices (noc) > number of dimenstions which means
    # find volume of the convex hull of the data
    hull_vol = fit_convhulln(data[data_dim, ], positions = FALSE)

    # Calcute the volume of non-simplex polytope (4 vertices in 2D, 5 in 3D or more)
    if(F){
      # find volume of the polytope fit, qhull requires at least 4 points in 2D
      # and more in higher dimensions so we have to add 20 points within
      # PCHA-fit polytope - more testing needed.
      # Currently program logic allows only for simplex polytopes
      archetypes = res$XC[data_dim, ]
      data_arc = cbind(archetypes, generate_data(archetypes, N_examples = 20,
                                                 jiiter = 0, size = 1))
      arc_vol = fit_convhulln(data_arc, positions = FALSE)
      res$arc_vol = arc_vol$vol
    } else {
      # Calculate the volume of simplex polytope
      archetypes = res$XC[data_dim, ]
      arch_red = archetypes - archetypes[, noc]
      res$arc_vol = abs(det(arch_red[, seq_len(noc - 1)]) /
                          factorial(noc - 1))
    }

    res$hull_vol = hull_vol$vol
    res$t_ratio = res$arc_vol / res$hull_vol

  } else if(volume_ratio == "variance_ratio") {

    # 2 variance ratio ---------------------------------------------------------
    # (variance of vertex positions within each dimension /
    #      variance of data in that dimension summed across dimensions)
    res = ParetoTI:::.cacl_var_in_dims(res, data,
                                       var_in_dims = T, normalise_var = T)

    dim_names = colnames(res$var_dim)[!colnames(res$var_dim) %in% "k"]
    res$t_ratio = mean(as.matrix(res$var_dim[, dim_names, with = FALSE]))

    res$hull_vol = NA
    res$arc_vol = NA

  } else {

    # 3 none

    if(volume_ratio == "t_ratio" & isTRUE(converge_else_fail)) message(paste0("Convex hull and t-ratio not computed for noc: ", noc," and nrow(data) = ", nrow(data),". fit_pch() can calculate volume or area of the polytope only when\nthe number of vertices (noc) > the number of dimensions (when polytope is convex):\ncheck that noc > nrow(data),\nselect only revelant dimensions or increase noc"))
    res$hull_vol = NA
    res$arc_vol = NA
    res$t_ratio = NA
  }
  res$var_vert = NA

  if(is.null(res$var_dim)){
    # calculate variance in positions only if it was not calculated before
    # (at variance_ratio step)
    res = ParetoTI:::.cacl_var_in_dims(res, data, var_in_dims, normalise_var)
  }

  res$total_var = NA
  # add summary table:
  res$summary = data.table::data.table(varexpl = as.numeric(res$varexpl),
                                       t_ratio = as.numeric(res$t_ratio),
                                       total_var = as.numeric(res$total_var),
                                       k = noc)
  res$call = match.call()
  res
}

##' @rdname fit_pch
##' @name k_fit_pch
##' @description \code{k_fit_pch()} finds polytopes of k dimensions in the data. This function applies \code{fit_pch()} to different k-s.
##' @param ks integer vector, dimensions of polytopes to be fit to data
##' @param bootstrap \code{k_fit_pch()}: use bootstrap to find average positions for each k? Also returns variability in vertex position to aid the selection of k. At excessive k position vary more.
##' @param simplex when testing multiple k using \code{k_fit_pch()} match dimensions to the number of archetypes? Use only on ordered principal components. If FALSE try all k for all dimensions in data. If simplex == TRUE test only simplex shapes (k=3 in 2D, k=4 in 3D, k=5 in 4D...). This assumes order of columns which is valid for principal components but may not be valid for untransformed data.
##' @return \code{k_fit_pch()}: object of class k_pch_fit (list) containing the same elements as pch_fit, but each is either a list of pch_fit elements (e.g. list of ks number of XC matrices) or a vector (which pch_fit element is one number). When length(ks) = 1 returns pch_fit.
##' @import clustermq
##' @export k_fit_pch
k_fit_pch = function(data, ks = 2:4, check_installed = TRUE,
                     bootstrap = FALSE, bootstrap_N = 10, sample_prop = 0.65,
                     bootstrap_type = c("s", "m", "cmq")[1],
                     seed = 345, simplex = FALSE,
                     var_in_dims = FALSE, normalise_var = TRUE, ...) {
  if(check_installed) .py_pcha_installed()
  # check that ks do not exceed dimensions when simplex is true
  if(nrow(data) < (max(ks) - 1) & isTRUE(simplex)) stop("simplex = TRUE but number of vertices (",
                                                        max(ks),") exceeds number of dimensions - 1 (", nrow(data),")")
  # run analysis for all k -----------------------------------------------------
  if(isTRUE(bootstrap)){
    res = lapply(ks, function(k) {
      # adjust number of dimensions to k
      if(isTRUE(simplex)) {
        data_dim = seq(1, k-1)
        if(k %in% c(1,2)) data_dim = seq(1, 2)
      } else data_dim = seq(1, nrow(data))
      # fit models
      fit_pch_bootstrap(data[data_dim,], noc = k, n = bootstrap_N,
                        check_installed = FALSE, type = bootstrap_type,
                        seed = seed, average = TRUE, sample_prop = sample_prop,
                        var_in_dims = var_in_dims,
                        normalise_var = normalise_var, ...)
    })
  } else {
    res = lapply(ks, function(k) {
      # adjust number of dimensions to k
      if(isTRUE(simplex)) {
        data_dim = seq(1, k-1)
        if(k %in% c(1,2)) data_dim = seq(1, 2)
      } else data_dim = seq(1, nrow(data))
      # fit models
      fit_pch(data = data[data_dim,], noc = k,
              ..., check_installed = FALSE,
              var_in_dims = var_in_dims,
              normalise_var = normalise_var)
    })
  }
  # combine results ------------------------------------------------------------
  if(length(ks) > 1){
    res = list(call = match.call(),
               pch_fits = .c_pch_fit_list(res))
    class(res) = "k_pch_fit"
    res$summary = res$pch_fits$summary
    res$pch_fits$summary = NULL
  } else {
    res = res[[1]]
  }
  res
}

##' @rdname fit_pch
##' @name fit_pch_bootstrap
##' @description \code{fit_pch_bootstrap()} uses resampling data with or without replacement (bootstrapping) to find robust positions of vertices of a polytope (Principal Convex Hull) describing the data. This function uses \code{fit_pch_resample()}.
##' @param n number of samples to be taken when bootstraping
##' @param sample_prop either NULL or the proportion of dataset that should be included in each sample. If NULL the polytope fitting algorithm is run n times on the same data which is useful for evaluating how often the algorithm gets stuck in local optima.
##' @param type one of s, m, cmq. s means single core processing using lapply. m means multi-core parallel procession using parLapply. cmq means multi-node parallel processing on a computing cluster using clustermq package. "See also" for details.
##' @param clust_options list of options for parallel processing. The default for "m" is list(cores = parallel::detectCores()-1, cluster_type = "PSOCK"). The default for "cmq" is list(memory = 2000, template = list(), n_jobs = 10, fail_on_error = FALSE). Change these options as required.
##' @param seed seed for reproducible random number generation. Works for all types of processing.
##' @param replace should resampling be with replacement? passed to \link[base]{sample.int}
##' @param average average archetype positions and varexpl? By default FALSE, return all fits to resampled data.
##' @param var_in_dims calculate variance in dimensions across vertices in addition to variance in vertices.
##' @param normalise_var normalise variance in position of vertices by variance in data in each dimention
##' @param reference align vertices (order_type="align") against reference polytope based on all data (TRUE), or align to the polytope from the first bootstrap iteration (FALSE). Second option can save time for large datasets
##' @return \code{fit_pch_bootstrap()} object of class b_pch_fit (list) containing the same elements as pch_fit, but each is either a list of pch_fit elements (e.g. list of n number of XC matrices) or a vector (which pch_fit element is one number).
##' @import clustermq
##' @export fit_pch_bootstrap
fit_pch_bootstrap = function(data, n = 3, sample_prop = NULL, check_installed = T,
                             type = c("s", "m", "cmq")[1], clust_options = list(),
                             seed = 235, replace = FALSE, average = FALSE,
                             order_type = c("cosine", "side", "align")[3],
                             var_in_dims = TRUE, normalise_var = TRUE,
                             reference = FALSE, ...) {
  if(check_installed) .py_pcha_installed()
  # single process -------------------------------------------------------------
  if(type == "s"){
    set.seed(seed)
    res = lapply(seq_len(n), fit_pch_resample, data, sample_prop,
                 replace = replace, order_type = order_type,
                 converge_else_fail = FALSE,
                 var_in_dims = var_in_dims,
                 normalise_var = normalise_var, ...)
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
    # set seed
    parallel::clusterSetRNGStream(cl, iseed = seed)
    res = parallel::parLapply(cl, seq_len(n), ParetoTI::fit_pch_resample, data,
                              sample_prop, replace = replace,
                              order_type = order_type,
                              converge_else_fail = FALSE,
                              var_in_dims = var_in_dims,
                              normalise_var = normalise_var, ...)
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
        res = clustermq::Q(fun = ParetoTI::fit_pch_resample, i = seq_len(n),
                           const = list(data = data, sample_prop = sample_prop,
                                        replace = replace,
                                        order_type = order_type,
                                        converge_else_fail = FALSE,
                                        var_in_dims = var_in_dims,
                                        normalise_var = normalise_var, ...),
                           seed = seed,
                           memory = options$memory, template = options$template,
                           n_jobs = options$n_jobs, rettype = "list",
                           fail_on_error = options$fail_on_error,
                           timeout = options$timeout)
      })
    })
  }
  # combine results ------------------------------------------------------------

  res = list(call = match.call(),
             pch_fits = .c_pch_fit_list(res))
  class(res) = "b_pch_fit"

  # match archetypes and reorder results ---------------------------------------

  if(isTRUE(order_type == "align")){
    # align to reference based on all data
    if(isTRUE(reference)){
      suppressMessages({
        ref = fit_pch(data = data, order_type = order_type,
                      converge_else_fail = FALSE, ...)
        if(is.null(ref)) return(NULL)
      })
      ref_XC = ref$XC
    } else {
      # align to reference - first polytope from bootstraping
      ref_XC = res$pch_fits$XC[[1]]
    }
    for (i in seq_len(length(res$pch_fits$XC))) {
      ind = align_arc(ref$XC, res$pch_fits$XC[[i]])$ind
      res$pch_fits$XC[[i]] = res$pch_fits$XC[[i]][, ind]
      res$pch_fits$S[[i]] = res$pch_fits$S[[i]][ind, ]
      res$pch_fits$C[[i]] = res$pch_fits$C[[i]][ind, ]
    }
  }

  # extract XC matrix for calculating averages and variance --------------------

  XC_array = simplify2array(res$pch_fits$XC)
  # average results ------------------------------------------------------------
  # update summary data.table
  res$summary = data.table::data.table(varexpl = mean(res$pch_fits$varexpl),
                                       t_ratio = mean(res$pch_fits$t_ratio),
                                       total_var = NA,
                                       k = res$pch_fits$summary$k[1])
  # average arhetype positions
  if(isTRUE(average)){
    res = average_pch_fits(res = res, XC_array = XC_array)
  }

  # calculate and include variance in positions --------------------------------

  dim. = c(dim(XC_array)[1] * dim(XC_array)[2], dim(XC_array)[3])
  res$var = matrix(matrixStats::rowVars(XC_array, dim. = dim.),
                   dim(XC_array)[1], dim(XC_array)[2])

  # calculate variance of data in each dimension
  if(isTRUE(normalise_var)) row_var_data = matrixStats::rowVars(data)

  # sum variance in position of each vertex across dimensions
  res$var_vert = colSums(res$var)
  # normalise variance in vertex positions by variance of data
  if(isTRUE(normalise_var)) res$var_vert = res$var_vert / sum(row_var_data)
  # find mean of variances of all vertices to get a single number
  res$total_var = mean(res$var_vert)

  # normalise variance in vertex positions by variance of data in each dimension
  if(isTRUE(normalise_var)) res$var = res$var / matrixStats::rowVars(data)

  # find mean variance in position for each dimension across vertices
  res = ParetoTI:::.cacl_var_in_dims(res, data, var_in_dims, normalise_var, XC_array)

  # update summary data.table  -------------------------------------------------
  res$summary[, total_var := res$total_var]
  res$summary = unique(res$summary)
  res
}

##' @rdname fit_pch
##' @name average_pch_fits
##' @description \code{average_pch_fits()} averages archetypes positions and relevant metrics across polytope fit obtained using bootstraping. Uses output of fit_pch_bootstrap() directly.
##' @param XC_array Used by average_pch_fits() inside fit_pch_bootstrap(). You should not to use it in most cases.
##' @return \code{average_pch_fits()} object of class pch_fit
##' @export average_pch_fits
average_pch_fits = function(res, XC_array = NULL){
  if(!is(res, "b_pch_fit")) stop("average_pch_fits(): res object provided is not b_pch_fit")
  if(is.null(XC_array)) XC_array = simplify2array(res$pch_fits$XC)
  res_aver = list(call = res$call)
  # calculate average XC matrix, set other S and C to NA
  res_aver$XC = rowMeans(XC_array, dims = 2)
  res_aver$S = NA
  res_aver$C = NA
  # calculate mean varexpl and t_ratio
  res_aver$SSE = mean(res$pch_fits$SSE)
  res_aver$varexpl = mean(res$pch_fits$varexpl)
  res_aver$hull_vol = mean(res$pch_fits$hull_vol)
  res_aver$arc_vol = mean(res$pch_fits$arc_vol)
  res_aver$t_ratio = mean(res$pch_fits$t_ratio)
  res_aver$var_vert = res$var_vert
  res_aver$var_dim = res$var_dim
  res_aver$total_var = res$total_var
  res_aver$summary = res$summary
  class(res_aver) = "pch_fit"
  res_aver
}

##' @rdname fit_pch
##' @name fit_pch_resample
##' @description \code{fit_pch_resample()} takes one sample of the data and fits a polytope (Principal Convex Hull) to data. This function uses \code{fit_pch()}.
##' @param i iteration number
##' @return \code{fit_pch_resample()} object of class pch_fit
##' @export fit_pch_resample
fit_pch_resample = function(i = 1, data, sample_prop = NULL, replace = FALSE,
                            var_in_dims = var_in_dims,
                            normalise_var = normalise_var, ...) {
  # do resampling of the data
  if(!is.null(sample_prop)){
    if(data.table::between(sample_prop[1], 0, 1)){
      col_ind = sample.int(ncol(data), round(ncol(data) * sample_prop[1], digits = 0),
                           replace = replace)
      data = data[, col_ind]
    } else stop("sample_prop should be NULL or a number between 0 and 1")
  }
  # fit polytope
  ParetoTI::fit_pch(data = data, ..., var_in_dims = var_in_dims,
                    normalise_var = normalise_var, check_installed = FALSE)
}

##' @rdname fit_pch
##' @name randomise_fit_pch
##' @description \code{randomise_fit_pch()} helps answer the question "how likely you are to observe data shaped as polytope given no relationship between variables?" This function disrupts the relationships between variables, find a polytope that best describes the data. It calculates variance explained, t-ratio of polytope volume to convex hull volume. Optionally, it uses bootstrapping of data points to measure variance in vertex positions. This function uses \code{randomise_fit_pch1()}. Number of vertices is determined automatically from \code{arc_data}.
##' @param arc_data observed polytope that describes the data. Class pch_fit, b_pch_fit and k_pch_fit objects generated by \code{fit_pch()} or \code{fit_pch_bootstrap()}
##' @param n_rand number of randomisation samples
##' @return \code{randomise_fit_pch()} S3 "r_pch_fit" object (list) containing: 1. data table with columns variance explained (rand_varexpl), t-ratio (rand_t_ratio), variance in positions of vertices (total_var) and number of vertices (k) for \code{n_rand} samples; 2. empirical p-value for those parameters for each number of vertices.
##' @import clustermq
##' @import data.table
##' @export randomise_fit_pch
randomise_fit_pch = function(data, arc_data, n_rand = 3, replace = FALSE,
                             bootstrap_N = c(NA, 50)[1], seed = 435,
                             volume_ratio = c("t_ratio", "variance_ratio", "none")[1],
                             maxiter = 1000, delta = 0,
                             order_type = "align", type = c("s", "m", "cmq")[1],
                             clust_options = list(), normalise_var = TRUE,
                             check_installed = T, ...) {
  if(check_installed) .py_pcha_installed()

  # extract some parameters of observed fit
  ks = sort(unique(arc_data$summary$k))
  var_in_dims = isTRUE(nrow(arc_data$var_dim) > 0) |
    isTRUE(nrow(arc_data$pch_fits$var_dim) > 0)
  # single process -------------------------------------------------------------

  if(type == "s"){
    set.seed(seed)
    res = lapply(seq_len(n_rand), randomise_fit_pch1, data = data, ks = ks,
                 replace = replace, bootstrap_N = bootstrap_N, seed = seed,
                 bootstrap_type = "s",
                 return_data = FALSE, return_arc = FALSE,
                 bootstrap_average = TRUE, volume_ratio = volume_ratio,
                 maxiter = maxiter, delta = delta,
                 order_type = order_type,
                 var_in_dims = var_in_dims,
                 normalise_var = normalise_var, ...)
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
    # set seed
    parallel::clusterSetRNGStream(cl, iseed = seed)
    res = parallel::parLapply(cl, seq_len(n_rand), ParetoTI::randomise_fit_pch1,
                              data = data, ks = ks,
                              replace = replace, bootstrap_N = bootstrap_N,
                              seed = seed, bootstrap_type = "s",
                              return_data = FALSE,
                              return_arc = FALSE, bootstrap_average = TRUE,
                              volume_ratio = volume_ratio,
                              maxiter = maxiter, delta = delta,
                              order_type = order_type,
                              var_in_dims = var_in_dims,
                              normalise_var = normalise_var, ...)
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
        res = clustermq::Q(fun = ParetoTI::randomise_fit_pch1, i = seq_len(n_rand),
                           const = list(data = data, ks = ks,
                                        replace = replace, bootstrap_N = bootstrap_N,
                                        seed = seed, bootstrap_type = "s",
                                        return_data = FALSE,
                                        return_arc = FALSE, bootstrap_average = TRUE,
                                        volume_ratio = volume_ratio,
                                        maxiter = maxiter, delta = delta,
                                        order_type = order_type,
                                        var_in_dims = var_in_dims,
                                        normalise_var = normalise_var, ...),
                           seed = seed,
                           memory = options$memory, template = options$template,
                           n_jobs = options$n_jobs, rettype = "list",
                           fail_on_error = options$fail_on_error,
                           timeout = options$timeout)
      })
    })
  }

  # combine results ------------------------------------------------------------

  # remove failed iterations: usual reason "PCHA failed to converge within max_iter"
  res = res[!vapply(res, is.null,
                    FUN.VALUE = logical(1))]

  rand_dist = lapply(res, function(x) x$summary)
  rand_dist = rbindlist(rand_dist)
  # melt this data.table
  setorder(rand_dist, k)
  rand_dist[, repl := seq_len(.N), by = .(k)]
  rand_dist = melt.data.table(rand_dist, id.vars = c("k", "repl"),
                              value.name = "var_rand", variable.name = "var_name")
  # add readable labels
  rand_dist[, var_label := .type_lab(unique(var_name), short = TRUE),
            by = .(var_name)]

  # match observed to randomised archetypes ------------------------------------

  obs_dist = copy(arc_data$summary)
  obs_dist[, k := as.character(k)]
  rand_dist[, k := as.character(k)]
  obs_dist = melt.data.table(obs_dist, id.vars = c("k"),
                             value.name = "var_obs", variable.name = "var_name")
  rand_dist = merge(rand_dist, obs_dist, by = c("k", "var_name"),
                    all = TRUE)
  # mark cases when observed value is at least as high in random data
  rand_dist[, obs_vs_rand := var_obs >= var_rand]
  # invert this for variance in positions
  rand_dist[var_name == "total_var", obs_vs_rand := var_obs <= var_rand]
  # substract proportion of those case from 1 to get p-value
  rand_dist[, p_value := 1 - mean(obs_vs_rand, na.rm = TRUE), by = .(var_name, k)]
  # Set cases with p-value=0 to 1 / number of non-na replicates, same for p-value=1
  rand_dist[p_value == 0, p_value := 1 / sum(!is.na(var_rand)),
            by = .(var_name, k)]
  rand_dist[p_value == 1, p_value := 1 - 1 / sum(!is.na(var_rand)),
            by = .(var_name, k)]

  obs_dist = unique(rand_dist[,.(k, var_name, var_obs, p_value)])

  # combine position variability in dimensions data  ---------------------------
  var_dim = lapply(res, function(x) x$var_dim)
  var_dim = rbindlist(var_dim)
  # if variability in dimensions was measured, calculate empirical pval --------
  if(nrow(var_dim) > 0){
    # melt observed data
    if(is(arc_data, "k_pch_fit")){
      var_dim_obs = arc_data$pch_fits$var_dim
    } else {
      var_dim_obs = arc_data$var_dim
    }
    var_dim_obs = melt.data.table(var_dim_obs, id.vars = c("k"),
                                  value.name = "var_obs", variable.name = "var_name")
    # add replicate labels and melt randomised data
    setorder(var_dim, k)
    var_dim[, repl := seq_len(.N),  by = .(k)]
    var_dim = melt.data.table(var_dim, id.vars = c("k", "repl"),
                              value.name = "var_rand", variable.name = "var_name")
    # merge observed and randomised
    var_dim = merge(var_dim, var_dim_obs, by = c("k", "var_name"),
                    all = TRUE)
    setorder(var_dim, k, -var_obs)

    # mark cases when observed variance in positions is at least as low in random data
    var_dim[, obs_vs_rand := var_obs < var_rand]
    # substract proportion of those case from 1 to get p-value
    var_dim[, p_value := 1 - mean(obs_vs_rand, na.rm = TRUE), by = .(var_name, k)]
    # Set cases with p-value=0 to 1 / number of non-na replicates, same for p-value=1
    var_dim[p_value == 0, p_value := 1 / sum(!is.na(var_rand)),
            by = .(var_name, k)]
    var_dim[p_value == 1, p_value := 1 - 1 / sum(!is.na(var_rand)),
            by = .(var_name, k)]
  } else var_dim = var_dim

  # return ---------------------------------------------------------------------

  res = list(call = match.call(),
             rand_dist = rand_dist,
             obs_dist = obs_dist,
             var_dim = var_dim)
  class(res) = "r_pch_fit"
  res
}

##' @rdname fit_pch
##' @name randomise_fit_pch1
##' @description \code{randomise_fit_pch1()} disrupts the relationships between variables (one sampling iteration), keeping the distribution of each variable constant, and fits a polytope (Principal Convex Hull) to data. This function uses \code{\link[ParetoTI]{rand_var}} and \code{fit_pch()}.
##' @param prob a vector of probability weights for obtaining the elements of the vector being sampled. Passed to \code{\link[base]{(sample.int}}.
##' @param bootstrap_N randomise_fit_pch1() and k_fit_pch(): integer, number of bootstrap samples on random data to measure variability in vertex positions. Set to NA fit once to all data instead of bootstraping. When this option is positive seed bootstrap_seed and sample_prop must be provided.
##' @param bootstrap_type \code{randomise_fit_pch1()} and \code{k_fit_pch()}: parallel processing type when bootstraping. Caution: avoid nested parallel processing, do not use "m" and "cmq" inside other parallel functions.
##' @param return_data return randomised data?
##' @param return_arc return archetype positions in randomised data?
##' @param bootstrap_average \code{randomise_fit_pch1()}: average positions and summary statistics when bootstraping? Passed to \code{average} argument of \code{fit_pch_bootstrap()}. When multiple ks this defaults to TRUE.
##' @return \code{randomise_fit_pch1()}: list containing function call, summary of the sample, optional data and optional position of vertices/archetypes.
##' @import data.table
##' @export randomise_fit_pch1
randomise_fit_pch1 = function(i = 1, data, ks = 2:4,
                              replace = FALSE, prob = NULL,
                              bootstrap_N = NA, seed = 435, sample_prop = 0.65,
                              bootstrap_type = c("s", "m", "cmq")[1],
                              return_data = FALSE, return_arc = FALSE,
                              bootstrap_average = FALSE,
                              volume_ratio = c("t_ratio", "variance_ratio", "none")[1],
                              var_in_dims = TRUE, normalise_var = TRUE, ...) {
  # randomise variables --------------------------------------------------------

  data = ParetoTI::rand_var(data, MARGIN = 1, replace = replace, prob = prob)

  # fit polytope ---------------------------------------------------------------

  if(is.na(bootstrap_N)) { # single fit ----------------------------------------
    # fitting many shapes (ks vertices) or fitting one shape
    arc_data = ParetoTI::k_fit_pch(data = data, ks = ks, check_installed = FALSE,
                                   bootstrap = FALSE, seed = seed,
                                   volume_ratio = volume_ratio,
                                   var_in_dims = var_in_dims,
                                   normalise_var = normalise_var,
                                   # prevent whole pipeline from failing
                                   # when PCHA doesn't converge:
                                   converge_else_fail = FALSE, ...)

  } else if(isTRUE(as.integer(bootstrap_N) >= 1)) { # bootstraping --------------
    if(length(ks) > 1){ # fitting many shapes (ks vertices)
      arc_data = ParetoTI::k_fit_pch(data = data, ks = ks, check_installed = TRUE,
                                     bootstrap = TRUE, bootstrap_N = bootstrap_N,
                                     sample_prop = sample_prop,
                                     bootstrap_type = bootstrap_type, seed = seed,
                                     replace = replace,
                                     var_in_dims = var_in_dims,
                                     normalise_var = normalise_var,
                                     ..., volume_ratio = volume_ratio)
    } else { # fitting many shapes (ks vertices)
      arc_data = ParetoTI::fit_pch_bootstrap(data, n = bootstrap_N, noc = ks,
                                             check_installed = FALSE,
                                             type = bootstrap_type, seed = seed,
                                             sample_prop = sample_prop,
                                             replace = replace, average = bootstrap_average,
                                             var_in_dims = var_in_dims,
                                             normalise_var = normalise_var,
                                             ..., volume_ratio = volume_ratio)
    }
  }

  ## combine results into summary ----------------------------------------------

  # escape if null
  if(is.null(arc_data)) return(NULL)

  res = list()
  res$call = match.call()
  res$summary = arc_data$summary
  ks_1 = length(ks) == 1
  if(ks_1){
    res$var_dim = arc_data$var_dim
  } else {
    res$var_dim = arc_data$pch_fits$var_dim
  }

  ## choose to return data and archetype positions -----------------------------

  if(isTRUE(return_data)) res$data = data else res$data = NULL
  if(isTRUE(return_arc)) res$arc_data = arc_data else res$arc_data = NULL
  res
}

##' @rdname fit_pch
##' @name fit_convhulln
##' @description \code{fit_convhulln()} computes smallest convex hull that encloses a set of points using \code{\link[geometry]{convhulln}} and returns positions of convex hull points (positions = TRUE).
##' @param positions return positions of convex hull points?
##' @return \code{fit_convhulln()} list with 3 elements: hull - matrix storing positions of points dim(points, dimensions)); area - surface area for 3D or higher and perimeter for 2D; vol - volume for 3D and surface area for 2D.
##' @export fit_convhulln
fit_convhulln = function(data, positions = TRUE) {
  # convhulln requires dim(examples, variables/dimentions)
  # in contrast to all other functions and PCHA
  hull = geometry::convhulln(t(data), options = "FA")
  if(isTRUE(positions)){
    hull$hull = vapply(seq(1, ncol(hull$hull)), function(i){
      data[hull$hull,i]
    }, FUN.VALUE = numeric(nrow(hull$hull) * 2))
    hull$hull = unique(hull$hull)
  } else hull$hull = NA
  hull
}

##' @rdname fit_pch
##' @name merge_arch_dist
##' @description \code{merge_arch_dist()} calculates distance to archtypes, and merges distance matrix to data that was used to identify archetypes, optionally adds other features of data points (through colData).
##' @param feature_data matrix with dim(dimensions, examples) where rownames are feature names and colnames are sample_id.
##' @param colData annotations of examples in feature_data - dim(examples, dimensions), e.g. colData in SingleCellExperiment object or output of \link[ParetoTI]{find_set_activity_AUCell}.
##' @param colData_id column in colData that contains values matching colnames of feature_data.
##' @param dist_metric how to describe distance to archetypes. Currently euclidian distance is implemented.
##' @param rank rank by distance metric (euclidian distance)?
##' @return \code{merge_arch_dist()} list: data.table with samples in columns and features in rows (speficied by sample_id column) and column names specifying archetypes
##' @export merge_arch_dist
##' @import data.table
merge_arch_dist = function(arch_data, data, feature_data,
                           colData = NULL, colData_id,
                           dist_metric = "euclidian", rank = TRUE){
  if(!is(arch_data, "pch_fit")) stop("arch_data should contain a single fit (pch_fit object): use fit_pch() or fit_pch_bootstrap() followed by average_pch_fits()")
  # find distance of data points to archetypes ---------------------------------
  dist = arch_dist(data, arch_data$XC, dist_metric = dist_metric)
  arc_col = colnames(dist)
  # if no rownames provided add V1, V2 ... Vn names
  if(is.null(rownames(dist))) rownames(dist) = paste0("V", seq_len(nrow(dist)))
  dist = as.data.table(dist, keep.rownames = "sample_id")
  # convert distances to ranks and scale between 0 and 1 -----------------------
  # (max rank for min distance)
  if(isTRUE(rank)){
    for (col in arc_col) {
      dist[, c(col) := frank(get(col), ties.method=c("average")) / .N]
    }
  }
  # merge gene data ------------------------------------------------------------
  # if no sample_id (colnames) provided add V1, V2 ... Vn names
  if(is.null(colnames(feature_data))) {
    colnames(feature_data) = paste0("V", seq_len(ncol(feature_data)))
  }
  features = as.data.table(t(feature_data), keep.rownames = "sample_id")
  features = merge(dist, features, by = "sample_id", all.x = T, all.y = F)
  # merge column annotations of data -------------------------------------------
  # (e.g. sample features in gene expression)
  if(!is.null(colData)) {
    features = merge(features, colData, by.x = "sample_id", by.y = colData_id,
                     all.x = T, all.y = F)
  }
  list(data = features,
       arc_col = arc_col)
}

.c_pch_fit_list = function(pch_fit_list){
  # remove failed fits
  pch_fit_list = pch_fit_list[!vapply(pch_fit_list, is.null,
                                      FUN.VALUE = logical(1))]
  # combine results
  list(XC = lapply(pch_fit_list, function(pch) pch$XC),
       S = lapply(pch_fit_list, function(pch) pch$S),
       C = lapply(pch_fit_list, function(pch) pch$C),
       SSE = vapply(pch_fit_list, function(pch) pch$SSE, FUN.VALUE = numeric(1L)),
       varexpl = vapply(pch_fit_list, function(pch) pch$varexpl, FUN.VALUE = numeric(1L)),
       hull_vol = vapply(pch_fit_list, function(pch) pch$hull_vol, FUN.VALUE = numeric(1L)),
       arc_vol = vapply(pch_fit_list, function(pch) pch$arc_vol, FUN.VALUE = numeric(1L)),
       t_ratio = vapply(pch_fit_list, function(pch) pch$t_ratio, FUN.VALUE = numeric(1L)),
       var_vert = lapply(pch_fit_list, function(pch) pch$var_vert),
       var_dim = data.table::rbindlist(lapply(pch_fit_list, function(pch) pch$var_dim)),
       total_var = vapply(pch_fit_list, function(pch) pch$total_var, FUN.VALUE = numeric(1L)),
       summary = data.table::rbindlist(lapply(pch_fit_list, function(pch) pch$summary)))
}

## .find_vertex_order orders vertices by cosine relative to the unit vector c(1, 1)
## Cosine is negative when vertex vector is to the left (counter-clockwise) of the unit vector
## Solution taken from here: https://stackoverflow.com/questions/13221873/determining-if-one-2d-vector-is-to-the-right-or-left-of-another
## XC2 is a matrix of dim(dimensions, archetype) that has only 2 dimentions
.find_vertex_order = function(XC2, noc, order_type = c("cosine", "side", "align")[3]){
  if(isTRUE(order_type == "side")){
    order_by = vapply(seq(1, noc), function(i) -XC2[1,i] + XC2[2,i],
                      FUN.VALUE = numeric(1L))
  } else if(isTRUE(order_type %in% c("cosine", "align"))) {
    order_by = vapply(seq(1, noc), function(i){
      sum(XC2[,i]) / sqrt(2 * sum(XC2[,i]^2))
    }, FUN.VALUE = numeric(1L))
  } else stop("order_type should be one of c(\"cosine\", \"side\", \"align\")")
  order(order_by)
}

##' @rdname fit_pch
##' @name print.pch_fit
##' @export print.pch_fit
print.pch_fit = function(res){
  cat("Minimal polytope with k vertices describing the data (S3 class pch_fit)\n\n")
  cat("Summary:  varexpl = variance within polytope\n          t_ratio = volume of polytope / volume of convex hull\n          total_var = total variance in positions of vertices\n            (by bootstraping, mean acrooss vertices)\n\n")
  print(res$summary)
  cat("\n vertex positions stored in 'object$XC' \n dim(variables/dimentions, vertices/archetypes)")
}

##' @rdname fit_pch
##' @name print.k_pch_fit
##' @export print.k_pch_fit
print.k_pch_fit = function(res){
  cat("Minimal polytopes with k vertices describing the data (S3 class k_pch_fit)\n\n")
  cat("Summary:  varexpl = variance within polytope\n          t_ratio = volume of polytope / volume of convex hull\n          total_var = total variance in positions of vertices\n            (by bootstraping, mean acrooss vertices)\n\n")
  print(res$summary)
  cat("\n vertex positions stored in 'object$pch_fits$XC' \n dim(variables/dimentions, vertices/archetypes)")
}

##' @rdname fit_pch
##' @name print.b_pch_fit
##' @export print.b_pch_fit
print.b_pch_fit = function(res){
  cat("Minimal polytopes describing the data \n calculated by bootstraping data points (cells) (S3 class b_pch_fit)\n\n")
  cat("Average summary:  varexpl = variance within polytope\n          t_ratio = volume of polytope / volume of convex hull\n          total_var = total variance in positions of vertices\n            (by bootstraping, mean acrooss vertices)\n\n")
  print(res$summary)
  cat("\nAll polytopes (fit to bootstrapped data), summary:\n")
  print(res$pch_fits$summary)
  cat("\n vertex positions stored in 'object$pch_fits$XC' \n dim(variables/dimentions, vertices/archetypes)")
}

##' @rdname fit_pch
##' @name print.r_pch_fit
##' @export print.r_pch_fit
print.r_pch_fit = function(res){
  cat("Background distribution of minimal polytopes fits \nto data with no relationships between variables (S3 class r_pch_fit)\n\n")
  cat("N randomisation trials: ",res$rand_dist[, max(repl)],"\n\n")
  cat("Summary of best-fit polytope to observed data (including p-value):\n\n")
  print(res$obs_dist)
  cat("\n          varexpl = variance within polytope\n          t_ratio = volume of polytope / volume of convex hull\n          total_var = total variance in positions of vertices\n            (by bootstraping, mean acrooss vertices)\n\n")
  cat("Function call:\n")
  print(res$call)
}

##' @rdname fit_pch
##' @name plot.r_pch_fit
##' @export plot.r_pch_fit
##' @description plot.r_pch_fit() plot distribution of t-ratio, total variance in vertex positions and variance explained used to calculate empirical p-value
##' @import ggplot2
plot.r_pch_fit = function(rand_arch,
                          ks = unique(rand_arch$rand_dist$k),
                          type = c("total_var", "varexpl", "t_ratio"),
                          colors = c("#1F77B4", "#D62728", "#2CA02C", "#17BED0", "#006400", "#FF7E0F"),
                          nudge_y = 0.8, nudge_x = 0.1, text_lab_size = 4, line_size = 0.5){

  # filter data
  rand_dist = copy(rand_arch$rand_dist[k %in% ks & var_name %in% type])
  rand_dist[, k := as.character(k)]
  # add p-value text only to the first row of each var_name and k
  rand_dist[, p_value_text := c(paste0("p-val:\n", signif(p_value[1], 2), ""),
                                rep("", .N-1)),
            by = .(var_name, k)]

  # plot
  ggplot(rand_dist, aes(x = var_rand, group = k, y = k,
                        color = k, fill = k)) +
    ggridges::geom_density_ridges(alpha = 0.5) +
    ylab("Density") +
    xlab("Metric value") +
    ggridges::theme_ridges() +
    theme(strip.text.y = element_text(angle = 0),
          plot.title = element_text(size=14,face="bold"),
          axis.title = element_text(size=14,face="bold"),
          axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 0),
          legend.position = "none",
          panel.grid.minor = element_line(size = 0),
          panel.grid.major.x = element_line(size = 0)) +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    geom_segment(aes(x = var_obs, xend = var_obs,
                     y = as.numeric(as.factor(k)), yend = as.numeric(as.factor(k)) + 0.95,
                     color = k, group = k)) +
    facet_grid(. ~ var_label, scales = "free_x") +
    geom_text(aes(x = var_obs, y = k,
                  label = p_value_text,
                  color = k, group = k), size = text_lab_size,
              nudge_y = nudge_y, nudge_x = nudge_x, inherit.aes = F)
}


##' @rdname fit_pch
##' @name plot_dim_var
##' @export plot_dim_var
##' @description plot_dim_var() plot distribution of variance in vertex positions in each dimension and corresponding empirical p-values
##' @import ggplot2
plot_dim_var = function(rand_arch,
                        ks = unique(rand_arch$rand_dist$k),
                        dim_names = c("V1", "V2", "V3", "V4", "V5", "V6"),
                        colors = c("#1F77B4", "#D62728", "#2CA02C",
                                   "#17BED0", "#006400", "#FF7E0F"),
                        nudge_y = 0.5, nudge_x = 0.5,
                        text_lab_size = 4, line_size = 0.5){

  # filter dimnames  and ks
  var_dim = copy(rand_arch$var_dim[k %in% ks & var_name %in% dim_names])

  # add p-value text only to the first row of each var_name and k
  var_dim[, p_value_text := c(paste0(signif(p_value[1], 2)),
                              rep("", .N-1)),
          by = .(var_name, k)]
  var_dim[, k := as.character(k)]

  ggplot(var_dim, aes(x = var_rand, group = var_name, y = var_name,
                      color = var_name, fill = var_name)) +
    ggridges::geom_density_ridges(alpha = 0.5) +
    scale_color_viridis_d() + scale_fill_viridis_d() +
    scale_x_log10() +
    facet_grid(. ~ k, scales = "free_x") +
    ggridges::theme_ridges() +
    theme(strip.text.y = element_text(angle = 0),
          plot.title = element_text(size=14,face="bold"),
          axis.title = element_text(size=14,face="bold"),
          axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 0),
          legend.position = "none",
          panel.grid.minor = element_line(size = 0),
          panel.grid.major.x = element_line(size = 0)) +
    geom_segment(aes(x = var_obs, xend = var_obs,
                     y = as.numeric(var_name), yend = as.numeric(var_name) + 0.95,
                     color = var_name, group = var_name)) +
    geom_text(aes(x = var_obs, y = as.numeric(var_name),
                  label = p_value_text,
                  color = var_name, group = var_name), size = text_lab_size,
              nudge_y = nudge_y, nudge_x = nudge_x) +
    xlab("Variance in vertex positions") +
    ylab("Dimension where variance is measured")
}

.cacl_var_in_dims = function(res, data, var_in_dims, normalise_var,
                             XC_array = NULL) {
  if(var_in_dims){
    # calculate variance in vertex positions in each dimension
    if(is(res, "pch_fit")){
      res$var_dim = matrix(matrixStats::rowVars(res$XC), nrow = 1, ncol = nrow(data))
      noc = ncol(res$XC)
    } else if(is(res, "b_pch_fit")){
      res$var_dim = matrix(XC_array, dim(XC_array)[1],
                           dim(XC_array)[2] * dim(XC_array)[3])
      res$var_dim = matrix(matrixStats::rowVars(res$var_dim),
                           nrow = 1, ncol = nrow(data))
      noc = dim(XC_array)[2]
    }
    # normalise by variance of data in each dimension
    if(normalise_var) res$var_dim = res$var_dim / matrix(matrixStats::rowVars(data),
                                                         nrow = 1, ncol = nrow(data))
    res$var_dim = data.table::data.table(res$var_dim)
    if(!is.null(rownames(data))) {
      data.table::setnames(res$var_dim, colnames(res$var_dim), rownames(data))
    }
    res$var_dim$k = noc
  } else {
    # add empty results if not needed
    res$var_dim = data.table::data.table(matrix(nrow = 0, ncol = nrow(data)))
    if(!is.null(rownames(data))) {
      data.table::setnames(res$var_dim, colnames(res$var_dim), rownames(data))
    }
    res$var_dim$k = integer(0)
  }
  res
}
