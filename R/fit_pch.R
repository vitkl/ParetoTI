##' Fit a polytope (Principal Convex Hull) to data using PCHA algorithm
##' @rdname fit_pch
##' @name fit_pch
##' @author Vitalii Kleshchevnikov
##' @description \code{fit_pch()} fits a polytope (Principal Convex Hull) to data using PCHA algorithm. All of the listed functions take input matrix dim(variables/dimentions, examples) and return output (XC) of dim(variables/dimentions, archetypes).
##' @details \code{fit_pch()} provides an R interface to python implementation of PCHA algorithm (Principal Convex Hull Analysis) by Ulf Aslak (https://github.com/ulfaslak/py_pcha) which was originally developed for Archetypal Analysis by Mørup et. al.
##' @param data numeric matrix in which to find archetypes, variables/dimentions in rows, examples in columns
##' @param noc integer, number of archetypes to find
##' @param I vector, entries of data to use for dictionary in C (optional)
##' @param U vector, entries of data to model in S (optional)
##' @param delta parameter that inflates original polytope(simplex) fit such that it may contain more points of the dataset
##' @param verbose if TRUE display messages
##' @param conv_crit The convergence criteria (default: 10^-6 relative change in SSE)
##' @param maxiter maximum number of iterations (default: 500 iterations)
##' @param check_installed if TRUE, check if python module py_pcha is found. Useful to set to FALSE for running analysis or within other functions
##' @param order_by integer, dimensions to be used for ordering vertices/archetypes. Vertices are ordered by angle (cosine) between c(1, 1) vector and a vector pointing to that vertex. Additional step finds when vertex vector is to the left (counter-clockwise) of the c(1, 1) vector.
##' @param order_by_side used for ordering. If TRUE use 2 dimentions, cosine distance to c(1,1) and measure to which side of the c(1,1) vector each vertex vector is located. If FALSE use more than 2 dimentions (provided via order_by) and use cosine distance to c(1,1, ..., 1) to order vertices.
##' @param convex_hull find volume of the convex hull of the data and the t-ratio? Caution. Computation time and memory use increse very quickly with dimensions. Do not use for more than 7-8 dimentions. Geometric figure should be at least simplex: qhull algorhirm will fail to find convel hull of flat 2D shapes in 3D, 3D shapes in 4D and so on.
##' @return \code{fit_pch()} object of class pch_fit (list) containing the following elements:
##' XC - numeric matrix, dim(I, noc)/dim(dimensions, archetypes) feature matrix (i.e. XC=data[,I]*C forming the archetypes);
##' S - numeric matrix, dim(noc, length(U)) matrix, S>=0 |S_j|_1=1;
##' C - numeric matrix, dim(noc, length(U)) matrix, S>=0 |S_j|_1=1;
##' SSE - numeric vector (1L), Sum of Squared Errors;
##' varexpl - numeric vector (1L), Percent variation explained by the model.
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
##' arc = fit_pch(data, noc=as.integer(3), delta=0.1)
##' # Fit the same polytope 3 times without resampling to test convergence of the algorithm.
##' arc_rob = fit_pch_bootstrap(data, n = 3, sample_prop = NULL,
##'                          noc=as.integer(3), delta=0.1)
##' # Fit the 10 polytopes to resampled datasets each time looking at 70% of examples.
##' arc_data = fit_pch_bootstrap(data, n = 10, sample_prop = 0.7,
##'                          noc=as.integer(3), delta=0.1)
##'
##' # Use local parallel processing to fit the 10 polytopes to resampled datasets each time looking at 70% of examples.
##' arc_data = fit_pch_bootstrap(data, n = 10, sample_prop = 0.7,
##'                          noc=as.integer(3), delta=0.1, type = "m")
##'
##' # Fit polytopes with 2-4 vertices
##' arc_ks = k_fit_pch(data, ks = 2:4, check_installed = T, delta=0.1)
fit_pch = function(data, noc = as.integer(3), I = NULL, U = NULL,
                   delta = 0, verbose = FALSE, conv_crit = 1e-6,
                   maxiter = 500, check_installed = T,
                   order_by = seq(1, nrow(data)), order_by_side = TRUE,
                   convex_hull = FALSE) {
  if(check_installed) .py_pcha_installed()
  # run PCHA
  res = py_PCHA$PCHA(X = data, noc = as.integer(noc), I = I, U = U,
                     delta = delta, verbose = verbose,
                     conv_crit = conv_crit, maxiter = maxiter)
  names(res) = c("XC", "S", "C", "SSE", "varexpl")
  # step sorting archetypes to improve reproducibility
  # (archetypes are inherently exchangeable)
  if(noc > 1) {
    XC2 = res$XC[order_by,]
    arch_order = .find_vertex_order(XC2, noc, order_by_side = order_by_side)
    res$XC = res$XC[, arch_order]
    res$S = res$S[arch_order, ]
    res$C = res$C[arch_order, ]
  } else {
    # when only one vertex make sure data is still in the matrix form
    res$XC = matrix(res$XC, length(res$XC), 1)
    res$S = matrix(res$S, 1, length(res$S))
    res$C = matrix(res$C, 1, 1)
  }

  if(isTRUE(convex_hull)){
    # find volume of the convex hull of the data
    hull_vol = fit_convhulln(data, positions = FALSE)
    # find volume of the polytope fit, qhull requires at least 4 points in 2D
    # and more in higher dimensions so we have to add 10 points within
    # PCHA-fit polytope
    archetypes = res$XC
    data_arc = cbind(archetypes, generate_data(archetypes, N_examples = 10,
                                               jiiter = 0, size = 1))
    arc_vol = fit_convhulln(t(data_arc), positions = FALSE)
    res$hull_vol = hull_vol$vol
    res$arc_vol = arc_vol$vol
    res$t_ratio = res$arc_vol / res$hull_vol
  } else {
    res$hull_vol = NA
    res$arc_vol = NA
    res$t_ratio = NA
  }
  res$call = match.call()
  class(res) = "pch_fit"
  res
}

##' @rdname fit_pch
##' @name k_fit_pch
##' @description \code{k_fit_pch()} finds polytopes of k dimensions in the data. This function applies \code{fit_pch()} to different k-s.
##' @param ks integer vector, dimensions of polytopes to be fit to data
##' @return \code{k_fit_pch()} object of class k_pch_fit (list) containing XC (list of length(ks), individual XC matrices), S (list of length(ks), individual S matrices),  C (list of length(ks), individual C matrices), SSE (vector, length(ks)); varexpl - (vector, length(ks)). When length(ks) = 1 returns pch_fit.
##' @import clustermq
##' @export k_fit_pch
k_fit_pch = function(data, ks = 2:4, check_installed = T, ...) {
  if(check_installed) .py_pcha_installed()
  # run analysis for all k -----------------------------------------------------
  res = lapply(ks, function(k) fit_pch(data = data, noc = k,
                                       ..., check_installed = F))
  # combine results ------------------------------------------------------------
  if(length(ks) > 1){
    res = list(call = match.call(),
               pch_fits = .c_pch_fit_list(res))
    class(res) = "k_pch_fit"
  } else {
    res = res[[1]]
  }
  res
}

##' @rdname fit_pch
##' @name fit_pch_bootstrap
##' @description \code{fit_pch_bootstrap()} Uses bootstrapping (resampling with) the data to find robust positions of vertices of a polytope (Principal Convex Hull) to data. This function uses \code{fit_pch_resample()}.
##' @param n number of samples to be taken when bootstraping
##' @param sample_prop either NULL or the proportion of dataset that should be included in each sample. If NULL the polytope fitting algorithm is run n times on the same data which is useful for evaluating how often the algorithm gets stuck in local optima.
##' @param type one of s, m, cmq. s means single core processing using lapply. m means multi-core parallel procession using parLapply. cmq means multi-node parallel processing on a computing cluster using clustermq package. "See also" for details.
##' @param clust_options list of options for parallel processing. The default for "m" is list(cores = parallel::detectCores()-1, cluster_type = "PSOCK"). The default for "cmq" is list(memory = 2000, template = list(), n_jobs = 10, fail_on_error = FALSE). Change these options as required.
##' @param seed seed for reproducible random number generation. Works for all types of processing.
##' @param replace should resampling be with replacement? passed to \link[base]{sample.int}
##' @param average average archetype positions and varexpl? By default NULL, return all fits to resampled data. Provide funtion such as mean() and median() to average those.
##' @return \code{fit_pch_bootstrap()} object of class b_pch_fit (list) containing XC (list of length n, individual XC matrices), S (list of length n, individual S matrices),  C (list of length n, individual C matrices), SSE (vector, length n); varexpl - (vector, length n).
##' @import clustermq
##' @export fit_pch_bootstrap
fit_pch_bootstrap = function(data, n = 3, sample_prop = NULL, check_installed = T,
                             type = c("s", "m", "cmq")[1], clust_options = list(),
                             seed = 235, replace = FALSE, average = NULL,  ...) {
  if(check_installed) .py_pcha_installed()
  # single process -------------------------------------------------------------
  if(type == "s"){
    set.seed(seed)
    res = lapply(seq_len(n), fit_pch_resample, data, sample_prop,
                 replace = replace, ...)
  }
  # multi-process --------------------------------------------------------------
  if(type == "m"){
    # set defaults or replace them with provided options
    default = list(cores = parallel::detectCores()-1, cluster_type = "PSOCK")
    default_retain = !names(default) %in% names(clust_options)
    options = c(default[default_retain], clust_options)
    # create cluster
    cl <- parallel::makeCluster(options$cores, type = options$cluster_type)
    # get library support needed to run the code
    parallel::clusterEvalQ(cl, {library(ParetoTI)})
    # put objects in place that might be needed for the code
    # parallel::clusterExport(cl, ls(envir = environment()), envir=environment())
    # set seed
    parallel::clusterSetRNGStream(cl, iseed = seed)
    res = parallel::parLapply(cl, seq_len(n), ParetoTI::fit_pch_resample, data,
                              sample_prop, replace = replace, ...)
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
    res = clustermq::Q(fun = ParetoTI::fit_pch_resample, i = seq_len(n),
                       const = list(data = data, sample_prop = sample_prop,
                                    replace = replace, ...),
                       seed = seed,
                       memory = options$memory, template = options$template,
                       n_jobs = options$n_jobs, rettype = "list",
                       fail_on_error = options$fail_on_error,
                       timeout = options$timeout)
  }
  # combine results ------------------------------------------------------------
  res = list(call = match.call(),
             pch_fits = .c_pch_fit_list(res))
  # average results ------------------------------------------------------------
  if(!is.null(average)){
    if(is.function(average)){
      # calculate average XC matrix, set other S and C to NA
      XC_array = simplify2array(arc_data4$pch_fits$XC)
      arc_data3$pch_fits$XC = list(rowMeans(XC_array, dims = 2))
      # calculate variance in positions - add option that returns mean positions and put that inside k_fit (this lets to output variances and varexpl in the same object)
      dim. = c(dim(XC_array)[1] * dim(XC_array)[2], dim(XC_array)[3])
      var = matrix(matrixStats::rowVars(XC_array,dim. = dim.),
                   dim(XC_array)[1], dim(XC_array)[2])
      matrix(matrixStats::rowMeans2(XC_array, dim. = dim.),
             dim(XC_array)[1], dim(XC_array)[2])
      # sum up variances for all vertices to get a single number - put this into fit_pch_bootstrap
      # calculate mean varexpl and t_ratio
    } else stop("average neither function nor NULL")
  }
  class(res) = "b_pch_fit"
  res
}

##' @rdname fit_pch
##' @name fit_pch_resample
##' @description \code{fit_pch_resample()} takes one sample of the data and fits a polytope (Principal Convex Hull) to data. This function uses \code{fit_pch()}.
##' @param i iteration number
##' @return \code{fit_pch_resample()} object of class pch_fit (list)
##' @export fit_pch_resample
fit_pch_resample = function(i = 1, data, sample_prop = NULL, replace = FALSE, ...) {
  # do resampling of the data
  if(!is.null(sample_prop)){
    if(data.table::between(sample_prop[1], 0, 1)){
      col_ind = sample.int(ncol(data), round(ncol(data) * sample_prop[1], digits = 0),
                           replace = replace)
      data = data[, col_ind]
    } else stop("sample_prop should be NULL or a number between 0 and 1")
  }
  # fit polytope
  ParetoTI::fit_pch(data = data, ..., check_installed = F)
}

##' @rdname fit_pch
##' @name randomise_fit_pch1
##' @description \code{randomise_fit_pch1()} helps answer the question "how likely you are to obtain the observed shape of the data given no relationship between variables?" disrupts the relationships between variables (one sample of the data), keeping the distribution of each variable constant, and fits a polytope (Principal Convex Hull) to data. This function uses \code{\link[ParetoTI]{rand_var}} and \code{fit_pch()}.
##' @param prob a vector of probability weights for obtaining the elements of the vector being sampled. Passed to \code{\link[base]{(sample.int}}.
##' @param bootstrap_N integer, number of bootstrap samples on random data to measure variability in vertex positions. When this option is chosen bootstrap_seed and sample_prop must be provided
##' @param bootstrap_type \code{randomise_fit_pch1()}: parallel processing type when bootstraping. Caution: avoid nested parallel processing, do not use "m" and "cmq" inside other parallel functions
##' @param return_data return randomised data?
##' @param return_arc return archetype positions in randomised data?
##' @return \code{randomise_fit_pch1()} object of class pch_fit (list)
##' @export randomise_fit_pch1
randomise_fit_pch1 = function(i = 1, data, true_fit = NULL,
                              replace = FALSE, prob = NULL,
                              bootstrap_N = NA, seed = 435,
                              bootstrap_type = c("s", "m", "cmq")[1],
                              return_data = FALSE, return_arc = FALSE, ...) {
  # randomise variables
  set.seed(seed)
  data = ParetoTI::rand_var(data, MARGIN = 1, replace = replace, prob = prob)
  # fit polytope
  if(is.na(bootstrap_N)) { # single
    arc_data = ParetoTI::fit_pch(data = data, ..., check_installed = F)
  } else if(isTRUE(as.integer(bootstrap_N) > 1)) { # with bootstrap
    arc_data = ParetoTI::fit_pch_bootstrap(data, n = bootstrap_N,
                                           check_installed = F, type = bootstrap_type,
                                           seed = seed,
                                           replace = replace, ...)
  }
  # compare to true_fit
  res = list()
  res$rand_varexpl = arc_data$varexpl
  res$rand_t_ratio = arc_data$t_ratio
  # return
  res$call = match.call()
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

.c_pch_fit_list = function(pch_fit_list){
  list(XC = lapply(pch_fit_list, function(pch) pch$XC),
       S = lapply(pch_fit_list, function(pch) pch$S),
       C = lapply(pch_fit_list, function(pch) pch$C),
       SSE = vapply(pch_fit_list, function(pch) pch$SSE, FUN.VALUE = numeric(1L)),
       varexpl = vapply(pch_fit_list, function(pch) pch$varexpl, FUN.VALUE = numeric(1L)),
       hull_vol = vapply(pch_fit_list, function(pch) pch$hull_vol, FUN.VALUE = numeric(1L)),
       arc_vol = vapply(pch_fit_list, function(pch) pch$arc_vol, FUN.VALUE = numeric(1L)),
       t_ratio = vapply(pch_fit_list, function(pch) pch$t_ratio, FUN.VALUE = numeric(1L)))
}

## .find_vertex_order orders vertices by cosine relative to the unit vector c(1, 1)
## Cosine is negative when vertex vector is to the left (counter-clockwise) of the unit vector
## Solution taken from here: https://stackoverflow.com/questions/13221873/determining-if-one-2d-vector-is-to-the-right-or-left-of-another
## XC2 is a matrix of dim(dimensions, archetype) that has only 2 dimentions
.find_vertex_order = function(XC2, noc, order_by_side = TRUE){
  if(isTRUE(order_by_side)){
    order_by_side = vapply(seq(1, noc), function(i) -XC2[1,i] + XC2[2,i],
                           FUN.VALUE = numeric(1L))
    cosine = order_by_side
  } else {
    cosine = vapply(seq(1, noc), function(i){
      sum(XC2[,i]) / sqrt(2 * sum(XC2[,i]^2))
    }, FUN.VALUE = numeric(1L))
  }
  order(cosine)
}
