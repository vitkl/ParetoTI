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
##' @param delta TO DO: find definition in matlab code
##' @param verbose TO DO: find definition in matlab code, probably "display some messages"
##' @param conv_crit TO DO: find definition in matlab code, probably "stop whe improvement in each iteraction is < conv_crit"
##' @param maxiter TO DO: find definition in matlab code, probably "max number of iterations"
##' @param check_installed if TRUE, check if python module py_pcha is found. Useful to set to FALSE for running analysis or within other functions
##' @return \code{fit_pch()} object of class pch_fit (list) containing the following elements:
##' XC - numeric matrix, I x noc feature matrix (i.e. XC=data[,I]*C forming the archetypes);
##' S - numeric matrix, noc x length(U) matrix, S>=0 |S_j|_1=1;
##' C - numeric matrix, noc x length(U) matrix, S>=0 |S_j|_1=1;
##' SSE - numeric vector (1L), Sum of Squared Errors;
##' varexlp - numeric vector (1L), Percent variation explained by the model.
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
##' # Fit the same polytope 3 times without subsampling to test convergence of the algorithm.
##' arc_rob = fit_pch_robust(data, n = 3, subsample = NULL,
##'                          noc=as.integer(3), delta=0.1)
##' # Fit the 10 polytopes to subsampled datasets each time looking at 70% of examples.
##' arc_data = fit_pch_robust(data, n = 10, subsample = 0.7,
##'                          noc=as.integer(3), delta=0.1)
##'
##' # Use local parallel processing to fit the 10 polytopes to subsampled datasets each time looking at 70% of examples.
##' arc_data = fit_pch_robust(data, n = 10, subsample = 0.7,
##'                          noc=as.integer(3), delta=0.1, type = "m")
fit_pch = function(data, noc = as.integer(3), I = NULL, U = NULL,
                   delta = 0, verbose = FALSE, conv_crit = 1e-6,
                   maxiter = 500, check_installed = T) {
  if(check_installed) .py_pcha_installed()
  # run PCHA
  res = py_PCHA$PCHA(X = data, noc = as.integer(noc), I = I, U = U,
                     delta = delta, verbose = verbose,
                     conv_crit = conv_crit, maxiter = maxiter)
  names(res) = c("XC", "S", "C", "SSE", "varexlp")
  # add step sorting archetypes to improve reproducibility (archetypes are inherently exchangeable)

  res$call = match.call()
  class(res) = "pch_fit"
  res
}

##' @rdname fit_pch
##' @name fit_pch_robust
##' @description \code{fit_pch_robust()} subsamples the data to find robust positions of vertices of a polytope (Principal Convex Hull) to data. This function uses \code{fit_pch_subsample()}.
##' @param n number of samples that should be taken
##' @param subsample either NULL or the proportion of dataset that should be included in each sample. If NULL the polytope fitting algorithm is run n times which is usefult for evaluating how often the algorithm gets stuck in local optima
##' @param type one of s, m, cmq. s means single core processing using lapply. m means multi-core parallel procession using parLapply. cmq means multi-node parallel processing on a computing cluster using clustermq package. "See also" for details.
##' @param clust_options list of options for parallel processing. The default for "m" is list(cores = parallel::detectCores()-1, cluster_type = "PSOCK"). The default for "cmq" is list(memory = 2000, template = list(), n_jobs = 10, fail_on_error = FALSE). Change these options as required.
##' @param seed seed for reproducible random number generation. Works for all types of processing.
##' @return \code{fit_pch_robust()} object of class r_pch_fit (list) containing XC (list of length n, individual XC matrices), S (list of length n, individual S matrices),  C (list of length n, individual C matrices), SSE (vector, length n); varexlp - (vector, length n).
##' @import clustermq
##' @export fit_pch_robust
fit_pch_robust = function(data, n = 3, subsample = NULL, check_installed = T,
                          type = c("s", "m", "cmq")[1], clust_options = list(),
                          seed = 235, ...) {
  if(check_installed) .py_pcha_installed()
  # single process -------------------------------------------------------------
  if(type == "s"){
    set.seed(seed)
    res = lapply(seq_len(n), fit_pch_subsample, data, subsample, ...)
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
    res = parallel::parLapply(cl, seq_len(n), ParetoTI::fit_pch_subsample, data,
                              subsample, ...)
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
    res = clustermq::Q(fun = ParetoTI::fit_pch_subsample, i = seq_len(n),
                       const = list(data = data, subsample = subsample, ...), # figure out how to supply this
                       seed = seed,
                       memory = options$memory, template = options$template,
                       n_jobs = options$n_jobs, rettype = "list",
                       fail_on_error = options$fail_on_error,
                       timeout = options$timeout)
  }
  # combine results ------------------------------------------------------------
  res = list(call = match.call(),
             pch_fits = .c_pch_fit_list(res))
  class(res) = "r_pch_fit"
  res
}

##' @rdname fit_pch
##' @name fit_pch_subsample
##' @description \code{fit_pch_subsample()} takes one subsample of the data and fits a polytope (Principal Convex Hull) to data. This function uses \code{fit_pch()}.
##' @param i iteration number
##' @return \code{fit_pch_subsample()} object of class pch_fit (list)
##' @export fit_pch_subsample
fit_pch_subsample = function(i = 1, data, subsample = NULL, ...) {
  if(!is.null(subsample)){
    if(data.table::between(subsample[1], 0, 1)){
      col_ind = sample.int(ncol(data), round(ncol(data) * subsample[1], digits = 0))
      data = data[, col_ind]
    } else stop("subsample should be NULL or a number between 0 and 1")
  }
  ParetoTI::fit_pch(data = data, ..., check_installed = F)
}

.c_pch_fit_list = function(pch_fit_list){
  list(XC = lapply(pch_fit_list, function(pch) pch$XC),
       S = lapply(pch_fit_list, function(pch) pch$S),
       C = lapply(pch_fit_list, function(pch) pch$C),
       SSE = vapply(pch_fit_list, function(pch) pch$SSE, FUN.VALUE = numeric(1L)),
       varexlp = vapply(pch_fit_list, function(pch) pch$varexlp, FUN.VALUE = numeric(1L)))
}
