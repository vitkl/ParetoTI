##' Find the smallest polytope (Principal Convex Hull) that contains most of the data
##' @rdname fit_pch
##' @name fit_pch
##' @author Vitalii Kleshchevnikov
##' @description \code{fit_pch()} fits a polytope (Principal Convex Hull) to data using PCHA algorithm. All of the listed functions take input matrix dim(variables/dimentions, examples) and return archetype positions (XC) of dim(variables/dimentions, archetypes).
##' @details \code{fit_pch()} provides an R interface to python implementation of PCHA algorithm (Principal Convex Hull Analysis) by Ulf Aslak (https://github.com/ulfaslak/py_pcha) which was originally developed for Archetypal Analysis by Mørup et. al.
##' @param data numeric matrix or object coercible to matrix in which to find archetypes, dim(variables/dimentions, examples)
##' @param noc integer, number of archetypes to find
##' @param I vector, entries of data to use for dictionary in C (optional)
##' @param U vector, entries of data to model in S (optional)
##' @param delta parameter that inflates original polytope(simplex) fit such that it may contain more points of the dataset
##' @param verbose if TRUE display messages
##' @param conv_crit The convergence criteria (default: 10^-4 relative change in SSE). Reduce to 1e-3 for reduced computation time but more approximate results. Increase to 1e-6 for improved accuracy (python and Matlab default). 1e-4 gives very similar results to 1e-5 or 1e-6 on datasets I looked at.
##' @param maxiter maximum number of iterations (default: 500 iterations)
##' @param check_installed if TRUE, check if python module py_pcha is found. Useful to set to FALSE for running analysis or within other functions
##' @param order_by integer, dimensions to be used for ordering archetypes. Archetypes are ordered by angle (cosine) between c(1, 1) vector and a vector pointing to that archetype. Additional step finds when archetype vector is to the left (counter-clockwise) of the c(1, 1) vector. When bootstraping archetypes can be aligned to reference and ordered (order_type == "align") by these dimensions.
##' @param order_type order archetypes by: cosine distance from c(1,1, ..., 1) vector ("cosine"), dot product that measures position to each side of the c(1,1) vector ("side"), align positions to reference when bootstraping(fit using all data, "align"). See \link[ParetoTI]{align_arc} When order_type is "align" archetypes are ordered by "cosine" first.
##' @param volume_ratio find volume of the convex hull of the data and the t-ratio ("t_ratio") or ratio of variance of archetype positions to variance of data ("variance_ratio") or do nothing ("none")? Caution! Geometric figure should be at least simplex: qhull algorhirm will fail to find convex hull of flat 2D shapes in 3D, 3D shapes in 4D and so on. So, for calculating this dimentions are selected based order of rows in data. Makes sense for principal components but not for original data. Caution 2! Computation time and memory use increse very quickly with dimensions. Do not use for more than 7-8 dimentions.
##' @param converge_else_fail throw an error and stop execution if PCHA did not converge in \code{maxiter} steps.
##' @param method method for archetypal analysis: PCHA (default), \link[stats]{kmeans}, \link[Seurat]{FindClusters} (louvain). Non-default methods use \code{data} and \code{noc}for data matrix and the number of archetypes, but require to pass additional arguments via method_options.
##' @param method_options additional arguments for non-default method, named list: \code{list(arg = "value")}.
##' @return \code{fit_pch()}: object of class pch_fit (list) containing the following elements:
##' XC - numeric matrix, dim(I, noc)/dim(dimensions, archetypes) feature matrix (i.e. XC=data[,I]*C forming the archetypes), or matrix with cluster centers (method = "kmeans");
##' S - numeric matrix of archetype values per example (e.g. cell) or cluster membership (method = "kmeans"), dim(noc, length(U)) matrix, S>=0 |S_j|_1=1;
##' C - numeric matrix, dim(noc, length(U)) matrix, S>=0 |S_j|_1=1;
##' SSE - numeric vector (1L), Sum of Squared Errors;
##' varexpl - numeric vector (1L), Percent variation explained by the model.
##' hull_vol - numeric vector (1L), Volume of convex hull of the data.
##' arc_vol - numeric vector (1L), Volume of polytope enclosed by archetypes.
##' t_ratio - numeric vector (1L), Ratio of \code{arc_vol} to \code{hull_vol}
##' var_vert - numeric vector (noc L), Variance in position of each archetype.
##' var_dim - numeric vector (noc L), Variance in positions of archetypes in each dimension.
##' total_var - numeric vector (1L), Mean variance in position of all archetypes.
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
##' # Find 3 archetypes in this data
##' arc = fit_pch(data, noc=as.integer(3), delta=0)
##' # Fit the model 3 times without resampling to test convergence of the algorithm.
##' arc_rob = fit_pch_bootstrap(data, n = 3, sample_prop = NULL,
##'                          noc=as.integer(3), delta=0)
##' # Fit 10 models to resampled datasets each time looking at 70% of examples.
##' arc_data = fit_pch_bootstrap(data, n = 10, sample_prop = 0.7,
##'                          noc=as.integer(3), delta=0)
##'
##' # Use local parallel processing to fit 10 models to resampled datasets each time looking at 70% of examples.
##' arc_data = fit_pch_bootstrap(data, n = 10, sample_prop = 0.7,
##'                          noc=as.integer(3), delta=0, type = "m")
##'
##' # Fit models with 2-4 archetypes
##' arc_ks = k_fit_pch(data, ks = 2:4, check_installed = T, delta=0)
##'
##' # Evaluate how much archetypes vary in randomised data, (variable shuffled
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
                   delta = 0, verbose = FALSE, conv_crit = 1e-4,
                   maxiter = 500, check_installed = T,
                   order_by = seq(1, nrow(data)),
                   order_type = c("align", "cosine", "side"),
                   volume_ratio = c("t_ratio", "variance_ratio", "none"),
                   converge_else_fail = TRUE,
                   var_in_dims = FALSE, normalise_var = TRUE,
                   method = c("pcha", "kmeans", "louvain", "poisson_aa", "aanet"),
                   method_options = list()) {

  # check arguments
  method = match.arg(method)
  volume_ratio = match.arg(volume_ratio)
  order_type = match.arg(order_type)


  if(method == "pcha"){

    if(check_installed) .py_pcha_installed()

    # run PCHA -------------------------------------------------------------------

    # coerce to matrix
    if(!is.matrix(data)) data = as.matrix(data)

    res = tryCatch({
      res = py_PCHA$PCHA(X = data, noc = as.integer(noc), I = I, U = U,
                         delta = delta, verbose = verbose,
                         conv_crit = conv_crit, maxiter = maxiter)
      names(res) = c("XC", "S", "C", "SSE", "varexpl")
      if(!is.null(rownames(data))) rownames(res$XC) = rownames(data)
      if(!is.null(colnames(data))) {
        colnames(res$S) = colnames(data)
      }

      class(res) = "pch_fit"
      res
    }, error = function(err) {
      if(isTRUE(converge_else_fail)) stop(paste0("fit_pch(noc = ",noc,") error: ", err))
      return(NULL)
    })

    method_res = NA

    #---------------------------------------------------------------------------

  } else if(method == "louvain") {

    # run louvain clustering using Seurat --------------------------------------

    # set defaults or replace them with provided options
    default = list(k.param = 20, prune.SNN = 1/15,
                   resolution = 1, n.start = 10, n.iter = 10,
                   graph = NA, lr = 0.2, noc_optim_iter = 200)
    default_retain = !names(default) %in% names(method_options)
    options = c(default[default_retain], method_options)

    # compute nearest neibours SNN graph
    if(is.na(options$graph)) {
      options$graph = Seurat::FindNeighbors(t(data),
                                            k.param = options$k.param,
                                            prune.SNN = options$prune.SNN,
                                            verbose = verbose)
    }


    # run clustering by optimising resolution to reach desired noc
    iter = 1
    clusters = NULL
    while((!isTRUE(uniqueN(clusters) == noc) | iter == 1) & # iterate until desired noc is reached
          iter <= options$noc_optim_iter ) { # or until noc_optim_iter exhausted

      # compute gradient of resolution
      if(iter == 1) grad = 0 else {
        grad = options$resolution * options$lr * log(noc / uniqueN(clusters))
      }
      options$resolution = abs(options$resolution + grad)
      # increment iterations
      iter = iter + 1

      # run clustering
      clusters = Seurat::FindClusters(options$graph$snn,
                                      resolution = options$resolution,
                                      n.start = options$n.start,
                                      n.iter = options$n.iter,
                                      verbose = verbose)
      clusters = as.integer(clusters[, 1])

    }

    # Create S as binary cluster membership matrix
    S = Matrix::sparseMatrix(i = clusters,
                             j = seq_len(length(clusters)), x = 1)
    # compute C so that X %*% C gives cluster averages
    C = S / matrix(Matrix::rowSums(S), nrow = nrow(S), ncol = ncol(S), byrow = FALSE)
    C = Matrix::t(C)

    # compute SSE (within cluster variance) and total variance
    ss = function(x) sum(scale(t(x), scale = FALSE)^2)
    totss = ss(data)
    tot.withinss = sum(vapply(unique(clusters), function(x){
      ss(data[, clusters == x])
    }, FUN.VALUE = numeric(1L)))

    # create pch_fit object to be returned
    res = list(XC = as.matrix(Matrix::crossprod(Matrix::t(data), C)),
               S = S, C = C,
               SSE = tot.withinss, varexpl = (totss - tot.withinss) / totss)
    if(!is.null(rownames(data))) rownames(res$XC) = rownames(data)
    colnames(res$XC) = NULL
    class(res) = "pch_fit"

    # change noc according to what louvain found
    noc = length(unique(clusters))

    method_res = list(resolution = options$resolution)

    #---------------------------------------------------------------------------

  } else if(method == "kmeans") {

    # run k-means --------------------------------------------------------------

    # set defaults or replace them with provided options
    default = list(iter.max = 10, nstart = 1,
                   algorithm = c("Hartigan-Wong", "Lloyd", "Forgy",
                                 "MacQueen"), trace=FALSE)
    default_retain = !names(default) %in% names(method_options)
    options = c(default[default_retain], method_options)

    res = kmeans(Matrix::t(data), centers = as.integer(noc),
                 iter.max = options$iter.max, nstart = options$nstart,
                 algorithm = options$algorithm, trace = options$trace)

    # Create S as binary cluster membership matrix
    S = Matrix::sparseMatrix(i = res$cluster,
                             j = seq_len(length(res$cluster)), x = 1)
    # compute C so that X %*% C gives cluster averages
    C = S / matrix(Matrix::rowSums(S), nrow = nrow(S), ncol = ncol(S), byrow = FALSE)
    C = Matrix::t(C)

    # create pch_fit object to be returned
    res = list(XC = t(res$centers), S = S, C = C,
               SSE = res$tot.withinss, varexpl = res$betweenss / res$totss)
    if(!is.null(rownames(data))) rownames(res$XC) = rownames(data)
    colnames(res$XC) = NULL
    class(res) = "pch_fit"

    method_res = NA

    #---------------------------------------------------------------------------

  }  else if(method == "poisson_aa"){

    # run probabilistic poisson archetypal analysis  ---------------------------

    # set defaults or replace them with provided options ========
    default = list(model_fun = ParetoTI::paa_poisson,
                   weight_alpha_prior = 0.8, c_alpha_prior = 0.001,
                   covar = NULL, precision = c("double", "single"),
                   optimiser = greta::adam(learning_rate = 0.3),
                   maxiter = maxiter, conv_crit = conv_crit,
                   initial_values = greta::initials(),
                   resolution = 1) # louvain initial values resolution
    default_retain = !names(default) %in% names(method_options)
    options = c(default[default_retain], method_options)

    # optionally: use louvain cluster centers as initial values
    if(isTRUE(options$initial_values == "louvain_centers")){
      # initial values ===============
      # compute Louvain clustering with Seurat to use as initial value
      clust = fit_pch(data, noc = noc, method = "louvain",
                      method_options = list(resolution = options$resolution),
                      volume_ratio = "none")
      c_init = Matrix::t(clust$C) + 1e-7 # add small random value to initialise correctly
      c_init = c_init / Matrix::rowSums(c_init)
      weights_init = Matrix::t(clust$S) + 1e-5 # add small random value to initialise correctly
      weights_init = weights_init / Matrix::rowSums(weights_init)
      arch_init = t(clust$XC)

      if(isTRUE(all.equal(options$model_fun, ParetoTI::paa_poisson))){
        options$initial_values = greta::initials(c_var = as.matrix(c_init),
                                                 weights_var = as.matrix(weights_init))
      } else if(isTRUE(all.equal(options$model_fun, ParetoTI::paa_poisson_free))) {
        options$initial_values = greta::initials(archetypes = as.matrix(arch_init),
                                                 weights_var = as.matrix(weights_init))
      }
    }


    # create greta / tensorflow model ========
    m = options$model_fun(t(data),                   # data: data points * dimensions
                          n_arc = noc,              # number of achetypes
                          weight_alpha_prior = options$weight_alpha_prior,
                          c_alpha_prior = options$c_alpha_prior,
                          covar = options$covar,            # covariates affecting the mean in addition to archetypes: data points * n_covariates
                          precision = options$precision
    )
    m$options = options

    # visualise computation graph ========
    if(verbose) plot(m$model)

    # add initial values if any provided ========
    # (this has to be done in the same environment as the components of the model)
    if(!isTRUE(all.equal(options$initial_values, greta::initials()))){

      evalq({
        init = options$initial_values

        vals = paste0(names(init), " = init[[",
                      seq_along(init), "]]", collapse = ", ")
        vals = paste0("initial_values = greta::initials(", vals, ")")
        vals = str2expression(vals)

        eval(expr = vals)
      }, envir = m)

    } else {

      evalq({initial_values = greta::initials()}, envir = m)

    }

    # solve model with optimisation ========
    evalq({
      tryCatch({
        opt_res = opt(model,
                      optimiser = options$optimiser,   # optimisation method used to find prior-adjusted maximum likelihood estimate
                      max_iterations = options$maxiter,
                      initial_values = initial_values,
                      tolerance = options$conv_crit, adjust = TRUE,
                      hessian = FALSE)
      }, error = function(err) {
        if(grepl("Error in while \\(self\\$it", err)) {
          stop(paste0(err, "\n\n try reducing learning rate with `method_options = list(optimiser = optimiser(learning_rate = 0.3))`\n
                         e.g. `method_options = list(optimiser = greta::adam(learning_rate = 0.3))`"))
        } else stop(err)
        return(NULL)
      })
    }, envir = m)

    opt_res = m$opt_res

    # extract results depending on the model ========
    if(isTRUE(all.equal(options$model_fun, ParetoTI::paa_poisson))){
      XC = t(opt_res$par$c %*% t(data))
      C = t(opt_res$par$c)
    } else if(isTRUE(all.equal(options$model_fun, ParetoTI::paa_poisson_free))) {
      XC = t(opt_res$par$archetypes)
      C = t(opt_res$par$weights)
    }
    # create pch_fit object ========
    res = list(XC = XC, S = t(opt_res$par$weights), C = C,
               SSE = opt_res$iterations, # number of iterations (e.i. did it converge?)
               varexpl = opt_res$value) # negative log-likelihood

    if(!is.null(rownames(data))) {
      rownames(res$XC) = rownames(data)
      colnames(res$S) = colnames(data)
      rownames(res$C) = colnames(data)
    }
    class(res) = "pch_fit"

    method_res = opt_res

  } else if(method == "aanet"){

    # run AAnet -------------------------------------------------------------------

    # coerce to matrix
    if(!is.matrix(data)) data = as.matrix(data)
    data = t(data) # transpose

    # set defaults or replace them with provided options ========
    default = list(noise_z_std = 0.0,
                   z_dim = list(256L, 64L), act_out = tensorflow::tf$nn$tanh,
                   learning_rate = 1e-3, gpu_mem = 0.4,
                   batch_size=128L, num_batches=5000L,
                   input_dim = ncol(data),
                   AAnet = reticulate::import_from_path("AAnet", path = "../AAnet/", convert = TRUE),
                   network = reticulate::import_from_path("network", path = "../AAnet/", convert = TRUE),
                   AAtools = reticulate::import_from_path("AAtools", path = "../AAnet/", convert = TRUE))

    default_retain = !names(default) %in% names(method_options)
    options = c(default[default_retain], method_options)

    # construct network  ========
    enc_net = network$Encoder(num_at=as.integer(noc), z_dim=options$z_dim)
    dec_net = network$Decoder(x_dim=options$input_dim, noise_z_std=options$noise_z_std, z_dim=options$z_dim, act_out=options$act_out)
    model = AAnet$AAnet(enc_net, dec_net, learning_rate=options$learning_rate, gpu_mem=options$gpu_mem)

    # train model ========
    model$train(aanet_data, batch_size=options$batch_size, num_batches=options$num_batches)

    # export results  ========
    # get archetypes in input space
    XC = t(model$get_ats_x())
    # compute loss function
    mse = model$compute_mse_loss(data = aanet_data)
    # compute S
    S = t(model$data2at(data = aanet_data))

    res = list(XC = XC, S = S, C = NULL, SSE = mse, varexpl = mse)
    if(!is.null(colnames(data))) rownames(res$XC) = colnames(data)
    if(!is.null(rownames(data))) {
      colnames(res$S) = rownames(data)
    }

    class(res) = "pch_fit"

    method_res = model

    #---------------------------------------------------------------------------

  } else stop("method should be pcha, poisson_aa, aanet, louvain or kmeans")

  if(is.null(res)) return(NULL)

  # step sorting archetypes to improve reproducibility
  # (archetypes are inherently exchangeable)
  if(noc > 1) {
    XC2 = res$XC[order_by,]
    arch_order = .find_archetype_order(XC2, noc, order_type = order_type)
    res$XC = res$XC[, arch_order]
    res$S = res$S[arch_order, ]
    res$C = res$C[, arch_order]
  } else {
    # when only one archetype make sure data is still in the matrix form
    res$XC = matrix(res$XC, length(res$XC), 1)
    res$S = matrix(res$S, 1, length(res$S))
    res$C = matrix(res$C, 1, 1)
  }

  # Calculate ratio of volumes to measure goodness of fit ----------------------
  # 1 t_ratio  -----------------------------------------------------------------
  # when calculating convex hull and volume of the polytope
  # adjust number of dimensions to noc
  data_dim = seq(1, noc-1)
  if(volume_ratio == "t_ratio" & nrow(data) >= length(data_dim) & noc > 2 & method != "poisson_aa"){
    # calculate volume or area of the polytope only when number of archetypes (noc) > number of dimenstions which means
    # find volume of the convex hull of the data
    hull_vol = fit_convhulln(data[data_dim, ], positions = FALSE)

    # Calcute the volume of non-simplex polytope (4 archetypes in 2D, 5 in 3D or more)
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
      res$arc_vol = abs(Matrix::det(arch_red[, seq_len(noc - 1)]) /
                          factorial(noc - 1))
    }

    res$hull_vol = hull_vol$vol
    res$t_ratio = res$arc_vol / res$hull_vol

  } else if(volume_ratio == "variance_ratio") {

    # 2 variance ratio ---------------------------------------------------------
    # (variance of archetype positions within each dimension /
    #      variance of data in that dimension summed across dimensions)
    res = ParetoTI:::.cacl_var_in_dims(res, data,
                                       var_in_dims = T, normalise_var = T)

    dim_names = colnames(res$var_dim)[!colnames(res$var_dim) %in% "k"]
    res$t_ratio = mean(as.matrix(res$var_dim[, dim_names, with = FALSE]))

    res$hull_vol = NA
    res$arc_vol = NA

  } else {

    # 3 NA

    if(volume_ratio == "t_ratio" & isTRUE(converge_else_fail)) message(paste0("Convex hull and t-ratio not computed for noc: ", noc," and nrow(data) = ", nrow(data),". fit_pch() can calculate volume or area of the polytope only when\nthe number of archetypes (noc) > the number of dimensions (when polytope is convex):\ncheck that noc > nrow(data),\nselect only revelant dimensions or increase noc"))
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

  # add method-speficic results
  res$method_res = method_res

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
##' @param bootstrap \code{k_fit_pch()}: use bootstrap to find average positions for each k? Also returns variability in archetype position to aid the selection of k. At excessive k position vary more.
##' @param simplex when testing multiple k using \code{k_fit_pch()} match dimensions to the number of archetypes? Use only on ordered principal components. If FALSE try all k for all dimensions in data. If simplex == TRUE test only simplex shapes (k=3 in 2D, k=4 in 3D, k=5 in 4D...). This assumes order of columns which is valid for principal components but may not be valid for untransformed data.
##' @return \code{k_fit_pch()}: object of class k_pch_fit (list) containing the same elements as pch_fit, but each is either a list of pch_fit elements (e.g. list of ks number of XC matrices) or a vector (which pch_fit element is one number). When length(ks) = 1 returns pch_fit.
##' @export k_fit_pch
k_fit_pch = function(data, ks = 2:4, check_installed = TRUE,
                     bootstrap = FALSE, bootstrap_N = 10, sample_prop = 0.65,
                     bootstrap_type = c("s", "m", "cmq")[1],
                     seed = 345, simplex = FALSE,
                     var_in_dims = FALSE, normalise_var = TRUE, ...) {

  if(check_installed) .py_pcha_installed()
  # coerce to matrix
  if(!is.matrix(data)) data = as.matrix(data)

  # check that ks do not exceed dimensions when simplex is true
  if(nrow(data) < (max(ks) - 1) & isTRUE(simplex)) stop("simplex = TRUE but number of archetypes (",
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
##' @description \code{fit_pch_bootstrap()} uses resampling data with or without replacement (bootstrapping) to find robust positions of archetypes of a polytope (Principal Convex Hull) describing the data. This function uses \code{fit_pch_resample()}.
##' @param n number of samples to be taken when bootstraping
##' @param sample_prop either NULL or the proportion of dataset that should be included in each sample. If NULL the polytope fitting algorithm is run n times on the same data which is useful for evaluating how often the algorithm gets stuck in local optima.
##' @param type one of s, m, cmq. s means single core processing using lapply. m means multi-core parallel procession using parLapply. cmq means multi-node parallel processing on a computing cluster using clustermq package. "See also" for details.
##' @param clust_options list of options for parallel processing. The default for "m" is list(cores = parallel::detectCores()-1, cluster_type = "PSOCK"). The default for "cmq" is list(memory = 2000, template = list(), n_jobs = 10, fail_on_error = FALSE). Change these options as required.
##' @param seed seed for reproducible random number generation. Works for all types of processing.
##' @param replace should resampling be with replacement? passed to \link[base]{sample.int}
##' @param average average archetype positions and varexpl? By default FALSE, return all fits to resampled data.
##' @param var_in_dims calculate variance in dimensions across archetypes in addition to variance in archetypes.
##' @param normalise_var normalise variance in position of archetypes by variance in data in each dimention
##' @param reference align archetypes (order_type="align") against reference polytope based on all data (TRUE), or align to the polytope from the first bootstrap iteration (FALSE). Second option can save time for large datasets
##' @return \code{fit_pch_bootstrap()} object of class b_pch_fit (list) containing the same elements as pch_fit, but each is either a list of pch_fit elements (e.g. list of n number of XC matrices) or a vector (which pch_fit element is one number).
##' @export fit_pch_bootstrap
fit_pch_bootstrap = function(data, n = 3, sample_prop = NULL, check_installed = T,
                             type = c("s", "m", "cmq")[1], clust_options = list(),
                             seed = 235, replace = FALSE, average = FALSE,
                             order_type = c("cosine", "side", "align")[3],
                             var_in_dims = TRUE, normalise_var = TRUE,
                             reference = FALSE, ...) {
  if(check_installed) .py_pcha_installed()
  # coerce to matrix
  if(!is.matrix(data)) data = as.matrix(data)

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
      ind = align_arc(ref_XC, res$pch_fits$XC[[i]])$ind
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

  # sum variance in position of each archetype across dimensions
  res$var_vert = colSums(res$var)
  # normalise variance in archetype positions by variance of data
  if(isTRUE(normalise_var)) res$var_vert = res$var_vert / sum(row_var_data)
  # find mean of variances of all archetypes to get a single number
  res$total_var = mean(res$var_vert)

  # normalise variance in archetype positions by variance of data in each dimension
  if(isTRUE(normalise_var)) res$var = res$var / matrixStats::rowVars(data)

  # find mean variance in position for each dimension across archetypes
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
  # calculate average XC matrix, set other C and S to NA because cells come from different samples
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
  res_aver$method_res = NA
  res_aver$total_var = res$total_var
  res_aver$summary = res$summary
  class(res_aver) = "pch_fit"
  res_aver
}

##' @rdname fit_pch
##' @name fit_pch_resample
##' @description \code{fit_pch_resample()} takes one sample of the data and fits a polytope (Principal Convex Hull) to data. This function uses \code{fit_pch()}.
##' @param i iteration number
##' @param replace fit_pch_resample/fit_pch_bootstrap: TRUE and FALSE are passed to \code{\link[base]{sample.int}} for density-based sampling, "geo_sketch" tells to sample examples using geosketch method preserving geometry of the data. See \link[ParetoTI]{geo_sketch}.
##' @return \code{fit_pch_resample()} object of class pch_fit
##' @export fit_pch_resample
fit_pch_resample = function(i = 1, data, sample_prop = NULL, replace = FALSE,
                            var_in_dims = var_in_dims,
                            normalise_var = normalise_var, ...) {

  # do resampling of the data -------------------------------------------------
  if(!is.null(sample_prop)){
    if(data.table::between(sample_prop[1], 0, 1)){

      if(replace != "geo_sketch" & replace %in% c(TRUE, FALSE)){

        # sample random examples with sample.int, density-based sampling -------
        col_ind = sample.int(ncol(data), round(ncol(data) * sample_prop[1], digits = 0),
                             replace = replace)

      } else if(replace == "geo_sketch") {

        # sample random examples with sample.int, geometry-based sampling ------
        col_ind = geo_sketch(data, N = as.integer(round(ncol(data) * sample_prop[1], digits = 0)),
                             use_PCs = FALSE, k = "auto", seed = NULL,
                             alpha = 0.1, max_iter = 200, verbose = 0,
                             check_installed = FALSE)

      } else stop("replace should be TRUE, FALSE or geo_sketch")

      data = data[, col_ind]

    } else stop("sample_prop should be NULL or a number between 0 and 1")
  }

  # fit polytope ---------------------------------------------------------------
  ParetoTI::fit_pch(data = data, ..., var_in_dims = var_in_dims,
                    normalise_var = normalise_var, check_installed = FALSE)

}

##' @rdname fit_pch
##' @name randomise_fit_pch
##' @description \code{randomise_fit_pch()} helps answer the question "how likely you are to observe data shaped as polytope given no relationship between variables?" This function disrupts the relationships between variables, find a polytope that best describes the data. It calculates variance explained, t-ratio of polytope volume to convex hull volume. Optionally, it uses bootstrapping of data points to measure variance in archetype positions. This function uses \code{randomise_fit_pch1()}. Number of archetypes is determined automatically from \code{arc_data}.
##' @param arc_data observed polytope that describes the data. Class pch_fit, b_pch_fit and k_pch_fit objects generated by \code{fit_pch()} or \code{fit_pch_bootstrap()}
##' @param n_rand number of randomisation samples
##' @return \code{randomise_fit_pch()} S3 "r_pch_fit" object (list) containing: 1. data table with columns variance explained (rand_varexpl), t-ratio (rand_t_ratio), variance in positions of archetypes (total_var) and number of archetypes (k) for \code{n_rand} samples; 2. empirical p-value for those parameters for each number of archetypes.
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
  # coerce to matrix
  if(!is.matrix(data)) data = as.matrix(data)

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
##' @param bootstrap_N randomise_fit_pch1() and k_fit_pch(): integer, number of bootstrap samples on random data to measure variability in archetype positions. Set to NA fit once to all data instead of bootstraping. When this option is positive seed bootstrap_seed and sample_prop must be provided.
##' @param bootstrap_type \code{randomise_fit_pch1()} and \code{k_fit_pch()}: parallel processing type when bootstraping. Caution: avoid nested parallel processing, do not use "m" and "cmq" inside other parallel functions.
##' @param return_data return randomised data?
##' @param return_arc return archetype positions in randomised data?
##' @param bootstrap_average \code{randomise_fit_pch1()}: average positions and summary statistics when bootstraping? Passed to \code{average} argument of \code{fit_pch_bootstrap()}. When multiple ks this defaults to TRUE.
##' @return \code{randomise_fit_pch1()}: list containing function call, summary of the sample, optional data and optional position of archetypes.
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
    # fitting many shapes (ks archetypes) or fitting one shape
    arc_data = ParetoTI::k_fit_pch(data = data, ks = ks, check_installed = FALSE,
                                   bootstrap = FALSE, seed = seed,
                                   volume_ratio = volume_ratio,
                                   var_in_dims = var_in_dims,
                                   normalise_var = normalise_var,
                                   # prevent whole pipeline from failing
                                   # when PCHA doesn't converge:
                                   converge_else_fail = FALSE, ...)

  } else if(isTRUE(as.integer(bootstrap_N) >= 1)) { # bootstraping --------------
    if(length(ks) > 1){ # fitting many shapes (ks archetypes)
      arc_data = ParetoTI::k_fit_pch(data = data, ks = ks, check_installed = TRUE,
                                     bootstrap = TRUE, bootstrap_N = bootstrap_N,
                                     sample_prop = sample_prop,
                                     bootstrap_type = bootstrap_type, seed = seed,
                                     replace = replace,
                                     var_in_dims = var_in_dims,
                                     normalise_var = normalise_var,
                                     ..., volume_ratio = volume_ratio)
    } else { # fitting many shapes (ks archetypes)
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
##' @param dist_metric how to describe distance to archetypes. Currently euclidean distance and archetype weights are implemented. Archetypes weights come from \code{arc_data$S} matrix so bootstrapped archetype positions cannot be used - only single fit to all data.
##' @param rank rank by distance metric (euclidean distance)?
##' @return \code{merge_arch_dist()} list: data.table with samples in rows (speficied by sample_id column) and features in columns (including distance to archetypes); and character vectors listing names of archetypes, feature columns and colData columns.
##' @export merge_arch_dist
##' @import data.table
merge_arch_dist = function(arc_data, data, feature_data,
                           colData = NULL, colData_id,
                           dist_metric = c("euclidean", "arch_weights")[1],
                           rank = TRUE){

  if(!is(arc_data, "pch_fit")) stop("arc_data should contain a single fit (pch_fit object): use fit_pch() or fit_pch_bootstrap() followed by average_pch_fits()")

  # find distance of data points to archetypes ---------------------------------
  if(dist_metric == "arch_weights"){

    dist = 1 - t(arc_data$S)
    if(is.null(colnames(dist))) {
      colnames(dist) = paste0("archetype_", seq_len(ncol(dist))) # archetype names
    }
    if(is.null(colnames(data))) {
      rownames(dist) = colnames(data) # data point names
    }

  } else if(dist_metric == "euclidean") {

    dist = arch_dist(data, arc_data$XC, dist_metric = dist_metric)

  }

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
    colData_col = colnames(colData)[!colnames(colData) %in% colData_id]
  } else {
    colData_col = NULL
  }
  list(data = features, arc_col = arc_col,
       features_col = rownames(feature_data),
       colData_col = colData_col)
}

.c_pch_fit_list = function(pch_fit_list){
  # remove failed fits
  pch_fit_list = pch_fit_list[!vapply(pch_fit_list, is.null,
                                      FUN.VALUE = logical(1))]

  if(is(pch_fit_list[[1]], "k_pch_fit")){
    pch_fit_list = lapply(pch_fit_list, function(pch) pch$pch_fits)
  }

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
       method_res = lapply(pch_fit_list, function(pch) pch$method_res),
       total_var = vapply(pch_fit_list, function(pch) pch$total_var, FUN.VALUE = numeric(1L)),
       summary = data.table::rbindlist(lapply(pch_fit_list, function(pch) pch$summary)))
}

##' @rdname fit_pch
##' @name c.k_pch_fit
##' @export c.k_pch_fit
c.k_pch_fit = function(...) {

  pch_fit_list = list(...)

  res = list()
  # combine function calls
  res$call = lapply(pch_fit_list, function(pch) pch$call)

  # combine k_pch_fit objects
  pch_list = lapply(pch_fit_list, function(pch) pch$pch_fits)
  res$pch_fits = list(XC = unlist(lapply(pch_list, function(pch) pch$XC), recursive = F),
                      S = unlist(lapply(pch_list, function(pch) pch$S), recursive = F),
                      C = unlist(lapply(pch_list, function(pch) pch$C), recursive = F),
                      SSE = unlist(lapply(pch_list, function(pch) pch$SSE)),
                      varexpl = unlist(lapply(pch_list, function(pch) pch$varexpl)),
                      hull_vol = unlist(lapply(pch_list, function(pch) pch$hull_vol)),
                      arc_vol = unlist(lapply(pch_list, function(pch) pch$arc_vol)),
                      t_ratio = unlist(lapply(pch_list, function(pch) pch$t_ratio)),
                      var_vert = unlist(lapply(pch_list, function(pch) pch$var_vert), recursive = F),
                      var_dim = data.table::rbindlist(lapply(pch_list, function(pch) pch$var_dim)),
                      total_var = unlist(lapply(pch_list, function(pch) pch$total_var)),
                      summary = data.table::rbindlist(lapply(pch_list, function(pch) pch$summary)))
  res$pch_fits$summary = NULL

  # combine summaries
  res$summary = rbindlist(lapply(pch_fit_list, function(pch) pch$summary))

  class(res) = "k_pch_fit"

  res
}

## .find_archetype_order orders archetypes by cosine relative to the unit vector c(1, 1)
## Cosine is negative when archetype vector is to the left (counter-clockwise) of the unit vector
## Solution taken from here: https://stackoverflow.com/questions/13221873/determining-if-one-2d-vector-is-to-the-right-or-left-of-another
## XC2 is a matrix of dim(dimensions, archetype) that has only 2 dimentions
.find_archetype_order = function(XC2, noc, order_type = c("cosine", "side", "align")[3]){
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
  cat("k representative archetypes (S3 class pch_fit)\n\n")
  cat("Summary:  varexpl = variance explained by data as weighted sum of archetypes\n          t_ratio = volume of polytope formed by archetypes / volume of convex hull\n          total_var = total variance in positions of archetypes\n            (by bootstraping, mean acrooss archetypes)\n\n")
  print(res$summary)
  cat("\n archetype positions stored in 'object$XC' \n dim(variables/dimentions, archetypes)")
}

##' @rdname fit_pch
##' @name print.k_pch_fit
##' @export print.k_pch_fit
print.k_pch_fit = function(res){
  cat("k representative archetypes (S3 class k_pch_fit)\n\n")
  cat("Summary:  varexpl = variance explained by data as weighted sum of archetypes\n          t_ratio = volume of polytope formed by archetypes / volume of convex hull\n          total_var = total variance in positions of archetypes\n            (by bootstraping, mean acrooss archetypes)\n\n")
  print(res$summary)
  cat("\n archetype positions stored in 'object$pch_fits$XC' \n dim(variables/dimentions, archetypes)")
}

##' @rdname fit_pch
##' @name print.b_pch_fit
##' @export print.b_pch_fit
print.b_pch_fit = function(res){
  cat("k representative archetypes \n calculated by bootstraping data points (cells) (S3 class b_pch_fit)\n\n")
  cat("Average summary:  varexpl = variance explained by data as weighted sum of archetypes\n          t_ratio = volume of polytope formed by archetypes / volume of convex hull\n          total_var = total variance in positions of archetypes\n            (by bootstraping, mean acrooss archetypes)\n\n")
  print(res$summary)
  cat("\nAll polytopes (fit to bootstrapped data), summary:\n")
  print(res$pch_fits$summary)
  cat("\n archetype positions stored in 'object$pch_fits$XC' \n dim(variables/dimentions, archetypes)")
}

##' @rdname fit_pch
##' @name print.r_pch_fit
##' @export print.r_pch_fit
print.r_pch_fit = function(res){
  cat("Background distribution of k representative archetypes \nin data with no relationships between variables (S3 class r_pch_fit)\n\n")
  cat("N randomisation trials: ",res$rand_dist[, max(repl)],"\n\n")
  cat("Summary of best-fit polytope to observed data (including p-value):\n\n")
  print(res$obs_dist)
  cat("\n          varexpl = variance explained by data as weighted sum of archetypes\n          t_ratio = volume of polytope formed by archetypes / volume of convex hull\n          total_var = total variance in positions of archetypes\n            (by bootstraping, mean acrooss archetypes)\n\n")
  cat("Function call:\n")
  print(res$call)
}

##' @rdname fit_pch
##' @name plot.r_pch_fit
##' @export plot.r_pch_fit
##' @description plot.r_pch_fit() plot distribution of t-ratio, total variance in archetype positions and variance explained used to calculate empirical p-value
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
##' @description plot_dim_var() plot distribution of variance in archetype positions in each dimension and corresponding empirical p-values
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
    xlab("Variance in archetype positions") +
    ylab("Dimension where variance is measured")
}

.cacl_var_in_dims = function(res, data, var_in_dims, normalise_var,
                             XC_array = NULL, most_extreme = NULL) {
  if(var_in_dims){
    # calculate variance in archetype positions in each dimension
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
    if(normalise_var) {
      # union of 100 most extreme points (20%, n%/ n dimensions) in each dimention
      if(is.null(most_extreme)){
        data_var = matrix(matrixStats::rowVars(data),
                          nrow = 1, ncol = nrow(data))
      } else {
        n_examples = ncol(data)
        n_dim = nrow(data)
        n_top_points = floor(n_examples * most_extreme * (1 / n_dim) / 2)
        data_var = apply(data, 1, function(x) {
          top_points = order(x)[c(seq(1, n_top_points),
                                  seq(n_examples - n_top_points + 1,
                                      n_examples))]
          top_points = unique(as.integer(top_points))
          var(x[top_points])
        })
        data_var = matrix(data_var, nrow = 1, ncol = nrow(data))
      }
      res$var_dim = res$var_dim / data_var
    }
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

##' @rdname fit_pch
##' @name pch_fit
##' @description \code{pch_fit()} a constructor function for the "pch_fit" class
##' @export pch_fit
pch_fit = function(XC, S, C, SSE, varexpl, arc_vol, hull_vol, t_ratio, var_vert,
                   var_dim, total_var, summary, call) {

  # add integrity checks

  # create object
  structure(list(XC = XC, S = S, C = C, SSE = SSE, varexpl = varexpl,
                 arc_vol = arc_vol, hull_vol = hull_vol, t_ratio = t_ratio,
                 var_vert = var_vert, var_dim = var_dim, total_var = total_var,
                 summary = summary, call = call),
            class = "pch_fit")
}

##' @rdname fit_pch
##' @name subset.k_pch_fit
##' @description \code{subset.k_pch_fit()} subsetting method for k_pch_fit object, allows selecting specific number of archetypes from k_fit_pch function output. E.g. useful to do PCA and UMAP projection for selected number of archetypes to aviod recalculation.
##' @export subset.k_pch_fit
subset.k_pch_fit = function(arc_data, ks) {

  # find the index of elements to be subset
  subset_ind = arc_data$summary[, k %in% ks]

  # subset elements
  arc_data$summary = arc_data$summary[subset_ind,]
  # list to subset:
  sub_names = names(arc_data$pch_fits)
  sub_names = sub_names[sub_names != "var_dim"]

  for (sub_name in sub_names) {
    arc_data$pch_fits[[sub_name]] = arc_data$pch_fits[[sub_name]][subset_ind]
  }

  if(nrow(arc_data$pch_fits$var_dim) == length(subset_ind)){
    arc_data$pch_fits$var_dim = arc_data$pch_fits$var_dim[subset_ind,]
  }

  # if subsetting only one element simplify k_pch_fit into pch_fit
  if(sum(subset_ind) == 1) {
    summary = arc_data$summary
    call = arc_data$call
    arc_data = arc_data$pch_fits
    arc_data$XC = arc_data$XC[[1]]
    arc_data$S = arc_data$S[[1]]
    arc_data$C = arc_data$C[[1]]
    arc_data$summary = summary
    arc_data$call = call
    class(arc_data) = "pch_fit"
  }

  arc_data

}
