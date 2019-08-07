##' Create probabilistic archetypal analysis model in greta (tensorflow backend)
##' @rdname paa_poisson
##' @name paa_poisson
##' @description \code{paa_poisson()} creates Archetypal Analysis greta model with Poisson likelihood function:
##' @description data = Poisson(mu)
##' @description mu = weights * archetypes
##' @description archetypes = c * data
##' @description Covariates can be added (experimental):
##' @description mu = exp(log(weights * archetypes) + beta * covar)
##' @param data matrix of positive count data: data points * dimensions. At the moment, greta does not support sparse matrix.
##' @param n_arc number of archetypes
##' @param weight_alpha_prior dirichet distrubution on archetype weights (sum to 1 for each data point). Value of 1 means flat distribution (all weights are equally likely to occur), less than 1 means sparse distribution (weights are likely to be close to 0 or close to 1), greater than 1 means the weights are less likely to be close to 0 or 1. To get a better idea on the meaning of this prior try sampling with different alpha values \code{alpha = c(0.8, 0.8); extraDistr::rdirichlet(1e5, alpha)}, where length(alpha) is the number of archetypes.
##' @param c_alpha_prior archetype construction weights c. Analogous to \code{weight_alpha_prior}, but set very close to 0 to enforce sparsity.
##' @param covar matrix of covariates affecting the expression (mu) in addition to archetypes: data points * n_covariates
##' @param precision argument for \link[greta]{model}. Use "single" for large datasets to reduce memory footprint
##' @return \code{poisson_regression()}: R environment containing the model and parameters as greta arrays
##' @export paa_poisson
##' @import greta
##' @examples
##' # create greta / tensorflow model
##' m = paa_poisson(data,              # data: data points * dimensions
##'                 n_arc = 7,              # number of achetypes
##'                 weight_alpha_prior = 0.8,
##'                 c_alpha_prior = 0.001,
##'                 precision = options$precision
##' )
##'
##' # visualise computation graph
##' plot(m$model)
##'
##' # solve model with optimisation
##' opt_res = opt(m$model,
##'                 optimiser = adam(learning_rate = 0.3),   # optimisation method used to find prior-adjusted maximum likelihood estimate: adam, l_bfgs_b,
##'                 max_iterations = 500,
##'                 tolerance = 1e-4, adjust = TRUE,
##'                 hessian = FALSE)
##'
##' # calculate archetypes using c matrix
##' archetypes = opt_res$par$c %*% data
paa_poisson = function(data,
                       n_arc = 7,              # number of achetypes
                       weight_alpha_prior = 0.8,
                       c_alpha_prior = 0.001,
                       covar = NULL,
                       precision = c("double", "single")
) {

  # define derived parameters ================
  n_points = nrow(data)        # number of data points (e.g. cells)
  n_dim = ncol(data)              # number of dimensions (e.g. genes)

  # make this a part of greta calculation ================
  weight_alpha = as_data(rep(weight_alpha_prior, n_arc)) # dirichlet alpha prior on archetype weights of data points
  c_alpha = as_data(rep(c_alpha_prior, n_points))   # dirichlet alpha prior on data point weights for archetype construction

  # convert data and covariates to greta arrays ================
  data = as_data(data)
  if(!is.null(covar)) covar = as_data(covar)

  # define priors on parameters ================
  # archetype weights of data points: dim(n_points, n_arc)
  weights = dirichlet(alpha = matrix(weight_alpha, nrow = 1, ncol = n_arc),
                      n_realisations = n_points, dimension = n_arc)
  # data point weights for archetype construction: dim(n_arc, n_points)
  c = dirichlet(alpha = matrix(c_alpha, nrow = 1, ncol = n_points),
                n_realisations = n_arc, dimension = n_points)

  # extract variable to enable setting initial values on dirichet parameters  ================
  # internal greta functions
  get_node = .internals$greta_arrays$get_node
  as.greta_array = .internals$greta_arrays$as.greta_array

  # dig out the greta array corresponding to the unconstrained version of weights
  weights_node = get_node(weights)
  weights_var_node = weights_node$children[[1]]
  weights_var = as.greta_array(weights_var_node)
  # dig out the greta array corresponding to the unconstrained version of c
  c_node = get_node(c)
  c_var_node = c_node$children[[1]]
  c_var = as.greta_array(c_var_node)


  # compute archetypes and expected value ======================================
  # expected value = (mu = dim(data points, dimensions))

  if(is.null(covar)) { # without cell covariates  ================

    # compute archetypes: dim(n_arc, n_dim) ================
    archetypes = c %*% data

    # compute averages ================
    mu = weights %*% archetypes

    # define the distribution over data
    distribution(data) = poisson(mu)

    # create model
    model = model(weights, c, precision = precision)

  } else { # with cell covariates / gene beta  ================

    # beta weights for covariates like nUMI or batch ================
    beta = normal(0, 1, dim = c(ncol(covar), n_dim))
    mu_covar = covar %*% beta

    # compute archetypes ================
    archetypes = c %*% data

    # compute averages ================
    mu_arch = weights %*% log(archetypes)
    mu = exp(mu_arch + mu_covar)

    # define the distribution over data ================
    distribution(data) = poisson(mu)

    # create model ================
    model = model(weights, c, beta, precision = precision)
  }

  # return environment of this function storing the greta arrays and the model
  environment()
}


##' @rdname paa_poisson
##' @name paa_poisson_free
##' @description \code{paa_poisson_free()} creates Archetypal Analysis greta model with Poisson likelihood function:
##' @description data = Poisson(mu)
##' @description mu = exp(weights * archetypes)
##' @description archetypes = Normal(mean = mean(data), sd = sd(data))
##' @description Covariates can be added (experimental):
##' @description mu = exp(weights * archetypes + beta * covar)
##' @param scale_data_sd by what factor to scale sd(data)
##' @param ... unused - temporary hack to make paa models easier to exchange
##' @export paa_poisson_free
paa_poisson_free = function(data,
                       n_arc = 7,              # number of achetypes
                       weight_alpha_prior = 0.8,
                       scale_data_sd = 0.2,
                       covar = NULL,
                       precision = c("double", "single"),
                       ...
) {

  # define derived parameters ================
  n_points = nrow(data)        # number of data points (e.g. cells)
  n_dim = ncol(data)           # number of dimensions (e.g. genes)
  # data average in each dimension
  data_mean = matrix(Matrix::colMeans(data), nrow = n_arc, ncol = n_dim, byrow = TRUE)
  data_mean = log(data_mean)
  # data sd in each dimension
  data_sd = matrix(DelayedMatrixStats::colSds(DelayedArray::DelayedArray(data)),
                   nrow = n_arc, ncol = n_dim, byrow = TRUE) * scale_data_sd
  data_sd = sqrt(data_sd) # if this evaluates to < 0 everything breaks

  # make this a part of greta calculation ================
  weight_alpha = as_data(rep(weight_alpha_prior, n_arc)) # dirichlet alpha prior on archetype weights of data points

  # convert covariates to greta arrays ================
  data = as_data(data)
  if(!is.null(covar)) covar = as_data(covar)

  # define priors on parameters ================
  # archetype weights of data points: dim(n_points, n_arc)
  weights = dirichlet(alpha = matrix(weight_alpha, nrow = 1, ncol = n_arc),
                      n_realisations = n_points, dimension = n_arc)
  # archetypes: dim(n_arc, n_dim)
  archetypes = normal(data_mean, data_sd)

  # extract variable to enable setting initial values on dirichet parameters  ================
  # internal greta functions
  get_node = .internals$greta_arrays$get_node
  as.greta_array = .internals$greta_arrays$as.greta_array

  # dig out the greta array corresponding to the unconstrained version of weights
  weights_node = get_node(weights)
  weights_var_node = weights_node$children[[1]]
  weights_var = as.greta_array(weights_var_node)

  # compute archetypes and expected value ======================================
  # expected value = (mu = dim(data points, dimensions))

  if(is.null(covar)) { # without cell covariates  ================

    # compute averages ================
    mu = weights %*% archetypes

    # define the distribution over data
    distribution(data) = poisson(exp(mu))

    # create model
    model = model(weights, archetypes, precision = precision)

  } else { # with cell covariates / gene beta  ================

    # beta weights for covariates like nUMI or batch ================
    beta = normal(0, 1, dim = c(ncol(covar), n_dim))
    mu_covar = covar %*% beta

    # compute averages ================
    mu_arch = weights %*% archetypes
    mu = exp(mu_arch + mu_covar)

    # define the distribution over data ================
    distribution(data) = poisson(mu)

    # create model ================
    model = model(weights, archetypes, beta, precision = precision)
  }

  # return environment of this function storing the greta arrays and the model
  environment()
}
