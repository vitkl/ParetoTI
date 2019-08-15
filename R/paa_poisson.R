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
  #data_mean = log(data_mean)
  # data sd in each dimension
  data_sd = matrix(DelayedMatrixStats::colSds(DelayedArray::DelayedArray(data)),
                   nrow = n_arc, ncol = n_dim, byrow = TRUE) * scale_data_sd
  #data_sd = sqrt(data_sd) # if this evaluates to < 0 everything breaks

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
  archetypes = normal(data_mean, data_sd, truncation = c(0, Inf))

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
    distribution(data) = poisson(mu)

    # create model
    model = model(weights, archetypes, precision = precision)

  } else { # with cell covariates / gene beta  ================
  }

  # return environment of this function storing the greta arrays and the model
  environment()
}

##' @rdname paa_poisson
##' @name paa_normal_free
##' @description \code{paa_normal_free()} creates Archetypal Analysis greta model with Normal likelihood function:
##' @description data = Normal(mean = mu, sd)
##' @description mu = weights * archetypes_mean
##' @description sd = weights * archetypes_sd
##' @description archetypes_mean = Normal(mean = mean(data), sd = sd(data))
##' @description archetypes_sd = Exponential(rate = 1 / sd(data))
##' @export paa_normal_free
paa_normal_free = function(data,
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
  #data_mean = log(data_mean)
  # data sd in each dimension
  data_sd = matrix(DelayedMatrixStats::colSds(DelayedArray::DelayedArray(data)),
                   nrow = n_arc, ncol = n_dim, byrow = TRUE) * scale_data_sd
  #data_sd = sqrt(data_sd) # if this evaluates to < 0 everything breaks

  # make this a part of greta calculation ================
  weight_alpha = as_data(matrix(rep(weight_alpha_prior, n_arc),
                                nrow = 1, ncol = n_arc)) # dirichlet alpha prior on archetype weights of data points

  # convert covariates to greta arrays ================
  data = as_data(data)
  if(!is.null(covar)) covar = as_data(covar)

  # define priors on parameters ================
  # archetype weights of data points: dim(n_points, n_arc)
  weights = dirichlet(alpha = weight_alpha,
                      n_realisations = n_points, dimension = n_arc)
  # archetypes: dim(n_arc, n_dim)
  archetypes_mean = normal(data_mean, data_sd)
  archetypes_sd = exponential(1 / data_sd)

  # extract variable to enable setting initial values on dirichet parameters  ================
  # internal greta functions
  get_node = .internals$greta_arrays$get_node
  as.greta_array = .internals$greta_arrays$as.greta_array

  # dig out the greta array corresponding to the unconstrained version of weights
  weights_node = get_node(weights)
  weights_var_node = weights_node$children[[1]]
  weights_var = as.greta_array(weights_var_node)

  # compute expected value ======================================
  # expected value = (mu = dim(data points, dimensions))

  if(is.null(covar)) { # without cell covariates  ================

    # compute averages ================
    mu = weights %*% archetypes_mean
    sd = weights %*% archetypes_sd

    # define the distribution over data
    distribution(data) = normal(mu, sd)

    # create model
    model = model(weights, archetypes_mean, archetypes_sd, precision = precision)

  } else { # with cell covariates / gene beta  ================
  }

  # return environment of this function storing the greta arrays and the model
  environment()
}

##' @rdname paa_poisson
##' @name paa_nb_free
##' @description \code{paa_nb_free()} creates Archetypal Analysis greta model with Negative binomial likelihood function:
##' @description data = NegativeBinomial(prob = prob, size = size)
##' @description prob = size / (size + mu)
##' @description mu = weights * archetypes_mean
##' @description size = weights * archetypes_size
##' @description archetypes_mean = Normal(mean = mean(data), sd = sd(data), truncation = c(0, Inf))
##' @description archetypes_size = Exponential(rate = (mean(data) + mean(data)^2) / var(data))
##' @export paa_nb_free
paa_nb_free = function(data,
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
  #data_mean = log(data_mean)
  # data sd in each dimension
  data_sd = matrix(DelayedMatrixStats::colSds(DelayedArray::DelayedArray(data)),
                   nrow = n_arc, ncol = n_dim, byrow = TRUE) * scale_data_sd
  #data_sd = sqrt(data_sd) # if this evaluates to < 0 everything breaks

  # make this a part of greta calculation ================
  weight_alpha = as_data(matrix(rep(weight_alpha_prior, n_arc),
                                nrow = 1, ncol = n_arc)) # dirichlet alpha prior on archetype weights of data points

  # convert covariates to greta arrays ================
  data = as_data(data)
  if(!is.null(covar)) covar = as_data(covar)

  # define priors on parameters ================
  # archetype weights of data points: dim(n_points, n_arc)
  weights = dirichlet(alpha = weight_alpha,
                      n_realisations = n_points, dimension = n_arc)
  # archetypes: dim(n_arc, n_dim)
  archetypes_mean = normal(data_mean, data_sd, truncation = c(0, Inf))
  archetypes_size = exponential(1 / data_sd) # exponential((data_mean + data_mean ^ 2) / data_sd ^ 2)

  # extract variable to enable setting initial values on dirichet parameters  ================
  # internal greta functions
  get_node = .internals$greta_arrays$get_node
  as.greta_array = .internals$greta_arrays$as.greta_array

  # dig out the greta array corresponding to the unconstrained version of weights
  weights_node = get_node(weights)
  weights_var_node = weights_node$children[[1]]
  weights_var = as.greta_array(weights_var_node)

  # compute expected value ======================================
  # expected value = (mu = dim(data points, dimensions))

  if(is.null(covar)) { # without cell covariates  ================

    # compute averages ================
    size = weights %*% archetypes_size
    mu = weights %*% archetypes_mean
    prob = size / (size + mu)

    # define the distribution over data
    distribution(data) = negative_binomial(size = size, prob = prob)

    # create model
    model = model(weights, archetypes_mean, archetypes_size, precision = precision)

  } else { # with cell covariates / gene beta  ================
  }

  # return environment of this function storing the greta arrays and the model
  environment()
}

##' @rdname paa_poisson
##' @name poisson_mixture
##' @description \code{poisson_mixture()} creates Archetypal Analysis greta model with Poisson likelihood function:
##' @description data = Poisson(mu)
##' @description mu = exp(weights * archetypes)
##' @description archetypes = Normal(mean = mean(data), sd = sd(data))
##' @param scale_data_sd by what factor to scale sd(data)
##' @param ... unused - temporary hack to make paa models easier to exchange
##' @export poisson_mixture
poisson_mixture = function(data,
                           n_arc = 7,              # number of achetypes
                           weight_alpha_prior = 0.8,
                           scale_data_sd = 0.2,
                           covar = NULL,
                           cell_weights = FALSE,
                           precision = c("double", "single"),
                           ...
) {

  # define derived parameters ================
  n_points = nrow(data)        # number of data points (e.g. cells)
  n_dim = ncol(data)           # number of dimensions (e.g. genes)
  # data average in each dimension
  data_mean = matrix(Matrix::colMeans(data), nrow = n_arc, ncol = n_dim, byrow = TRUE)
  #data_mean = log(data_mean)
  # data sd in each dimension
  data_sd = matrix(DelayedMatrixStats::colSds(DelayedArray::DelayedArray(data)),
                   nrow = n_arc, ncol = n_dim, byrow = TRUE) * scale_data_sd
  #data_sd = sqrt(data_sd) # if this evaluates to < 0 everything breaks

  # make this a part of greta calculation ================
  weight_alpha = as_data(matrix(weight_alpha_prior, nrow = 1, ncol = n_arc)) # dirichlet alpha prior on archetype weights of data points

  # convert covariates to greta arrays ================
  data = as_data(data)
  if(!is.null(covar)) covar = as_data(covar)

  # define priors on parameters ================
  # archetype weights of data points: dim(n_points, n_arc)
  if(cell_weights){

    weights = dirichlet(alpha = weight_alpha,
                        n_realisations = n_points, dimension = n_arc)
    weights2 = greta_array(weights, dim = c(n_arc * n_points, 1)) %*% ones(1, n_dim)
    weights2 = greta_array(weights2, dim = c(n_arc, n_points, n_dim))

  } else {

    weights = dirichlet(alpha = weight_alpha,
                        n_realisations = 1, dimension = n_arc)
    weights2 = t(weights)

  }
  # archetypes: dim(n_arc, n_dim)
  archetypes = normal(data_mean, data_sd, truncation = c(0, Inf))

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

    # compute mixture components ================
    distr = list()
    for (arc in seq_len(n_arc)) {
      distr = c(distr, list(poisson(ones(n_points, 1) %*% archetypes[arc,])))
    }

    # define the mixture distribution over data
    mixture_distribution = greta:::mixture_distribution
    distribution(data) <- greta:::distrib("mixture", distr, weights2, dim = NULL)

    # create model
    model = model(weights, archetypes, precision = precision)

  } else { # with cell covariates / gene beta  ================
  }

  # return environment of this function storing the greta arrays and the model
  environment()
}

# # simulate a mixture of poisson random variables and try to recover the
# # parameters with a Bayesian model
# set.seed(1)
# x <- cbind(c(rpois(800, 30),
#              rpois(200, 100),
#              rpois(400, 50)),
#            c(rpois(800, 200),
#              rpois(200, 150),
#              rpois(400, 10)))
# n_points = nrow(x)       # number of data points (e.g. cells)
# n_dim = ncol(x)
# n_arc = 3
# weight_alpha = 0.5
#
# plot_arc(data = t(x))
#
# weights = dirichlet(alpha = matrix(weight_alpha, nrow = 1, ncol = n_arc),
#                                n_realisations = 1, dimension = n_arc)
# weights = t(weights)
# #weights = greta_array(weights, dim = c(n_arc * n_points, 1)) %*% ones(1, n_dim)
# #weights = greta_array(weights, dim = c(n_arc, n_points, n_dim))
#
# # weights <- t(extraDistr::rdirichlet(n_points, rep(weight_alpha, n_arc)))
# # weights = array(weights, dim = c(n_arc * n_points, 1)) %*% array(1, c(1, n_dim))
# # weights = array(weights, dim = c(n_arc, n_points, n_dim))
#
# rates <- normal(0, 10, truncation = c(0, Inf), dim = c(n_arc, n_dim))
# distr = list()
# for (arc in seq_len(n_arc)) {
#   distr = c(distr, list(poisson(ones(n_points, 1) %*% rates[arc,])))
# }
# mixture_distribution = greta:::mixture_distribution
# distribution(x) <- greta:::distrib("mixture", distr, weights, dim = NULL)
#
# m <- model(rates)
# plot(m)
# #draws_rates <- mcmc(m, n_samples = 500)
# opt_result = opt(m, adam(), max_iterations = 5000)
# opt_result$iterations
#
# # check the mixing probabilities after fitting using calculate()
# # (you could also do this within the model)
# draws_weights <- calculate(weights, draws_rates)
#
# # get the posterior means
# summary(draws_rates)$statistics[, "Mean"]
# summary(draws_weights)$statistics[, "Mean"]
#
# # PCHA
# arc_data = fit_pch(t(x), noc = 3L)
# plot_arc(arc_data, t(x))
#
# # Mixture plot
# mixt_data = copy(arc_data)
# mixt_data$XC = t(array(summary(draws_rates)$statistics[, "Mean"], c(3, 2)))
# mixt_data$XC = t(opt_result$par$rates)
# plot_arc(mixt_data, t(x))
