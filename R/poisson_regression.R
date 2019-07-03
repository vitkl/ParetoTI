##' Poisson regression model to find effects of sample covariates on variables
##' @name poisson_regression
##' @rdname poisson_regression
##' @author Vitalii Kleshchevnikov
##' @description Poisson regression models gene expression (Y) as a function of gene mean and sample covariates (X)
##' mu = beta * X
##' Y ~ Poisson(mu)
##' @param data data matrix (data points * dimensions), can be dense and sparse matrix, SummarizedExperiment/SingleCellExperiment, Seurat (counts slot is used). \code{poisson_regression()} accepts only a dense matrix at the moment (limitation of greta).
##' @param family character naming the data distribution
##' @param covar matrix (data points * covariates) or vector of column names (for compute_gene_deviance() and SingleCellExperiment, Seurat) containing covariates affecting expression in addition to gene mean (coverage, batch). Adding this will find genes whose deviance (residuals) is unexplained both by these covariates and Poisson noise (covar = NULL tests Poisson noise alone).
##' @param precision argument for \link[greta]{model}. Use "single" for large datasets to reduce memory footprint
##' @param verbose logical, plot greta model structure and print messages?
##' @param optimiser method to use for finding regression coefficients and deviance when adding covariates see \link[greta]{opt}
##' @param max_iterations number of iterations to run optimiser for.
##' @param tolerance the numerical tolerance for the solution, the optimiser stops when the (absolute) difference in the joint density between successive iterations drops below this level.
##' @return \code{compute_gene_deviance()}: list containing the deviance vector with dimension names (genes) as names, beta coefficient matrix (dimensions * coeffs) and greta model used to compute those. For SingleCellExperiment the same object with beta coeffecients and deviance as rowData is returned. For Seurat the same object is returned updated with beta coeffecients and deviance in \code{Seurat::GetAssay(obj, "RNA")@meta.features}.
##' @export compute_gene_deviance
##' @examples
##' # Use fake data as example
##' # Random data that fits into the triangle
##' set.seed(4355)
##' arc_data = generate_arc(arc_coord = list(c(7, 3, 10), c(12, 17, 11), c(30, 20, 9)),
##'                         mean = 0, sd = 1)
##' data = generate_data(arc_data$XC, N_examples = 1e4, jiiter = 0.04, size = 0.9)
##' # Take Poisson sample with the mean defined by each entry of the data matrix
##' # (this create Poisson-distributed positive integer data)
##' data = matrix(rpois(length(data), (data)), nrow(data), ncol(data))
##'
##' # Compute deviance from the mean (residuals for Poisson data)
##' dev = compute_gene_deviance(t(data))
##' # As you can see, the third dimension has lowest deviance
##' dev
##' # because the vertices of the triangle have almost identical position in third dimension.
##' plot_arc(arch_data = arc_data, data = data,
##'          which_dimensions = 1:3, data_alpha = 0.5)
##' # You can use deviance to find which dimension have variability to be explained with Archetypal Analysis
##'
##' # Create a probabilistic Poisson regression model with greta
##' # to study effects of covariates on Poisson data (requires greta installed)
##' \dontrun{
##' model = poisson_regression(t(data),
##'                      covar = matrix(rnorm(ncol(data)), ncol(data), 1))
##' # plot the structure of tensorflow computation graph
##' plot(model$model)
##'
##' # find parameters using adam optimiser
##' res = greta::opt(model$model, optimiser = greta::adam(), max_iterations = 500)
##' # did the model converge before 500 iterations?
##' res$iterations
##' # Value of Poisson negative log likelihood (see greta documentation for details)
##' res$value
##' # View beta parameters for each dimension (columns), log(mean) in the first row,
##' # covariate coefficients in the subsequent rows
##' res$par$beta
##' }
compute_gene_deviance = function(data,
                                 family = "poisson",
                                 covar = NULL,
                                 precision = c("double", "single"),
                                 verbose = FALSE,
                                 optimiser = greta::adam(),
                                 max_iterations = 5000,
                                 tolerance = 1e-06
) {

  # extract data =================================================================== #
  if (is(data, "SummarizedExperiment")) {

    y = Matrix::t(assay(data, "counts"))

  } else if(is(data, "Seurat")) {

    y = Matrix::t(data[["RNA"]]@counts)

  } else if(is(data, "matrix") | is(data, "dgCMatrix")) {
    y = data
  } else stop("data should be one of matrix, dgCMatrix, SummarizedExperiment, SingleCellExperiment, Seurat")

  if(!is.null(covar)){ # use greta model with covariates (slow)

    # extract covariates from SummarizedExperiment, SingleCellExperiment, Seurat
    if (is.character(covar) & is(data, "SummarizedExperiment")) {

      covar = as.matrix(colData(data)[, covar])
      # add bias term to find averages
      covar = cbind(matrix(log(mu), nrow = ncol(y), ncol = 1), covar)

    } else if(is.character(covar) & is(data, "Seurat")) {

      covar = as.matrix(seu@meta.data[, covar])
      # add bias term to find averages
      covar = cbind(matrix(log(mu), nrow = ncol(y), ncol = 1), covar)

    }

    # greta does not understand sparse matrix yet so convert to matrix
    y = as.matrix(y)

    # create poisson regression model & find parameters ============================== #
    m = poisson_regression(y, covar, precision)

    # visualise computation graph
    if(verbose) print(plot(m$model))

    # find parameters by optimisation
    opt_res = opt(m$model, optimiser = optimiser, max_iterations = max_iterations,
                  tolerance = tolerance, adjust = FALSE,
                  hessian = FALSE)

    # computer deviance ============================================================== #
    # Extract expected reads and beta
    mu = opt_res$par$mu
    beta = t(opt_res$par$beta)
    colnames(beta) = c("mean", colnames(covar))

  } else { # use analytical solution

    mu = Matrix::colMeans(y) # compute mean
    beta = matrix(log(mu), nrow = ncol(y), ncol = 1)
    mu = matrix(rep(mu, nrow(y)), nrow = nrow(y), byrow = TRUE)

    m = NULL
    colnames(beta) = "mean"
  }

  # name beta matrix
  rownames(beta) = colnames(y)
  colnames(beta) = paste0("beta_", colnames(beta))

  # Compute deviance
  term1 = y / mu
  term1[term1 == 0] = 1
  term1 = y * log(term1) # in the limit y*log(y) for y=0 -> 0
  deviance = 2 * Matrix::colSums(term1 - y + mu)
  names(deviance) = colnames(y)

  # return depending on input ====================================================== #
  if (is(data, "SummarizedExperiment")) {

    rowData(data)[, "deviance"] = deviance
    for (n in colnames(beta)) {
      rowData(data)[, n] = beta[rownames(data), n]
    }

    return(data) # sce with beta and deviance

  } else if(is(data, "Seurat")) {

    data[["RNA"]] = Seurat::AddMetaData(object = data[["RNA"]],
                                        metadata = deviance[rownames(data)],
                                        col.name = "deviance")
    for (n in colnames(beta)) {
      data[["RNA"]] = Seurat::AddMetaData(object = data[["RNA"]],
                                          metadata = beta[rownames(data), n],
                                          col.name = n)
    }

    return(data)

  } else if(is(data, "matrix") | is(data, "dgCMatrix")) {

    return(list(deviance = deviance, beta = beta, model = m))

  }

}
##' @return \code{poisson_regression()}: R environment containing the model and parameters as greta arrays
##' @param beta_mean prior mean for coefficients
##' @param beta_sd prior sd for coefficients, use small values to regularise (e.g. penalise coefficients that deviate too much from 0)
##' @export poisson_regression
poisson_regression = function(data, covar = NULL,
                              beta_mean = 0,
                              beta_sd = 3,
                              precision = c("double", "single")){

  # define parameters
  n_beta = 1 # number of beta coefficients in the linear model
  if(!is.null(covar)) n_beta = n_beta + ncol(covar)
  beta = greta::normal(beta_mean, beta_sd, dim = c(n_beta, ncol(data))) # beta coefficients

  # define predictors
  x = matrix(rep(1, nrow(data)), nrow = nrow(data))
  if(!is.null(covar)) x = cbind(x, covar)
  x = greta::as_data(x)

  data = greta::as_data(data)

  # define distribution over data
  mu = exp(greta::`%*%`(x, beta))
  greta::distribution(data) = greta::poisson(mu)

  # construct model
  model = greta::model(beta, mu, precision = precision)

  environment()
}
