##' Project polytope fit and data in low dimentions (PCA)
##' @param arc_data objects of class "pch_fit", "b_pch_fit", "k_pch_fit" storing the position of archetypes and other data produced by \code{\link[ParetoTI]{fit_pch}}(). arc_data$XC is matrix of dim(dimensions, archetypes) or list where each element is XC matrix from an independent run of the polytope fitting algorithm.
##' @param data matrix of data in which archetypes/polytope were found, dim(variables/dimentions, examples)
##' @param n_dim number of principal component dimensions
##' @param proj_matrix crossprod(proj_matrix, data) projects data/archetypes into principal component space
##' @examples
##' # Random data that fits into the triangle
##' set.seed(4355)
##' arc_data = generate_arc(arc_coord = list(c(5, 0), c(-10, 15), c(-30, -20)),
##'                           mean = 0, sd = 1)
##' data = generate_data(arc_data$XC, N_examples = 1e4, jiiter = 0.04, size = 0.9)
project = function(arc_data, data, n_dim = 3, proj_matrix = NULL){
  s = irlba::irlba(data, 50)
  cds@reducedDimK = tcrossprod(diag(s$d), s$v)
  s = svd::propack.svd(X = data, neig = n_dim)

  crossprod(s$u, logcounts(data_run))
  tcrossprod(diag(s$d), s$v)
  all.equal(tcrossprod(diag(s$d), s$v), crossprod(diag(s$d), t(s$v)))

  crossprod(svd_pc$u, counts(sc))[,1:5];
  tcrossprod(diag(svd_pc$d[1:3]), svd_pc$v)[,1:5];
  tcrossprod(diag(base_pc$d[1:3]), base_pc$v)[,1:5];
  t(reducedDim(scater_pc, 1))[1:3, 1:5]
}
