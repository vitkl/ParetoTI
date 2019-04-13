##' Generate random data within specified archetypes
##' @rdname generate_arc
##' @name generate_arc
##' @author Vitalii Kleshchevnikov
##' @description \code{generate_arc()} generates the matrix of archetypes using specified coordinates (some noise added).
##' @param arc_coord list of archetype coordinates, one numeric vector (length of N dimensions) per each archetype.
##' @param mean mean of random distribution added to arc_coord
##' @param sd standard deviationn of random distribution added to arc_coord
##' @return \code{generate_arc()} object of class "random_arc" (similar to "pch_fit"), element XC is a matrix of archetypes of dim(dimensions, archetypes)
##' @export generate_arc
##' @seealso \code{\link[ParetoTI]{fit_pch}}, \code{\link[ParetoTI]{arch_dist}}
##' @examples
##' # Random data that fits into the triangle
##' set.seed(4355)
##' archetypes = generate_arc(arc_coord = list(c(5, 0), c(-10, 15), c(-30, -20)),
##'                           mean = 0, sd = 1)
##' data = generate_data(archetypes$XC, N_examples = 1e4, jiiter = 0.04, size = 0.9)
##' # Find Euclidian distance between data points and archetypes
##' distance = sqrt(arch_dist(data, archetypes))
##' # Find Euclidian distance between archetypes
##' arc_distance = sqrt(arch_dist(archetypes, archetypes))
generate_arc = function(arc_coord = list(c(5, 0), c(-10, 15), c(-30, -20)),
                        mean = 0, sd = 1){
  n_arc = length(arc_coord)
  N_dim = length(arc_coord[[1]])
  archetypes = matrix(rnorm(n_arc * N_dim, mean, sd), n_arc, N_dim)
  for (i in seq_len(n_arc)) {
    archetypes[i,] = archetypes[i,] + arc_coord[[i]]
  }

  archetypes = list(XC = t(archetypes),
                    S = NA, C = NA, SSE = NA,
                    varexpl = NA, call = match.call())
  class(archetypes) = "random_arc"
  archetypes
}

##' @rdname generate_arc
##' @name generate_data
##' @param archetypes matrix of archetypes of dim(dimensions, archetypes)
##' @param N_examples number of examples to be generated
##' @param jiiter add noise to weigth so that data is not a perfect polytope (e.g. triangle, see examples)
##' @param size scale the data within a polytope
##' @description \code{generate_data()} produces matrix of random data that fits a polytope defined by archetypes by multiplying position of archetypes by random weigths (that sum to 1)
##' @return \code{generate_data()} matrix of archetypes of dim(dimensions, examples)
##' @export generate_data
generate_data = function(archetypes, N_examples = 1e4, jiiter = 0.1, size = 1) {
  n_arc = ncol(archetypes)
  weights = matrix(runif(N_examples * n_arc, 0, 1), N_examples, n_arc)
  weights = (weights / rowSums(weights)) * size
  noise = matrix(rnorm(N_examples * n_arc, 0, jiiter), N_examples, n_arc)
  weights = weights + noise
  t(weights %*% t(archetypes))
}
