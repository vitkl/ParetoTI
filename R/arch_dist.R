##' Compute distance to archetypes
##' @rdname arch_dist
##' @name arch_dist
##' @author Vitalii Kleshchevnikov
##' @description \code{arch_dist()} calculates distance from every point to every archetype give matrices that contain this data.
##' @param data matrix of dim(examples, dimensions)
##' @param archetypes matrix of dim(archetypes, dimensions)
##' @return matrix of distances to archetypes of dim(examples, archetypes)
##' @export arch_dist
##' @seealso \code{\link[ParetoTI]{fit_pch}}, \code{\link[ParetoTI]{generate_arc}}, \code{\link[ParetoTI]{generate_data}}
##' @examples
##' # Triangle with sides of 3,4,5 units length
##' arch_dist(matrix(c(0,0),1,2), matrix(c(3,4),1,2))
##' # Random data that fits into the triangle
##' set.seed(4355)
##' archetypes = generate_arc(arc_coord = list(c(5, 0), c(-10, 15), c(-30, -20)),
##'                           mean = 0, sd = 1, N_dim = 2)
##' data = generate_data(archetypes, N_examples = 1e4, jiiter = 0.04, size = 0.9)
##' # Find Euclidian distance between data points and archetypes
##' distance = sqrt(arch_dist(data, archetypes))
##' # Find Euclidian distance between archetypes
##' arc_distance = sqrt(arch_dist(archetypes, archetypes))
arch_dist = function(data, archetypes){
  data = t(data)
  archetypes = lapply(seq_len(nrow(archetypes)), function(i) archetypes[i,])
  vapply(archetypes, function(arc, data){
    diff = arc - data
    colSums(diff^2)
  }, FUN.VALUE = numeric(ncol(data)), data)
}
