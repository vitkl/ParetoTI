##' Match/align arhetypes
##' @rdname align_arc
##' @name align_arc
##' @description \code{align_arc()} matches archetypes in arc2 to arc1 by exhaustively testing all possible pairings.
##' @details Exhaustive search works better than optimisation using SANN method for k=2:9 (163.0747ms vs 354.6311ms for 3 dimensions, k=2:5 are faster than 2ms). This optimisation method essentially samples the space rather than doing gradient descent. This is because total distance as a function of pairings of archetypes is not a differentiable function. For this reason, optimisation will often fail to find correct pairing.
##' @param arc1 reference matrix of archetype positions dim(dimensions, archetypes)
##' @param arc2 matrix of archetype positions dim(dimensions, archetypes) to be aligned with arc1
##' @return \code{align_arc()} list containing: total distance between archetypes in arc1 and arc2 (dist), integer vector specifying indices (columns) of arc2 that match archetypes in arc1 (ind)
##' @export align_arc
align_arc = function(arc1, arc2) {
  if(!isTRUE(all.equal(dim(arc1), dim(arc2)))) stop("align_arc() trying to match different number of archetypes")
  arc_distance = arch_dist(arc1, arc2)
  arc2_combs = gen_permut(nrow(arc_distance))
  distances = vapply(seq(1, nrow(arc_distance)), function(i, arc2_combs){
    arc_distance[i, arc2_combs[,i]]
  }, FUN.VALUE = numeric(nrow(arc2_combs)), arc2_combs)
  distances = rowSums(distances)
  list(dist = min(distances), ind = arc2_combs[which.min(distances),])
}

##' @rdname align_arc
##' @name gen_permut
##' @description \code{gen_permut()} used for exhaustive search generates a matrix of all possible permutations of n elements. Each row is a different permutation.
##' @param n number of element to permute
##' @details \code{gen_permut()} function is taken from https://stackoverflow.com/questions/11095992/generating-all-distinct-permutations-of-a-list-in-r
##' @return \code{gen_permut()} a matrix of all possible gen_permut (in rows)
##' @export gen_permut
gen_permut = function(n){
  if(n==1){
    return(matrix(1))
  } else {
    sp = gen_permut(n-1)
    p = nrow(sp)
    A = matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] = cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}

##' @rdname align_arc
##' @name align_arc_opt
##' @description \code{align_arc_opt()} does the same as \code{align_arc()} but using optimisation (method = "SANN") instead of exhaustively trying all combinations.
##' @param control options for SANN methods. See \link[stats]{optim}
##' @return same as \code{align_arc()}
##' @export align_arc_opt
##' @examples
##' # Generate data
##' set.seed(4355)
##' archetypes = generate_arc(arc_coord = list(c(5, 0), c(-10, 15), c(-30, -20)),
##'                           mean = 0, sd = 1)
##' data = generate_data(archetypes, N_examples = 1e4, jiiter = 0.04, size = 0.9)
##' dim(data)
##' # fit polytopes to 2 subsamples of the data
##' arc_data = fit_pch_bootstrap(data, n = 2, sample_prop = 0.65, seed = 2543,
##'                         order_type = "align", noc = as.integer(6),
##'                         delta = 0, type = "s")
##' # align archetypes by exhausive search
##' align_arc(arc_data$pch_fits$XC[[1]], arc_data$pch_fits$XC[[2]])
##' # align archetypes by optimisation
##' set.seed(123)
##' align_arc_opt(arc_data$pch_fits$XC[[1]], arc_data$pch_fits$XC[[2]])
##' # restart multiple times to see if it finds real match
##' res = lapply(seq(1, 100), function(i) align_arc_opt(arc_data$pch_fits$XC[[1]], arc_data$pch_fits$XC[[2]]))
##' # show best pairing
##' res[[which.min(sapply(res, function(r) r$dist))]]
##'
##' # Measure which approach is faster depending on k (uncomment)
##' #for (i in 2:12) {
##' #    # fit polytopes to 2 subsamples
##' #    arc_data = fit_pch_bootstrap(data, n = 2, sample_prop = 0.65, seed = 2543,
##' #                             order_type = "align", noc = as.integer(i),
##' #                             delta = 0, type = "s")
##' #    print(i)
##' #    # measure how much time each method takes
##' #    print(microbenchmark::microbenchmark(align_arc(arc_data$pch_fits$XC[[1]], arc_data$pch_fits$XC[[2]]),
##' #                                         align_arc_opt(arc_data$pch_fits$XC[[1]], arc_data$pch_fits$XC[[2]]),
##' #                                         times = 10))
##' #}
align_arc_opt = function(arc1, arc2, control = list(maxit = 10000, temp = 2000,
                                                    REPORT = 500)) {
  arc_distance = arch_dist(arc1, arc2)
  res = optim(par = c(1:nrow(arc_distance)),  # Initial sequence
              fn = function(sq) {  # Target function / distance function
                sum(vapply(seq(1, nrow(arc_distance)), function(i, sq){
                  arc_distance[i, sq[i]]
                }, FUN.VALUE = numeric(1), sq))
              },
              gr = function(sq) {  # Generate new candidate sequence
                idx = seq(1, NROW(arc_distance))
                for(i in seq(1, 4)){
                  changepoints = sample(idx, size = 2, replace = FALSE)
                  tmp = sq[changepoints[1]]
                  sq[changepoints[1]] = sq[changepoints[2]]
                  sq[changepoints[2]] = tmp
                  sq
                }
                sq
                #sample(sq, replace = FALSE)
              },
              method = "SANN",
              control = control)
  list(dist = res$value, ind = res$par)
}
