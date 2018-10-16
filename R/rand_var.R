##' Permute each column of a matrix or sample from it's ECDF (with replacement)
##' @rdname rand_var
##' @name rand_var
##' @author Vitalii Kleshchevnikov
##' @description This function uses permutations or sampling with replacement to disrupt relationships between variables but keep the distribution of each variable constant. Variables should be in columns of a matrix.
##' @param x matrix to be permuted
##' @param MARGIN which dimension of a matrix to permute? Passed to \code{\link[base]{apply}}.
##' @param replace should sampling be with replacement? Passed to \code{\link[base]{sample.int}}.
##' @param prob a vector of probability weights for obtaining the elements of the vector being sampled. Passed to \code{\link[base]{(sample.int}}.
##' @return Matrix of the same dimention as the original matrix but with values in each column permuted.
##' @export rand_var
##' @examples
##' set.seed(4354)
##' N = 1000
##' x = matrix(rnorm(N * 10 * N), N * 10, N)
##' x_rand = rand_var(x, MARGIN = 2, replace = FALSE, prob = NULL)
rand_var = function(x, MARGIN = 2, replace = FALSE, prob = NULL){
  nrows = nrow(x)
  ncols = ncol(x)
  coln = colnames(x)
  rown = rownames(x)
  # x_rand = vapply(x, FUN = function(x_i, nrows) x_i[sample.int(nrows)], nrows,
  #                FUN.VALUE = numeric(nrows)) # often gives
  # "Error: vector memory exhausted (limit reached?)",
  # it is not the slowest part of PCHA computation so that's OK
  # apply sample.int to each column
  if(MARGIN == 1) n = ncols else if(MARGIN == 2) n = nrows else stop("MARGIN not 1 or 2")
  x_rand = apply(x, MARGIN = MARGIN, FUN = function(x_i, n) {
    x_i[sample.int(n, n, replace = replace, prob = prob)]
  }, n) # median time (ms) mat of dim(10000, 1000) 742.5874
  x_rand = t(x_rand)

  colnames(x_rand) = coln
  rownames(x_rand) = rown
  x_rand
}
