##' Permute each column of a matrix or sample from it's ECDF (with replacement)
##' @rdname rand_var
##' @name rand_var
##' @author Vitalii Kleshchevnikov
##' @description This function uses permutations or sampling with replacement to disrupt relationships between variables but keep the distribution of each variable constant. Variables should be in columns of a matrix.
##' @param x matrix to be permuted
##' @param replace should sampling be with replacement? Passed to [base](sample.int).
##' @param prob a vector of probability weights for obtaining the elements of the vector being sampled. Passed to [base](sample.int).
##' @return Matrix of the same dimention as the original matrix but with values in each column permuted.
##' @export rand_var
##' @seealso \code{\link{}}, \code{\link{}}
##' @examples
##' set.seed(4354)
##' N = 1000
##' x = matrix(rnorm(N * 10 * N), N * 10, N)
##'
rand_var = function(x, replace = FALSE, prob = NULL){
  nrows = nrow(x)
  coln = colnames(x)
  rown = rownames(x)
  # x_rand = vapply(x, FUN = function(x_i, nrows) x_i[sample.int(nrows)], nrows,
  #                FUN.VALUE = numeric(nrows)) # often gives
  # "Error: vector memory exhausted (limit reached?)",
  # it is not the slowest part of PCHA computation so that's OK
  # apply sample.int to each column
  x_rand = apply(x, 2, FUN = function(x_i, nrows) {
    x_i[sample.int(nrows, nrows, replace = replace, prob = prob)]
  }, nrows) # median time (ms) mat of dim(10000, 1000) 742.5874

  colnames(x_rand) = coln
  rownames(x_rand) = rown
  x_rand
}
