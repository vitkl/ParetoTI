##' Permute each column of a matrix
##' @rdname rand_var
##' @name rand_var
##' @author Vitalii Kleshchevnikov
##' @description This function uses permutations to disrupts relationships between variables but keep the distribution of each variable constant. Variables should be in columns.
##' @param X matrix to be permuted
##' @return matrix of the same dimention as the original matrix but with values in each column permuted
##' @export rand_var
##' @seealso \code{\link{}}, \code{\link{}}
##' @examples
##' set.seed(4354)
##' N = 1000
##' X = matrix(rnorm(N * 10 * N), N * 10, N)
rand_var = function(X){
  # apply sample.int to each column
  nrows = nrow(X)
  coln = colnames(X)
  rown = rownames(X)
  # X_rand = vapply(X, FUN = function(X_i, nrows) X_i[sample.int(nrows)], nrows,
  #                FUN.VALUE = numeric(nrows)) # often gives
  # "Error: vector memory exhausted (limit reached?)",
  # it is not the slowest part of PCHA computation so that's OK
  X_rand = apply(X, 2, FUN = function(X_i, nrows) X_i[sample.int(nrows)], nrows) # median time (ms) mat of dim(10000,1000) 742.5874
  colnames(X_rand) = coln
  rownames(X_rand) = rown
  X_rand
}
