library(tensorflow)
tfe_enable_eager_execution()
library(greta)

variable()
# generate sparse predictors ====================
set.seed(1)
N_points = 500
N_dims = 50
sparse_factor = 4
dense_prop = N_dims / sparse_factor
sparse_prop = N_dims - dense_prop
X = matrix(sample(c(rpois(dense_prop * N_points, 6), # generate observations
                    rep(0, sparse_prop * N_points))), # add 10s
           N_dims, N_points)
X_sparse = as(X, "dgCMatrix")

# generate dense response ====================
y = matrix(rnorm(N_points), N_points, 1)

# greta regression model ====================
# prior on coefficients
beta = normal(0, 10, dim = c(N_dims, 1))
beta_mean = normal(0, 10, dim = c(1, 1))
sd = exponential(0.1)

# compute the mean using standard operations
#mu = beta_mean + t(X) %*% beta

# compute mean using the custom sparse operation
mu = sparse_matmul(Matrix::t(X_sparse), beta)

# define likelihood
distribution(y) = normal(mu, sd = sd)

m = model(beta, beta_mean)

plot(m)

# fit the model with opt() ====================
res = opt(m, optimiser = adam())

# fit the model with mcmc() ====================
draws = mcmc(m, chains = 4)
bayesplot::mcmc_trace(draws)



sparse_matmul = function(x_sparse, y) {

  if(!is(x_sparse, "dgTMatrix")) x_sparse = as(x_sparse, "dgTMatrix")

  # extract indices of non-sparse values
  ind = cbind(x_sparse@i, x_sparse@j)
  ind = split(ind, seq_len(nrow(ind)))
  names(ind) = NULL
  shape = c(dim(x_sparse))
  out_dim = c(nrow(x_sparse), ncol(y))
  x_sparse = as.numeric(x_sparse@x)

  greta:::op("sparse_matmul_op", x_sparse, y, operation_args = list(ind = ind,
                                                                    shape = shape),
             tf_operation = tf_sparse_matmul, dim = out_dim)

}

tf_sparse_matmul = function(x_sparse, y, ind, shape) {

  # covert dgCMatrix into sparse tensor
  sp_a = tf$sparse$SparseTensor(indices=ind, values=x_sparse, dense_shape=shape)

  #sparse_dense_matmult_batch(sp_a, y)
  tf$sparse$matmul(sp_a, y)

}

sparse_dense_matmult_batch = function(sp_a, b) {

  map_function = function(x){

    sparse_slice = tf$sparse$reshape(tf$sparse$slice(
      sp_a, c(x[1], 0L, 0L), c(1L, sp_a$dense_shape[1], sp_a$dense_shape[2])),
      c(sp_a$dense_shape[1], sp_a$dense_shape[2]))

    dense_slice = tf$reshape(tf$slice(
      b, c(x[2], 0L, 0L), c(1L, b$shape[1], b$shape[2])),
      c(b$shape[1], b$shape[2]))

    mult_slice = tf$sparse$matmul(sparse_slice, dense_slice)

    mult_slice
  }


  elems = list(tf$range(0L, sp_a$dense_shape[0], delta=1L, dtype=tf$int64),
               tf$range(0L, b$shape[0], delta=1L, dtype=tf$int64))
  print(elems)

  tf$map_fn(map_function, elems, dtype=tf$float64, back_prop=TRUE)
}

# load tensorflow and enable immediate execution
library(tensorflow)
tfe_enable_eager_execution()

set.seed(3248)
# create sparse tensor
dim = c(2, 4, 2) # dim
ind = sapply(dim, function(x) sample.int(x, 16, replace = T)) # random indices
sp_a = tf$sparse$SparseTensor(indices=ind, values=as.numeric(rep(2, nrow(ind))),
                              dense_shape=as.integer(dim))
sp_a = tf$dtypes$cast(sp_a, dtype = tf$float64) # change dtype (R converts to double by default)
# create dense tensor
b = tf$constant(array(rnorm(2 * 2 * 3), dim = c(2, 2, 3)))
# matrix multiply
sparse_dense_matmult_batch(sp_a, b)
