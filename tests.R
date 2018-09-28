# Update man pages and install the package
library(roxygen2)
library(devtools)
document()
install()

devtools::load_all()

devtools::install_github("vitkl/ParetoTI", dependencies = T)
library(ParetoTI)

# test fitPCH
set.seed(4354)
N = 100
data = matrix(rnorm(N * 10 * N), N * 10, N)
dim(data)
archetypes = fit_pch(data, noc = as.integer(3), delta = 0.1)

# install python and / or py_pcha module
install_py_pcha(method = "conda")
ParetoTI::install_py_pcha(method = "virtualenv")
reticulate::py_discover_config("py_pcha")

set.seed(4354)
N = 500
data = matrix(rnorm(N * 10 * N), N * 10, N)
dim(data)

microbenchmark::microbenchmark({
  # Fit a polytope with 3 vertices to data matrix
  arc = fit_pch(data, noc=as.integer(3), delta=0.1)
}, {
  # Fit the same polytope 3 times without subsampling to test convergence of the algorithm.
  arc_rob = fit_pch_robust(data, n = 3, subsample = NULL,
                           noc=as.integer(3), delta=0.1)
}, {
  # Fit the 10 polytopes to subsampled datasets each time looking at 70% of examples.
  arc_data = fit_pch_robust(data, n = 10, subsample = 0.7, seed = 2543,
                            noc=as.integer(3), delta=0.1)
}, {
  # Use local parallel processing to fit the 10 polytopes to subsampled datasets each time looking at 70% of examples.
  arc_data = fit_pch_robust(data, n = 10, subsample = 0.7, seed = 2543,
                            noc=as.integer(3), delta=0.1, type = "m")
}, times = 3)#, {
#  # Use parallel processing on a computing cluster with clustermq to fit the 10 polytopes to subsampled datasets each time looking at 70% of examples.
arc_data_cmq = fit_pch_robust(data, n = 20, subsample = 0.95, seed = 2543,
                          noc=as.integer(3), delta=0.1, type = "cmq")
#})


# comparing RPCHA and R-python inferface
devtools::install_github("gokceneraslan/RPCHA", dependencies = T)
devtools::install_github("vitkl/ParetoTI", dependencies = T)
library(PCHA)
library(ParetoTI)
##' @param distance matrix of dim(examples, archetypes)
perform = function(distance, slope = rep(-0.003, 3), y_intercept = 1){
  rep(slope, nrow(slope))
  performance = matrix(NA, nrow = nrow(distance), ncol = ncol(distance))
  for(i in seq_len(ncol(distance))){
    distance[,i] = distance[,i] * slope[i] + y_intercept
  }
  performance[performance < 0] = 0
  performance
}
comp_fitness = function(performance) {
  fitness = rowSums(performance)
  fitness[rowSums(performance == 0) > 0] = 0
  fitness = fitness / max(fitness)
}
# test fitPCH
set.seed(4355)
archetypes = generate_arc(arc_coord = list(c(5, 0), c(-10, 15), c(-30, -20)),
                          mean = 0, sd = 1, N_dim = 2)
data = generate_data(archetypes, N_examples = 1e4, jiiter = 0.04, size = 0.9)

distance = sqrt(arch_dist(data, archetypes))
arc_distance = sqrt(arch_dist(archetypes, archetypes))
#performance = perform(distance, slope = c(-0.3, -0.2, -0.4), y_intercept = 10)
#plot(distance[,1], performance[,1])
#fitness = comp_fitness(performance)
#plot(distance[,1], fitness)

#within_polytope = rowSums(distance <= max(arc_distance)) == ncol(distance) &
#  rowSums(distance) <= sum(arc_distance[upper.tri(arc_distance)])
#sum(within_polytope)
#data = data[fitness > 0.9,]

s = svd(t(data))
PCs_data = t(diag(s$d) %*% t(s$v))
proj_arc = archetypes %*% t(s$u)
PCs = rbind(PCs_data, proj_arc)
data_w_arc = rbind(data, archetypes)
color = c(rep("black", nrow(PCs_data)), rep("red", nrow(proj_arc)))

plot(PCs[,1], PCs[,2], col = color)
plot(data_w_arc[,1], data_w_arc[,2], col = color)
plot(data[,1], data[,2])

plot(archetypes[,1], archetypes[,2])
plot(archetypes[,2], archetypes[,3])
plot(archetypes[,3], archetypes[,4])
s$u %*% t(archetypes)
data = matrix(rnorm(N * 10 * N), N * 10, N)
microbenchmark::microbenchmark({
  # Fit a polytope with 3 vertices to data matrix
  arc = ParetoTI::fit_pch(data, noc=as.integer(3), delta=0.1, conv_crit = 1e-06, maxiter = 500, verbose = FALSE)
}, {
  #RPCHA
  res = PCHA::PCHA(data, noc=as.integer(3), delta=0.1, conv_crit = 1e-06, maxiter = 500, verbose = FALSE)
}, times = 10)

devtools::load_all()

# set directory for user libraries, update pip, setuptools, wheel in that environment
#export PYTHONUSERBASE=$vk7/software/python_libs/
#python -m pip install --user -i https://pypi.python.org/simple -U pip distribute
#python -m pip install --user -i https://pypi.python.org/simple --upgrade pip setuptools wheel virtualenv -U pip --user distribute
