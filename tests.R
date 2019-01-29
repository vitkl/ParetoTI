# Update man pages and install the package
library(roxygen2)
library(devtools)
document()
install()

devtools::load_all("../ParetoTI/")

#install.packages("BiocManager") # for installing BioConductor dependencies
BiocManager::install("vitkl/ParetoTI", dependencies = T)
BiocManager::install("vitkl/ParetoTI", dependencies = c("Depends", "Imports", "LinkingTo", "Suggests"))

library(ParetoTI)

devtools::install_url("http://spams-devel.gforge.inria.fr/hitcounter2.php?file=file/36615/spams-R-v2.6-2017-03-22.tar.gz")

# Let's try SPAMS package method (written in c++)
library(ParetoTI)
library(spams)
# Random data that fits into the triangle (2D)
set.seed(4355)
archetypes = generate_arc(arc_coord = list(c(5, 0), c(-10, 15), c(-30, -20)),
                          mean = 0, sd = 1)
data = generate_data(archetypes$XC, N_examples = 1e4, jiiter = 0.04, size = 0.99)
microbenchmark::microbenchmark({
  s_a <- spams.archetypalAnalysis(X = data, p = 3, Z0 = NULL, returnAB = FALSE,
                                  robust=FALSE, epsilon=1e-3, computeXtX=TRUE,
                                  stepsFISTA=0, stepsAS=100, randominit=TRUE,
                                  numThreads=-1)},{
                                    # Compare to PCHA
  s_p <- fit_pch(data, noc = as.integer(3), delta = 0, conv_crit = 1e-03)},
  times = 3)
archetypes$XC
s_a
s_p$XC
align_arc(s_a, archetypes$XC)
align_arc(s_p$XC, archetypes$XC)

res = lapply(seq_len(10), function(i) fit_pch(data, noc = as.integer(3),
                                        delta = 0, conv_crit = 0.3*1e-03))
sapply(res, function(res_1) align_arc(res_1$XC, archetypes$XC)$dist)
sapply(res, function(res_1) align_arc(res_1$XC, archetypes$XC)$ind)

grid = expand.grid(c(1, 3, 5, 8), c(1e-02, 1e-03, 1e-04, 1e-05, 1e-06))
grid = sort(grid[, 1] * grid[, 2])
res = lapply(grid, function(i) fit_pch(data, noc = as.integer(3),
                                              delta = 0, conv_crit = i))
plot(log10(grid), sapply(res, function(res_1) align_arc(res_1$XC, archetypes$XC)$dist))
grid
sapply(res, function(res_1) align_arc(res_1$XC, archetypes$XC)$dist)
sapply(res, function(res_1) align_arc(res_1$XC, archetypes$XC)$ind)

rob = lapply(grid, function(i) fit_pch_bootstrap(data, n = 20, type = "m", sample_prop = 0.65, seed = 2543, noc=as.integer(3), delta=0, conv_crit = i))

plot(log10(grid), sapply(res, function(res_1) rob$total_var))
grid
sapply(res, function(res_1) rob$total_var)

# create python environment and install py_pcha module
install_py_pcha(method = "conda")
ParetoTI::install_py_pcha(method = "virtualenv")
reticulate::py_discover_config("py_pcha")

R.utils::setOption("ParetoTI_envname", "py_pcha_cython")
getOption("ParetoTI_envname")
library(ParetoTI)
library(ggplot2)
# Random data that fits into the triangle (2D)
set.seed(4355)
archetypes = generate_arc(arc_coord = list(c(5, 0), c(-10, 15), c(-30, -20)),
                          mean = 0, sd = 1)
data = generate_data(archetypes$XC, N_examples = 1e4, jiiter = 0.04, size = 0.99)
plot_arc(arch_data = archetypes, data = data,
         which_dimensions = 1:2) +
  theme_bw()
# Plot data as 2D density rather than points
plot_arc(arch_data = archetypes, data = data,
         which_dimensions = 1:2, geom = ggplot2::geom_bin2d) +
  theme_bw()

# Random data that fits into the triangle (3D)
set.seed(4355)
archetypes = generate_arc(arc_coord = list(c(5, 0, 4), c(-10, 15, 0), c(-30, -20, -5)),
                          mean = 0, sd = 1)
data = generate_data(archetypes$XC, N_examples = 1e4, jiiter = 0.04, size = 0.99)
plot_arc(arch_data = archetypes, data = data,
         which_dimensions = 1:3)

# test fitPCH
arc_data = fit_pch(data, noc = as.integer(3), delta = 0)
plot_arc(arch_data = arc_data, data = data,
         which_dimensions = 1:3)
plot_arc(arch_data = arc_data, data = data,
         which_dimensions = 1:2) +
  theme_bw()


speed_test = microbenchmark::microbenchmark({
  # Fit a polytope with 3 vertices to data matrix
  arc = fit_pch(data, noc=as.integer(3), delta=0)
}, {
  # Fit the same polytope 3 times without subsampling to test convergence of the algorithm.
  arc_rob_conv = fit_pch_bootstrap(data, n = 3, sample_prop = NULL,
                                   noc=as.integer(3), delta=0)
}, {
  # Fit the 20 polytopes to subsampled datasets each time looking at 65% of examples.
  arc_data_rob = fit_pch_bootstrap(data, n = 20, sample_prop = 0.65, seed = 2543,
                                   noc=as.integer(3), delta=0)
}, {
  # Use local parallel processing to fit the 20 polytopes to subsampled datasets each time looking at 65% of examples.
  arc_data_rob_m = fit_pch_bootstrap(data, n = 20, sample_prop = 0.65, seed = 2543,
                                     noc=as.integer(3), delta=0, type = "m")
}, times = 5)
speed_test_cmq = microbenchmark::microbenchmark({
  # Use parallel processing on a computing cluster with clustermq to fit the 20 polytopes to subsampled datasets each time looking at 65% of examples.
  arc_data_rob_cmq = fit_pch_bootstrap(data, n = 200, sample_prop = 0.65, seed = 2543,
                                       noc = as.integer(3),
                                       delta = 0, type = "cmq",
                                       clust_options = list(memory = 1000, n_jobs = 10))
}, times = 5)

align_arc(arc$XC, archetypes$XC)
align_arc(average_pch_fits(arc_data_rob_m)$XC, archetypes$XC)

## Does bootstrap average give a better approximation of true vertex positions?
pcha_bench = function(conv_crit) {
  library(ParetoTI)
  res = sapply(1:50, function(i){
    archetypes = generate_arc(arc_coord = list(c(5, 0), c(-10, 15), c(-30, -20)),
                              mean = 0, sd = 1)
    data = generate_data(archetypes$XC, N_examples = 1e4, jiiter = 0.04, size = 0.99)
    arc = fit_pch(data, noc=as.integer(3), delta=0)
    arc_data_rob_m = fit_pch_bootstrap(data, n = 20, sample_prop = 0.65, seed = NULL,
                                       noc=as.integer(3), delta=0, type = "s",
                                       conv_crit = conv_crit)
    vect = c(align_arc(arc$XC, archetypes$XC)$dist,
             align_arc(average_pch_fits(arc_data_rob_m)$XC, archetypes$XC)$dist)
    names(vect) = c("all", "bootstrap")
    vect
  })
  res = as.data.table(res, keep.rownames = "type")
  res = melt.data.table(res, id.vars = "type")
  res[, conv_crit := conv_crit]
  res
}
# generate convergence values
conv_vals = expand.grid(c(1, 3, 5, 8), c(1e-02, 1e-03, 1e-04, 1e-05, 1e-06))
conv_vals = sort(conv_vals[, 1] * conv_vals[, 2])
res = clustermq::Q(pcha_bench, c(1e-03, 1e-04, 1e-05, 1e-06),
                   memory = 2000, n_jobs = 4, seed = 4534)
res = rbindlist(res)
saveRDS(res, "../../PCHA_accuracy_bootstrap_vs_all_data/conv_crit_res.rds")
# yes!
ggplot(res, aes(x = value, fill = type)) +
  geom_density(aes(y=..count.. + 1), alpha = 0.5) +
  facet_wrap( ~ conv_crit) +
  scale_y_log10() +
  theme_bw()

plot_arc(arch_data = arc_data_rob_cmq, data = data,
         which_dimensions = 1:3, line_size = 1.5)
plot_arc(arch_data = arc_data_rob_cmq, data = data,
         which_dimensions = 1:2, line_size = 1) +
  theme_bw()

# test function for different k
arc_ks = k_fit_pch(data, ks = 1:4, check_installed = T, delta=0)
plot_arc(arch_data = arc_ks, data = data,
         which_dimensions = 1:3, type = "all", arch_size = 2,
         colors = c("#D62728", "#1F77B4", "#2CA02C", "#17BED0", "grey"))
plot_arc(arch_data = arc_ks, data = data,
         which_dimensions = 1:2, type = "all", arch_size = 2,
         colors = c("#D62728", "#1F77B4", "#2CA02C", "#17BED0", "grey")) +
  theme_bw()

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
                          mean = 0, sd = 1)
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
  arc = ParetoTI::fit_pch(data, noc=as.integer(3), delta=0, conv_crit = 1e-06, maxiter = 500, verbose = FALSE)
}, {
  #RPCHA
  res = PCHA::PCHA(data, noc=as.integer(3), delta=0, conv_crit = 1e-06, maxiter = 500, verbose = FALSE)
}, times = 10)
res$XC = as.matrix(res$XC)
class(res) = "pch_fit"
plot_arc(arch_data = res, data = data,
         which_dimensions = 1:3)

# set directory for user libraries, update pip, setuptools, wheel in that environment
#export PYTHONUSERBASE=$vk7/software/python_libs/
#python -m pip install --user -i https://pypi.python.org/simple -U pip distribute
#python -m pip install --user -i https://pypi.python.org/simple --upgrade pip setuptools wheel virtualenv -U pip --user distribute

library(AnnotationHub)
hub = AnnotationHub()
selected_name = names(query(hub, "OrgDb"))[mcols(query(hub, "OrgDb"))$taxonomyid == 9606]
org_db = hub[[selected_name]]
org_db

pkgdown::build_site()


## illustrations for rotation report


plot_arc(arch_data = arc_data, data = data,
         which_dimensions = 1:2,
         nudge = c(0, 0.1),
         colors = c("#FFC003", "#D62728")) +
  theme_bw() +
  xlab("PC1") + ylab("PC2")

archetypes = generate_arc(arc_coord = list(c(5, 0, 4), c(-20, 10, 0), c(-10, -10, -5)),
                          mean = 0, sd = 1)
data = generate_data(archetypes$XC, N_examples = 5*1e2, jiiter = 0.04, size = 0.9)
arc_data = fit_pch(data, noc = as.integer(3), # number of vertices = 3
                   delta = 0)
plot_arc(arch_data = arc_data, data = data,
         which_dimensions = 1:2,
         nudge = c(0, 0.1),
         colors = c("#747171", "#D62728")) +
  theme_bw() +
  xlab("PC1") + ylab("PC2")


## testing cython PCHA
python -m timeit -n 1 -r 100 --verbose -s 'import numpy as np; from py_pcha.PCHA import PCHA; dimensions = 15; examples = 1000; X = np.random.random((dimensions, examples))' 'XC, S, C, SSE, varexpl = PCHA(X, noc=3, delta=0, maxiter = 2000)'



# profile.py

import pstats, cProfile

import numpy as np; from py_pcha.PCHA import PCHA;
dimensions = 15; examples = 10000; X = np.random.random((dimensions, examples))

cProfile.runctx("PCHA(X, noc=3, delta=0, maxiter = 2000)", globals(), locals(), "Profile.prof")

s = pstats.Stats("Profile.prof")
s.strip_dirs().sort_stats("time").print_stats()


X; noc = 3; I=None; U=None; delta=0; verbose=False; conv_crit=1E-6; maxiter=2000

export CFLAGS=-I/Users/vk7/anaconda3/envs/py_pcha_cython/lib/python2.7/site-packages/numpy/core/include
python setup.py build_ext --inplace

import numpy as np
from py_pcha.PCHA import PCHA

dimensions = 15
examples = 100
X = np.random.random((dimensions, examples))

XC, S, C, SSE, varexpl = PCHA(X, noc=3, delta=0)


X_array = X
X = np.asmatrix(X)
N, M = X.shape
I = range(M)
U = range(M)

SST = np.sum(np.diag(X[:, U] * X[:, U].T))
