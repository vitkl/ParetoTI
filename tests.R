# Update man pages and install the package
library(roxygen2)
library(devtools)
document()
install()

unloadNamespace("ParetoTI")
devtools::load_all("../ParetoTI/")

#install.packages("BiocManager") # for installing BioConductor dependencies
BiocManager::install("vitkl/ParetoTI", dependencies = T)
BiocManager::install("vitkl/ParetoTI", dependencies = c("Depends", "Imports", "LinkingTo"))

library(ParetoTI)

# Random data that fits into the triangle
set.seed(4355)
arc_data = generate_arc(arc_coord = list(c(5, 1, 4), c(10, 15, 1), c(30, 20, 5)),
                          mean = 0, sd = 1)
data = generate_data(arc_data$XC, N_examples = 1e4, jiiter = 0.04, size = 0.9)
#arc_data = fit_pch_bootstrap(data, n = 10, sample_prop = 0.5, noc = as.integer(3))
arc_data = fit_pch(data, noc = as.integer(3))
# Plot
plot_arc(arc_data = arc_data, data = data,
         which_dimensions = 1:2, data_alpha = 0.5) +
  ggplot2::theme_bw()

# Project to PCs (in this case just rotate to align x-axis with
# the axis of most variation because the data is already 2D)
pcs = project_to_pcs(arc_data, data, n_dim = 3,
                     pc_method = c("svd", "irlba")[1],
                     zscore = F, log2 = F, offset = 2)
# Plot in PC coordinates
plot_arc(arc_data = pcs$arc_data, data = pcs$data,
         which_dimensions = 1:2, data_alpha = 0.5) +
  ggplot2::theme_bw()

# Project from PCs back to expression
projected = project_from_pc(pcs$arc_data, pcs$s,
                            undo_zscore = F, undo_log2 = F, offset = 2)

# Plot plot in projected coordinates
plot_arc(arc_data = projected, data = data,
         which_dimensions = 1:2, data_alpha = 0.5) +
  ggplot2::theme_bw()

devtools::install_url("http://spams-devel.gforge.inria.fr/hitcounter2.php?file=file/36615/spams-R-v2.6-2017-03-22.tar.gz")

# Let's try SPAMS package method (written in c++)
library(ParetoTI)
library(spams)
# Random data that fits into the triangle (2D)
set.seed(4355)
archetypes = generate_arc(arc_coord = list(c(5, 0), c(-10, 15), c(-30, -20)),
                          mean = 0, sd = 1)
data = generate_data(archetypes$XC, N_examples = 1e4, jiiter = 0.04, size = 0.99)
data = matrix(rnorm(2*1e4), 2, 1e4)
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
plot_arc(arc_data = archetypes, data = data,
         which_dimensions = 1:2) +
  ggplot2::theme_bw()
# Plot data as 2D density rather than points
plot_arc(arc_data = archetypes, data = data,
         which_dimensions = 1:2, geom = ggplot2::geom_bin2d) +
  ggplot2::theme_bw()

# Random data that fits into the triangle (3D)
set.seed(4355)
archetypes = generate_arc(arc_coord = list(c(5, 0, 4), c(-10, 15, 0), c(-30, -20, -5)),
                          mean = 0, sd = 1)
data = generate_data(archetypes$XC, N_examples = 1e3, jiiter = 0.04, size = 0.99)
plot_arc(arc_data = archetypes, data = data,
         which_dimensions = 1:3, data_alpha = 0.5)

# test fitPCH
arc_data = fit_pch(data, noc = as.integer(3), delta = 0)
plot_arc(arc_data = arc_data, data = data,
         which_dimensions = 1:3)
plot_arc(arc_data = arc_data, data = data,
         which_dimensions = 1:2, data_alpha = 0.5) +
  ggplot2::theme_bw()

colnames( arc_data$XC) = paste0("A_", 1:3)
rownames( arc_data$XC) = paste0("R_", 1:3)
rownames( data) = paste0("R_", 1:3)

# test projection to PCs
pcs = project_to_pcs(arc_data, data, n_dim = 3, pc_method = c("svd", "irlba")[1])
plot_arc(arc_data = pcs$arc_data, data = pcs$data,
         which_dimensions = 1:2, data_alpha = 0.5) +
  ggplot2::theme_bw()

# test projection from PCs
 pcs$s$u %*% pcs$arc_data$XC
 arc_data$XC

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
  arc_data_rob_m = fit_pch_bootstrap(data, n = 200, sample_prop = 0.65, seed = 2543,
                                     noc=as.integer(3), delta=0, type = "m")
}, times = 5)
speed_test_cmq = microbenchmark::microbenchmark({
  # Use parallel processing on a computing cluster with clustermq to fit the 20 polytopes to subsampled datasets each time looking at 65% of examples.
  arc_data_rob_cmq = fit_pch_bootstrap(data, n = 200, sample_prop = 0.65, seed = 2543,
                                       noc = as.integer(3),
                                       delta = 0, type = "cmq",
                                       clust_options = list(memory = 1000, n_jobs = 10))
}, times = 5)

library(ggplot2)
library(ParetoTI)
set.seed(4355)
archetypes = generate_arc(arc_coord = list(c(5, 0, 4, 1, 0, 6), c(-10, 15, 0, 0, 1, 4), c(-30, -20, -5, 1, 0, 5)),
                          mean = 0, sd = 1)
data = generate_data(archetypes$XC, N_examples = 1e3, jiiter = 0.04, size = 0.99)
arc_ks = k_fit_pch(data, ks = 2:6, check_installed = T,
                   bootstrap = T, bootstrap_N = 50, maxiter = 500,
                   bootstrap_type = "s", clust_options = list(cores = 3),
                   seed = 2543, replace = "geo_sketch",
                   volume_ratio = "none", # set to "none" if too slow
                   delta=0, conv_crit = 1e-04, order_type = "align",
                   sample_prop = 0.1)
# Show variance explained by k-vertex model on top of k-1 model (each k separately)
plot_arc_var(arc_ks, type = "res_varexpl", point_size = 2, line_size = 1.5) + theme_bw()

# Show variance in position of vertices obtained using bootstraping
# - use this to find largest k that has low variance
plot_arc_var(arc_ks, type = "total_var", point_size = 2, line_size = 1.5) +
  theme_bw() +
  ylab("Mean variance in position of vertices")

align_arc(arc$XC, archetypes$XC)
align_arc(average_pch_fits(arc_data_rob_m)$XC, archetypes$XC)

arc_data_rob_avg = average_pch_fits(arc_data_rob_m)
weights = solve.qr(qr(arc_data_rob_avg$XC), data)
hist(weights)

#------------------------------------------------------------------------------=
# compare to k-means
library(ParetoTI)
devtools::load_all("../ParetoTI/")
library(ggplot2)
library(cowplot)
set.seed(4355)
# generate data
archetypes = generate_arc(arc_coord = list(c(5, 0, 4), c(-10, 15, 0), c(-30, -20, -5)),
                          mean = 0, sd = 1)
data = generate_data(archetypes$XC, N_examples = 1e3, jiiter = 0.04, size = 0.99)

# find archetypes
arc = fit_pch(data, noc = 4)
# find clusters
clusters = fit_pch(data, noc = 4, method = "kmeans")

plot_grid(plot_arc(arc_data = arc, data = data,
                   which_dimensions = 1:2) + ylim(-18, 17),
          plot_arc(arc_data = clusters, data = data,
                   which_dimensions = 1:2,
                   data_lab = as.character(apply(clusters$S, 2, which.max))) + ylim(-18, 17),
          align = "vh")

# bootstrap with kmeans
clusters_rob = fit_pch_bootstrap(data, n = 200, sample_prop = 0.65, seed = 2543,
                                 noc=4, method = "kmeans")
plot_arc(arc_data = clusters_rob, data = data,
         which_dimensions = 1:2,
         data_lab = as.character(apply(clusters$S, 2, which.max))) + ylim(-18, 17)

# trying different number of clusters
arc_ks = k_fit_pch(data, ks = 2:5,
                   bootstrap = T, bootstrap_N = 200, maxiter = 500,
                   bootstrap_type = "s", clust_options = list(cores = 3),
                   seed = 2543, replace = FALSE,
                   volume_ratio = "none", # set to "none" if too slow
                   order_type = "align", sample_prop = 0.65, reference = T, method = "kmeans")

# Show variance explained by k-vertex model on top of k-1 model (each k separately)
plot_arc_var(arc_ks, type = "res_varexpl", point_size = 2, line_size = 1.5) + theme_bw()

# Show variance in position of vertices obtained using bootstraping
# - use this to find largest k that has low variance
plot_arc_var(arc_ks, type = "total_var", point_size = 2, line_size = 1.5) +
  theme_bw() +
  ylab("Mean variance in position of vertices")
#------------------------------------------------------------------------------=

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

plot_arc(arc_data = arc_data_rob_cmq, data = data,
         which_dimensions = 1:3, line_size = 1.5)
plot_arc(arc_data = arc_data_rob_cmq, data = data,
         which_dimensions = 1:2, line_size = 1) +
  theme_bw()

# test function for different k
arc_ks = k_fit_pch(data, ks = 1:4, check_installed = T, delta=0)
plot_arc(arc_data = arc_ks, data = data,
         which_dimensions = 1:3, type = "all", arch_size = 2,
         colors = c("#D62728", "#1F77B4", "#2CA02C", "#17BED0", "grey"))
plot_arc(arc_data = arc_ks, data = data,
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
plot_arc(arc_data = res, data = data,
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


plot_arc(arc_data = arc_data, data = data,
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
plot_arc(arc_data = arc_data, data = data,
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

library(ParetoTI)
# Random data that fits into the triangle (2D)
set.seed(4355)
archetypes3 = generate_arc(arc_coord = list(c(5, 0, 17, 4), c(-10, 15, -10, 8), c(-20, -20, 5, 1)),
                          mean = 0, sd = 1)
data3 = generate_data(archetypes3$XC, N_examples = 1e4, jiiter = 0.04, size = 0.99)
p1 = plot_arc(arc_data = archetypes3, data = data3,
         which_dimensions = 1:2) +
  ggplot2::theme_bw()

# Random data that fits into the triangle (2D)
archetypes2 = generate_arc(arc_coord = list(c(5, -2, 3, -1)*3, c(-2, 1.5, 2, 3)*3),
                          mean = 0, sd = 1)
data2 = generate_data(archetypes2$XC, N_examples = 1e4, jiiter = 0.04, size = 0.99)
p2 = plot_arc(arc_data = archetypes2, data = data2,
         which_dimensions = 1:2) +
  ggplot2::theme_bw()

data_mix = rbind(data3, data2)
arc_data = k_fit_pch(data_mix, k = 1:8, delta = 0, volume_ratio = "none")

u_data = arch_to_umap(arc_data, data_mix, method = "umap-learn", metric = "euclidean")
p3 = plot_arc(arc_data = u_data$arc_data, data = u_data$data,
         which_dimensions = 1:2, line_size = 0) +
  ggplot2::theme_bw()
p3

cowplot::plot_grid(p2, p1)
