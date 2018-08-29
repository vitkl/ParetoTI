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
#  # Use local parallel processing to fit the 10 polytopes to subsampled datasets each time looking at 70% of examples.
arc_data = fit_pch_robust(data, n = 100, subsample = 0.7, seed = 2543,
                          noc=as.integer(3), delta=0.1, type = "cmq")
#})

devtools::load_all()

# set directory for user libraries, update pip, setuptools, wheel in that environment
#export PYTHONUSERBASE=$vk7/software/python_libs/
#python -m pip install --user -i https://pypi.python.org/simple -U pip distribute
#python -m pip install --user -i https://pypi.python.org/simple --upgrade pip setuptools wheel virtualenv -U pip --user distribute
