## ParetoTI R package 
### Pareto Task Inference In R (based on ParTI)

This package allows you to perform Pareto Task Inference in R as described in Yuval Hart & Uri Alon paper in Nature Methods (2015):
    [Inferring biological tasks using Pareto analysis of high-dimensional data.](https://www.nature.com/articles/nmeth.3254)
    
This package realises polytope dimentionality and polytope position fitting using a python 2.7 implementation of PCHA algorithm. This algorithm was originally implemented in [Matlab by Morten Mørup](http://www.mortenmorup.dk/MMhomepageUpdated_files/Page327.htm) and re-implemented in python by Ulf Aslak Jensen.    

Statistical significance of this fit is determined using permutations of the dataset that disrupt relationships between variables but keep the distribution of each variable constant.

### Using the package

**Installation**  

You need to install development version of this package from github.com and install a python dependency (PCHA algorithm implementation). If you want to take advantage of high-performance computing cluster you need to do extra set up, please follow the instructions here: https://github.com/mschubert/clustermq.

```r
# Install ParetoTI package, this should also install reticulate package, if not - install manually.
devtools::install_github("vitkl/ParetoTI", dependencies = T)
# Load package
library(ParetoTI)
# Install python dependency into conda python environment and install py_pcha module
ParetoTI::install_py_pcha(method = "conda")
# If no conda manager installed on your machine, try this (uncomment):
## ParetoTI::install_py_pcha(method = "virtualenv")
# If this fails, install python 2.7 Anaconda distribution, then use 'method = "conda"'.
# Finally, check that py_pcha library is successfully installed and discoverable
reticulate::py_discover_config("py_pcha")
```

**Example: Fitting polytope to random data (finding Pareto front)**  

```r
library(ParetoTI)
library(ggplot2)

# Generate random data that fits into the triangle (3D)
set.seed(4355)
archetypes = generate_arc(arc_coord = list(c(5, 0, 4), c(-10, 15, 0), c(-30, -20, -5)),
                          mean = 0, sd = 1)
data = generate_data(archetypes$XC, N_examples = 1e4, jiiter = 0.04, size = 0.9)

# Fit polytope to those data
arc_data = fit_pch(data, noc = as.integer(3), delta = 0)

# Show results as interactive 3D scatterplot using plotly
plot_arc(arch_data = arc_data, data = data,
         which_dimensions = 1:3)
# Plot static 2D scatterplot using ggplot2
plot_arc(arch_data = arc_data, data = data,
         which_dimensions = 1:2) +
  theme_bw()
# Plot data as 2D density rather than scatterplot
plot_arc(arch_data = arc_data, data = data,
    which_dimensions = 1:2, geom = ggplot2::geom_bin2d) +
  theme_bw()
```

### Development and further updates

It is currently under development and enables only polytope fitting and statistical significance tests.
Other analyses, in particular, models of feature enrichment at archetype will be included later.

Alternative methods for polytope fitting will be considered such as SPAMS implemented in c++ used via R interface: https://www.stat.berkeley.edu/~yuansi.chen/demo/demo.html, https://arxiv.org/pdf/1405.6472.pdf, http://spams-devel.gforge.inria.fr/downloads.html
