## ParetoTI R package 
### Pareto Task Inference In R (based on ParTI)

  This package allows to find tasks that cells need to perform and trade-offs between them. 
  
  Caution: the package in currently in development and testing. If you want to use it please contact Vitalii Kleshchevnikov (vk7 at sanger.ac.uk)
  
  Need to perform multiple tasks and natural selection put cells on a Pareto front, a narrow subspace where performance at those tasks is optimal. How important the tasks are in the environment puts cells at different locations along Pareto front. This reflects trade-off in performance at those tasks. Pareto front in the performance space translates into simple shapes gene expression of cell population. By finding minimal simplex polytope (triangle in 2D, tetrahedron in 3D, 5-vertex polytope in 4D) that encloses most of the data you can describe within cell-type heterogeniety. Cells near each vertex are specialists at one tasks, cells withing the shape perform a weighted combination of tasks. You can indentify the cellular tasks by finding what is special about cells closest to each vertex. This relies on recent work by Uri Alon group that showed that Pareto front is equal to minimal polytope defined by specialist phenotypes (convex hull defined by archetypes) and developed a matlab package ParTI for performing this analysis.
  
  See [Manual](https://vitkl.github.io/ParetoTI/articles/introduction.html).
  
  ParTI matlab package is described in more detail in Yuval Hart & Uri Alon paper in Nature Methods (2015):
    [Inferring biological tasks using Pareto analysis of high-dimensional data.](https://www.nature.com/articles/nmeth.3254)
    
  This ParetoTI R package realises very similar procedure with a few differences. 
  
  Polytope dimentionality and polytope positions are found using a python 2.7 implementation of PCHA algorithm. This algorithm was originally implemented in [Matlab by Morten Mørup](http://www.mortenmorup.dk/MMhomepageUpdated_files/Page327.htm) and re-implemented in python by Ulf Aslak Jensen.    

  As in original package, statistical significance of this fit is determined using permutations of the dataset that disrupt relationships between variables but keep the distribution of each variable constant. We implement bootstraping method to measure variability in vertex position. This adds additional selection criteria for best fit polytopes. Excessive number of vertices will lead to higher variance in positions. 

  Features, genes whose expression decreases with distance from each vertex are identified using Generalised Additive Models (cubic splines) that define a smooth function of gene expression as a function of distance from vertices. We define p-value as the probability that the first derivative of this function is below zero (function is decreasing).

### Using the package

**Installation**  

You need to install development version of this package from github.com and install a python dependency (PCHA algorithm implementation). If you want to take advantage of high-performance computing cluster you need to do extra set up, please follow the instructions here: https://github.com/mschubert/clustermq.

```r
# Install ParetoTI package, this should also install reticulate package, if not - install manually.
install.packages("BiocManager") # for installing BioConductor dependencies
BiocManager::install("vitkl/ParetoTI", dependencies = T)

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

**Example: Fitting polytope to random triangle-shaped data (finding Pareto front)**  

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

### Development and further improvements

It is currently under development and enables polytope fitting, statistical significance test by permutation, evaluating variance in vertex position by bootstraping, feature enrichment at archetype using the first derivative of Generalised Additive Model, measuring gene set activities in each cell with subsequent enrichment at archetype.

Key improvement that can be made is an implementation of PCHA algorhitm in c++ and for sparse matrices. Alternative methods for polytope fitting will be considered such as SPAMS implemented in c++ used via R interface: https://www.stat.berkeley.edu/~yuansi.chen/demo/demo.html, https://arxiv.org/pdf/1405.6472.pdf, http://spams-devel.gforge.inria.fr/downloads.html
