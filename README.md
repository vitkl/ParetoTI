## ParetoTI R package 
### R toolbox for Archetypal Analysis and Pareto Task Inference for unsupervised characterisation of single cell data

  This package allows describing a continuum of cell states by identifying representative cells with most distinctive expression (archetypes). Detecting which genes and functional annotations (e.g. GO) make those cells distinctive gives marker genes and predicted functions for each state.   
  
  Archetypal Analysis avoids drawing arbitrary cluster boundaries when a continuum of cell states is apparent (e.g. within-cell-type heterogeniety, multi-step differentiation trajectories, [this artificial example](https://vitkl.github.io/ParetoTI/articles/Comparison_to_kmeans.html#find-archetypes-with-pcha-and-cluster-centers-with-k-means)). Yet the method provides a similar and easy to understand conceptual framework for unsupervised discovery of cell subpopulations.    
  
  In theory (see [citation](https://www.nature.com/articles/nmeth.3254)), archetypes are cells specialised at a combination of tasks, and the continuum of transcriptional states between archetypes mirrors the trade-offs between tasks. In practice, many factors lead could to a distinct expression profile (archetype), including technical bias, cell doublets, stress response induced by cell dissociation and other factors. So it is necessary to keep [the small-world large-world distinction](https://youtu.be/4WVelCswXo4?t=2256) in mind: specific statistical/computational method will find the best archetypes (specialised cells) it could find conditional on that method's assumptions about the data - rather than the true archetypes you are trying to discover. So you have to supervise and rule out spurious causes the same way you would do for annotating cell clusters.    
  
  Caution: the package in currently in development and testing. If you encounter problems please contact Vitalii Kleshchevnikov (vk7 at sanger.ac.uk) or report bugs on [github](https://github.com/vitkl/ParetoTI/issues).
  
  [![DOI](https://zenodo.org/badge/146516634.svg)](https://zenodo.org/badge/latestdoi/146516634)
  
#### Background
  
  Need to perform multiple tasks and natural selection put cells on a Pareto front, a narrow subspace where performance at those tasks is optimal. How important the tasks are in the environment puts cells at different locations along Pareto front. This reflects trade-off in performance at those tasks. Pareto front in the performance space translates into simple shapes gene expression of cell population. By finding minimal simplex polytope (triangle in 2D, tetrahedron in 3D, 5-vertex polytope in 4D) that encloses most of the data you can describe within cell-type heterogeniety. This is done with archetypal analysis method (matrix factorisation) that finds the most distictive representative cells at "the corners of the data". This makes it similar to clustering methods and allows interpreting archetypes as cells (rather than dimensions). The theory suggests that cells near each archetype/vertex are specialists at one tasks, whicle cells between archetypes perform a weighted combination of tasks. You can indentify the cellular tasks by finding what is special about cells closest to each archetype/vertex. This package is inspired and partially relies on recent work by Uri Alon's group that showed that Pareto front is equal to minimal polytope defined by specialist phenotypes (archetypes) and their matlab package ParTI for performing this analysis.
  
  ParTI matlab package is described in more detail in Yuval Hart & Uri Alon paper in Nature Methods (2015):
  [Inferring biological tasks using Pareto analysis of high-dimensional data.](https://www.nature.com/articles/nmeth.3254)

#### Vignette

  Please follow this [example](https://vitkl.github.io/ParetoTI/articles/Hepatocyte_example.html) reproducing archetypes of liver cells from:
  [Continuum of Gene-Expression Profiles Provides Spatial Division of Labor within a Differentiated Cell Type. Miri Adler et al.](https://www.sciencedirect.com/science/article/pii/S2405471218304824)

### Short guide: installation and basic functionality

#### Installation 

You need to install development version of this package from github.com and install a python dependencies (PCHA algorithm implementation, and (optionally) keras, geosketch, fbpca, umap-learn). If you want to take advantage of high-performance computing cluster you need to do extra set up, please follow the instructions here: https://github.com/mschubert/clustermq.

```r
# Install ParetoTI package, this should also install reticulate package, if not - install manually.
install.packages("BiocManager") # for installing BioConductor dependencies
BiocManager::install("vitkl/ParetoTI", dependencies = c("Depends", "Imports", "LinkingTo", "Suggests"))
```
```r
# Load package
library(ParetoTI)
# If package does not load because "py_pcha is not found" make sure you do not have
# a python environment already loaded in this R session (e.g. restart R and try loading again).

# Install python dependencies (like py_pcha) into python conda environment,
# and (optionally) install *extra_packages*.
ParetoTI::install_py_pcha(method = "conda", 
                          extra_packages = c("tensorflow", "pandas", "keras", "h5py",
                                        "geosketch", "pydot", "sklearn", "umap-learn"))
```
```r
# If no conda manager installed on your machine, install python 2.7 miniconda distribution:
# https://conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation
# or try virtualenv (uncomment):
## ParetoTI::install_py_pcha(method = "virtualenv")
# You can also install these python modules directly in terminal,
# you just ensure that the conda environment is named "reticulate_PCHA" (uncomment):
## conda create -n reticulate_PCHA python=2.7.13 pip
## source activate reticulate_PCHA && pip install --upgrade py_pcha numpy scipy datetime tensorflow pandas keras h5py geosketch pydot sklearn umap-learn
```
```r
# Finally, check that py_pcha library is successfully installed and discoverable
reticulate::py_discover_config("py_pcha")
```

**Example: Finding archetypes in data created by simulating a triangle  **  

```r
library(ParetoTI)
library(ggplot2)

# Generate random data that fits into the triangle (3D)
set.seed(4355)
archetypes = generate_arc(arc_coord = list(c(5, 0, 4), c(-10, 15, 0), c(-30, -20, -5)),
                          mean = 0, sd = 1)
data = generate_data(archetypes$XC, N_examples = 1e4, jiiter = 0.04, size = 0.9)
```
```r
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
```r
# Project to UMAP coordinates (3D -> 2D)
arc_umap = arch_to_umap(arc_data, data, which_dimensions = 1:2,
                        method = c("naive", # implemented in R and slow
                                   "umap-learn")) # requires python module
plot_arc(arch_data = arc_umap$arch_data, data = arc_umap$data,
    which_dimensions = 1:2) +
    theme_bw()
```
```r
# Project to tSNE coordinates (3D -> 2D, requires Rtsne package)
arc_tsne = arch_to_tsne(arc_data, data, which_dimensions = 1:2)
plot_arc(arch_data = arc_tsne$arch_data, data = arc_tsne$data,
    which_dimensions = 1:2) +
    theme_bw()
```

### Details and comparison to ParTI in Matlab
    
  This ParetoTI R package realises very similar procedure with a few differences. 
  
  Archetypal analysis is performed using PCHA algorithm. For that our R package interfaces with a python 2.7 implementation of PCHA. This algorithm was originally implemented in [Matlab by Morten Mørup](http://www.mortenmorup.dk/MMhomepageUpdated_files/Page327.htm) and re-implemented in python by Ulf Aslak Jensen.    

  As in the original package, statistical significance of how well polytope defined by archetypes describes the data is determined using permutations of the dataset that disrupt relationships between variables but keep the distribution of each variable constant. We implement bootstraping method to measure variability in vertex position. This adds additional selection criteria for best fit polytopes. Excessive number of vertices will lead to higher variance in positions. 
  
  Archetypes can be visualised in 2D and 3D.

  Features, genes whose expression decreases with distance from each vertex are identified using Wilcox test as those highly expressed near archetypes. As more direct but slower approach to answering this question we also implemented Generalised Additive Model (cubic splines) that describe gene expression as a a smooth function of distance from archetypes. In that case we define p-value as the probability that the first derivative of this function is below zero (function is decreasing). This approach can be used for plotting the expression of selected genes.
  
  We also strongly recommend to remove doublet cells using scrublet (https://github.com/AllonKleinLab/scrublet, or other doublet prediction methods) because doublets result in spurious archetypes and make selecting the number of archetypes harder. You can vary the scrublet score threshold until ParetoTI no longer finds arhetyopes with high scrublet score. Also, for convenience, we provide an interface to logistic regression model (keras) for classifying cells with logistic regression, and Geometric Sketch method (geosketch) for reducing the size of the data while preseving rare populations. Both can be done to aid Pareto Task Inference analysis. 

### Development and further improvements

It is currently under development and enables polytope fitting, statistical significance test by permutation, evaluating variance in vertex position by bootstraping, feature enrichment at archetype using the Wilcox test and the first derivative of Generalised Additive Model, measuring gene set activities in each cell with subsequent enrichment at archetypes.

Alternative methods for archetypal analysis are being considered: AANet (https://github.com/KrishnaswamyLab/AAnet/) and SPAMS implemented in c++ via R interface (https://arxiv.org/pdf/1405.6472.pdf, http://spams-devel.gforge.inria.fr/downloads.html)
