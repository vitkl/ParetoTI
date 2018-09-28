## ParetoTI R package 
### Pareto Task Inference In R (based on ParTI)

This package allows you to perform Pareto Task Inference in R as described in Yuval Hart & Uri Alon paper in Nature Methods (2015):
    [Inferring biological tasks using Pareto analysis of high-dimensional data.](https://www.nature.com/articles/nmeth.3254)
    
This package realises polytope dimentionality and polytope position fitting using a python 2.7 implementation of PCHA algorithm. This algorithm was originally implemented in [Matlab by Morten Mørup](http://www.mortenmorup.dk/MMhomepageUpdated_files/Page327.htm) and re-implemented in python by Ulf Aslak Jensen.    

Statistical significance of this fit is determined using permutations of the dataset that disrupt relationships between variables but keep the distribution of each variable constant.

It is currently under development and enables only polytope fitting and statistical significance tests.
Other analyses, in particular, models of feature enrichment at archetype will be included later.

Alternative methods for polytope fitting will be considered such as SPAMS implemented in c++ and using R interface: https://www.stat.berkeley.edu/~yuansi.chen/demo/demo.html, https://arxiv.org/pdf/1405.6472.pdf, http://spams-devel.gforge.inria.fr/downloads.html
