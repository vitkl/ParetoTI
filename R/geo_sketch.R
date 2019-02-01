##' Create downsampled representation of the data with Geometric Sketch
##' @rdname geo_sketch
##' @name geo_sketch
##' @description \code{geo_sketch()} Create a downsampled representation of the data while preserving the structure and rare populations with Geometric Sketch method (https://github.com/brianhie/geosketch, https://www.biorxiv.org/content/10.1101/536730v1). This method splits gene expression space into hypercubes (covering boxes) of equal volume and samples uniformly from those cubes. This reduces the density of the data and cell numbers while preserving global structure. Useful for fitting polytopes to large datasets. First step - find PCs using Facebook PCA method, second step - define hypercubes and sample cells.
##' @param data matrix, sparse matrix or SingleCellExperiment (SummarizedExperiment): dim(dimensions * examples)
##' @param N number of cells to sample
##' @param assay_slot slot in data (SingleCellExperiment, SummarizedExperiment) containing matrix to use for PCA and geometric sampling
##' @param use_PCs logical, use PCs (TRUE) or data directly (FALSE)?
##' @param PCs number of PCs to use. Identified using Facebook implementation of PCA (fbpca).
##' @param k Number of covering boxes. When `'auto'` and replace is `True`, draws sqrt(X.shape[0]) covering boxes. When `'auto'` and replace is `False`, draws N covering boxes.
##' @param seed Random number generation seed passed to numpy
##' @param alpha Binary search halts when it obtains between `k * (1 - alpha)` and `k * (1 + alpha)` covering boxes.
##' @param max_iter Maximum iterations at which to terminate binary seach in rare case of non-monotonicity of covering boxes with box side length.
##' @param verbose report progress
##' @return \code{geo_sketch()} integer vector of cell indices
##' @import reticulate
##' @export geo_sketch
geo_sketch = function(data, N, assay_slot = "logcounts", use_PCs = TRUE, PCs = 100, k = "auto", seed = 5232, replace = FALSE, alpha = 0.1, max_iter = 200, verbose = 0) {

  # check if python modules are installed
  .py_fbpca_installed()

  if(is(data, "SingleCellExperiment") || is(data, "SummarizedExperiment")) {
    # extract assay and convert to dgCMatrix
    data = as(SummarizedExperiment::assay(data, assay_slot), "dgCMatrix")
  } else {
    # convert any other object to dgCMatrix
    data = as(data, "dgCMatrix")
  }

  if(use_PCs){
    # find PCs
    s = fbpca$pca(Matrix::t(data), k = as.integer(PCs))
    names(s) = c("U", "s", "Vt")

    # construct reduced dimension space
    X_dimred = Matrix::tcrossprod(s$U[, seq_len(PCs)], diag(s$s[seq_len(PCs)]))
  } else {
    X_dimred = Matrix::t(data)
  }

  # segment space into hypercubes and find representative sample
  unlist(geosketch$gs(X_dimred, N = as.integer(N),
                      k=k, seed=as.integer(seed), replace=replace,
                      alpha=alpha, max_iter=as.integer(max_iter), verbose=verbose))

}

.py_fbpca_installed = function() {
  # check if python package is availlable and give a helpful error message
  err = tryCatch(fbpca$pca(matrix(1:60, 10, 6), k =2), error = function(e) e)
  if(!is.null(err$message)) {
    if(grepl("Python", err$message)){
      stop(paste0(err$message,
                  ", fbpca and geosketch missing, please use install_py_pcha() to install"))
    }
  }
}
