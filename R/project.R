##' Project archetypes and data in PCA dimentions
##' @name project_to_pcs
##' @rdname project_to_pcs
##' @description project_to_pcs() projects archetypes (\code{arc_data}) and data points (\code{arc_data}) to PC space. Archetypes are projected into PC space of data, e.i. archetypes do not affect PCA and are just projected afterwards.
##' @param arc_data objects of class "pch_fit", "b_pch_fit" storing the position of archetypes and other data produced by \code{\link[ParetoTI]{fit_pch}}(). arc_data$XC is matrix of dim(dimensions, archetypes) or list where each element is XC matrix from an independent run of the archetypal analysis.
##' @param data matrix of data used in archetypal analysis, dim(variables/dimentions, examples)
##' @param n_dim number of principal component dimensions
##' @param s list 's' containing SVD decomposition results (U, d, Vt), standard deviation and mean of genes used for decomposition (sd, means)
##' @param pc_method method to use for finding PCs: \code{\link[base]{svd}} or \code{\link[irlba]{irlba}}
##' @param log2 log2-transform before to z-scoring and PC-projection
##' @param offset log2 transformation offset (e.g. \code{log2(x + offset)})
##' @param zscore standardise (substract the mean and divide by standard deviation) prior to PC-projection
##' @return project_to_pcs(): list with projected $data, archetypes ($arc_data) and $s list of decomposition matrices, sds and means
##' @export project_to_pcs
##' @examples
##' # Random data that fits into the triangle
##' set.seed(4355)
##' arc_data = generate_arc(arc_coord = list(c(5, 0), c(-10, 15), c(-30, -20)),
##'                           mean = 0, sd = 1)
##' data = generate_data(arc_data$XC, N_examples = 1e4, jiiter = 0.04, size = 0.9)
##' # Plot
##' plot_arc(arch_data = arc_data, data = data,
##'          which_dimensions = 1:2, data_alpha = 0.5) +
##'          ggplot2::theme_bw()
##'
##' # Project to PCs (in this case just rotate to align x-axis with
##' #                  the axis of most variation because the data is already 2D)
##' pcs = project_to_pcs(arc_data, data, n_dim = 2, pc_method = c("svd", "irlba")[1])
##' # Plot in PC coordinates
##' plot_arc(arch_data = pcs$arc_data, data = pcs$data,
##'          which_dimensions = 1:2, data_alpha = 0.5) +
##'          ggplot2::theme_bw()
##'
##' # Project from PCs back to expression
##' projected = project_from_pc(pcs$arc_data, pcs$s,
##'                             undo_zscore = FALSE, undo_log2 = FALSE)
##'
##' # Plot plot in projected coordinates
##' plot_arc(arch_data = projected, data = data,
##'          which_dimensions = 1:2, data_alpha = 0.5) +
##'          ggplot2::theme_bw()
project_to_pcs = function(arc_data = NULL, data, n_dim = nrow(data), s = NULL,
                          pc_method = c("svd", "irlba"),
                          log2 = FALSE, offset = 1, zscore = FALSE){

  pc_method = match.arg(pc_method)

  # check object class and extract archetypes  --------------------------------
  if(!(is(arc_data, "b_pch_fit") | is(arc_data, "pch_fit") |
       is(arc_data, "random_arc") | is.null(arc_data))) {
    stop("arc_data should be of class pch_fit, b_pch_fit or random_arc")
  }

  # log2-transform --------------------------------
  if(log2) {

    # transform data
    data = log2(data + offset)

    # transform archetypes
    if(is(arc_data, "pch_fit") | is(arc_data, "random_arc")) {
      arc_data$XC = log2(arc_data$XC + offset)
    } else if(is(arc_data, "b_pch_fit")) {
      arc_data$pch_fits$XC = lapply(arc_data$pch_fits$XC, function(XC) {
        log2(XC + offset)
      })
    }
  }

  # zscore --------------------------------
  if(zscore) {

    # find the mean and sd
    mean = Matrix::rowMeans(data)
    if(is(data, "sparseMatrix")){
      sd = DelayedMatrixStats::rowSds(DelayedArray::DelayedArray(data))
    } else {
      sd = matrixStats::rowSds(data)
    }
    names(mean) = rownames(data)
    names(sd) = rownames(data)


    # transform data
    data = (data - matrix(mean, nrow = nrow(data), ncol = ncol(data), byrow = F)) /
      matrix(sd, nrow = nrow(data), ncol = ncol(data), byrow = F)

    # transform archetypes
    if(is(arc_data, "pch_fit") | is(arc_data, "random_arc")) {

      arc_data$XC = (arc_data$XC - matrix(mean, nrow = nrow(data), ncol = ncol(arc_data$XC), byrow = F)) /
        matrix(sd, nrow = nrow(data), ncol = ncol(arc_data$XC), byrow = F)

    } else if(is(arc_data, "b_pch_fit")) {
      arc_data$pch_fits$XC = lapply(arc_data$pch_fits$XC, function(XC) {

        (XC - matrix(mean, nrow = nrow(data), ncol = ncol(XC), byrow = F)) /
          matrix(sd, nrow = nrow(data), ncol = ncol(XC), byrow = F)

      })
    }
  } else {
    mean = NULL
    sd = NULL
  }

  # find PCs by svd-decomposition  --------------------------------
  if(is.null(s)){
    if(pc_method == "irlba") {
      s = irlba::irlba(data, nv = n_dim)
    } else if(pc_method == "svd"){
      s = svd(data, nu = n_dim, nv = n_dim)
    }
  }
  # name dimensions
  rownames(s$u) = rownames(data)
  colnames(s$u) = paste0("PC", seq_len(n_dim))
  # add mean and sd to s
  s$means = mean
  s$sds = sd

  # project data  --------------------------------
  colnames = colnames(data) # save data point names
  data = crossprod(s$u, data) # equivalent to t(s$u) %*% data
  colnames(data) = colnames # reassign data point names
  rownames(data) = colnames(s$u)


  # project archetypes  --------------------------------
  if(is(arc_data, "pch_fit") | is(arc_data, "random_arc")) { # single fit

    colnames = colnames(arc_data$XC) # save names
    arc_data$XC = crossprod(s$u, arc_data$XC) # project
    rownames(arc_data$XC) = colnames(s$u) # reassign names
    colnames(arc_data$XC) = colnames

  } else if(is(arc_data, "b_pch_fit")) { # multiple fits

    arc_data$pch_fits$XC = lapply(arc_data$pch_fits$XC, function(XC) {
      colnames = colnames(XC) # save names
      XC = crossprod(s$u, XC) # project
      rownames(XC) = colnames(s$u) # reassign names
      colnames(XC) = colnames
      XC
    })

  }

  # return projected data, archetypes and s list of decomposition matrices  ----
  list(data = data, arc_data = arc_data, s = s)
}


##' @rdname project_to_pcs
##' @name project_from_pc
##' @description project_from_pc() projects archetypes and data points to original space provided SVD decomposition results. Optionally do the reverse of log2 transformation to obtain normalised expression space.
##' @param undo_zscore undo z-scoring by multiplying by standard deviation and adding the mean? Undo z-scoring precedes exponentiation.
##' @param undo_log2 undo log2-transformation by exponentiating and substracting pseudocount?
##' @return project_from_pc(): archetypes projected to data space
##' @export project_from_pc
project_from_pc = function(arc_data, s, undo_zscore = FALSE,
                           undo_log2 = FALSE, offset = 1){

  # check object class ---------------------------------------------
  if(!(is(arc_data, "b_pch_fit") | is(arc_data, "pch_fit") |
       is(arc_data, "random_arc"))) {
    stop("arc_data should be of class pch_fit or b_pch_fit")
  }

  # check elements of s ---------------------------------------------
  if(mean(c("u", "d", "v") %in% names(s)) != 1) {
    stop("list 's' should contain SVD decomposition results (u, d, v)")
  } else if(undo_zscore & mean(c("sds", "means") %in% names(s)) != 1) {
    stop("list 's' should contain SVD decomposition results (u, d, v), AND
         standard deviation and mean of genes used for z-scoring (sds, means)")
  }

  # project to original space ---------------------------------------------
  if(is(arc_data, "pch_fit") | is(arc_data, "random_arc")) { # single fit

    colnames = colnames(arc_data$XC) # save names
    arc_data$XC = s$u %*% arc_data$XC # project
    rownames(arc_data$XC) = rownames(s$u) # reassign dimension names
    colnames(arc_data$XC) = colnames

  } else if(is(arc_data, "b_pch_fit")) { # multiple fits

    arc_data$pch_fits$XC = lapply(arc_data$pch_fits$XC, function(XC) {
      colnames = colnames(XC) # save names
      XC = s$u %*% XC # project
      rownames(XC) = rownames(s$u) # reassign names
      colnames(XC) = colnames
      XC
    })

  }


  # rescale log space by the mean and sd ---------------------------------------------
  if(undo_zscore) {
    if(is(arc_data, "pch_fit") | is(arc_data, "random_arc")) { # single fit

      if(!is.null(rownames(arc_data$XC))) { # use dimension names to filter sds and means
        dim_names = rownames(arc_data$XC)
      } else {
        dim_names = seq_len(nrow(arc_data$XC))
      }
      # rescale
      arc_data$XC = arc_data$XC *
        matrix(s$sds[dim_names], nrow = nrow(arc_data$XC),
               ncol = ncol(arc_data$XC), byrow = F) +
        matrix(s$means[dim_names], nrow = nrow(arc_data$XC),
               ncol = ncol(arc_data$XC), byrow = F)

    } else if(is(arc_data, "b_pch_fit")) { # multiple fits

      arc_data$pch_fits$XC = lapply(arc_data$pch_fits$XC, function(XC) {

        if(!is.null(rownames(XC))) { # use dimension names to filter sds and means
          dim_names = rownames(XC)
        } else {
          dim_names = seq_len(nrow(XC))
        }
        # rescale
        XC * matrix(s$sds[dim_names], nrow = nrow(XC),
                    ncol = ncol(XC), byrow = F) +
          matrix(s$means[dim_names], nrow = nrow(XC),
                 ncol = ncol(XC), byrow = F)

      })
    }
  }

  # exponentiate and substract pseudocounts ---------------------------------------------
  if(undo_log2) {

    if(is(arc_data, "pch_fit") | is(arc_data, "random_arc")) { # single fit

      # exponentiate
      arc_data$XC = 2^arc_data$XC - 1

    } else if(is(arc_data, "b_pch_fit")) { # multiple fits

      arc_data$pch_fits$XC = lapply(arc_data$pch_fits$XC, function(XC) {
        # exponentiate
        2^XC - 1
      })
    }
  }

  arc_data
}
