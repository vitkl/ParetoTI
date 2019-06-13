##' Find and plot similarity structure between archetypes
##' @rdname plot_arch_struct
##' @name plot_arch_struct
##' @author Vitalii Kleshchevnikov
##' @description Find similarity structure between archetypes by comparing archetype weights for cells, measuring similarity of archetypes in expression space and in marker gene. genesXarch option shows
##' @param arc object containing archetypal analysis results (class "pch_fit")
##' @param type one of c("cells", "space", "marker_genes", "genesXarch"). "cells": distances are computed over archetype weights of cells (S matrix). "space": measures distances between archetypes in space, PCs or gene expression, depending on which space was used to find archetypes (XC matrix). "marker_genes": distance is computed base on lists marker genes. "genesXarch": show both archetypes and dimensions, both clustered with hierarchical clustering. When \code{marker_genes} is provided "genesXarch" shows marker gene membership across archetypes rather than expression values of those genes (PCs) at archetypes.
##' @param dist_fun function use to compute distance, should take one argument and compute distances between rows.
##' @param marker_genes filtered data.table containing markers for each archetype, normally the output of \code{\link[ParetoTI]{get_top_decreasing}} stored in $enriched_genes. When \code{type = "genesXarch"} dimensions are filtered using gene names supplied in \code{marker_genes}.
##' @param marker_genes_mean show the strength of marker gene association with archetypes? (mean_diff in the output of \code{\link[ParetoTI]{get_top_decreasing}} stored in $enriched_genes). By default is FALSE - show only marker gene membership (0/1).
##' @param marker_genes_mean_col which column in \code{marker_genes} stores the strength of marker gene association with archetypes?
##' @return Matrix of the same dimention as the original matrix but with values in each column permuted.
##' @export plot_arch_struct
plot_arch_struct = function(arc, type = c("cells", "space", "marker_genes", "genesXarch")[1],
                            dist_fun = function(x) dist(x, method = "euclidean"),
                            marker_genes = NULL, marker_genes_mean = F,
                            marker_genes_mean_col = "mean_diff") {

  if(!is(arc, "pch_fit")) stop("arc should be of class pch_fit (single archetype fit)")

  # find distances
  if(type == "cells") {

    dist_res = arc$S
    if(is.null(colnames(arc$XC))) {
      colnames(dist_res) = paste0("archetype_", seq(1, ncol(dist_res)))
    } else rownames(dist_res) = colnames(arc$XC)
    dist_res = dist_fun(dist_res)

  } else if(type == "space" | type == "genesXarch") {

    dist_res = arc$XC
    if(is.null(colnames(arc$XC))) colnames(dist_res) = paste0("archetype_", seq(1, ncol(dist_res)))
    dist_res = dist_fun(t(dist_res))

  } else if(type == "marker_genes") {

    dist_res = marker_genes
    dist_res = dcast.data.table(dist_res, arch_name ~ genes, fill = 0,
                                value.var = "p", fun.aggregate = length)
    dist_res = as.matrix(dist_res, rownames = "arch_name")
    dist_res = dist_fun(dist_res)

  }

  # convert to matrix
  dist_dt = as.matrix(dist_res)

  # do hierarchical clustering to order archetypes
  dist_hc = hclust(dist_res)
  names_ord = colnames(dist_dt)[dist_hc$order]
  dist_dt = dist_dt[names_ord, names_ord]

  if(type == "genesXarch"){

    dist_dt = arc$XC # overwrite using the original archetype matrix
    if(is.null(colnames(arc$XC))) colnames(dist_dt) = paste0("archetype_", seq(1, ncol(dist_dt)))

    if(!is.null(marker_genes)) {
      if(marker_genes_mean){ # use matrix of marker gene effect-size
        dist_dt = dcast.data.table(marker_genes, genes ~ arch_name,
                                   value.var = marker_genes_mean_col, fill = 0)

      } else { # use binary matrix of marker genes instead
        dist_dt = dcast.data.table(marker_genes, genes ~ arch_name,
                                   value.var = marker_genes_mean_col, fill = 0,
                                   fun.aggregate = length)
      }

      dist_dt = as.matrix(dist_dt, rownames = "genes")
    }

    fact_dist = dist_fun(dist_dt)
    fact_hc = hclust(fact_dist)
    fact_names_ord = rownames(dist_dt)[fact_hc$order]
    dist_dt = dist_dt[fact_names_ord, names_ord]

  }

  # convert to data.table for plotting
  dist_dt = as.data.table(dist_dt, keep.rownames = "arch_col_1")
  dist_dt = melt.data.table(dist_dt, id.vars = "arch_col_1",
                            value.name = "dist", variable.name = "arch_col_2")

  if(type == "genesXarch"){ # use dimension names when looking at genesXarch plot
    dist_dt[, arch_col_1 := factor(arch_col_1, levels = fact_names_ord)]
    lab = ""
  } else {
    dist_dt[, arch_col_1 := factor(arch_col_1, levels = names_ord)]
    lab = "dist"
  }
  dist_dt[, arch_col_2 := factor(arch_col_2, levels = names_ord)]
  dist_dt[as.character(arch_col_1) == as.character(arch_col_2), dist := NA]

  # plot
  ggplot(dist_dt, aes(arch_col_2, arch_col_1, color = dist, fill = dist)) +
    geom_tile() +
    scale_color_viridis_c(direction = -1) + scale_fill_viridis_c(direction = -1) +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = -45, hjust = 0)) +
    labs(fill = lab, color = lab)
}
