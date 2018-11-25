##' Find activity of a gene set using AUCell model (continuos and binary)
##' @rdname find_set_activity_AUCell
##' @name find_set_activity_AUCell
##' @description \code{find_set_activity_AUCell()} finds activity of each gene set in each cell by combining AUCell functions into a pipeline. Aerts lab that developed AUCell recommends adjusting the threshold when binarising activities (Details: https://www.bioconductor.org/packages/devel/bioc/vignettes/AUCell/inst/doc/AUCell.html#determine-the-cells-with-the-given-gene-signatures-or-active-gene-sets). See documentations for details: \link[AUCell]{AUCell_buildRankings}, \link[AUCell]{AUCell_calcAUC}, \link[AUCell]{AUCell_exploreThresholds}.
##' @param expr_mat expression matrix (genes in rows, cells in columns) or one of: dgCMatrix, ExpressionSet, and SummarizedExperiment or SingleCellExperiment both of which require assay_name.
##' @param assay_name name of assay in SummarizedExperiment or SingleCellExperiment, normally counts or logcounts
##' @param aucMaxRank argument for \link[AUCell]{AUCell_calcAUC}. Threshold to calculate the AUC.In a simplified way, the AUC value represents the fraction of genes, within the top X genes in the ranking, that are included in the signature. The parameter 'aucMaxRank' allows to modify the number of genes (maximum ranking) that is used to perform this computation. By default, it is set to 0.05 of the total number of genes in the rankings. Common values may range from 0.01 to 0.3.
##' @param gene_sets data.table or coercible to data.table that contains gene set annotations.
##' @param gene_col column in gene_sets storing gene identifiers.
##' @param set_id_col column in gene_sets storing set identifiers.
##' @param set_name_col column in gene_sets storing readable set names.
##' @param binary binarise gene set activities using \link[AUCell]{AUCell_exploreThresholds}?
##' @param nCores number of cores for parallel processing. See AUCell docs for details.
##' @param plotHist plot the AUC histograms? \link[AUCell]{AUCell_exploreThresholds}.
##' @param ... other arguments passed to \link[AUCell]{AUCell_exploreThresholds}.
##' @return \code{find_set_activity_AUCell()} data.table of gene set activities with cell in rows and gene sets in columns. Column titled "cells" contains cell ids (column names of expr_mat).
##' @export find_set_activity_AUCell
##' @import data.table
find_set_activity_AUCell = function(expr_mat, assay_name = "logcounts",
                                    aucMaxRank = nrow(expr_mat) * 0.05,
                                    gene_sets, gene_col = "ALIAS",
                                    set_id_col = "GOALL", set_name_col = "TERM",
                                    binary = FALSE, nCores = 1,
                                    plotHist = FALSE, plotStats = TRUE,  ...){

  if(!is.data.table(try(as.data.table(gene_sets), silent = TRUE))) stop("gene_sets should be data.table or coercible to data.table: data.frame or matrix")
  gene_sets = as.data.table(gene_sets)

  if(is(expr_mat, "matrix") | is(expr_mat, "dgCMatrix") | is(expr_mat, "ExpressionSet")){
    cells_rankings = AUCell::AUCell_buildRankings(expr_mat, nCores = nCores,
                                                  plotStats = plotStats)
  } else if (is(expr_mat, "SummarizedExperiment")) {
    cells_rankings = AUCell::AUCell_buildRankings(expr_mat, nCores = nCores,
                                                  assayName = assay_name,
                                                  plotStats = plotStats)
  }

  gene_set_list = split(gene_sets[, get(gene_col)],
                        gene_sets[, get(set_id_col)])

  cells_AUC = AUCell::AUCell_calcAUC(gene_set_list, cells_rankings,
                                     aucMaxRank = aucMaxRank,
                                     nCores = nCores)
  if(binary){
    cells_assignment = AUCell::AUCell_exploreThresholds(cells_AUC,
                                                        plotHist = plotHist,
                                                        nCores = nCores, assign = TRUE, ...)
    # extract cell to go assignments
    cells_assignment = lapply(cells_assignment, function(x) x$assignment)
    cells_assignment = unlist2(cells_assignment)
    cells_assignment = data.table(cells = cells_assignment,
                                  sets = names(cells_assignment))
    # Map set names
    ids = unique(cells_assignment$sets)
    if(set_id_col == set_name_col) {
      GO_names = gene_sets[match(ids, get(set_id_col)),
                           .(get(set_name_col))]
      setnames(GO_names, colnames(GO_names), "sets")
    } else {
      GO_names = gene_sets[match(ids, get(set_id_col)),
                           .(get(set_id_col), get(set_name_col))]
      setnames(GO_names, colnames(GO_names), c("sets", set_name_col))
    }

    cells_assignment = merge(cells_assignment, GO_names,
                             by = "sets", all.x = T, all.y = F)
    # convert to cells * sets matrix
    cells_assignment = dcast.data.table(cells_assignment, cells ~ get(set_name_col),
                                        value.var = set_name_col,
                                        fun.aggregate = length, fill = 0)
  } else {
    cells_assignment = as.data.table(t(assay(cells_AUC)), keep.rownames = "cells")
    # Map set names
    ids = colnames(cells_assignment)
    ids = ids[ids != "cells"]
    terms = gene_sets[match(ids, get(set_id_col)), get(set_name_col)]
    setnames(cells_assignment, ids, terms)
  }
  cells_assignment
}

##' @rdname find_set_activity_AUCell
##' @name find_set_activity_pseudoinv
##' @description \code{find_set_activity_pseudoinv()} finds activity of each gene set in each cell by solving this matrix equation: expression = activities x gene_assignment_to_sets => activities = pseudoinverse(gene_assignment_to_sets) x expression. Based on code from Inferelator package by Richard Bonneau lab.
##' @param noself Remove self-interactions from set annotations (when sets are TF targets)
##' @export find_set_activity_pseudoinv
##' @import data.table
find_set_activity_pseudoinv = function(expr_mat, assay_name = "logcounts",
                                       gene_sets, gene_col = "ALIAS",
                                       set_id_col = "GOALL", set_name_col = "TERM",
                                       noself = FALSE) {
  # generate gene set assignment matrix from data.table
  gene_set_mat = dcast.data.table(gene_sets, get(gene_col) ~ get(set_id_col),
                                  value.var = set_id_col,
                                  fun.aggregate = length, fill = 0)
  gene_set_mat = as.matrix(gene_set_mat, rownames = "gene_col")

  # find non-empty gene sets
  non_empty_cols = colSums2(gene_set_mat != 0) > 0

  # remove self-interaction of TFs when gene sets are TFs
  tfs = colnames(gene_set_mat)
  tfs = intersect(tfs, rownames(gene_set_mat))
  if (noself) {
    diag(gene_set_mat[tfs, tfs]) = 0
  }

  # extract gene expression matrix from object provided
  if (is(expr_mat, "matrix")) NULL else if (is(expr_mat, "dgCMatrix")){
    expr_mat = as.matrix(expr_mat)
  } else if (is(expr_mat, "ExpressionSet")){
    expr_mat = as.matrix(Biobase::exprs(expr_mat))
  } else if (is(expr_mat, "SummarizedExperiment")) {
    expr_mat = as.matrix(SummarizedExperiment::assay(expr_mat, assay_name))
  }

  # Adjust the size of gene_set_mat to match dimensions of expr_mat
  gene_set_mat = gene_set_mat[match(rownames(expr_mat), rownames(gene_set_mat)),]
  rownames(gene_set_mat) = rownames(expr_mat)
  gene_set_mat[is.na(gene_set_mat)] = 0

  activities = matrix(0, ncol(gene_set_mat), ncol(expr_mat),
                      dimnames=list(colnames(gene_set_mat), colnames(expr_mat)))

  if (any(non_empty_cols)) {
    activities[non_empty_cols, ] = corpcor::pseudoinverse(gene_set_mat[, non_empty_cols, drop=FALSE]) %*% expr_mat
  }

  # when no targets are available for a TF use expression as activity
  use_exp = intersect(colnames(gene_set_mat)[!non_empty_cols], rownames(expr_mat))
  activities[use_exp, ] = expr_mat[use_exp, ]
  activities = as.data.table(t(activities), keep.rownames = "cells")

  # Map set names
  ids = colnames(activities)
  ids = ids[ids != "cells"]
  terms = gene_sets[match(ids, get(set_id_col)), get(set_name_col)]
  setnames(activities, ids, terms)
  activities
}
