##' Annotate archetypes with custom labels
##' @rdname annotate_archetypes
##' @name annotate_archetypes
##' @author Vitalii Kleshchevnikov
##' @description \code{annotate_archetypes()} adds custom labels to archetype indices: \code{paste0(label, "_", index)}
##' @param arc_data objects of class "pch_fit", "b_pch_fit", "k_pch_fit" storing the position of archetypes, and other data from \code{\link[ParetoTI]{fit_pch}}() run. arc_data$XC is matrix of dim(dimensions, archetypes) or list where each element is XC matrix from an independent run of the archetypal analysis.
##' @param ... named arguments specifying label for indicex, e.g. \code{Epithelial = c(1, 5, 9)}
##' @return \code{annotate_archetypes()} archetypal analysis result object, the same as input but with named archetypes
##' @export annotate_archetypes
##' @seealso \code{\link[ParetoTI]{fit_pch}}, \code{\link[ParetoTI]{arch_dist}}
##' @examples
##' # Random data that fits into the triangle
##' set.seed(4355)
##' archetypes = generate_arc(arc_coord = list(c(5, 0), c(-10, 15), c(-30, -20)),
##'                           mean = 0, sd = 1)
##' data = generate_data(archetypes$XC, N_examples = 1e4, jiiter = 0.04, size = 0.9)
##'
##' # Find 3 archetypes in this data
##' arc = fit_pch(data, noc=as.integer(3), delta=0)
##'
##' # Annotate
##' arc = annotate_archetypes(arc, top = 2, bottom = c(1, 3))
##' # Plot
##' plot_arc(arc, data, which_dimensions = 1:2, arc_names_num = F)
annotate_archetypes = function(arc, ...) {

  # check object class
  if(is(arc, "b_pch_fit")) {
    n_arch = ncol(arc$pch_fits$XC[[1]])
  } else if(is(arc, "pch_fit")) {
    n_arch = ncol(arc$XC)
  } else stop("arc should be pch_fit or b_pch_fit")

  vert = paste0("", seq(1, n_arch))

  annots = list(...)

  for (i in seq_len(length(annots))) {
    ind = vert %in% paste0("", annots[[i]])
    vert[ind] = paste0(names(annots)[i], "_", vert[ind])
  }

  # name archetypes and return pch_fit back
  if(is(arc, "b_pch_fit")) {
    arc$pch_fits$XC = lapply(arc$pch_fits$XC, function(XC){
      colnames(XC) = vert
      XC
    })
  } else if(is(arc, "pch_fit")) {
    colnames(arc$XC) = vert
  }

  arc
}
