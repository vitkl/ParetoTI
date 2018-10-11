##'
##'
##' Project in low dimentions (PCA) and process data and polytopes fits for plotting
##' @param arc_data list of matrices storing the position of archetypes, class "pch_fit", "r_pch_fit". Each element of a list represents an independent run of the polytope fitting algorithm
##' @param data matrix of data in which archetypes/polytope were found, dim(variables/dimentions, examples)
##' @examples
##' # Random data that fits into the triangle
##' set.seed(4355)
##' arc_data = generate_arc(arc_coord = list(c(5, 0), c(-10, 15), c(-30, -20)),
##'                           mean = 0, sd = 1, N_dim = 2)
##' data = generate_data(archetypes, N_examples = 1e4, jiiter = 0.04, size = 0.9)
project = function(arc_data, data){

}
process_for_plotly = function(arc_data, data){

}

.arc_data_table = function(arc_data, data, type = c("average", "all")[1], average_func = mean){
  if(is(arc_data, "pch_fit") | is(arc_data, "random_arc")){
    arc_data = as.data.table(t(arc_data$XC))
    arc_data$lab = "archetypes"
    data = as.data.table(t(data))
    data$lab = "data"
  }
  if(is(arc_data, "r_pch_fit")){
    arc_data$pch_fits$XC = arc_data$pch_fits$XC[!sapply(arc_data$pch_fits$XC, is.null)]
    arc_data = lapply(seq(1, length(arc_data$pch_fits$XC)), function(i){
      arc_data = as.data.table(t(arc_data$pch_fits$XC[[i]]))
      arc_data$lab = paste0("archetypes", i)
      arc_data
    })
    arc_data = rbindlist(arc_data)
    data = as.data.table(t(data))
    data$lab = "data"
  }
  rbind(arc_data, data)
}
##' add a line between all archetypes
##'@param arc_data data.table that contains positions of data points and archetypes
##'@param label character specifying which data point to connect (when multiple replicates of fitted polytopes)
.archLines = function(arc_data, label = "archetypes1", type = c("average", "all")[1], average_func = mean){
  arch_lines = arc_data[grepl(label, arc_data$lab)]
  if(type == "average"){
    arch_lines[, arch_id := seq(1, .N), by = lab]
    col_names = colnames(arch_lines)
    col_names = col_names[!col_names %in% c("lab", "arch_id")]
    for (col_name in col_names) {
      arch_lines[, c(col_name) := average_func(get(col_name)), by = "arch_id"]
    }
    arch_lines$lab = label
    arch_lines = unique(arch_lines)
  }
  aa = data.table(1:nrow(arch_lines))
  aa = aa[, .(V2 = 1:nrow(arch_lines)), by = V1]
  aa = aa[V1 != V2]
  aa[, V3 := paste0(sort(c(V1,V2)), collapse = ""), by = .(V1,V2)]
  aa = unique(aa[, .(c(V1, V2)), by = V3])
  arch_lines[aa$V1,]
}
