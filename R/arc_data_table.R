#use environment() to store these tables
##'
##'
##' Project in low dimentions (PCA) and process data and polytopes fits for plotting
##' @param arc_data list of matrices storing the position of archetypes, class "pch_fit", "r_pch_fit". Each element of a list represents an independent run of the polytope fitting algorithm
##' @param data matrix of data in which archetypes/polytope were found
project = function(arc_data, data){

}
process_for_plotly = function(arc_data, data){

}

.arc_data_table = function(arc_data, data){
  archetypes = archetypes[!sapply(archetypes, is.null)]
  if(is(archetypes, "pch_fit")){

  }
  if(is(archetypes, "r_pch_fit")){

  }
  # transpose matrices from pcha
  nrepl = length(archetypes)
  narch = nrow(archetypes[[1]])
  archetypes = Reduce(rbind, archetypes)
  as.data.table(archetypes)
  archetypes = data.table(PC1 = archetypes[,1],
                          PC2 = archetypes[,2],
                          PC3 = archetypes[,3],
                          lab = c(rep(paste0("archetypes", 1:nrepl),
                                      each = narch)))
  setorder(archetypes, lab, PC1, PC2, PC3)
  #archetypes[, archetype_num := c(rep(paste0("archetypes", 1:narch),
  #                                   times = nrepl))]
  as.data.table(mat)
  Arch_data = data.table(PC1 = mat[,1],
                         PC2 = mat[,2],
                         PC3 = mat[,3],
                         lab = rep("perturbations", nrow(mat)))
  rbind(archetypes, Arch_data)
}
##' add a line between all archetypes
##'@param Arch_data data.table that contains positions of data points and archetypes
##'@param label character specifying which data point to connect (when multiple replicates of fitted polytopes)
.archLines = function(Arch_data, label = "archetypes1"){
  arch_lines = Arch_data[lab == label]
  aa = data.table(1:nrow(arch_lines))
  aa = aa[, .(V2 = 1:nrow(arch_lines)), by = V1]
  aa = aa[V1 != V2]
  aa[, V3 := paste0(sort(c(V1,V2)), collapse = ""), by = .(V1,V2)]
  aa = unique(aa[, .(c(V1, V2)), by = V3])
  arch_lines[aa$V1,]
}

.archLinesDT = function(Arch_data, archetypes, line_obj_name = "arch_lines_", envir = .GlobalEnv){
  nrepl = length(archetypes)
  for (i in 1:nrepl) {
    var = paste0(line_obj_name, i)
    val = paste0("archetypes", i)
    assign(var, archLines(Arch_data, label = val), envir = envir)
  }
}
