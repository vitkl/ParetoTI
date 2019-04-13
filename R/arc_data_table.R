.arc_data_table = function(arc_data, data,
                           data_lab = "data", arc_lab = "archetypes",
                           which_dimensions = 1:2){
  # if single fit just add label column
  if(is(arc_data, "pch_fit") | is(arc_data, "random_arc")){

    arch_ids = colnames(arc_data$XC)

    arc_data = as.data.table(t(arc_data$XC[which_dimensions,]))
    # user label for archetypes is only used when looking at single fit
    arc_data$lab = arc_lab

    # add archetype id if archetypes have names
    arc_data[, arch_id := ifelse(rep(is.null(arch_ids), .N), seq(1, .N), arch_ids)]

  }
  # if multiple fits from bootstrap or different k - combine all matrices
  # in one data.table with each matrix getting unique label
  if(is(arc_data, "b_pch_fit")){

    arc_data$pch_fits$XC = arc_data$pch_fits$XC[!sapply(arc_data$pch_fits$XC, is.null)]
    arc_data = lapply(seq(1, length(arc_data$pch_fits$XC)), function(i){

      arch_ids = colnames(arc_data$pch_fits$XC[[i]])

      arc_data = as.data.table(t(arc_data$pch_fits$XC[[i]][which_dimensions,]))
      arc_data$lab = paste0("archetypes", i)

      # add archetype id column if archetypes have names
      arc_data[, arch_id := ifelse(rep(is.null(arch_ids), .N), seq(1, .N), arch_ids)]
      arc_data

    })
    arc_data = rbindlist(arc_data)
  }
  if(is(arc_data, "k_pch_fit")){
    arc_data$pch_fits$XC = arc_data$pch_fits$XC[!sapply(arc_data$pch_fits$XC, is.null)]
    arc_data = lapply(seq(1, length(arc_data$pch_fits$XC)), function(i){

      if(nrow(arc_data$pch_fits$XC[[i]]) < length(which_dimensions)) return(NULL)
      arc_data = as.data.table(t(arc_data$pch_fits$XC[[i]][which_dimensions,]))

      # label by number of archetypes rather than iteration
      arc_data$lab = paste0("archetypes", nrow(arc_data))
      # add archetype id column if archetypes have names
      arc_data[, arch_id := seq(1, .N)]

      arc_data
    })
    arc_data = arc_data[!sapply(arc_data, is.null)]
    arc_data = rbindlist(arc_data)
  }
  data = as.data.table(t(data[which_dimensions,]))
  data$lab = data_lab
  list(arc_data = arc_data, # archetypes
       data = data) # data points
}
##' add a line between all archetypes
##'@param arc_data data.table that contains positions of archetypes
##'@param pch_fit not used
##'@param arc_lab labels for archetypes
##'@param type average across arc_lab (average)?
.archLines = function(arc_data, arc_lab = "archetypes", pch_fit = TRUE,
                      type = c("average", "all")[1], average_func = mean){

  arch_lines = copy(arc_data)

  if(type == "average"){

    col_names = colnames(arch_lines)
    col_names = col_names[!col_names %in% c("lab", "arch_id")]
    for (col_name in col_names) {
      arch_lines[, c(col_name) := average_func(get(col_name)), by = "arch_id"]
    }

    arch_lines$lab = arc_lab
    arch_lines = unique(arch_lines)

  }

  aa = data.table(1:nrow(arch_lines))
  aa = aa[, .(V2 = 1:nrow(arch_lines)), by = V1]
  aa = aa[V1 != V2]
  aa[, V3 := paste0(sort(c(V1,V2)), collapse = ""), by = .(V1,V2)]
  aa = unique(aa[, .(c(V1, V2)), by = V3])
  arch_lines[aa$V1,]
}
