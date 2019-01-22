##' Logistic regression model to classify vertices
##' @rdname predict_vert
##' @name build_logit_vert
##' @description build_logit_vert() defines keras model
##' @export build_logit_vert
##' @import data.table
build_logit_vert = function(features_list, penalty = 0.7,
                            regularizer = keras::regularizer_l2,
                            initializer = "random_uniform",
                            activation = "softmax",
                            loss = "categorical_crossentropy",
                            metrics = list("accuracy")) {
  model = keras::keras_model_sequential() %>%
    keras::layer_dense(units = ncol(features_list$y), activation = activation,
                       kernel_regularizer = regularizer(l = penalty),
                       bias_regularizer = regularizer(l = penalty),
                       kernel_initializer = initializer,
                       bias_initializer = initializer,
                       input_shape = ncol(features_list$X))
  model = keras::compile(model, loss = loss,
                         optimizer = keras::optimizer_sgd(nesterov = TRUE),
                         metrics = metrics)
  model
}

##' @rdname predict_vert
##' @name make_features
##' @description make_features() generate features for predicting vertex by expression profile of several genes with distance from that vertex
##' @export make_features
##' @import data.table
make_features = function(data_attr, features_list, n_bins = 10, n_samples = 10,
                         select_bins = seq_len(n_bins),
                         combine_bins = seq_len(n_bins)){
  # define predictors  ---------------------------------------------
  X_dt = data_attr$data[, c(data_attr$arc_col, features_list), with = F]


  # melt archetypes: example = example of each archetype (n_samples * n_vert total)
  X_dt = melt.data.table(X_dt, id.vars = c(features_list),
                         variable.name = "label", value.name = "distance")
  # order gene expression by distance from archetype (new var = order)
  X_dt[, order := frank(distance, ties.method = "average"), by = .(label)]
  setorder(X_dt, label, order)

  # put cells into bins
  X_dt[, bin := {
    max_ord = .N
    bins = rep(seq_len(n_bins), each = round(max_ord / n_bins))
    if(length(bins) > max_ord) bins = bins[seq_len(max_ord)]
    c(bins, rep(n_bins, max_ord - length(bins)))
  }, by = .(label)]
  # select a subset of bins
  X_dt = X_dt[bin %in% select_bins]
  # combine bins, e.g. 1->1, 2-9->2, 10->10
  X_dt[, bin := combine_bins[bin]]

  # create examples * features matrix ---- deterministic
  X_dt[, example := sample.int(n_samples, .N, replace = T), by = .(bin, label)]
  X_dt[, example := paste0(label, "_",example)]
  # create examples * features matrix ---- bootstraping
  #to do

  # melt genes (keeping order within genes - name preditors by pasting gene_bin)
  X_dt = melt.data.table(X_dt, id.vars = c("example", "label", "distance", "order", "bin"),
                         variable.name = "gene", value.name = "expression")
  X_dt[, predictor := paste0(gene, "_", bin)]
  setorder(X_dt, example, gene, bin)

  # calculate mean expression in each bin (predictor = gene _ bin)
  X_dt[, pred_val := mean(expression), by = .(predictor, example)]
  X_dt[, mean_dist := mean(distance), by = .(predictor, example)]


  X_dt = unique(X_dt[, .(example, label, predictor, pred_val, mean_dist, gene)])
  # Scale (z-score) expression in bins between 0 and 1 to ensure other datasets work better
  X_dt[, pred_val := scale(pred_val, TRUE, TRUE), by = .(gene, example)]

  # dcast examples and transpose to produce examples * features matrix
  X = dcast.data.table(X_dt, example ~ predictor, value.var = "pred_val")
  X = as.matrix(X, rownames = "example")

  # define labels ---------------------------------------------
  # produce matrix with vertex classes
  y = dcast.data.table(unique(X_dt[, .(example, label)]),
                       example ~ label, fun.aggregate = length,
                       value.var = "label", fill = 0)
  y = as.matrix(y, rownames = "example")
  list(X = X, y = y, X_dt = X_dt,
       settings = list(features_list = features_list, n_bins = n_bins,
                       n_samples = n_samples, select_bins = select_bins,
                       combine_bins = combine_bins))
}

##' @rdname predict_vert
##' @name predict_vert
##' @description predict_vert() predict vertices using a model and attribute data (distance to vertex + feature of cells)
##' @export predict_vert
##' @import data.table
predict_vert = function(vert_model, data_attr, features_data) {
  # generate features the same way as feature data
  all_X = make_features(data_attr,
                        features_list = features_data$settings$features_list,
                        n_bins = features_data$settings$n_bins, n_samples = 1,
                        select_bins = features_data$settings$select_bins,
                        combine_bins = features_data$settings$combine_bins)

  # predict classes
  classes = predict_classes(vert_model, all_X$X, batch_size = NULL, verbose = 0,
                            steps = NULL)
  names(classes) = colnames(features_data$y)[classes + 1]

  probabilities = predict(vert_model, all_X$X,
                          batch_size = NULL, verbose = 0, steps = NULL)
  rownames(probabilities) = rownames(all_X$X)
  colnames(probabilities) = colnames(features_data$y)
  list(probabilities = probabilities, classes = classes)
}

##' @rdname predict_vert
##' @name split_train_val
##' @description split_train_val() split data into training and validation. Shuffle examples to ensure that model sees all vertices in training and validation
##' @param val_prop proportion of samples used for validation.
##' @param per_class Use proportion of samples within each class for validation? By default the function takes class into account when splitting data (TRUE)
##' @export split_train_val
##' @import data.table
split_train_val = function(features_data, val_prop = 0.3, per_class = TRUE) {

  if(isTRUE(per_class)) {
    val_ind = as.integer(unlist(apply(features_data$y, 2, FUN = function(x){
      class_ind = which(x == 1)
      sample(class_ind, round(length(class_ind) * val_prop, digits = 0))
    })))
  } else {
    val_ind = sample.int(nrow(features_data$y),
                         round(nrow(features_data$y) * val_prop, digits = 0))
  }

  list(val_data = list(features_data$X[val_ind,],
                       features_data$y[val_ind,]),
       train_data = list(X = features_data$X[-val_ind,],
                         y = features_data$y[-val_ind,]))
}

##' @rdname predict_vert
##' @name plot_confusion_vert
##' @description plot_confusion_vert() plot vertex probabilities as confusion matrix
##' @export plot_confusion_vert
##' @import data.table
plot_confusion_vert = function(predicted) {
  predicted = as.data.table(predicted, keep.rownames = "observed_vertex")
  predicted = melt.data.table(predicted, id.vars = "observed_vertex",
                              variable.name = "predicted_vertex",
                              value.name = "probability")
  predicted[, observed_vertex := gsub("_1$", "", observed_vertex)]
  ggplot2::ggplot(predicted, ggplot2::aes(predicted_vertex, observed_vertex,
                                          color = probability, fill = probability))+
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c() + ggplot2::scale_color_viridis_c() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = -25))
}

##' @rdname predict_vert
##' @name get_feature_weights
##' @description get_feature_weights() get gene weights and ranked features from keras model
##' @export get_feature_weights
##' @import data.table
get_feature_weights = function(vert_model, features_data){
  # extract weights ---------------------------------------------
  bias = get_weights(vert_model)[[2]]
  gene_weights = get_weights(vert_model)[[1]]
  rownames(gene_weights) = colnames(features_data$X)
  colnames(gene_weights) = colnames(features_data$y)

  # top-10 genes for each cluster
  top = data.table(ind = seq_len(nrow(gene_weights)))
  for (col in colnames(gene_weights)) {
    ind = order(gene_weights[,col], decreasing = T)
    top = cbind(top, rownames(gene_weights)[ind])
    setnames(top, "V2", col)
  }
  list(top = top, gene_weights = gene_weights)
}
