##' Classify cells with TensorFlow (logistic regression)
##' @rdname fit_logistic_model
##' @name fit_logistic_model
##' @description fit_logistic_model() Uses TensorFlow to fit logistic regression model for classifying cells by attribute in colData slot of \code{sce}. Most arguments have sensible defaults. Use history plot to determine if model is learning and that it performs equally well on training and validation sets.
##' @param sce \link[SingleCellExperiment]{SingleCellExperiment} object.
##' @param y Optionally, you can provide your own cell labels as \code{y} matrix (dim = cells * classes) where each row should sum to 1. This can be continuous cell labels such as those from archetypal analysis (\link[ParetoTI]{fit_pch}) or NMF (non-negativa matrix factorisation). In such case, kullback_leibler_divergence is a more suitable cost function.
##' @param assay_slot Assay slot in \code{sce} containing gene expression matrix (default is logcounts).
##' @param y_col Column in\code{colData(sce)} containing cell labels.
##' @param activation Activation function. Using "softmax" gives logistic regression.
##' @param loss Loss or Cost function that evaluates the difference between preditions and true labels and is minimised during model training. Use \code{"categorical_crossentropy"} for discrete class labels and \code{"kullback_leibler_divergence"} for continuous labels that sum to 1.
##' @param regularizer Function to penalise high values of model parameters/weights and reduce overfitting (memorisation of the dataset). Details: \link[keras]{regularizer_l1}. L1 regularisation tends to push most weights to 0 (thus acting as feature selection method) and enforce sparse weights. L2 regularisation also reduces overfitting by keeping most weights small but does not shrink them to 0. Set to \code{NULL} to not regularise. Both weights and bias are regularised the same way.
##' @param penalty Regularisation penalty between 0 and 1. The higher the penalty the more stringent is the regularisation. Very high values can lead to poor model performance due to high bias (limited flexibility). Sensible values: 0.01 for regularizer_l1 and 0.5 for regularizer_l2. Change this parameter based history plot to make sure the model performs equally well on training and validation sets.
##' @param initializer Method of initialising weights and bias. You do not normally need to change this. See https://keras.io/initializers/ for details.
##' @param optimizer Which optimiser should be used to fit the model. You do not normally need to change this. See https://keras.io/optimizers/ for details.
##' @param metrics Metrics that evaluate performance of the model. Usually this is accuracy and loss function. See https://keras.io/metrics/ for details
##' @param epochs Number of training epochs. You do not normally need to change this.
##' @param validation_split What proportion of cells should be used for validation.
##' @param validation_split_per_class Do the validation split within each class to maintains proportion of classes in training and validation sets (TRUE)?
##' @param callback_early_stopping_patience Number of epochs to wait for improvement before stopping early.
##' @param batch_size Look at data in batches of \code{batch_size} cells. All batches will be seen in each epoch.
##' @param shuffle Logical, whether to shuffle the training data before each epoch? Details: \link[keras]{fit.keras.engine.training.Model}
##' @param verbose Logical. Show and plot diagnosic output (TRUE)?
##' @param model Provide your own keras/tensorflow model. Output units must be equal to the number of classes (columns of y), input_shape must be equal to nrow(sce). This can be used to extend logistic regression model by adding hidden layers
##' @export fit_logistic_model
##' @import data.table
##' @import SingleCellExperiment
##' @examples
##' # download PBMC data as SingleCellExperiment object
##'
##' # split in 2 parts
##'
##' # fit logistic regression model to 1st part
##'
##' # use this model to predict cell types in the second part
##'
fit_logistic_model = function(sce, y = NULL,
                              assay_slot = "logcounts", y_col = "Cell_class",
                              activation = "softmax",
                              loss = c("categorical_crossentropy", "kullback_leibler_divergence")[1],
                              regularizer = keras::regularizer_l1,
                              penalty = 0.01, initializer = "random_uniform",
                              optimizer = keras::optimizer_sgd(lr = 0.01, nesterov = TRUE),
                              metrics = list("accuracy"), epochs = 100,
                              validation_split = 0.3, validation_split_per_class = TRUE,
                              callback_early_stopping_patience = 50,
                              batch_size = 1000, shuffle = TRUE,
                              verbose = TRUE, model = NULL){

  if(is.null(y)) { # if label matrix is not provided - create
    # generate label matrix
    y = data.table(ident = colData(sce)[, y_col], cell_id = colnames(sce))
    y = dcast.data.table(y, cell_id ~ ident, fun.aggregate = length)
    y = as.matrix(y, rownames = "cell_id")
    y = as(y, "dgCMatrix")
  }

  features_data = list(y = y, X = Matrix::t(as(assay(sce, assay_slot), "dgCMatrix")))

  if(!is.null(regularizer)) regularizer = regularizer(l = penalty)

  if(is.null(model)) {
    # define logistic regression keras model
    build_log_reg = function() {
      model = keras::keras_model_sequential()
      model = keras::layer_dense(object = model, units = ncol(y), activation = activation,
                                 kernel_regularizer = regularizer,
                                 bias_regularizer = regularizer,
                                 kernel_initializer = initializer,
                                 bias_initializer = initializer,
                                 input_shape = nrow(sce))
      model = keras::compile(object = model, loss = loss,
                             optimizer = optimizer,
                             metrics = metrics)
      model
    }

    model = build_log_reg()
  }

  # The patience parameter is the amount of epochs to check for improvement.
  early_stop = keras::callback_early_stopping(monitor = "val_loss",
                                              patience = callback_early_stopping_patience)

  # produce training and validation sets
  t_v_set = split_train_val(features_data, val_prop = validation_split,
                            per_class = validation_split_per_class)

  # Fit the model and store training stats
  history = keras::fit(model, x = t_v_set$train_data$X,
                       y = t_v_set$train_data$y,
                       epochs = epochs,
                       batch_size = batch_size,
                       validation_data = t_v_set$val_data,
                       verbose = 0, shuffle = shuffle,
                       callbacks = list(early_stop))

  if(verbose) print(plot(history))

  # extract weights
  bias = keras::get_weights(model)[[2]]
  gene_weights = keras::get_weights(model)[[1]]
  rownames(gene_weights) = rownames(sce)
  colnames(gene_weights) = colnames(y)

  # top genes for each cluster
  top = data.table(ind = seq_len(nrow(gene_weights)))
  for (col in colnames(gene_weights)) {
    ind = order(gene_weights[,col], decreasing = T)
    gene_names = rownames(gene_weights)
    gene_names[gene_weights[,col] < 1e-3 & gene_weights[,col] > -1e-3] = NA
    top = cbind(top, gene_names[ind])
    setnames(top, "V2", col)
  }

  if(verbose) print(top[1:10])

  # Caclculate probability for each cluster
  probabilities = keras::predict_proba(model, Matrix::t(as(assay(sce, assay_slot), "dgCMatrix")),
                                       batch_size = batch_size,
                                       verbose = 0, steps = NULL)
  colnames(probabilities) = colnames(y)
  rownames(probabilities) = colnames(sce)

  if(verbose){
    # cells with strong prediction for several classes
    p1 = ggplot2::qplot(rowSums(probabilities > 0.3), geom = "histogram") +
      ggtitle("N cells with strong prediction for one or several classes (p > 0.3)")
    print(p1)
    # and max prediction for any class
    p2 = ggplot2::qplot(rowMax(probabilities), geom = "histogram") +
      ggtitle("Max probability for any class")
    print(p2)
  }

  # use entropy to measure strength in predictions
  class_entropy = apply(probabilities, 1, FUN = entropy::entropy.plugin)

  # compute class labels by selecting a label with max probability
  pred_class = colnames(probabilities)[apply(probabilities, 1, which.max)]

  # calculate and plot confusion matrix for classes
  confusion = table(colData(sce)[, y_col], pred_class)
  if(verbose) print(plot_confusion(confusion) + ggtitle("Confusion between classes"))

  structure(list(probabilities = probabilities,
                 pred_class = pred_class, class_entropy = class_entropy,
                 gene_weights = gene_weights, top = top, model = model,
                 history = history, call = match.call(), y = y,
                 confusion = confusion), class = "logistic_model_fit_TF")
}

##' @rdname fit_logistic_model
##' @name plot_confusion
##' @description plot_confusion() plot confusion matrix: number of cells assigned to observed class and predicted class
##' @param confusion The confusion table generated by table(), OR the output of ParetoTI::fit_logistic_model(), OR the output of ParetoTI::predict_logistic_prob()
##' @param normalize Normalise so that cell of each observed class sum to 1?
##' @param text_color Color of on-plot text showing absolute numbers of cells.
##' @export plot_confusion
##' @import data.table
plot_confusion = function(confusion, normalize = FALSE, text_color = "grey60") {

  if(is(confusion, "logistic_model_fit_TF")) {
    confusion = confusion$confusion
  } else if(is(confusion, "logistic_model_prediction")) {
    confusion = confusion$confusion
  } else if(is(confusion, "table")) {
    confusion = confusion
  } else stop("confusion argument should be the confusion table generated by table(), OR the output of ParetoTI::fit_logistic_model(), OR output of ParetoTI::predict_logistic_prob()")


  confusion_dt = as.data.table(confusion)
  setnames(confusion_dt, colnames(confusion_dt), c("observed_class", "predicted_class", "N"))

  if(normalize){
    confusion_dt[, proportion := N / sum(N), by = "observed_class"]
    g = ggplot2::ggplot(confusion_dt,
                        ggplot2::aes(predicted_class, observed_class,
                                     color = proportion, fill = proportion))
  } else {
    g = ggplot2::ggplot(confusion_dt, ggplot2::aes(predicted_class, observed_class,
                                                   color = N, fill = N))
  }

  g + ggplot2::geom_tile() +
    ggplot2::geom_text(aes(label = N, x = predicted_class, y = observed_class),
                       color = text_color) +
    ggplot2::scale_fill_viridis_c() + ggplot2::scale_color_viridis_c() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = -25))
}

##' @rdname fit_logistic_model
##' @name predict_logistic_prob
##' @description predict_logistic_prob() predict class assignments of cells in SingleCellExperiment object. Genes used to train the model are selected automatically.
##' @param model_res Output of ParetoTI::fit_logistic_model(), class "logistic_model_fit_TF"
##' @param ref_y Reference cell labels as \code{ref_y} matrix (dim = cells * classes) where each row should sum to 1. Optional, it can be used to match continuous cell labels across datasets.
##' @param ref_y_col Reference class column in \code{colData(sce)}. For example, not annotated clusters assigned by a clustering algorhitm. Optional, it can be used to match discrete cell labels across datasets.
##' @export predict_logistic_prob
##' @import data.table
predict_logistic_prob = function(sce, model_res, assay_slot = "logcounts",
                                 ref_y = NULL, ref_y_col = NULL,
                                 batch_size = 1000, verbose = TRUE){

  if(!is(model_res, "logistic_model_fit_TF")) stop("model_res should be output of ParetoTI::fit_logistic_model()")

  gene_names = rownames(model_res$gene_weights)

  # Caclculate probability for each cluster
  probabilities = keras::predict_proba(model_res$model,
                                       Matrix::t(as(assay(sce[gene_names,], assay_slot), "dgCMatrix")),
                                       batch_size = batch_size,
                                       verbose = 0, steps = NULL)
  colnames(probabilities) = colnames(model_res$y)
  rownames(probabilities) = colnames(sce)

  if(verbose){
    # cells with strong prediction for several classes
    p1 = ggplot2::qplot(rowSums(probabilities > 0.3), geom = "histogram") +
      ggtitle("N cells with strong prediction for one or several classes (p > 0.3)")
    print(p1)
    # and max prediction for any class
    p2 = ggplot2::qplot(rowMax(probabilities), geom = "histogram") +
      ggtitle("Max probability for any class")
    print(p2)
  }

  # use entropy to measure strength in predictions
  class_entropy = apply(probabilities, 1, FUN = entropy::entropy.plugin)

  # compute class labels by selecting a label with max probability
  pred_class = colnames(probabilities)[apply(probabilities, 1, which.max)]

  # calculate and plot confusion matrix for classes
  if(!is.null(ref_y_col)) {
    # discrete labels
    confusion = table(colData(sce)[, ref_y_col], pred_class)
    if(verbose) print(plot_confusion(confusion) + ggtitle("Confusion between classes"))

  } else if(!is.null(ref_y)) {

  } else confusion = NULL

  structure(list(probabilities = probabilities,
                 pred_class = pred_class, class_entropy = class_entropy,
                 call = match.call(), confusion = confusion),
            class = "logistic_model_prediction")
}
