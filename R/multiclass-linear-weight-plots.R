sb_multiclass_linear_weight_plots <- function(fit) {
  if (fit$method == "xgboost_linear") {
    booster <- sb_load_xgboost_state_model(list(model_raw = fit$engine,
      params = fit$parameters, features = fit$features), "xgboost_binary_linear")$booster
    model <- jsonlite::fromJSON(rawToChar(xgboost::xgb.save.raw(booster, raw_format = "json")))
    weights <- model$learner$gradient_booster$model$weights
    coefficients <- matrix(weights, nrow = length(fit$classes))
    classes <- fit$classes
    features <- c(fit$features, "Bias")
    note <- "Original full-precision gblinear class-margin weights and learned biases from the saved model. The global base score is separate. These are not tree gain or reference-class logistic coefficients."
  } else {
    coefficients <- fit$engine$W
    features <- colnames(coefficients)
    at <- match(features, fit$feature_codes)
    features[!is.na(at)] <- fit$features[at[!is.na(at)]]
    classes <- fit$classes[match(as.character(fit$engine$ClassNames), fit$codes)]
    note <- "Original LIBLINEAR one-versus-rest weights on standardized inputs, including the bias weight. Not multinomial reference-class coefficients; class probabilities are subsequently normalized by the native predictor."
  }
  d <- data.frame(feature = rep(features, each = nrow(coefficients)),
    class = rep(classes, ncol(coefficients)), weight = as.vector(coefficients))
  d$feature <- factor(d$feature, levels = rev(features))
  feature <- weight <- NULL
  p <- ggplot2::ggplot(d, ggplot2::aes(weight, feature)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
    ggplot2::geom_point(colour = "#2878A0") + ggplot2::facet_wrap(~class) +
    sb_viz_theme() + ggplot2::labs(title = "Class-specific linear weights", x = "Weight",
      y = NULL, caption = note)
  list(plots = list(linear_weights = p), tables = list(linear_weights = d))
}
