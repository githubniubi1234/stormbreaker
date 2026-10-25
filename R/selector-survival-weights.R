sb_survival_weight_selector_ids <- function() c("ridge_cox_weights", "plsrcox_weights", "svm_survival_weights")

sb_survival_weight_selector_params <- function(id) {
  model <- switch(id,
    ridge_cox_weights = list(alpha = 0, inner_folds = 5L, seed = 1L, lambda_rule = "lambda.min"),
    plsrcox_weights = sb_plsrcox_params(),
    svm_survival_weights = sb_survivalsvm_params())
  list(model = model, fraction = .5, minimum = 2L, rule = "absolute_coefficient_times_training_sd")
}

sb_survival_weight_selector_adapter <- function(id, params, capture_history = FALSE) {
  list(id = id, params = c(params, if (capture_history) list(capture_history = TRUE)), needs_outcome = TRUE,
    fit = function(x, time, event, group_ids = NULL) {
      x <- sb_validate_predictor_frame(x, "x")
      evidence <- list()
      if (id == "ridge_cox_weights") {
        tuned <- sb_fit_lasso_nonzero(x, time, event, params$model$inner_folds,
          params$model$seed, group_ids, alpha = 0, id = id, capture_history = capture_history)
        beta <- tuned$scores$coefficient[match(names(x), tuned$scores$feature)]
        evidence <- list(lambda = tuned$lambda, tuning_curve = tuned$tuning_curve,
          fold_id = tuned$fold_id, group_ids = group_ids, coefficient_path = tuned$coefficient_path)
      } else if (id == "plsrcox_weights") {
        controls <- params$model
        controls$components <- min(controls$components, ncol(x))
        model <- sb_fit_plsrcox(x, time, event, controls)
        beta <- model$coefficient
        evidence$components <- controls$components
        if (capture_history) evidence$loadings <- model$engine$pp
      } else {
        model <- sb_fit_survivalsvm(x, time, event, params$model)
        beta <- model$coefficient
      }
      sd <- vapply(x, stats::sd, numeric(1))
      score <- abs(beta * sd)
      if (length(score) != ncol(x) || any(!is.finite(score))) stop("Native screening weights are unavailable.")
      eligible <- which(score > 0)
      if (!length(eligible)) stop("No nonzero screening weights; no fallback genes were added.")
      n_keep <- min(length(eligible), max(params$minimum, ceiling(ncol(x) * params$fraction)))
      chosen <- eligible[order(-score[eligible], eligible)][seq_len(n_keep)]
      selected <- seq_len(ncol(x)) %in% chosen
      structure(list(id = id, params = params, source_features = names(x),
        selected_features = names(x)[selected],
        scores = data.frame(feature = names(x), coefficient = as.numeric(beta),
          training_sd = unname(sd), score = unname(score), selected = selected),
        evidence = evidence), class = "stormbreaker_selector_fit")
    },
    apply = function(fitted, x) {
      x <- sb_validate_predictor_frame(x, "x")
      if (!identical(names(x), fitted$source_features)) stop("Predictors differ from the fitted screening model.")
      x[fitted$selected_features]
    })
}
