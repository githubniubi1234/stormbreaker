sb_binary_model_selector_ids <- function() {
  c("lasso_binary_nonzero", "stepglm_binary_aic",
    "glmboost_binary_nonzero", "rf_binary_importance",
    "enet_binary_nonzero", "gbm_binary_importance", "svm_rfe_binary", "gamboost_binary_nonzero", "xgboost_binary_importance", "plsrglm_binary_weights")
}

sb_binary_model_selector_options <- function(id, options = NULL) {
  defaults <- switch(id,
    plsrglm_binary_weights = list(),
    xgboost_binary_importance = list(seed = 1L),
    gamboost_binary_nonzero = list(seed = 1L),
    svm_rfe_binary = list(inner_folds = 5L, seed = 1L),
    lasso_binary_nonzero = list(inner_folds = 5L, seed = 1L, lambda_rule = "lambda.1se"),
    enet_binary_nonzero = list(inner_folds = 5L, seed = 1L, lambda_rule = "lambda.1se", alpha = .5),
    gbm_binary_importance = list(seed = 1L),
    stepglm_binary_aic = list(direction = "both"),
    glmboost_binary_nonzero = list(iterations = 100L, shrinkage = .1),
    rf_binary_importance = list(trees = 500L, seed = 1L, importance = "permutation", threshold = 0))
  if (is.null(defaults)) sb_abort("Unknown binary model selector.", "stormbreaker_unsupported_selector")
  if (!is.null(options) && (!is.list(options) || (length(options) &&
      (is.null(names(options)) || anyNA(names(options)) || anyDuplicated(names(options)) ||
       !all(names(options) %in% names(defaults)))))) sb_abort(
    "Unknown or duplicated selector_options; see sb_pipeline help.", "stormbreaker_bad_selector_options")
  defaults[names(options)] <- options
  for (n in intersect(names(defaults), c("inner_folds", "seed", "iterations", "trees"))) {
    defaults[[n]] <- sb_check_positive_integer(defaults[[n]], n, minimum = if (n == "inner_folds") 3L else 1L)
  }
  if (id %in% c("lasso_binary_nonzero", "enet_binary_nonzero")) defaults$lambda_rule <- sb_validate_glmnet_lambda_rule(defaults$lambda_rule)
  if (id == "enet_binary_nonzero" && (!is.numeric(defaults$alpha) || length(defaults$alpha) != 1L ||
      !is.finite(defaults$alpha) || defaults$alpha <= 0 || defaults$alpha > 1)) {
    sb_abort("Enet selection requires alpha in (0,1].", "stormbreaker_bad_selector_options")
  }
  if (id == "stepglm_binary_aic") defaults$direction <- sb_validate_stepglm_direction(defaults$direction)
  if (id == "glmboost_binary_nonzero" && (!is.numeric(defaults$shrinkage) ||
      length(defaults$shrinkage) != 1L || !is.finite(defaults$shrinkage) ||
      defaults$shrinkage <= 0 || defaults$shrinkage > 1)) sb_abort(
    "shrinkage must be in (0,1].", "stormbreaker_bad_selector_options")
  if (id == "rf_binary_importance") {
    if (!is.character(defaults$importance) || length(defaults$importance) != 1L ||
        is.na(defaults$importance) || !defaults$importance %in% c("permutation", "gini") ||
        !is.numeric(defaults$threshold) || length(defaults$threshold) != 1L ||
        !is.finite(defaults$threshold)) sb_abort(
      "RF needs permutation/gini importance and a finite threshold.", "stormbreaker_bad_selector_options")
  }
  defaults
}

sb_fit_binary_model_selector <- function(x, outcome, id, params, group_ids = NULL, capture_history = FALSE) {
  if (id == "svm_rfe_binary") {
    fit <- sb_fit_svm_rfe(x, outcome, group_ids, params$inner_folds, params$seed)
    fit$params <- params
    return(fit)
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (length(outcome) != nrow(x) || anyNA(outcome) || !all(outcome %in% c(0, 1)) ||
      length(unique(outcome)) != 2L) sb_abort("Binary selection needs both 0/1 classes.", "stormbreaker_bad_binary_outcome")
  original_names <- names(x)
  work <- x
  names(work) <- paste0("feature_", seq_len(ncol(x)))
  evidence <- list()
  if (id %in% c("lasso_binary_nonzero", "enet_binary_nonzero")) {
    alpha <- if (id == "enet_binary_nonzero") params$alpha else 1
    tuned <- sb_select_glmnet_logistic_lambda(work, outcome,
      inner_folds = params$inner_folds, seed = params$seed,
      selection_rule = params$lambda_rule, group_ids = group_ids, alpha = alpha,
      capture_path = capture_history)
    model <- glmnet::glmnet(as.matrix(work), outcome, family = "binomial", alpha = alpha,
      lambda = tuned$lambda, standardize = TRUE)
    score <- as.numeric(stats::coef(model))[-1L]
    selected <- score != 0
    evidence <- list(lambda = tuned$lambda, tuning_curve = tuned$curve,
      fold_id = tuned$fold_id, group_ids = group_ids)
    if (capture_history) {
      evidence$coefficient_path <- tuned$coefficient_path
      rownames(evidence$coefficient_path$beta) <- original_names
    }
  } else if (id == "plsrglm_binary_weights") {
    controls <- sb_specialized_binary_linear_params("plsrglm_binary")
    controls$components <- min(controls$components, ncol(work))
    model <- sb_fit_specialized_binary_linear(work, outcome, "plsrglm_binary", controls)
    score <- abs(model$coefficient * vapply(work, stats::sd, numeric(1)))
    eligible <- which(score > 0 & is.finite(score))
    n_keep <- min(length(eligible), max(2L, ceiling(ncol(work) / 2)))
    chosen <- eligible[order(-score[eligible], eligible)][seq_len(n_keep)]
    selected <- seq_len(ncol(work)) %in% chosen
    evidence <- list(rule = "absolute_coefficient_times_training_sd_top_half_minimum_two",
      components = model$effective$components)
    if (capture_history) evidence$loadings <- model$diagnostic_components$loadings
  } else if (id == "xgboost_binary_importance") {
    native_params <- sb_xgboost_params("xgboost_binary")
    native_params$seed <- params$seed
    model <- sb_fit_xgboost(work, outcome, "xgboost_binary", native_params)
    score <- sb_xgboost_selector_gain(model)
    selected <- score > 0
    evidence <- list(rule = "positive_native_gain", xgboost_version = as.character(utils::packageVersion("xgboost")))
  } else if (id == "gamboost_binary_nonzero") {
    native_params <- sb_classical_binary_params("gamboost_binary")
    native_params$seed <- params$seed
    model <- sb_fit_classical_binary(work, outcome, "gamboost_binary", native_params)$engine
    score <- sb_mboost_spline_selector_scores(model, names(work))
    selected <- score != 0
    loss <- as.numeric(mboost::risk(model))
    evidence <- list(rule = "active_bbs_spline_terms", score_definition = "spline_coefficient_L2_norm",
      training_loss = data.frame(step = seq_along(loss) - 1L, loss = loss))
  } else if (id == "gbm_binary_importance") {
    if (!requireNamespace("gbm", quietly = TRUE)) sb_abort("Install gbm for GBM selection.", "stormbreaker_missing_dependency")
    model <- withr::with_seed(params$seed, gbm::gbm.fit(x = as.matrix(work), y = outcome, distribution = "bernoulli",
      n.trees = 100L, interaction.depth = 1L, n.minobsinnode = 10L,
      shrinkage = .1, bag.fraction = 1, nTrain = nrow(work), keep.data = FALSE,
      verbose = FALSE, var.names = names(work)))
    importance <- summary(model, n.trees = 100L, plotit = FALSE)
    score <- importance$rel.inf[match(names(work), importance$var)]
    selected <- score > 0
    evidence <- list(rule = "positive_relative_influence", trees = 100L,
      training_loss = data.frame(trees = seq_along(model$train.error), loss = model$train.error))
  } else if (id %in% c("stepglm_binary_aic", "glmboost_binary_nonzero")) {
    if (id == "stepglm_binary_aic") {
      native <- sb_fit_stepglm_binary(list(x = work, outcome = outcome),
        list(direction = params$direction, k = 2, max_steps = 1000L), capture_history)
    } else {
      if (!requireNamespace("mboost", quietly = TRUE)) sb_abort("Install mboost for glmBoost selection.", "stormbreaker_missing_dependency")
      native <- sb_fit_glmboost_binary(list(x = work, outcome = outcome),
        list(iterations = params$iterations, shrinkage = params$shrinkage, center = FALSE), capture_history)
    }
    score <- as.numeric(native$coefficient[names(work)])
    selected <- if (id == "stepglm_binary_aic") names(work) %in% native$effective$selected_features else score != 0
    evidence <- native$effective
    evidence$selected_features <- original_names[match(evidence$selected_features, names(work))]
  } else {
    if (!requireNamespace("randomForest", quietly = TRUE)) sb_abort("Install randomForest for RF selection.", "stormbreaker_missing_dependency")
    model <- withr::with_seed(params$seed, randomForest::randomForest(x = work, y = factor(outcome, levels = c(0, 1)),
      ntree = params$trees, importance = TRUE, keep.forest = FALSE))
    importance <- randomForest::importance(model, type = if (params$importance == "permutation") 1 else 2, scale = FALSE)
    score <- as.numeric(importance[names(work), 1L])
    selected <- score > params$threshold
    if (capture_history) evidence$oob_error <- data.frame(
      trees = seq_len(nrow(model$err.rate)), model$err.rate, check.names = FALSE)
  }
  if (length(score) != ncol(x) || any(!is.finite(score))) sb_abort(
    "Selector returned invalid feature scores.", "stormbreaker_invalid_selector_scores")
  if (!any(selected)) sb_abort("No features met the declared selection rule; no fallback genes were added.", "stormbreaker_empty_selection")
  structure(list(id = id, source_features = original_names,
    selected_features = original_names[selected], params = params,
    scores = data.frame(feature = original_names, score = score, selected = selected),
    evidence = evidence), class = "stormbreaker_selector_fit")
}

sb_binary_model_selector_adapter <- function(id, params, capture_history = FALSE) {
  list(id = id, params = params, capture_history = capture_history, needs_outcome = TRUE,
    fit = function(x, outcome, group_ids = NULL) {
      if (capture_history) sb_fit_binary_model_selector(x, outcome, id, params, group_ids, TRUE)
      else sb_fit_binary_model_selector(x, outcome, id, params, group_ids)
    },
    apply = function(fitted, x) {
      x <- sb_validate_predictor_frame(x, "x")
      if (!inherits(fitted, "stormbreaker_selector_fit") || fitted$id != id ||
          !identical(names(x), fitted$source_features)) sb_abort(
        "Predictor columns do not match the fitted selector.", "stormbreaker_selector_schema_mismatch")
      x[fitted$selected_features]
    })
}
