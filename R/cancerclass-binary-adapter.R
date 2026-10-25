sb_cancerclass_available <- function() {
  requireNamespace("cancerclass", quietly = TRUE) &&
    requireNamespace("Biobase", quietly = TRUE)
}

sb_cancerclass_binary_params <- function() {
  list(
    feature_method = "welch.test",
    hparam = 0.75,
    distance = "cor",
    feature_count = "all",
    score = "z",
    calibration = "oof_logistic",
    inner_folds = 5L,
    inner_seed = 1L
  )
}

sb_cancerclass_expression_set <- function(x, outcome) {
  sample_names <- paste0("sample_", seq_len(nrow(x)))
  expression <- t(data.matrix(x))
  colnames(expression) <- sample_names
  phenotype <- data.frame(
    class = factor(outcome, levels = c(0L, 1L)),
    row.names = sample_names
  )
  Biobase::ExpressionSet(
    assayData = expression,
    phenoData = Biobase::AnnotatedDataFrame(phenotype)
  )
}

sb_fit_cancerclass_native_predictor <- function(x, outcome, params) {
  search_before <- search()
  withr::local_package("Biobase", quietly = TRUE)
  attached_dependencies <- setdiff(
    search(),
    c(search_before, "package:Biobase")
  )
  withr::defer({
    for (package in attached_dependencies) {
      if (package %in% search()) {
        base::detach(package, character.only = TRUE)
      }
    }
  }, priority = "last")
  fitted <- suppressMessages(cancerclass::fit(
    sb_cancerclass_expression_set(x, outcome),
    class = "class",
    method = params$feature_method,
    hparam = params$hparam
  ))
  predictor <- if (methods::is(fitted, "predictor")) {
    methods::slot(fitted, "predictor")
  } else {
    NULL
  }
  valid <- is.matrix(predictor) && nrow(predictor) == ncol(x) &&
    setequal(rownames(predictor), names(x)) &&
    all(c(params$feature_method, "0", "1") %in% colnames(predictor))
  if (!valid) {
    sb_abort(
      "`cancerclass` did not return one usable centroid row per feature.",
      "stormbreaker_invalid_cancerclass_fit"
    )
  }
  predictor
}

sb_cancerclass_score_from_predictor <- function(
    predictor, x, features, distance) {
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), features)) {
    sb_abort(
      "Prediction columns do not match the fixed cancerclass schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  data <- t(data.matrix(x))
  data <- data[rownames(predictor), , drop = FALSE]
  distance_zero <- cancerclass::get.d(
    data,
    predictor[, "0"],
    method = distance
  )
  distance_one <- cancerclass::get.d(
    data,
    predictor[, "1"],
    method = distance
  )
  class_distance <- cancerclass::get.d(
    predictor[, "0"],
    predictor[, "1"],
    method = distance
  )
  score <- (distance_zero - distance_one) / class_distance
  valid <- length(class_distance) == 1L && !is.na(class_distance) &&
    is.finite(class_distance) && class_distance > 0 &&
    length(score) == nrow(x) && !anyNA(score) && all(is.finite(score))
  if (!valid) {
    sb_abort(
      "`cancerclass` produced an invalid native correlation-distance score.",
      "stormbreaker_invalid_cancerclass_score"
    )
  }
  unname(score)
}

sb_fit_cancerclass_binary <- function(x, outcome, params) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      "A fixed cancerclass analysis fold must contain both classes.",
      "stormbreaker_cancerclass_single_class"
    )
  }
  if (!sb_cancerclass_available()) {
    sb_abort(
      paste0(
        "Packages `cancerclass` and `Biobase` are required for the fixed ",
        "cancerclass learner."
      ),
      "stormbreaker_missing_cancerclass"
    )
  }

  fold_id <- sb_make_stratified_fold_id(
    outcome,
    folds = params$inner_folds,
    seed = params$inner_seed,
    stratum_label = "class"
  )
  oof_score <- numeric(nrow(x))
  for (fold in seq_len(params$inner_folds)) {
    analysis_rows <- which(fold_id != fold)
    assessment_rows <- which(fold_id == fold)
    predictor <- sb_fit_cancerclass_native_predictor(
      x[analysis_rows, , drop = FALSE],
      outcome[analysis_rows],
      params
    )
    oof_score[assessment_rows] <- sb_cancerclass_score_from_predictor(
      predictor,
      x[assessment_rows, , drop = FALSE],
      names(x),
      params$distance
    )
  }

  calibration_warnings <- character()
  calibration <- withCallingHandlers(
    stats::glm.fit(
      x = cbind("(Intercept)" = 1, native_z_score = oof_score),
      y = outcome,
      family = stats::binomial()
    ),
    warning = function(warning) {
      calibration_warnings <<- c(
        calibration_warnings,
        conditionMessage(warning)
      )
      invokeRestart("muffleWarning")
    }
  )
  coefficients <- unname(as.numeric(calibration$coefficients))
  valid_calibration <- isTRUE(calibration$converged) &&
    length(coefficients) == 2L && !anyNA(coefficients) &&
    all(is.finite(coefficients))
  if (!valid_calibration) {
    sb_abort(
      "The inner out-of-fold cancerclass probability calibration failed.",
      "stormbreaker_cancerclass_calibration_failure"
    )
  }
  names(coefficients) <- c("intercept", "native_z_score")
  final_predictor <- sb_fit_cancerclass_native_predictor(x, outcome, params)

  structure(
    list(
      features = names(x),
      params = params,
      predictor = final_predictor,
      calibration_coefficients = coefficients,
      calibration_fold_id = fold_id,
      calibration_oof_score = oof_score,
      calibration_warnings = unique(calibration_warnings)
    ),
    class = "stormbreaker_cancerclass_binary"
  )
}

sb_cancerclass_raw_score <- function(model, x) {
  if (!inherits(model, "stormbreaker_cancerclass_binary")) {
    sb_abort(
      "The supplied model is not a fitted fixed cancerclass learner.",
      "stormbreaker_bad_cancerclass_fit"
    )
  }
  if (!sb_cancerclass_available()) {
    sb_abort(
      paste0(
        "Packages `cancerclass` and `Biobase` are required for fixed ",
        "cancerclass prediction."
      ),
      "stormbreaker_missing_cancerclass"
    )
  }
  sb_cancerclass_score_from_predictor(
    model$predictor,
    x,
    model$features,
    model$params$distance
  )
}

sb_predict_cancerclass_binary <- function(model, x) {
  score <- sb_cancerclass_raw_score(model, x)
  probability <- stats::plogis(
    model$calibration_coefficients[[1L]] +
      model$calibration_coefficients[[2L]] * score
  )
  if (length(probability) != length(score) || anyNA(probability) ||
    any(!is.finite(probability)) ||
    any(probability < 0 | probability > 1)) {
    sb_abort(
      "The fixed cancerclass learner produced invalid probabilities.",
      "stormbreaker_invalid_cancerclass_probability"
    )
  }
  unname(probability)
}

sb_extract_cancerclass_calibration_coefficients <- function(model) {
  data.frame(
    term = c("(Intercept)", "native_z_score"),
    coefficient = unname(model$calibration_coefficients),
    stringsAsFactors = FALSE
  )
}

sb_extract_cancerclass_binary_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_cancerclass_binary")) {
    sb_bad_learner_state("The fitted cancerclass learner state is invalid.")
  }
  structure(
    list(
      learner = "cancerclass_binary",
      features = model$features,
      params = model$params,
      predictor = model$predictor,
      calibration_coefficients = model$calibration_coefficients,
      score_reference = "cancerclass_native_welch_correlation_z",
      probability_formula = "plogis(intercept + slope * native_z_score)",
      training_data_retention = "none"
    ),
    class = c(
      "stormbreaker_cancerclass_binary_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_cancerclass_binary_learner_state <- function(state, x) {
  x <- sb_validate_learner_state_x(state, x, "cancerclass_binary")
  valid <- is.matrix(state[["predictor"]]) &&
    is.numeric(state[["calibration_coefficients"]]) &&
    length(state[["calibration_coefficients"]]) == 2L &&
    !anyNA(state[["calibration_coefficients"]]) &&
    all(is.finite(state[["calibration_coefficients"]]))
  if (!valid) {
    sb_bad_learner_state("The saved cancerclass learner state is invalid.")
  }
  model <- structure(
    list(
      features = state[["features"]],
      params = state[["params"]],
      predictor = state[["predictor"]],
      calibration_coefficients = state[["calibration_coefficients"]]
    ),
    class = "stormbreaker_cancerclass_binary"
  )
  list(probability = sb_predict_cancerclass_binary(model, x))
}

sb_cancerclass_binary_learner_adapter <- function() {
  params <- sb_cancerclass_binary_params()
  list(
    id = "cancerclass_binary",
    params = params,
    fit = function(x, outcome) {
      sb_fit_cancerclass_binary(x, outcome, params)
    },
    predict = sb_predict_cancerclass_binary,
    extract_coefficients = sb_extract_cancerclass_calibration_coefficients,
    extract_state = sb_extract_cancerclass_binary_state,
    predict_state = sb_predict_cancerclass_binary_learner_state
  )
}
