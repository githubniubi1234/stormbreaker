sb_validate_naive_bayes_var_smoothing <- function(var_smoothing) {
  valid <- is.numeric(var_smoothing) && length(var_smoothing) == 1L &&
    !is.na(var_smoothing) && is.finite(var_smoothing) &&
    var_smoothing >= 0
  if (!valid) {
    sb_abort(
      "`var_smoothing` must be one finite, non-negative number.",
      "stormbreaker_bad_naive_bayes_var_smoothing"
    )
  }
  as.numeric(var_smoothing)
}

sb_population_column_variance <- function(x) {
  centered <- sweep(x, 2L, colMeans(x), "-")
  colMeans(centered^2)
}

sb_fit_gaussian_naive_bayes <- function(x, outcome, var_smoothing) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      "A Gaussian Naive Bayes analysis fold must contain both classes.",
      "stormbreaker_naive_bayes_single_class"
    )
  }
  var_smoothing <- sb_validate_naive_bayes_var_smoothing(var_smoothing)
  matrix_x <- as.matrix(x)
  global_variance <- sb_population_column_variance(matrix_x)
  epsilon <- unname(var_smoothing * max(global_variance))
  classes <- c(0L, 1L)
  counts <- vapply(classes, function(class) {
    sum(outcome == class)
  }, integer(1))
  means <- do.call(rbind, lapply(classes, function(class) {
    colMeans(matrix_x[outcome == class, , drop = FALSE])
  }))
  variances <- do.call(rbind, lapply(seq_along(classes), function(index) {
    class_x <- matrix_x[outcome == classes[[index]], , drop = FALSE]
    sb_population_column_variance(class_x) + epsilon
  }))
  rownames(means) <- rownames(variances) <- as.character(classes)
  colnames(means) <- colnames(variances) <- names(x)
  valid_variance <- is.finite(epsilon) &&
    all(is.finite(variances)) && all(variances > 0)
  if (!valid_variance) {
    sb_abort(
      paste0(
        "Gaussian Naive Bayes did not produce one positive, finite ",
        "variance per class and predictor."
      ),
      "stormbreaker_invalid_naive_bayes_variance"
    )
  }
  priors <- counts / length(outcome)
  names(counts) <- names(priors) <- as.character(classes)

  structure(
    list(
      features = names(x),
      classes = classes,
      var_smoothing = var_smoothing,
      epsilon = epsilon,
      counts = counts,
      priors = priors,
      means = means,
      variances = variances
    ),
    class = "stormbreaker_gaussian_naive_bayes"
  )
}

sb_predict_gaussian_naive_bayes <- function(model, x) {
  if (!inherits(model, "stormbreaker_gaussian_naive_bayes")) {
    sb_abort(
      "The supplied model is not a fitted Gaussian Naive Bayes learner.",
      "stormbreaker_bad_naive_bayes_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the Gaussian Naive Bayes schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  matrix_x <- as.matrix(x)
  class_names <- as.character(model$classes)
  log_joint <- vapply(class_names, function(class) {
    centered <- sweep(matrix_x, 2L, model$means[class, ], "-")
    scaled_squared_distance <- sweep(
      centered^2,
      2L,
      model$variances[class, ],
      "/"
    )
    log(model$priors[[class]]) - 0.5 * (
      sum(log(2 * pi * model$variances[class, ])) +
        rowSums(scaled_squared_distance)
    )
  }, numeric(nrow(matrix_x)))
  log_joint <- matrix(
    log_joint,
    nrow = nrow(matrix_x),
    ncol = length(class_names),
    dimnames = list(NULL, class_names)
  )
  probability <- unname(stats::plogis(
    log_joint[, "1"] - log_joint[, "0"]
  ))
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      "Gaussian Naive Bayes produced invalid class-one probabilities.",
      "stormbreaker_invalid_naive_bayes_prediction"
    )
  }
  probability
}

sb_extract_gaussian_naive_bayes_parameters <- function(model) {
  if (!inherits(model, "stormbreaker_gaussian_naive_bayes")) {
    sb_abort(
      "The supplied model is not a fitted Gaussian Naive Bayes learner.",
      "stormbreaker_bad_naive_bayes_fit"
    )
  }
  values <- do.call(rbind, lapply(model$classes, function(class) {
    class_name <- as.character(class)
    data.frame(
      feature = model$features,
      class = as.integer(class),
      mean = unname(model$means[class_name, ]),
      variance = unname(model$variances[class_name, ]),
      standard_deviation = unname(sqrt(model$variances[class_name, ])),
      stringsAsFactors = FALSE
    )
  }))
  rownames(values) <- NULL
  list(
    var_smoothing = model$var_smoothing,
    epsilon = model$epsilon,
    priors = data.frame(
      class = model$classes,
      n = as.integer(model$counts[as.character(model$classes)]),
      prior = unname(model$priors[as.character(model$classes)]),
      stringsAsFactors = FALSE
    ),
    values = values
  )
}

sb_gaussian_naive_bayes_learner_adapter <- function(var_smoothing = 1e-9) {
  var_smoothing <- sb_validate_naive_bayes_var_smoothing(var_smoothing)
  params <- list(
    var_smoothing = var_smoothing,
    variance_estimator = "population",
    prior = "empirical"
  )
  list(
    id = "gaussian_naive_bayes",
    params = params,
    fit = function(x, outcome) {
      sb_fit_gaussian_naive_bayes(x, outcome, var_smoothing)
    },
    predict = sb_predict_gaussian_naive_bayes,
    extract_parameters = sb_extract_gaussian_naive_bayes_parameters,
    extract_state = function(model, context = list()) {
      sb_extract_naive_bayes_learner_state(model)
    },
    predict_state = sb_predict_naive_bayes_learner_state
  )
}
