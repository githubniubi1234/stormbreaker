sb_empty_binary_moderate_calibration_knots <- function() {
  data.frame(
    role = character(),
    quantile_probability = numeric(),
    logit_probability = numeric(),
    probability = numeric(),
    stringsAsFactors = FALSE
  )
}

sb_empty_binary_moderate_calibration_patients <- function() {
  data.frame(
    row_instance = integer(),
    patient_id = character(),
    fold = integer(),
    outcome = integer(),
    probability = numeric(),
    logit_probability = numeric(),
    calibrated_probability = numeric(),
    absolute_error = numeric(),
    stringsAsFactors = FALSE
  )
}

sb_empty_binary_moderate_calibration_curve <- function() {
  data.frame(
    probability = numeric(),
    logit_probability = numeric(),
    calibrated_probability = numeric(),
    n_row_instances = integer(),
    stringsAsFactors = FALSE
  )
}

sb_empty_binary_moderate_calibration_folds <- function() {
  data.frame(
    fold = integer(),
    status = character(),
    n = integer(),
    n_scored = integer(),
    n_positive = integer(),
    n_negative = integer(),
    n_boundary_probability = integer(),
    error_class = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  )
}

sb_binary_moderate_calibration_empty_fit <- function() {
  terms <- c("(Intercept)", "calibration_ns1", "calibration_ns2")
  coefficients <- data.frame(
    term = terms,
    estimate = rep(NA_real_, length(terms)),
    stringsAsFactors = FALSE
  )
  names(coefficients$estimate) <- terms
  list(
    status = "not_fitted",
    method = "binomial_glm_three_knot_restricted_natural_cubic_spline",
    family = "binomial",
    link = "logit",
    basis = "splines_ns_intercept_false",
    terms = terms,
    coefficients = coefficients,
    rank = NA_integer_,
    converged = NA,
    boundary = NA,
    iterations = NA_integer_,
    deviance = NA_real_,
    null_deviance = NA_real_,
    separation_method = "bounded_linear_program_signed_margin",
    separation_status = "not_checked",
    separation_detected = NA,
    separation_solver_status = NA_integer_,
    separation_objective = NA_real_,
    separation_tolerance = sqrt(.Machine$double.eps),
    separation_minimum_signed_margin = NA_real_,
    separation_maximum_signed_margin = NA_real_,
    warnings = character(),
    error_class = NA_character_,
    error_message = NA_character_
  )
}

sb_binary_moderate_calibration_empty_result <- function(scope, n = 0L) {
  safe_scope <- if (is.character(scope) && length(scope) >= 1L) {
    as.character(scope[[1L]])
  } else {
    NA_character_
  }
  list(
    metric = "binary_calibration_curve",
    calibration_level = "moderate",
    estimand = "smoothed_observed_class_1_probability",
    method = "logistic_three_knot_restricted_natural_cubic_spline",
    scope = safe_scope,
    purpose = "descriptive_evaluation_only",
    prediction_scale = "class_1_probability",
    prediction_transform = "logit",
    probability_clipping = "none",
    positive_class = 1L,
    spline = list(
      type = "restricted_natural_cubic_spline",
      n_knots = 3L,
      knot_probabilities = c(0.10, 0.50, 0.90),
      quantile_type = 7L,
      outcome_tuning = "none",
      basis = "splines_ns",
      intercept = FALSE,
      boundary_knots = c(NA_real_, NA_real_),
      interior_knots = NA_real_
    ),
    status = "incomplete",
    input_status = "not_checked",
    ici = NA_real_,
    e50 = NA_real_,
    e90 = NA_real_,
    n = as.integer(n),
    n_scored = 0L,
    n_positive = 0L,
    n_negative = 0L,
    n_boundary_probability = 0L,
    n_unique_probability = 0L,
    n_folds = 0L,
    n_evaluable_folds = 0L,
    knots = sb_empty_binary_moderate_calibration_knots(),
    fit = sb_binary_moderate_calibration_empty_fit(),
    curve = sb_empty_binary_moderate_calibration_curve(),
    patient_calibration = sb_empty_binary_moderate_calibration_patients(),
    per_fold = sb_empty_binary_moderate_calibration_folds(),
    error_class = NA_character_,
    error_message = NA_character_
  )
}

sb_binary_moderate_calibration_fail <- function(
    result,
    status,
    input_status,
    error_class,
    error_message) {
  result$status <- status
  result$input_status <- input_status
  result$error_class <- error_class
  result$error_message <- error_message
  result
}

sb_binary_moderate_calibration_error_metric <- function(error, scope, n = 0L) {
  result <- sb_binary_moderate_calibration_empty_result(scope, n)
  result$status <- "metric_failed"
  result$input_status <- "unexpected_error"
  result$error_class <- class(error)[[1L]]
  result$error_message <- conditionMessage(error)
  result
}

sb_binary_moderate_calibration_fold_audit <- function(patients, status) {
  if (all(is.na(patients$fold))) {
    return(sb_empty_binary_moderate_calibration_folds())
  }
  fold_ids <- sort(unique(patients$fold))
  rows <- lapply(fold_ids, function(fold_id) {
    index <- patients$fold == fold_id
    successful <- all(status[index] == "success")
    data.frame(
      fold = as.integer(fold_id),
      status = if (successful) "success" else "prediction_failed",
      n = as.integer(sum(index)),
      n_scored = if (successful) as.integer(sum(index)) else 0L,
      n_positive = as.integer(sum(patients$outcome[index] == 1L)),
      n_negative = as.integer(sum(patients$outcome[index] == 0L)),
      n_boundary_probability = as.integer(sum(
        patients$probability[index] %in% c(0, 1),
        na.rm = TRUE
      )),
      error_class = if (successful) {
        NA_character_
      } else {
        "stormbreaker_invalid_binary_moderate_calibration_prediction"
      },
      error_message = if (successful) {
        NA_character_
      } else {
        "One or more saved predictions in this fold failed."
      },
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

sb_binary_moderate_calibration_knots <- function(logit_probability) {
  values <- unname(stats::quantile(
    logit_probability,
    probs = c(0.10, 0.50, 0.90),
    names = FALSE,
    type = 7
  ))
  data.frame(
    role = c("lower_boundary", "interior", "upper_boundary"),
    quantile_probability = c(0.10, 0.50, 0.90),
    logit_probability = values,
    probability = stats::plogis(values),
    stringsAsFactors = FALSE
  )
}

sb_binary_moderate_calibration_basis <- function(logit_probability, knots) {
  basis <- splines::ns(
    logit_probability,
    knots = knots[[2L]],
    Boundary.knots = knots[c(1L, 3L)],
    intercept = FALSE
  )
  colnames(basis) <- c("calibration_ns1", "calibration_ns2")
  basis
}

sb_binary_moderate_calibration_separation <- function(design, outcome) {
  column_scale <- apply(abs(design), 2L, max)
  normalized_design <- sweep(design, 2L, column_scale, "/")
  outcome_sign <- ifelse(outcome == 1L, 1, -1)
  signed_design <- sweep(normalized_design, 1L, outcome_sign, "*")
  n_terms <- ncol(signed_design)
  objective <- c(colSums(signed_design), -colSums(signed_design))
  constraints <- rbind(
    cbind(signed_design, -signed_design),
    rep(1, 2L * n_terms)
  )
  solver <- tryCatch(
    lpSolve::lp(
      direction = "max",
      objective.in = objective,
      const.mat = constraints,
      const.dir = c(rep(">=", nrow(signed_design)), "<="),
      const.rhs = c(rep(0, nrow(signed_design)), 1)
    ),
    error = function(error) error
  )
  tolerance <- sqrt(.Machine$double.eps)
  if (inherits(solver, "condition")) {
    return(list(
      status = "failed",
      detected = NA,
      solver_status = NA_integer_,
      objective = NA_real_,
      tolerance = tolerance,
      minimum_signed_margin = NA_real_,
      maximum_signed_margin = NA_real_,
      error_class = class(solver)[[1L]],
      error_message = conditionMessage(solver)
    ))
  }
  if (!identical(as.integer(solver$status), 0L)) {
    return(list(
      status = "failed",
      detected = NA,
      solver_status = as.integer(solver$status),
      objective = unname(solver$objval),
      tolerance = tolerance,
      minimum_signed_margin = NA_real_,
      maximum_signed_margin = NA_real_,
      error_class =
        "stormbreaker_binary_moderate_calibration_separation_check_failed",
      error_message = paste0(
        "The bounded separation linear program returned solver status ",
        as.integer(solver$status), "."
      )
    ))
  }

  positive <- solver$solution[seq_len(n_terms)]
  negative <- solver$solution[n_terms + seq_len(n_terms)]
  separator <- positive - negative
  signed_margin <- drop(signed_design %*% separator)
  objective_value <- unname(solver$objval)
  minimum_margin <- unname(min(signed_margin))
  maximum_margin <- unname(max(signed_margin))
  feasible <- minimum_margin >= -tolerance
  detected <- objective_value > tolerance && feasible &&
    maximum_margin > tolerance
  if (objective_value > tolerance && !feasible) {
    return(list(
      status = "failed",
      detected = NA,
      solver_status = as.integer(solver$status),
      objective = objective_value,
      tolerance = tolerance,
      minimum_signed_margin = minimum_margin,
      maximum_signed_margin = maximum_margin,
      error_class =
        "stormbreaker_binary_moderate_calibration_separation_check_failed",
      error_message = paste0(
        "The bounded separation linear program returned a solution outside ",
        "its signed-margin tolerance."
      )
    ))
  }
  list(
    status = if (detected) "detected" else "not_detected",
    detected = detected,
    solver_status = as.integer(solver$status),
    objective = objective_value,
    tolerance = tolerance,
    minimum_signed_margin = minimum_margin,
    maximum_signed_margin = maximum_margin,
    error_class = NA_character_,
    error_message = NA_character_
  )
}

sb_binary_moderate_calibration_curve_table <- function(patients) {
  probability <- sort(unique(patients$probability))
  first <- match(probability, patients$probability)
  data.frame(
    probability = probability,
    logit_probability = patients$logit_probability[first],
    calibrated_probability = patients$calibrated_probability[first],
    n_row_instances = as.integer(tabulate(
      match(patients$probability, probability),
      nbins = length(probability)
    )),
    stringsAsFactors = FALSE
  )
}

sb_binary_moderate_calibration_metric <- function(scored, scope) {
  n <- if (is.data.frame(scored)) nrow(scored) else 0L
  result <- sb_binary_moderate_calibration_empty_result(scope, n)
  required <- c("patient_id", "outcome", "probability", "status")
  required_once <- is.data.frame(scored) && n >= 1L &&
    all(vapply(required, function(name) {
      sum(names(scored) == name) == 1L
    }, logical(1)))
  fold_count <- if (is.data.frame(scored)) sum(names(scored) == "fold") else 0L
  valid_scope <- is.character(scope) && length(scope) == 1L &&
    !is.na(scope) && nzchar(scope)
  if (!required_once || fold_count > 1L || !valid_scope) {
    return(sb_binary_moderate_calibration_fail(
      result,
      "incomplete",
      "input_failed",
      "stormbreaker_bad_binary_moderate_calibration_input",
      paste0(
        "Binary moderate calibration requires one non-empty scored table, ",
        "one valid scope, required columns, and at most one fold column."
      )
    ))
  }

  patient_id <- scored$patient_id
  outcome <- scored$outcome
  probability <- scored$probability
  prediction_status <- scored$status
  fold <- if (fold_count == 1L) scored$fold else rep(NA_integer_, n)
  valid_patient_id <- is.character(patient_id) && length(patient_id) == n &&
    !anyNA(patient_id) && all(nzchar(trimws(patient_id)))
  valid_outcome <- (is.numeric(outcome) || is.logical(outcome)) &&
    length(outcome) == n && !anyNA(outcome) &&
    (!is.numeric(outcome) || all(is.finite(outcome))) &&
    all(outcome %in% c(0, 1))
  valid_probability_type <- is.numeric(probability) &&
    length(probability) == n
  valid_status <- is.character(prediction_status) &&
    length(prediction_status) == n && !anyNA(prediction_status)
  valid_fold <- fold_count == 0L || (
    is.numeric(fold) && length(fold) == n && !anyNA(fold) &&
      all(is.finite(fold)) && all(fold == trunc(fold))
  )
  if (!valid_patient_id || !valid_outcome || !valid_probability_type ||
      !valid_status || !valid_fold) {
    return(sb_binary_moderate_calibration_fail(
      result,
      "incomplete",
      "input_failed",
      "stormbreaker_bad_binary_moderate_calibration_input",
      paste0(
        "Binary moderate calibration requires valid aligned patient IDs, ",
        "0/1 outcomes, probabilities, statuses, and optional integer folds."
      )
    ))
  }

  valid_probability <- !anyNA(probability) && all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  boundary <- !is.na(probability) & is.finite(probability) &
    probability %in% c(0, 1)
  interior <- !is.na(probability) & is.finite(probability) &
    probability > 0 & probability < 1
  logit_probability <- rep(NA_real_, n)
  logit_probability[interior] <- stats::qlogis(probability[interior])
  patients <- data.frame(
    row_instance = seq_len(n),
    patient_id = patient_id,
    fold = if (fold_count == 1L) as.integer(fold) else rep(NA_integer_, n),
    outcome = as.integer(outcome),
    probability = unname(probability),
    logit_probability = logit_probability,
    calibrated_probability = rep(NA_real_, n),
    absolute_error = rep(NA_real_, n),
    stringsAsFactors = FALSE
  )
  result$patient_calibration <- patients
  result$n_positive <- as.integer(sum(patients$outcome == 1L))
  result$n_negative <- as.integer(sum(patients$outcome == 0L))
  result$n_boundary_probability <- as.integer(sum(boundary))
  result$n_unique_probability <- as.integer(length(unique(
    probability[!is.na(probability) & is.finite(probability)]
  )))
  result$per_fold <- sb_binary_moderate_calibration_fold_audit(
    patients,
    prediction_status
  )
  result$n_folds <- as.integer(nrow(result$per_fold))
  result$n_evaluable_folds <- as.integer(sum(
    result$per_fold$status == "success"
  ))

  if (any(prediction_status != "success")) {
    return(sb_binary_moderate_calibration_fail(
      result,
      "incomplete",
      "prediction_failed",
      "stormbreaker_invalid_binary_moderate_calibration_prediction",
      "Complete successful probabilities are required for calibration."
    ))
  }
  if (!valid_probability) {
    return(sb_binary_moderate_calibration_fail(
      result,
      "incomplete",
      "input_failed",
      "stormbreaker_bad_binary_moderate_calibration_input",
      "Successful binary predictions must be finite probabilities in [0, 1]."
    ))
  }
  result$n_scored <- as.integer(n)
  if (result$n_positive == 0L || result$n_negative == 0L) {
    return(sb_binary_moderate_calibration_fail(
      result,
      "incomplete",
      "not_evaluable",
      "stormbreaker_binary_moderate_calibration_one_class",
      "Both outcome classes are required for binary moderate calibration."
    ))
  }
  if (result$n_boundary_probability > 0L) {
    return(sb_binary_moderate_calibration_fail(
      result,
      "incomplete",
      "not_evaluable",
      "stormbreaker_binary_moderate_calibration_boundary_probability",
      paste0(
        "Exact probabilities zero or one cannot enter a finite-logit ",
        "calibration fit without clipping, which is not performed."
      )
    ))
  }
  if (any(!is.finite(logit_probability))) {
    return(sb_binary_moderate_calibration_fail(
      result,
      "incomplete",
      "not_evaluable",
      "stormbreaker_binary_moderate_calibration_nonfinite_logit",
      "The unmodified probabilities did not produce finite logits."
    ))
  }

  knots <- sb_binary_moderate_calibration_knots(logit_probability)
  result$knots <- knots
  knot_values <- knots$logit_probability
  result$spline$boundary_knots <- unname(knot_values[c(1L, 3L)])
  result$spline$interior_knots <- unname(knot_values[[2L]])
  if (anyNA(knot_values) || any(!is.finite(knot_values)) ||
      any(diff(knot_values) <= 0)) {
    return(sb_binary_moderate_calibration_fail(
      result,
      "incomplete",
      "not_evaluable",
      "stormbreaker_binary_moderate_calibration_duplicate_knots",
      paste0(
        "The fixed 10th, 50th, and 90th percentile logit knots are not ",
        "finite and strictly increasing."
      )
    ))
  }

  basis <- tryCatch(
    sb_binary_moderate_calibration_basis(logit_probability, knot_values),
    error = function(error) error
  )
  if (inherits(basis, "condition")) {
    return(sb_binary_moderate_calibration_fail(
      result,
      "fit_failed",
      "success",
      class(basis)[[1L]],
      conditionMessage(basis)
    ))
  }
  design <- cbind("(Intercept)" = 1, basis)
  captured_warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      stats::glm.fit(
        x = design,
        y = patients$outcome,
        family = stats::binomial(),
        control = stats::glm.control()
      ),
      warning = function(warning) {
        captured_warnings <<- c(captured_warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error) error
  )
  result$fit$warnings <- unique(captured_warnings)
  if (inherits(fit, "condition")) {
    result$fit$status <- "failed"
    result$fit$error_class <- class(fit)[[1L]]
    result$fit$error_message <- conditionMessage(fit)
    return(sb_binary_moderate_calibration_fail(
      result,
      "fit_failed",
      "success",
      result$fit$error_class,
      result$fit$error_message
    ))
  }

  terms <- result$fit$terms
  coefficient_values <- unname(fit$coefficients[terms])
  names(coefficient_values) <- terms
  result$fit$coefficients$estimate <- coefficient_values
  result$fit$rank <- as.integer(fit$rank)
  result$fit$converged <- isTRUE(fit$converged)
  result$fit$boundary <- isTRUE(fit$boundary)
  result$fit$iterations <- as.integer(fit$iter)
  result$fit$deviance <- unname(fit$deviance)
  result$fit$null_deviance <- unname(fit$null.deviance)
  # Saturated fitted probabilities alone do not establish separation. Keep
  # this known warning and still require every fit and separation check below.
  unexpected_warnings <- setdiff(
    captured_warnings,
    "glm.fit: fitted probabilities numerically 0 or 1 occurred"
  )
  if (length(unexpected_warnings) > 0L) {
    result$fit$status <- "failed"
    result$fit$error_class <-
      "stormbreaker_binary_moderate_calibration_glm_warning"
    result$fit$error_message <- paste(unique(captured_warnings), collapse = " ")
    return(sb_binary_moderate_calibration_fail(
      result,
      "fit_failed",
      "success",
      result$fit$error_class,
      result$fit$error_message
    ))
  }
  if (!result$fit$converged) {
    result$fit$status <- "failed"
    result$fit$error_class <-
      "stormbreaker_binary_moderate_calibration_nonconvergence"
    result$fit$error_message <- "The binary calibration GLM did not converge."
    return(sb_binary_moderate_calibration_fail(
      result,
      "fit_failed",
      "success",
      result$fit$error_class,
      result$fit$error_message
    ))
  }
  if (result$fit$boundary) {
    result$fit$status <- "failed"
    result$fit$error_class <-
      "stormbreaker_binary_moderate_calibration_boundary_fit"
    result$fit$error_message <- "The binary calibration GLM stopped at a boundary."
    return(sb_binary_moderate_calibration_fail(
      result,
      "fit_failed",
      "success",
      result$fit$error_class,
      result$fit$error_message
    ))
  }
  valid_coefficients <- identical(result$fit$rank, 3L) &&
    length(coefficient_values) == 3L && !anyNA(coefficient_values) &&
    all(is.finite(coefficient_values))
  if (!valid_coefficients) {
    result$fit$status <- "failed"
    result$fit$error_class <-
      "stormbreaker_binary_moderate_calibration_singular_fit"
    result$fit$error_message <- paste0(
      "The binary calibration GLM returned missing, singular, or non-finite ",
      "basis coefficients."
    )
    return(sb_binary_moderate_calibration_fail(
      result,
      "fit_failed",
      "success",
      result$fit$error_class,
      result$fit$error_message
    ))
  }

  # IRLS can report convergence with finite coefficients even when separation
  # makes the logistic maximum-likelihood estimate nonexistent. This post-fit
  # linear program checks that data-dependent failure without using a
  # coefficient cutoff or replacing the requested calibration model.
  separation <- sb_binary_moderate_calibration_separation(
    design,
    patients$outcome
  )
  result$fit$separation_status <- separation$status
  result$fit$separation_detected <- separation$detected
  result$fit$separation_solver_status <- separation$solver_status
  result$fit$separation_objective <- separation$objective
  result$fit$separation_tolerance <- separation$tolerance
  result$fit$separation_minimum_signed_margin <-
    separation$minimum_signed_margin
  result$fit$separation_maximum_signed_margin <-
    separation$maximum_signed_margin
  if (identical(separation$status, "failed")) {
    result$fit$status <- "failed"
    result$fit$error_class <- separation$error_class
    result$fit$error_message <- separation$error_message
    return(sb_binary_moderate_calibration_fail(
      result,
      "fit_failed",
      "success",
      result$fit$error_class,
      result$fit$error_message
    ))
  }
  if (isTRUE(separation$detected)) {
    result$fit$status <- "failed"
    result$fit$error_class <-
      "stormbreaker_binary_moderate_calibration_separation"
    result$fit$error_message <- paste0(
      "The binary calibration design has complete or quasi-complete outcome ",
      "separation, so no finite logistic maximum-likelihood estimate is ",
      "reported."
    )
    return(sb_binary_moderate_calibration_fail(
      result,
      "fit_failed",
      "success",
      result$fit$error_class,
      result$fit$error_message
    ))
  }

  calibrated_probability <- stats::plogis(drop(
    design %*% coefficient_values
  ))
  valid_calibrated <- !anyNA(calibrated_probability) &&
    all(is.finite(calibrated_probability)) &&
    all(calibrated_probability >= 0 & calibrated_probability <= 1)
  if (!valid_calibrated) {
    result$fit$status <- "failed"
    result$fit$error_class <-
      "stormbreaker_binary_moderate_calibration_invalid_probability"
    result$fit$error_message <- paste0(
      "The binary calibration GLM produced invalid probabilities."
    )
    return(sb_binary_moderate_calibration_fail(
      result,
      "fit_failed",
      "success",
      result$fit$error_class,
      result$fit$error_message
    ))
  }

  absolute_error <- abs(calibrated_probability - probability)
  result$patient_calibration$calibrated_probability <- calibrated_probability
  result$patient_calibration$absolute_error <- absolute_error
  result$curve <- sb_binary_moderate_calibration_curve_table(
    result$patient_calibration
  )
  result$ici <- unname(mean(absolute_error))
  result$e50 <- unname(stats::quantile(
    absolute_error,
    probs = 0.50,
    names = FALSE,
    type = 7
  ))
  result$e90 <- unname(stats::quantile(
    absolute_error,
    probs = 0.90,
    names = FALSE,
    type = 7
  ))
  result$fit$status <- "success"
  result$status <- "success"
  result$input_status <- "success"
  result
}
