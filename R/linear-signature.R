sb_signature_optional_text <- function(value, argument) {
  if (is.null(value)) {
    return(NULL)
  }
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(trimws(value))) {
    sb_abort(
      paste0(argument, " must be NULL or one non-empty string."),
      "stormbreaker_bad_signature_metadata"
    )
  }
  value
}

#' Define a locked linear prediction signature
#'
#' sb_linear_signature() records one explicit, fixed linear formula without
#' fitting, feature selection, tuning, or outcome access. Coefficient names are
#' the exact predictor column names required later by [sb_predict_signature()].
#' Their order is retained.
#'
#' For a binary signature, the only supported interpretation is a published
#' logistic formula: probability(class 1) equals
#' plogis(intercept + sum(x_j * coefficient_j)).
#' The intercept must be supplied explicitly, including when its stated value
#' is zero; an omitted binary intercept cannot define an absolute probability.
#'
#' For a survival signature, the linear score is returned as a risk value, with
#' larger values meaning higher predicted risk. If and only if a prediction
#' horizon, its time unit, and the baseline survival at zero linear score are
#' all supplied, survival probability is calculated as
#' S(h | x) = S0(h) ^ exp(intercept + sum(x_j * coefficient_j)).
#'
#' The function deliberately does not infer gene identifiers, map probes,
#' transform or standardize expression, impute values, orient coefficients,
#' recover an omitted intercept or baseline survival, or verify a citation.
#' A caller-reviewed one-to-one map can be applied separately with
#' [sb_apply_feature_map()]; signature prediction itself remains exact-name
#' only.
#' source and input_description are retained as free-text audit metadata only.
#' The caller remains responsible for establishing that predictor columns use
#' the exact identifiers and scale required by the original formula. A horizon
#' probability is not available when the publication gives only a relative-risk
#' score.
#'
#' This object is a reproducible formula carrier, not a literature-evidence
#' registry entry and not proof that a signature was independently validated.
#' Like ordinary R objects it is mutable; "locked" means that prediction does
#' not learn from evaluation data, not that the object is cryptographically
#' immutable.
#'
#' @param signature_id One non-empty identifier for this signature.
#' @param task_type Exactly "binary" or "survival".
#' @param coefficients A non-empty finite numeric vector named by unique,
#'   non-empty predictor column names.
#' @param intercept For a binary signature, one explicitly supplied finite
#'   numeric intercept. For a survival signature, one finite numeric intercept;
#'   `NULL` uses the conventional zero-reference linear predictor.
#' @param source Optional free-text citation or source identifier. It is stored
#'   but not verified.
#' @param input_description Optional free-text description of the required
#'   feature identifiers and preprocessing scale. It is stored but not applied.
#' @param prediction_horizon For survival only, an optional finite positive
#'   horizon. It must be supplied together with time_unit and
#'   baseline_survival.
#' @param time_unit For a survival horizon, exactly "day", "month", or "year".
#' @param baseline_survival For a survival horizon, the published baseline
#'   survival probability at zero linear score, between zero and one inclusive.
#'
#' @return A stormbreaker_linear_signature containing the exact formula,
#'   prediction semantics, optional source metadata, and required feature order.
#' @export
#'
#' @examples
#' signature <- sb_linear_signature(
#'   "example_binary_signature",
#'   "binary",
#'   coefficients = c(gene_a = 0.7, gene_b = -0.4),
#'   intercept = -0.2,
#'   input_description = "Values prepared on the formula's stated scale"
#' )
#' signature
sb_linear_signature <- function(
    signature_id,
    task_type,
    coefficients,
    intercept = NULL,
    source = NULL,
    input_description = NULL,
    prediction_horizon = NULL,
    time_unit = NULL,
    baseline_survival = NULL) {
  if (!is.character(signature_id) || length(signature_id) != 1L ||
      is.na(signature_id) || !nzchar(trimws(signature_id))) {
    sb_abort(
      "signature_id must be one non-empty string.",
      "stormbreaker_bad_signature_id"
    )
  }
  if (!is.character(task_type) || length(task_type) != 1L ||
      is.na(task_type) || !task_type %in% c("binary", "survival")) {
    sb_abort(
      "task_type must be exactly \"binary\" or \"survival\".",
      "stormbreaker_bad_signature_task"
    )
  }
  coefficient_names <- names(coefficients)
  valid_coefficients <- is.numeric(coefficients) &&
    is.null(dim(coefficients)) && length(coefficients) > 0L &&
    !anyNA(coefficients) && all(is.finite(coefficients)) &&
    !is.null(coefficient_names) && !anyNA(coefficient_names) &&
    all(nzchar(trimws(coefficient_names))) &&
    anyDuplicated(coefficient_names) == 0L
  if (!valid_coefficients) {
    sb_abort(
      paste0(
        "coefficients must be a non-empty finite numeric vector named by ",
        "unique, non-empty predictor columns."
      ),
      "stormbreaker_bad_signature_coefficients"
    )
  }
  if (is.null(intercept)) {
    if (identical(task_type, "binary")) {
      sb_abort(
        paste0(
          "A binary logistic signature requires an explicitly stated ",
          "intercept; use intercept = 0 only when zero is the stated value."
        ),
        "stormbreaker_missing_signature_intercept"
      )
    }
    intercept <- 0
  }
  if (!is.numeric(intercept) || length(intercept) != 1L ||
      is.na(intercept) || !is.finite(intercept)) {
    sb_abort(
      "intercept must be one finite number.",
      "stormbreaker_bad_signature_intercept"
    )
  }
  source <- sb_signature_optional_text(source, "source")
  input_description <- sb_signature_optional_text(
    input_description,
    "input_description"
  )

  horizon_fields <- c(
    !is.null(prediction_horizon),
    !is.null(time_unit),
    !is.null(baseline_survival)
  )
  if (identical(task_type, "binary") && any(horizon_fields)) {
    sb_abort(
      "Survival-horizon fields cannot be attached to a binary signature.",
      "stormbreaker_bad_signature_horizon"
    )
  }
  if (identical(task_type, "survival") &&
      any(horizon_fields) && !all(horizon_fields)) {
    sb_abort(
      paste0(
        "A survival horizon requires prediction_horizon, time_unit, ",
        "and baseline_survival together."
      ),
      "stormbreaker_bad_signature_horizon"
    )
  }
  if (all(horizon_fields)) {
    prediction_horizon <- sb_validate_endpoint_horizon(prediction_horizon)
    time_unit <- sb_validate_endpoint_time_unit(time_unit)
    valid_baseline <- is.numeric(baseline_survival) &&
      length(baseline_survival) == 1L && !is.na(baseline_survival) &&
      is.finite(baseline_survival) && baseline_survival >= 0 &&
      baseline_survival <= 1
    if (!valid_baseline) {
      sb_abort(
        "baseline_survival must be one finite number between 0 and 1.",
        "stormbreaker_bad_signature_baseline_survival"
      )
    }
    baseline_survival <- as.numeric(baseline_survival)
  }

  prediction <- if (identical(task_type, "binary")) {
    list(
      output = "probability",
      positive_class = 1L,
      link = "logit"
    )
  } else {
    list(
      output = "risk",
      direction = "higher_linear_score_means_higher_risk"
    )
  }
  if (all(horizon_fields)) {
    prediction <- list(
      output = c("risk", "survival_probability"),
      direction = "higher_linear_score_means_higher_risk",
      prediction_horizon = prediction_horizon,
      time_unit = time_unit,
      baseline_survival = baseline_survival,
      formula = "baseline_survival^exp(linear_score)"
    )
  }
  result <- list(
    signature_id = signature_id,
    pipeline_id = paste0("signature__", signature_id),
    task_type = task_type,
    learner = "locked_linear_signature",
    features = coefficient_names,
    coefficients = coefficients,
    intercept = as.numeric(intercept),
    prediction = prediction,
    source = source,
    input_description = input_description
  )
  if (all(horizon_fields)) {
    result <- append(
      result,
      list(
        prediction_horizon = prediction_horizon,
        time_unit = time_unit,
        baseline_survival = baseline_survival
      ),
      after = match("intercept", names(result))
    )
  }
  structure(result, class = "stormbreaker_linear_signature")
}

sb_validate_linear_signature <- function(signature) {
  if (!inherits(signature, "stormbreaker_linear_signature") ||
      !is.list(signature)) {
    sb_abort(
      "signature must be returned by sb_linear_signature().",
      "stormbreaker_bad_linear_signature"
    )
  }
  sb_linear_signature(
    signature_id = signature[["signature_id"]],
    task_type = signature[["task_type"]],
    coefficients = signature[["coefficients"]],
    intercept = signature[["intercept"]],
    source = signature[["source"]],
    input_description = signature[["input_description"]],
    prediction_horizon = signature[["prediction_horizon"]],
    time_unit = signature[["time_unit"]],
    baseline_survival = signature[["baseline_survival"]]
  )
}

#' @export
print.stormbreaker_linear_signature <- function(x, ...) {
  cat("<stormbreaker_linear_signature>\n")
  cat("Signature:", x$signature_id, "| Task:", x$task_type, "\n")
  cat("Required numeric features:", length(x$features), "\n")
  if ("prediction_horizon" %in% names(x)) {
    cat(
      "Prediction horizon:", format(x$prediction_horizon),
      "[", x$time_unit, "]\n"
    )
  }
  if (!is.null(x$source)) {
    cat("Source metadata:", x$source, "\n")
  }
  invisible(x)
}

sb_signature_safe_metadata <- function(signature) {
  signature_id <- if (is.list(signature) &&
      sb_external_scalar_character(signature[["signature_id"]])) {
    signature[["signature_id"]]
  } else {
    NA_character_
  }
  task_type <- if (is.list(signature) &&
      sb_external_scalar_character(signature[["task_type"]]) &&
      signature[["task_type"]] %in% c("binary", "survival")) {
    signature[["task_type"]]
  } else {
    NA_character_
  }
  list(
    signature_id = signature_id,
    pipeline_id = if (is.na(signature_id)) {
      NA_character_
    } else {
      paste0("signature__", signature_id)
    },
    task_type = task_type,
    learner = "locked_linear_signature",
    prediction_horizon = if (is.list(signature)) {
      signature[["prediction_horizon"]]
    } else {
      NULL
    },
    time_unit = if (is.list(signature)) signature[["time_unit"]] else NULL
  )
}

sb_signature_prediction_rows <- function(
    metadata,
    patient_ids,
    status,
    error_class,
    error_message,
    values = NULL) {
  n <- length(patient_ids)
  rows <- data.frame(
    patient_id = patient_ids,
    pipeline_id = rep(metadata$pipeline_id, n),
    learner = rep(metadata$learner, n),
    stringsAsFactors = FALSE
  )
  if (identical(metadata$task_type, "binary")) {
    rows$positive_class <- rep.int(1L, n)
    rows$probability <- if (is.null(values)) {
      rep(NA_real_, n)
    } else {
      values$probability
    }
  } else if (identical(metadata$task_type, "survival")) {
    rows$risk <- if (is.null(values)) rep(NA_real_, n) else values$risk
    if (!is.null(metadata$prediction_horizon)) {
      rows$prediction_horizon <- rep(metadata$prediction_horizon, n)
      rows$time_unit <- rep(metadata$time_unit, n)
      rows$survival_probability <- if (is.null(values)) {
        rep(NA_real_, n)
      } else {
        values$survival_probability
      }
    }
  }
  rows$status <- rep(status, n)
  rows$error_class <- rep(error_class, n)
  rows$error_message <- rep(error_message, n)
  rows
}

sb_signature_prediction_failure <- function(
    signature,
    id_col,
    patient_ids,
    stage,
    error) {
  metadata <- sb_signature_safe_metadata(signature)
  error_class <- class(error)[[1L]]
  error_message <- conditionMessage(error)
  result <- list(
    status = "failed",
    task_type = metadata$task_type,
    pipeline_id = metadata$pipeline_id,
    learner = metadata$learner,
    id_col = if (is.character(id_col) && length(id_col) == 1L &&
        !is.na(id_col)) {
      id_col
    } else {
      NA_character_
    },
    signature = signature,
    predictions = sb_signature_prediction_rows(
      metadata,
      patient_ids,
      "failed",
      error_class,
      error_message
    ),
    error_stage = stage,
    error_class = error_class,
    error_message = error_message
  )
  structure(result, class = "stormbreaker_final_prediction")
}

sb_select_signature_prediction_features <- function(
    new_data,
    id_col,
    features) {
  if (id_col %in% features) {
    sb_abort(
      "The prediction ID column cannot also be a signature feature.",
      "stormbreaker_overlapping_signature_roles"
    )
  }
  counts <- vapply(features, function(feature) {
    sum(names(new_data) == feature)
  }, integer(1L))
  if (any(counts == 0L)) {
    sb_abort(
      paste0(
        "Prediction data are missing signature features: ",
        paste(features[counts == 0L], collapse = ", "),
        "."
      ),
      "stormbreaker_missing_signature_features"
    )
  }
  if (any(counts > 1L)) {
    sb_abort(
      paste0(
        "Prediction data duplicate signature features: ",
        paste(features[counts > 1L], collapse = ", "),
        "."
      ),
      "stormbreaker_duplicate_signature_feature_columns"
    )
  }
  positions <- vapply(features, function(feature) {
    which(names(new_data) == feature)
  }, integer(1L))
  x <- new_data[positions]
  names(x) <- features
  sb_validate_predictor_frame(x, "new_data signature features")
}

sb_signature_linear_score <- function(signature, x) {
  score <- rep(signature$intercept, nrow(x))
  for (feature in signature$features) {
    score <- score + x[[feature]] * signature$coefficients[[feature]]
  }
  if (anyNA(score) || any(!is.finite(score))) {
    sb_abort(
      "The fixed signature formula produced a non-finite linear score.",
      "stormbreaker_nonfinite_signature_score"
    )
  }
  unname(score)
}

sb_signature_prediction_values <- function(signature, x) {
  score <- sb_signature_linear_score(signature, x)
  if (identical(signature$task_type, "binary")) {
    probability <- stats::plogis(score)
    return(list(probability = unname(probability)))
  }
  values <- list(risk = score)
  if ("prediction_horizon" %in% names(signature)) {
    baseline <- signature$baseline_survival
    probability <- if (identical(baseline, 0)) {
      rep(0, length(score))
    } else {
      baseline^exp(score)
    }
    valid <- !anyNA(probability) && all(is.finite(probability)) &&
      all(probability >= 0 & probability <= 1)
    if (!valid) {
      sb_abort(
        "The fixed signature formula produced an invalid survival probability.",
        "stormbreaker_bad_signature_probability"
      )
    }
    values$survival_probability <- unname(probability)
  }
  values
}

sb_signature_prediction_provenance <- function() {
  list(
    package = "stormbreaker",
    package_version = as.character(utils::packageVersion("stormbreaker")),
    r_version = R.version.string,
    method = "locked_linear_signature_without_fitting"
  )
}

#' Predict from a locked linear signature without outcomes
#'
#' sb_predict_signature() applies the exact formula stored by
#' [sb_linear_signature()] to a non-empty new-patient data frame. It reads only
#' the requested ID column and the signature's named numeric features. Unrelated
#' columns, including accidentally supplied endpoint columns, do not enter the
#' calculation. No coefficient, intercept, baseline survival, feature,
#' transformation, threshold, or model component is learned or changed.
#'
#' The result is a stormbreaker_final_prediction, so it can be passed directly
#' to [sb_evaluate_external()], [sb_compare_external()], and
#' [sb_bootstrap_external_comparison()]. Binary signatures return the
#' probability of class 1. Survival signatures always return a linear risk
#' score and return a horizon survival probability only when the signature
#' explicitly stored all required horizon information.
#'
#' Feature names are matched exactly and reordered to the coefficient order.
#' Missing, duplicated, non-numeric, missing-valued, or non-finite signature
#' inputs produce a stage-aware failed prediction. The function never drops a
#' feature, imputes a value, refits a formula, or substitutes another model.
#'
#' @param signature A stormbreaker_linear_signature returned by
#'   [sb_linear_signature()].
#' @param new_data A non-empty data frame containing one unique patient ID and
#'   every exactly named numeric signature feature.
#' @param id_col Name of the patient-ID column in new_data.
#'
#' @return A successful or failed stormbreaker_final_prediction compatible
#'   with the existing prediction-first external-evaluation workflow.
#' @export
#'
#' @examples
#' signature <- sb_linear_signature(
#'   "example_binary_signature",
#'   "binary",
#'   coefficients = c(gene_a = 0.7, gene_b = -0.4),
#'   intercept = -0.2
#' )
#' new_data <- data.frame(
#'   patient_id = c("p1", "p2"),
#'   gene_a = c(0.1, 0.5),
#'   gene_b = c(-0.2, 0.4)
#' )
#' sb_predict_signature(signature, new_data, "patient_id")
sb_predict_signature <- function(signature, new_data, id_col) {
  id_attempt <- sb_final_fit_attempt(
    sb_extract_final_prediction_ids(new_data, id_col)
  )
  if (!is.null(id_attempt$error)) {
    return(sb_signature_prediction_failure(
      signature,
      id_col,
      character(),
      "ids",
      id_attempt$error
    ))
  }
  patient_ids <- id_attempt$value

  signature_attempt <- sb_final_fit_attempt(
    sb_validate_linear_signature(signature)
  )
  if (!is.null(signature_attempt$error)) {
    return(sb_signature_prediction_failure(
      signature,
      id_col,
      patient_ids,
      "signature",
      signature_attempt$error
    ))
  }
  signature <- signature_attempt$value
  metadata <- sb_signature_safe_metadata(signature)

  feature_attempt <- sb_final_fit_attempt(
    sb_select_signature_prediction_features(
      new_data,
      id_col,
      signature$features
    )
  )
  if (!is.null(feature_attempt$error)) {
    return(sb_signature_prediction_failure(
      signature,
      id_col,
      patient_ids,
      "features",
      feature_attempt$error
    ))
  }
  value_attempt <- sb_final_fit_attempt(
    sb_signature_prediction_values(signature, feature_attempt$value)
  )
  if (!is.null(value_attempt$error)) {
    return(sb_signature_prediction_failure(
      signature,
      id_col,
      patient_ids,
      "formula",
      value_attempt$error
    ))
  }

  result <- list(
    status = "success",
    task_type = signature$task_type,
    pipeline_id = signature$pipeline_id,
    learner = signature$learner,
    id_col = id_col,
    signature = signature,
    predictions = sb_signature_prediction_rows(
      metadata,
      patient_ids,
      "success",
      NA_character_,
      NA_character_,
      value_attempt$value
    ),
    provenance = sb_signature_prediction_provenance(),
    error_stage = NA_character_,
    error_class = NA_character_,
    error_message = NA_character_
  )
  structure(result, class = "stormbreaker_final_prediction")
}
