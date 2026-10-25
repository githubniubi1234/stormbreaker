sb_validate_stepcox_direction <- function(direction) {
  valid <- is.character(direction) && length(direction) == 1L &&
    !is.na(direction) && direction %in% c("forward", "backward", "both")
  if (!valid) {
    sb_abort(
      paste0(
        "`stepcox_direction` must be one of `forward`, `backward`, or `both`."
      ),
      "stormbreaker_bad_stepcox_direction"
    )
  }
  unname(as.character(direction))
}

sb_stepcox_selection_path <- function(model, internal_features,
                                      source_features, final_aic) {
  if (is.null(model$anova)) {
    return(data.frame(
      Step = "",
      AIC = final_aic,
      stringsAsFactors = FALSE
    ))
  }
  path <- as.data.frame(model$anova, stringsAsFactors = FALSE)
  actions <- as.character(path$Step)
  actions[is.na(actions)] <- ""
  for (index in seq_along(actions)) {
    action <- trimws(actions[[index]])
    if (!nzchar(action)) {
      actions[[index]] <- ""
      next
    }
    operation <- substr(action, 1L, 1L)
    internal <- trimws(substring(action, 2L))
    position <- match(internal, internal_features)
    if (!is.na(position) && operation %in% c("+", "-")) {
      actions[[index]] <- paste(operation, source_features[[position]])
    }
  }
  path$Step <- actions
  rownames(path) <- NULL
  path
}

sb_fit_stepcox <- function(x, time, event, direction) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  direction <- sb_validate_stepcox_direction(direction)
  if (any(time < 0)) {
    sb_abort(
      "`time` must contain non-negative values.",
      "stormbreaker_bad_adapter_time"
    )
  }
  if (!all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only 0/1 values.",
      "stormbreaker_bad_adapter_event"
    )
  }

  source_features <- names(x)
  internal_features <- paste0(
    ".stormbreaker_stepcox_x",
    seq_along(source_features)
  )
  analysis <- x
  names(analysis) <- internal_features
  analysis$.stormbreaker_stepcox_time <- time
  analysis$.stormbreaker_stepcox_event <- event
  response <- paste0(
    "survival::Surv(",
    ".stormbreaker_stepcox_time, .stormbreaker_stepcox_event)"
  )
  full_formula <- stats::reformulate(
    internal_features,
    response = response,
    env = environment()
  )
  null_formula <- stats::as.formula(
    paste(response, "~ 1"),
    env = environment()
  )
  start_formula <- if (identical(direction, "forward")) {
    null_formula
  } else {
    full_formula
  }
  # Tighten likelihood convergence without changing the native warning or
  # Cholesky thresholds. Candidate-model warnings still stop selection.
  cox_control <- survival::coxph.control(
    eps = 1e-11,
    iter.max = 50L,
    toler.inf = survival::coxph.control()$toler.inf
  )

  model <- withCallingHandlers(
    {
      start <- survival::coxph(
        start_formula,
        data = analysis,
        ties = "efron",
        singular.ok = FALSE,
        x = TRUE,
        y = TRUE,
        model = TRUE,
        control = cox_control
      )
      stats::step(
        start,
        scope = list(lower = null_formula, upper = full_formula),
        direction = direction,
        trace = 0,
        k = 2
      )
    },
    warning = function(warning) {
      sb_abort(
        paste0(
          "The StepCox selector raised a warning: ",
          conditionMessage(warning)
        ),
        "stormbreaker_stepcox_warning",
        parent = warning
      )
    }
  )
  if (NROW(model$y) != nrow(analysis)) {
    sb_abort(
      "The StepCox selector dropped analysis patients.",
      "stormbreaker_selector_rows_dropped"
    )
  }

  coefficients <- stats::coef(model)
  if (length(coefficients) < 1L) {
    sb_abort(
      paste0(
        "StepCox selected no predictor with direction `", direction, "`."
      ),
      "stormbreaker_stepcox_empty_selection"
    )
  }
  selected_positions <- match(names(coefficients), internal_features)
  valid_coefficients <- !anyNA(selected_positions) &&
    !anyNA(coefficients) && all(is.finite(coefficients))
  if (!valid_coefficients) {
    sb_abort(
      "The StepCox selector produced an invalid fitted model.",
      "stormbreaker_invalid_stepcox_fit"
    )
  }
  selected_features <- source_features[selected_positions]
  final_aic <- unname(stats::extractAIC(model, k = 2)[[2L]])
  if (!is.finite(final_aic)) {
    sb_abort(
      "The StepCox selector produced a non-finite AIC.",
      "stormbreaker_invalid_stepcox_fit"
    )
  }

  scores <- data.frame(
    feature = source_features,
    original_order = as.integer(seq_along(source_features)),
    selected = source_features %in% selected_features,
    coefficient = NA_real_,
    stringsAsFactors = FALSE
  )
  scores$coefficient[selected_positions] <- unname(coefficients)

  structure(
    list(
      id = "stepcox",
      direction = direction,
      criterion = "AIC",
      penalty = 2,
      start_model = if (identical(direction, "forward")) "null" else "full",
      ties = "efron",
      cox_control = cox_control,
      source_features = source_features,
      selected_features = selected_features,
      scores = scores,
      final_aic = final_aic,
      selection_path = sb_stepcox_selection_path(
        model,
        internal_features,
        source_features,
        final_aic
      )
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_apply_stepcox <- function(fitted, x) {
  valid_fit <- inherits(fitted, "stormbreaker_selector_fit") &&
    identical(fitted$id, "stepcox") &&
    is.character(fitted$direction) && length(fitted$direction) == 1L &&
    fitted$direction %in% c("forward", "backward", "both") &&
    is.character(fitted$source_features) &&
    is.character(fitted$selected_features) &&
    length(fitted$selected_features) > 0L &&
    anyDuplicated(fitted$selected_features) == 0L &&
    all(fitted$selected_features %in% fitted$source_features)
  if (!valid_fit) {
    sb_abort(
      "The fitted selector is not a valid StepCox state.",
      "stormbreaker_bad_selector_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), fitted$source_features)) {
    sb_abort(
      paste0(
        "Predictor columns must exactly match the fitted selector schema: ",
        paste(fitted$source_features, collapse = ", "),
        "."
      ),
      "stormbreaker_selector_schema_mismatch"
    )
  }
  x[fitted$selected_features]
}

sb_stepcox_selector_adapter <- function(direction) {
  direction <- sb_validate_stepcox_direction(direction)
  list(
    id = "stepcox",
    params = list(direction = direction),
    needs_outcome = TRUE,
    fit = function(x, time, event) {
      sb_fit_stepcox(x, time, event, direction)
    },
    apply = sb_apply_stepcox
  )
}
