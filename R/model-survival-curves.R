#' Add a full baseline curve to a saved linear survival model
#'
#' @param fit Successful portable linear Cox final fit with a saved horizon.
#' @param development The original, unchanged development survival task.
#'   Test or validation outcomes must never be supplied. IDs, sample count,
#'   time units and the saved horizon hazard are checked, but those checks
#'   cannot prove that all supplied values are the original training values.
#' @return A copy of fit with `survival_curve` metadata and time/hazard table.
#'   Existing coefficients, learner state, and horizon predictions are unchanged.
#' @details Estimates only the missing baseline at fixed saved coefficients,
#'   using the saved Breslow or Efron tie method. Unsupported baseline methods
#'   are refused. No variable selection, tuning or coefficient fitting occurs.
#'   Store the returned model in a new file to preserve the original artifact.
#' @export
sb_add_survival_curve <- function(fit, development) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  state <- fit$learner_state
  baseline <- state$baseline_hazard
  if (!state$learner %in% sb_linear_survival_learner_ids() ||
    is.null(baseline) || !baseline$ties %in% c("breslow", "efron")) {
    sb_abort(
      "A linear Cox model with a saved Breslow/Efron baseline is required.",
      "stormbreaker_unsupported_survival_curve"
    )
  }
  development <- sb_validate_task_contract(development)
  ids <- as.character(development$data[[development$id_col]])
  if (!setequal(ids, fit$development$patient_ids) ||
    length(ids) != length(fit$development$patient_ids) ||
    !identical(development$time_unit, baseline$time_unit)) {
    sb_abort(
      "Supply the original development task with matching IDs and time units.",
      "stormbreaker_curve_development_mismatch"
    )
  }
  x <- sb_select_final_feature_frame(
    development$data, state$features,
    "Development data", "development fitted features"
  )
  risk <- drop(as.matrix(x) %*% state$coefficient)
  time <- development$data[[development$time_col]]
  event <- development$data[[development$event_col]]
  horizon <- baseline$prediction_horizon
  sb_validate_prediction_horizon(horizon, time)
  event_times <- sort(unique(time[event == 1L]))
  increments <- vapply(event_times, function(t) {
    at_risk <- time >= t
    deaths <- time == t & event == 1L
    n_deaths <- sum(deaths)
    maximum <- max(risk[at_risk])
    denominator <- sum(exp(risk[at_risk] - maximum))
    if (baseline$ties == "breslow") {
      return(exp(log(n_deaths) - maximum - log(denominator)))
    }
    death_risk <- sum(exp(risk[deaths] - maximum))
    sum(exp(-maximum - log(denominator - (seq_len(n_deaths) - 1) / n_deaths * death_risk)))
  }, numeric(1))
  hazard <- cumsum(increments)
  if (any(!is.finite(hazard)) || any(increments <= 0)) {
    sb_abort(
      "Recovered baseline hazard is not finite and nondecreasing.",
      "stormbreaker_invalid_baseline_hazard"
    )
  }
  at_horizon <- sum(increments[event_times <= horizon])
  if (!isTRUE(all.equal(at_horizon, baseline$cumulative_baseline_hazard,
    tolerance = 1e-8, check.attributes = FALSE
  ))) {
    sb_abort(
      "Recovered baseline does not match the saved horizon; check original training data and method.",
      "stormbreaker_curve_baseline_mismatch"
    )
  }
  grid <- sort(unique(c(0, event_times, horizon, max(time))))
  values <- c(0, hazard)[findInterval(grid, event_times) + 1L]
  fit$survival_curve <- list(
    table = data.frame(
      time = grid,
      cumulative_hazard = values, baseline_survival = exp(-values)
    ),
    time_unit = baseline$time_unit, prediction_horizon = horizon,
    max_followup = max(time), source = baseline$source, ties = baseline$ties,
    reference = "x_beta = 0; not an average patient",
    development_n = length(ids), method = "original_development_fixed_coefficients"
  )
  fit
}

#' Predict full survival step curves without reading patient outcomes
#'
#' @param fit A final fit augmented by [sb_add_survival_curve()].
#' @param new_data New feature data with unique identifiers; other columns ignored.
#' @param id_col Unique identifier column name.
#' @param times Optional finite nonnegative time points within training follow-up.
#'   Sorted unique points are returned. NULL uses the full saved step grid.
#' @return A `stormbreaker_survival_curves` list with baseline table and
#'   time-by-patient survival matrix, IDs, horizon and unchanged final prediction.
#'   Full step-grid values are retained separately for faithful plotting even
#'   when only a few prediction times were requested.
#'   No uncertainty intervals are inferred from point estimates.
#' @export
sb_survival_curves <- function(fit, new_data, id_col, times = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  curve <- fit$survival_curve
  if (is.null(curve)) {
    sb_abort(
      "Add the original-training baseline with sb_add_survival_curve first.",
      "stormbreaker_missing_survival_curve"
    )
  }
  ids <- sb_extract_final_prediction_ids(new_data, id_col)
  x <- sb_select_final_prediction_features(new_data, id_col, fit$required_features)
  if (is.null(times)) times <- curve$table$time
  if (!is.numeric(times) || !length(times) || any(!is.finite(times)) ||
    any(times < 0 | times > curve$max_followup)) {
    sb_abort(
      "Times must be finite, nonnegative and within original training follow-up.",
      "stormbreaker_curve_time_unsupported"
    )
  }
  times <- sort(unique(times))
  hazard <- curve$table$cumulative_hazard[findInterval(times, curve$table$time)]
  risk <- drop(as.matrix(x) %*% fit$learner_state$coefficient)
  probability <- vapply(
    hazard, function(h) sb_survival_probability_from_hazard(risk, h),
    numeric(length(ids))
  )
  probability <- matrix(probability, nrow = length(ids), ncol = length(times))
  full_probability <- vapply(
    curve$table$cumulative_hazard,
    function(h) sb_survival_probability_from_hazard(risk, h), numeric(length(ids))
  )
  full_probability <- matrix(full_probability, nrow = length(ids), ncol = nrow(curve$table))
  structure(
    list(
      baseline = data.frame(
        time = times, cumulative_hazard = hazard,
        survival = exp(-hazard)
      ), survival = t(probability), patient_ids = ids,
      plot_baseline = data.frame(
        time = curve$table$time,
        survival = curve$table$baseline_survival
      ), plot_survival = t(full_probability),
      time_unit = curve$time_unit, prediction_horizon = curve$prediction_horizon,
      max_followup = curve$max_followup, reference = curve$reference,
      prediction = sb_predict_final(fit, new_data, id_col)
    ),
    class = "stormbreaker_survival_curves"
  )
}

#' Plot a baseline and selected individual survival curves
#'
#' @param x A result from [sb_survival_curves()].
#' @param rows Patient row numbers, or integer(0) for baseline only.
#' @param baseline Whether to include the zero-score baseline.
#' @return A ggplot step curve with the saved-horizon probability markers.
#'   Patient row numbers replace original IDs, but predictions remain sensitive.
#' @export
sb_plot_survival_curves <- function(x, rows = 1L, baseline = TRUE) {
  sb_require_curve_plotting()
  if (!inherits(x, "stormbreaker_survival_curves")) {
    sb_abort(
      "Supply an sb_survival_curves result.", "stormbreaker_bad_survival_curves"
    )
  }
  if (!is.numeric(rows) || anyNA(rows) || any(rows != floor(rows)) ||
    any(rows < 1 | rows > ncol(x$survival)) || anyDuplicated(rows) ||
    !is.logical(baseline) || length(baseline) != 1L || is.na(baseline) ||
    (!length(rows) && !baseline)) {
    sb_abort(
      "Choose valid unique patient row numbers and/or the baseline.",
      "stormbreaker_bad_survival_curves"
    )
  }
  d <- data.frame(time = numeric(), survival = numeric(), series = character())
  if (baseline) {
    d <- data.frame(
      time = x$plot_baseline$time,
      survival = x$plot_baseline$survival, series = "Baseline (x beta = 0)"
    )
  }
  for (i in rows) {
    d <- rbind(d, data.frame(
      time = x$plot_baseline$time,
      survival = x$plot_survival[, i], series = paste("Observation", i)
    ))
  }
  mark <- d[d$time == x$prediction_horizon, , drop = FALSE]
  labels <- paste(mark$series, sprintf("%.1f%%", 100 * mark$survival), collapse = " | ")
  palette <- if (length(rows)) {
    stats::setNames(
      rep(c("#24547A", "#B86236", "#3A806E", "#805A8D"), length.out = length(rows)),
      paste("Observation", rows)
    )
  } else {
    character()
  }
  if (baseline) palette <- c("Baseline (x beta = 0)" = "#838B91", palette)
  time <- survival <- series <- NULL
  ggplot2::ggplot(d, ggplot2::aes(x = time, y = survival, colour = series)) +
    ggplot2::geom_step(direction = "hv", linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = x$prediction_horizon, linetype = 3, colour = "#7A8791") +
    ggplot2::geom_point(data = mark, size = 2) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(
      title = "Model-predicted survival", colour = NULL,
      x = paste("Time |", x$time_unit), y = "Survival probability",
      subtitle = sb_viz_wrap(paste("Horizon:", x$prediction_horizon, x$time_unit, "|", labels)),
      caption = sb_viz_wrap(paste(
        "Fixed model coefficients and original-training baseline.",
        "The zero-score baseline is not an average patient. No uncertainty intervals.",
        "No extrapolation beyond observed training follow-up."
      ))
    ) +
    sb_viz_theme()
}
