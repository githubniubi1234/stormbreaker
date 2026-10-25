sb_akritas_learner_ids <- function() "akritas_survival"

sb_akritas_params <- function() list(lambda = 0.5)

sb_require_akritas_dependencies <- function() {
  for (package in c("survivalmodels", "distr6")) if (!requireNamespace(package, quietly = TRUE))
    sb_abort(paste0("Install `", package, "` for native Akritas prediction."),
      "stormbreaker_missing_dependency")
}

sb_akritas_distribution <- function(basis) {
  if (ncol(basis) == 1L) distr6::Empirical$new(basis) else distr6::EmpiricalMV$new(basis)
}

sb_akritas_empirical_cdf <- function(basis, x) {
  if (ncol(basis) == 1L) {
    return(as.numeric(sb_akritas_distribution(basis)$cdf(data = x)))
  }
  # Same C routine as distr6::EmpiricalMV$.cdf, without reconstructing R6
  # support/parameter objects for already validated numeric matrices.
  native <- utils::getFromNamespace("C_EmpiricalMVCdf", "distr6")
  as.numeric(native(x = x, data = basis))
}

sb_fit_akritas_survival <- function(x, time, event, learner, params) {
  sb_require_akritas_dependencies()
  inputs <- sb_validate_classical_survival_inputs(x, time, event, learner)
  basis <- unname(as.matrix(inputs$x))
  # The native C engine expects decreasing follow-up time. Keep the empirical
  # distribution basis in input order, and reorder only the risk-set inputs.
  ord <- order(inputs$time, decreasing = TRUE)
  list(learner = learner, features = names(inputs$x), params = params,
    basis = basis, truth = unname(cbind(inputs$time, inputs$event)[ord, , drop = FALSE]),
    cdf_train = sb_akritas_empirical_cdf(basis, basis)[ord],
    grid = sort(unique(inputs$time)))
}

sb_akritas_curves <- function(model, x, times = model$grid) {
  sb_require_akritas_dependencies()
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) sb_abort(
    "Prediction columns differ from the fitted Akritas predictors.", "stormbreaker_learner_schema_mismatch")
  # survivalmodels 0.1.191 retains this native engine but removed the exported
  # R wrapper. Use the installed engine, not a replacement nearest-neighbour fit.
  native <- utils::getFromNamespace("C_Akritas", "survivalmodels")
  cdf <- sb_akritas_empirical_cdf(model$basis, unname(as.matrix(x)))
  distinct_cdf <- unique(cdf)
  # The native query input is this scalar CDF, not the predictor row. Compute
  # identical inputs once within this call, then restore patient order exactly.
  curve <- native(model$truth, as.numeric(times), model$grid,
    model$cdf_train, distinct_cdf,
    model$params$lambda)
  curve <- curve[match(cdf, distinct_cdf), , drop = FALSE]
  if (!is.matrix(curve) || !identical(dim(curve), c(nrow(x), length(times))) ||
      any(!is.finite(curve)) || any(curve < 0 | curve > 1))
    sb_abort("Native Akritas returned invalid survival curves.", "stormbreaker_invalid_akritas_prediction")
  curve
}

sb_predict_akritas_risk <- function(model, x) {
  curve <- sb_akritas_curves(model, x)
  sb_akritas_risk_from_curves(curve, model$grid)
}

sb_akritas_risk_from_curves <- function(curve, grid) {
  # Integral of the right-continuous survival step function up to max analysis
  # follow-up, not the historical wrapper's incomplete event-mass mean.
  before <- cbind(1, curve)[, seq_along(grid), drop = FALSE]
  -as.numeric(before %*% diff(c(0, grid)))
}

sb_akritas_direct_source <- function() "native_akritas_conditional_product_limit"

sb_akritas_survival_at_horizon <- function(model, x_analysis, time, event,
                                         x_assessment, prediction_horizon, ...) {
  horizon <- sb_validate_prediction_horizon(prediction_horizon, time)
  list(prediction_horizon = horizon, source = sb_akritas_direct_source(),
    ties = "native_grouped_event_times",
    probability = as.numeric(sb_akritas_curves(model, x_assessment, horizon)))
}

sb_extract_akritas_state <- function(model, context = list()) {
  state <- model
  state$training_data_retention <- "predictor_basis_and_paired_outcomes_required_for_conditional_product_limit"
  state$survival_probability <- if (!is.null(context$prediction_horizon)) list(
    prediction_horizon = sb_validate_prediction_horizon(context$prediction_horizon, context$time),
    time_unit = context$time_unit, source = sb_akritas_direct_source(), ties = "native_grouped_event_times") else NULL
  class(state) <- c("stormbreaker_akritas_state", "stormbreaker_learner_state")
  state
}

sb_predict_akritas_state <- function(state, x) {
  x <- sb_validate_learner_state_x(state, x, state$learner)
  times <- sort(unique(c(state$grid, state$survival_probability$prediction_horizon)))
  curves <- sb_akritas_curves(state, x, times)
  result <- list(risk = sb_akritas_risk_from_curves(
    curves[, match(state$grid, times), drop = FALSE], state$grid))
  if (!is.null(state$survival_probability)) {
    p <- state$survival_probability
    result$survival_probability <- as.numeric(curves[, match(p$prediction_horizon, times)])
    result$prediction_horizon <- p$prediction_horizon
    result$time_unit <- p$time_unit
  }
  result
}

sb_akritas_learner_adapter <- function(learner, params = sb_akritas_params()) {
  list(id = learner, params = params,
    fit = function(x, time, event) sb_fit_akritas_survival(x, time, event, learner, params),
    predict = sb_predict_akritas_risk, predict_survival = sb_akritas_survival_at_horizon,
    extract_state = sb_extract_akritas_state, predict_state = sb_predict_akritas_state)
}
