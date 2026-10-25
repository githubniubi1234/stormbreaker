sb_pamr_selector_ids <- function() c("pamr_binary_nonzero", "pamr_survival_nonzero")

sb_pamr_selector_params <- function(id) {
  if (id == "pamr_survival_nonzero") sb_pamr_survival_params() else
    sb_classical_binary_params("pamr_binary")
}

sb_pamr_selector_state <- function(engine, x, id, params) {
  selected_index <- pamr::pamr.predict(engine, t(as.matrix(x)),
    threshold = params$threshold, type = "nonzero")
  centroids <- pamr::pamr.predict(engine, t(as.matrix(x)),
    threshold = params$threshold, type = "centroid")
  score <- apply(abs((centroids - engine$centroid.overall) / engine$sd), 1L, max)
  if (length(score) != ncol(x) || any(!is.finite(score))) sb_abort(
    "PAMR returned invalid shrunken-centroid scores.", "stormbreaker_invalid_pamr_selector")
  if (!length(selected_index)) sb_abort(
    "PAMR selected no nonzero shrunken-centroid features; no genes were added.",
    "stormbreaker_empty_selection")
  selected <- seq_len(ncol(x)) %in% selected_index
  structure(list(id = id, params = params, source_features = names(x),
    selected_features = names(x)[selected],
    scores = data.frame(feature = names(x), score = unname(score), selected = selected),
    evidence = list(rule = "native_nonzero_shrunken_centroids",
      pamr_version = as.character(utils::packageVersion("pamr")))),
    class = "stormbreaker_selector_fit")
}

sb_pamr_selector_adapter <- function(id, params) {
  list(id = id, params = params, needs_outcome = TRUE,
    fit = if (id == "pamr_survival_nonzero") {
      function(x, time, event) {
        model <- sb_fit_pamr_survival(x, time, event, "pamr_survival", params)
        sb_pamr_selector_state(model$engine, x, id, params)
      }
    } else {
      function(x, outcome) {
        model <- sb_fit_classical_binary(x, outcome, "pamr_binary", params)
        sb_pamr_selector_state(model$engine, x, id, params)
      }
    },
    apply = function(fitted, x) {
      x <- sb_validate_predictor_frame(x, "x")
      if (!identical(names(x), fitted$source_features)) stop("Predictors differ from fitted PAMR selection.")
      x[fitted$selected_features]
    })
}
