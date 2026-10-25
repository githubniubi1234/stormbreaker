# Same training-only coordinates used by the formal-run adapters. This removes
# the known constant-sum direction, not low-variance principal components.
sb_constant_sum_projection <- function(x) {
  matrix_x <- as.matrix(x)
  p <- ncol(matrix_x)
  sums <- rowSums(matrix_x)
  projection <- list(features = names(x), basis = NULL, center = rep(0, p), row_sum = mean(sums))
  if (p > 1L && max(abs(sums - mean(sums))) < 1e-10) {
    basis <- stats::contr.helmert(p)
    projection$basis <- sweep(basis, 2L, sqrt(colSums(basis * basis)), "/")
    projection$center <- colMeans(matrix_x)
    reconstructed <- sweep(as.matrix(sb_constant_sum_project(x, projection)) %*%
      t(projection$basis), 2L, projection$center, "+")
    if (max(abs(reconstructed - matrix_x)) >= 1e-8) {
      sb_abort("Constant-sum coordinates did not reconstruct training inputs.",
        "stormbreaker_constant_sum_mismatch")
    }
  }
  projection
}

sb_constant_sum_project <- function(x, projection) {
  matrix_x <- as.matrix(x[, projection$features, drop = FALSE])
  if (is.null(projection$basis)) return(as.data.frame(matrix_x))
  if (max(abs(rowSums(matrix_x) - projection$row_sum)) > 1e-8) {
    sb_abort("Prediction input does not obey the fitted constant-sum representation.",
      "stormbreaker_constant_sum_mismatch")
  }
  transformed <- sweep(matrix_x, 2L, projection$center, "-") %*% projection$basis
  colnames(transformed) <- paste0("contrast_", seq_len(ncol(transformed)))
  as.data.frame(transformed)
}

sb_is_constant_sum_state <- function(state) {
  if (!is.list(state) || is.data.frame(state)) return(FALSE)
  identical(state$representation, "constant_sum_orthonormal_v1") ||
    identical(state$round2_representation, "constant_sum_orthonormal_v1") ||
    identical(state$variant, "brglm2_constant_sum_v1")
}

sb_constant_sum_learner_adapter <- function(adapter, constructor) {
  id <- adapter$id
  if (!id %in% c("brglm2_logistic", "lda_binary", "hda_binary", "fda_binary", "mda_binary", "loclda_binary")) {
    return(adapter)
  }
  native <- adapter
  make_native <- function(params) constructor(id, params)
  adapter$fit <- function(x, outcome) {
    projection <- sb_constant_sum_projection(x)
    if (is.null(projection$basis)) return(native$fit(x, outcome))
    inner_params <- native$params
    if (id == "loclda_binary") {
      inner_params$k_max <- sb_loclda_neighbors(inner_params, ncol(x), nrow(x), length(unique(outcome)))
    }
    inner <- make_native(inner_params)$fit(sb_constant_sum_project(x, projection), outcome)
    structure(list(learner = id, features = names(x), params = native$params,
      projection = projection, inner_params = inner_params, inner = inner),
      class = "stormbreaker_constant_sum_model")
  }
  adapter$predict <- function(model, x) {
    if (!inherits(model, "stormbreaker_constant_sum_model")) return(native$predict(model, x))
    make_native(model$inner_params)$predict(model$inner, sb_constant_sum_project(x, model$projection))
  }
  adapter$extract_state <- function(model, context = list()) {
    if (!inherits(model, "stormbreaker_constant_sum_model")) return(native$extract_state(model, context))
    structure(list(learner = id, features = model$features, params = model$params,
      representation = "constant_sum_orthonormal_v1", projection = model$projection,
      inner_params = model$inner_params,
      inner = make_native(model$inner_params)$extract_state(model$inner, context)),
      class = c("stormbreaker_constant_sum_state", "stormbreaker_learner_state"))
  }
  adapter$predict_state <- function(state, x) {
    if (!sb_is_constant_sum_state(state)) return(native$predict_state(state, x))
    params <- if (is.null(state$inner_params)) native$params else state$inner_params
    if (id == "hda_binary") loadNamespace("hda")
    make_native(params)$predict_state(sb_upgrade_legacy_names(state$inner),
      sb_constant_sum_project(x, state$projection))
  }
  if (!is.null(native$extract_coefficients)) {
    adapter$extract_coefficients <- function(model) {
      if (!inherits(model, "stormbreaker_constant_sum_model")) return(native$extract_coefficients(model))
      coefficients <- make_native(model$inner_params)$extract_coefficients(model$inner)
      beta <- as.vector(model$projection$basis %*% coefficients$coefficient[-1L])
      data.frame(term = c("(Intercept)", model$features),
        coefficient = c(coefficients$coefficient[[1L]] - sum(beta * model$projection$center), beta))
    }
  }
  adapter
}

sb_constant_sum_linear_parameters <- function(state) {
  inner <- sb_upgrade_legacy_names(state$inner)
  parameters <- if (inner$learner == "lda_binary") sb_discriminant_parameters(inner) else
    sb_extended_logistic_parameters(inner)
  if (is.null(state$projection$basis)) return(parameters)
  beta <- as.vector(state$projection$basis %*% parameters$coefficient)
  parameters$coefficient <- stats::setNames(beta, state$projection$features)
  parameters$intercept <- parameters$intercept - sum(beta * state$projection$center)
  if (!is.null(parameters$means)) {
    parameters$means <- sweep(parameters$means %*% t(state$projection$basis),
      2L, state$projection$center, "+")
    colnames(parameters$means) <- state$projection$features
  }
  if (!is.null(parameters$quadratic)) {
    parameters$quadratic <- matrix(0, length(beta), length(beta),
      dimnames = list(state$projection$features, state$projection$features))
  }
  parameters
}
