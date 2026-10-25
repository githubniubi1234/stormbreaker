sb_ncvreg_path_diagnostic <- function(fit, type) {
  model <- fit$diagnostic_artifacts$learner
  path <- model$diagnostic_path
  if (!is.list(path) || !is.numeric(path$lambda) || !length(path$lambda) ||
      any(!is.finite(path$lambda)) || any(path$lambda <= 0) ||
      !is.numeric(path$bic) || length(path$bic) != length(path$lambda) ||
      any(!is.finite(path$bic))) sb_abort(
    "The original MCP/SCAD path is missing or invalid. Older fits may lack path capture.",
    "stormbreaker_invalid_diagnostic")
  beta <- as.matrix(path$beta)
  expected <- if (model$learner %in% sb_ncvreg_binary_learner_ids()) c("(Intercept)", model$features) else model$features
  if (!is.numeric(beta) || ncol(beta) != length(path$lambda) ||
      !identical(rownames(beta), expected) || any(!is.finite(beta))) sb_abort(
    "Native MCP/SCAD coefficient path does not match the fitted features.",
    "stormbreaker_invalid_diagnostic")
  if (type == "bic_path") {
    tab <- data.frame(lambda = path$lambda, bic = path$bic)
  } else {
    tab <- data.frame(lambda = rep(path$lambda, each = nrow(beta)),
      feature = rep(rownames(beta), ncol(beta)), coefficient = as.vector(beta))
  }
  structure(list(type = type, table = tab, learner = fit$pipeline$learner,
    pipeline_id = fit$pipeline$id, engine_class = path$engine_class,
    selected_lambda = model$effective$selected_lambda,
    source = "native original-fit MCP/SCAD path; intercept explicitly labelled when present"),
    class = "stormbreaker_algorithm_diagnostic")
}
