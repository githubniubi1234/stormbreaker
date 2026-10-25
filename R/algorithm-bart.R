#' Inspect original BART posterior predictions
#' @param fit Successful binary, discrete-time survival or AFT BART final fit.
#' @param newdata Explicit query data. No training queries are implicit.
#' @return All retained posterior draws and their mean and pointwise 95 percent
#'   equal-tail credible limits (empirical quantiles, type 7). Binary values are
#'   class-1 probabilities; survival values are survival probabilities at native
#'   time points; AFT values are latent log-time locations, not realized times.
#'   These describe posterior parameter uncertainty, not prediction intervals
#'   for individual outcomes, validation performance or simultaneous bands.
#' @export
sb_bart_posterior <- function(fit, newdata = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  id <- s$learner
  if (!id %in% c("bart_binary", "bart_survival", "bart_aft_survival")) sb_abort(
    "This posterior view requires BART.", "stormbreaker_unsupported_diagnostic")
  if (is.null(newdata)) sb_abort("Supply explicit queries for posterior BART predictions.",
    "stormbreaker_unavailable_diagnostic")
  x <- sb_select_final_feature_frame(newdata, s$features, "Query data", "query fitted features")
  if (!nrow(x)) sb_abort("Queries must contain at least one row.", "stormbreaker_invalid_diagnostic")
  e <- if (is.null(s$engine_raw)) s$engine else unserialize(s$engine_raw)
  n <- nrow(x)
  time <- NA_real_
  if (id == "bart_aft_survival") {
    draws <- sb_bart_aft_log_time_draws(sb_restore_rule_aft_model(s), x)
    scale <- "Latent log-time location"
  } else {
    matrix_x <- as.matrix(x)
    if (id == "bart_survival") {
      time <- as.numeric(e$times)
      matrix_x <- cbind(t = rep(time, n), matrix_x[rep(seq_len(n), each = length(time)), , drop = FALSE])
      if (length(e$rm.const)) matrix_x <- matrix_x[, e$rm.const, drop = FALSE]
    } else matrix_x <- sb_bart_prediction_matrix(e, matrix_x)
    fn <- get(if (id == "bart_binary") "predict.pbart" else "predict.survbart", asNamespace("BART"))
    invisible(utils::capture.output(p <- fn(e, matrix_x, mc.cores = 1L, openmp = FALSE)))
    draws <- if (id == "bart_binary") p$prob.test else p$surv.test
    scale <- if (id == "bart_binary") "Class-1 probability" else "Survival probability"
  }
  if (!is.matrix(draws) || ncol(draws) != n * length(time) || nrow(draws) < 2L ||
      any(!is.finite(draws))) sb_abort("Original posterior draws are incomplete.", "stormbreaker_invalid_diagnostic")
  interval <- apply(draws, 2L, stats::quantile, probs = c(.025, .975), names = FALSE, type = 7)
  tab <- data.frame(query = rep(seq_len(n), each = length(time)), time = rep(time, n),
    mean = colMeans(draws), lower = interval[1L, ], upper = interval[2L, ])
  structure(list(type = "bart_posterior", table = tab, draws = draws, score_scale = scale,
    learner = id, pipeline_id = fit$pipeline$id, source = "Original saved BART posterior tree draws; native prediction",
    interpretation = paste("Pointwise 95% equal-tail posterior credible limits, not simultaneous bands or validation.",
      "AFT values are latent log-time locations, not realized outcome prediction intervals.",
      "No new MCMC, training, surrogate model or patient prediction card.")), class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_bart_posterior <- function(x) {
  query <- time <- mean <- lower <- upper <- NULL
  if (x$learner == "bart_survival") {
    ribbon <- do.call(rbind, lapply(split(x$table, x$table$query), function(d) {
      z <- d[rep(seq_len(nrow(d)), each = 2L)[- (2L * nrow(d))], ]
      z$time <- c(d$time[1L], rep(d$time[-1L], each = 2L))
      z
    }))
    p <- ggplot2::ggplot(x$table, ggplot2::aes(time, mean)) +
      ggplot2::geom_ribbon(data = ribbon, ggplot2::aes(ymin = lower, ymax = upper), fill = "#2878A0", alpha = .18) +
      ggplot2::geom_step(colour = "#2878A0") + ggplot2::facet_wrap(~query, labeller = ggplot2::label_both) +
      ggplot2::labs(x = "Time")
  } else p <- ggplot2::ggplot(x$table, ggplot2::aes(factor(query), mean, ymin = lower, ymax = upper)) +
    ggplot2::geom_pointrange(colour = "#2878A0") + ggplot2::labs(x = "Query")
  p + ggplot2::labs(title = "BART posterior predictions", y = x$score_scale) + sb_viz_theme()
}
