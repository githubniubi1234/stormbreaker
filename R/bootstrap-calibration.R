sb_calibration_smooth <- function(y, probability, grid) {
  if (length(unique(probability)) < 2L) stop("Calibration needs at least two distinct predictions.")
  curve <- stats::lowess(probability, y, f = 2/3, iter = 0)
  stats::approx(curve$x, curve$y, xout = grid, ties = mean, rule = 1)$y
}

sb_calibration_predict <- function(fit, task) {
  if (!identical(fit$status, "success")) stop(fit$error_message)
  prediction <- sb_predict_final(fit, task$data[c(task$id_col, task$predictors)], task$id_col)
  if (!identical(prediction$status, "success")) stop(prediction$error_message)
  prediction$predictions$probability
}

#' Bootstrap optimism-corrected binary calibration with full pipeline refitting
#' @param task Binary task with one independent patient per row. Grouped tasks
#'   are not supported by this ordinary-patient bootstrap interface.
#' @param pipeline A predeclared binary pipeline. Its selection and tuning are
#'   rerun within every bootstrap sample; selection between pipelines is outside
#'   this procedure and is not corrected by it.
#' @param reps Number of ordinary bootstrap samples, default 1000.
#' @param seed Random seed.
#' @param grid Optional increasing prediction grid inside (0,1). The default
#'   spans the original predictions' 5th to 95th percentiles.
#' @return Saved apparent fit, calibration table, pointwise effective replicate
#'   counts, bootstrap differences and failure reasons. Corrected calibration is
#'   apparent minus mean bootstrap training-minus-original calibration.
#'   LOWESS uses span 2/3, zero robustness iterations, linear interpolation and
#'   no extrapolation. Values are not clipped to the unit interval. This is internal
#'   optimism correction, not external validation or a confidence band.
#' @export
sb_bootstrap_calibration <- function(task, pipeline, reps = 1000L, seed = 1L, grid = NULL) {
  task <- sb_validate_binary_task_contract(task)
  if (!is.null(task$group_col)) stop("Use one independent patient per row; grouped bootstrap is not implemented.")
  reps <- sb_viz_positive_integer(reps, "reps")
  if (reps < 2L) stop("Use at least two bootstrap replicates.")
  withr::local_seed(seed)
  fit <- sb_fit_final(task, pipeline)
  p <- sb_calibration_predict(fit, task)
  y <- task$data[[task$outcome_col]]
  if (is.null(grid)) grid <- seq(stats::quantile(p, .05), stats::quantile(p, .95), length.out = 80)
  if (!is.numeric(grid) || length(grid) < 2L || any(!is.finite(grid)) ||
      any(grid <= 0 | grid >= 1) || any(diff(grid) <= 0)) stop("grid must increase strictly inside (0,1).")
  apparent <- sb_calibration_smooth(y, p, grid)
  differences <- matrix(NA_real_, reps, length(grid))
  status <- data.frame(replicate = seq_len(reps), status = "failed", reason = "")
  group_col <- make.unique(c(names(task$data), ".bootstrap_patient"))[ncol(task$data) + 1L]
  for (b in seq_len(reps)) {
    result <- tryCatch({
      i <- sample.int(task$n, task$n, replace = TRUE)
      d <- task$data[i, , drop = FALSE]
      d[[group_col]] <- as.character(d[[task$id_col]])
      d[[task$id_col]] <- paste0("boot_", seq_len(nrow(d)))
      boot <- sb_binary_task(d, task$id_col, task$outcome_col, task$predictors,
        decision_thresholds = task$decision_thresholds, group_col = group_col)
      f <- sb_fit_final(boot, pipeline)
      train <- sb_calibration_predict(f, boot)
      original <- sb_calibration_predict(f, task)
      sb_calibration_smooth(d[[task$outcome_col]], train, grid) -
        sb_calibration_smooth(y, original, grid)
    }, error = identity)
    if (inherits(result, "error")) status$reason[b] <- conditionMessage(result) else {
      differences[b, ] <- result
      status$status[b] <- if (any(is.finite(result))) "success" else "failed"
      if (status$status[b] == "failed") status$reason[b] <- "No overlapping calibration support."
    }
  }
  effective <- colSums(is.finite(differences))
  optimism <- colMeans(differences, na.rm = TRUE)
  optimism[effective < 2L] <- NA_real_
  table <- data.frame(probability = grid, apparent = apparent,
    optimism = optimism, corrected = apparent - optimism, effective_reps = effective)
  structure(list(fit = fit, curve = table, replicates = differences, status = status,
    requested_reps = reps, seed = seed, n = task$n,
    smoother = "LOWESS f=2/3 iter=0; no extrapolation"),
    class = "stormbreaker_bootstrap_calibration")
}

#' Plot saved apparent and optimism-corrected calibration
#' @param x A result from [sb_bootstrap_calibration()].
#' @return A ggplot. Gaps in correction support remain gaps; no confidence band.
#' @export
sb_plot_bootstrap_calibration <- function(x) {
  if (!inherits(x, "stormbreaker_bootstrap_calibration")) stop("Supply sb_bootstrap_calibration output.")
  sb_require_curve_plotting()
  d <- rbind(data.frame(x = x$curve$probability, y = x$curve$apparent, curve = "Apparent"),
    data.frame(x = x$curve$probability, y = x$curve$corrected, curve = "Bias-corrected"))
  y <- curve <- NULL
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y, linetype = curve)) +
    ggplot2::geom_line(linewidth = .7, na.rm = TRUE) +
    ggplot2::geom_abline(intercept = 0, slope = 1, colour = "grey60", linetype = 3) +
    ggplot2::scale_linetype_manual(values = c(Apparent = 2, `Bias-corrected` = 1)) +
    ggplot2::labs(title = "Calibration", x = "Predicted probability", y = "Observed proportion", linetype = NULL) +
    sb_viz_theme() + ggplot2::theme(legend.position = "bottom")
  attr(p, "calibration_diagnostics") <- x[c("status", "requested_reps", "smoother")]
  p
}
