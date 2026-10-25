#' Plot cumulative/dynamic survival ROC at multiple time points
#' @param task Survival task with one independent patient per row.
#' @param risk Finite numeric prognostic scores, higher means worse prognosis.
#'   Unnamed scores follow task row order; named scores are aligned by patient
#'   IDs. Use one saved model's score, not separately selected models per time.
#' @param times Positive time points in the task's declared time unit.
#' @return List with plot, coordinates and per-time AUC/status table.
#' @details Cases have observed events at or before t; controls survive beyond
#'   t. Earlier censored patients are neither cases nor controls. IPCW uses
#'   the supplied cohort's marginal reverse Kaplan-Meier censoring estimate:
#'   cases receive 1/G(T-), controls 1/G(t). Tied scores enter together and AUC
#'   uses trapezoidal integration (half credit for ties). This assumes marginal
#'   independent censoring. No confidence intervals or automatic model ranking
#'   are produced; unsupported horizons remain in the status table. Do not pool
#'   incompatible fold-specific score scales. The score is not a probability.
#' @export
surv_plot_roc <- function(task, risk, times) {
  sb_require_curve_plotting()
  task <- sb_risk_plot_task(task)
  risk <- sb_risk_plot_scores(risk, task)
  if (!is.numeric(times) || !length(times) || any(!is.finite(times)) || any(times <= 0) || anyDuplicated(times)) stop("Supply unique positive times in the task time unit.")
  time <- task$data[[task$time_col]]
  event <- task$data[[task$event_col]]
  censoring <- sb_fit_censoring_km(time, event)
  curves <- list()
  statistics <- lapply(sort(times), function(t) {
    cases <- event == 1 & time <= t
    controls <- time > t
    row <- data.frame(time = t, auc = NA_real_, n_cases = sum(cases), n_controls = sum(controls),
      n_early_censored = sum(event == 0 & time <= t), status = "not_evaluable", reason = "")
    if (!any(cases) || !any(controls)) {
      row$reason <- "No cases or no controls at this time."
      return(row)
    }
    gcase <- vapply(time[cases], function(s) sb_censoring_survival_at(censoring, s, left_limit = TRUE), numeric(1))
    gcontrol <- sb_censoring_survival_at(censoring, t)
    if (any(!is.finite(c(gcase, gcontrol))) || any(c(gcase, gcontrol) <= 0)) {
      row$reason <- "Censoring survival is zero or unavailable."
      return(row)
    }
    weights_case <- weights_control <- numeric(length(risk))
    weights_case[cases] <- 1 / gcase
    weights_control[controls] <- 1 / gcontrol
    ord <- order(risk, decreasing = TRUE)
    last <- c(which(diff(risk[ord]) != 0), length(risk))
    sensitivity <- c(0, cumsum(weights_case[ord])[last] / sum(weights_case))
    fpr <- c(0, cumsum(weights_control[ord])[last] / sum(weights_control))
    row$auc <- sum(diff(fpr) * (utils::head(sensitivity, -1) + utils::tail(sensitivity, -1)) / 2)
    row$status <- "success"
    curves[[as.character(t)]] <<- data.frame(time = t, threshold = c(Inf, risk[ord][last]),
      fpr = fpr, sensitivity = sensitivity, label = sprintf("%g %s: AUC %.3f", t, task$time_unit, row$auc))
    row
  })
  statistics <- do.call(rbind, statistics)
  coordinates <- if (length(curves)) do.call(rbind, curves) else
    data.frame(time = numeric(), threshold = numeric(), fpr = numeric(), sensitivity = numeric(), label = character())
  fpr <- sensitivity <- label <- NULL
  p <- ggplot2::ggplot(coordinates, ggplot2::aes(fpr, sensitivity, colour = label)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70") +
    ggplot2::geom_path(linewidth = .7) + ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "1 - Specificity", y = "Sensitivity", colour = NULL) + sb_viz_theme()
  list(plot = p, coordinates = coordinates, statistics = statistics,
    method = "cumulative_dynamic_ipcw_marginal_reverse_km", time_unit = task$time_unit)
}
