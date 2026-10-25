sb_risk_cox_estimate <- function(d, type = "group") {
  tryCatch({
    value <- if (type == "group") as.numeric(d$group == "High") else d$risk
    if (length(unique(value)) < 2L || !any(d$event == 1)) stop("Insufficient groups, variation or events.")
    model <- withCallingHandlers(survival::coxph(survival::Surv(time, event) ~ value, data = d),
      warning = function(w) stop(conditionMessage(w)))
    te <- unname(stats::coef(model)[1]); se <- sqrt(stats::vcov(model)[1, 1])
    if (!is.finite(te) || !is.finite(se) || se <= 0) stop("Cox estimate is not finite.")
    data.frame(log_hr = te, se = se, hr = exp(te), lower = exp(te - stats::qnorm(.975) * se),
      upper = exp(te + stats::qnorm(.975) * se), p = 2 * stats::pnorm(-abs(te / se)),
      status = "success", reason = "")
  }, error = function(e) data.frame(log_hr = NA_real_, se = NA_real_, hr = NA_real_,
    lower = NA_real_, upper = NA_real_, p = NA_real_, status = "not_evaluable", reason = conditionMessage(e)))
}

#' Plot one model's performance across cohorts
#' @param performance Data frame with unique cohort names and a numeric column
#'   named c_index or auc. Supply one model and one time horizon at a time.
#' @param metric Metric column to display.
#' @param title Optional short title, for example '1-year AUC'.
#' @return A ggplot. Estimates are supplied, not recalculated or ranked.
#' @export
sb_plot_cohort_metrics <- function(performance, metric = c("c_index", "auc"), title = NULL) {
  metric <- match.arg(metric)
  sb_require_curve_plotting()
  if (!is.data.frame(performance) || !all(c("cohort", metric) %in% names(performance)) ||
      !nrow(performance) || anyNA(performance$cohort) || anyDuplicated(performance$cohort) ||
      !is.numeric(performance[[metric]]) || any(!is.finite(performance[[metric]])) ||
      any(performance[[metric]] < 0 | performance[[metric]] > 1)) stop("Supply unique cohorts and finite metric values in [0,1].")
  d <- data.frame(cohort = factor(performance$cohort, levels = rev(performance$cohort)), value = performance[[metric]])
  cohort <- value <- NULL
  ggplot2::ggplot(d, ggplot2::aes(value, cohort, fill = cohort)) +
    ggplot2::geom_col(width = .72, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", value)), hjust = -.1, size = 3.5) +
    ggplot2::scale_x_continuous(limits = c(0, 1.08), breaks = seq(0, 1, .2)) +
    ggplot2::labs(x = if (metric == "c_index") "C-index" else "AUC", y = NULL, title = title) + sb_viz_theme()
}

#' Pool independent cohorts' risk-score Cox hazard ratios
#' @param cohorts Named list of surv_risk_plots results for independent cohorts.
#'   Use the same saved model and training cutoff in every cohort. Do not include
#'   a pooled cohort alongside its constituent cohorts; IDs must not overlap.
#' @param type group compares High with Low. continuous reports HR per one unit
#'   of the saved model score, without per-cohort standardization.
#' @return Cohort estimates, fixed and random effects, heterogeneity, and a
#'   forest plot. Failed cohort estimates are returned and not silently pooled.
#' @details Univariate Cox models use Efron ties and Wald intervals. Meta-analysis
#'   uses metafor inverse-variance fixed effects and DerSimonian-Laird random
#'   effects with normal intervals. At least two evaluable independent cohorts
#'   are required. Matching IDs can detect overlap, not prove independence.
#'   Training-cohort associations are apparent, not independent validation.
#' @export
surv_meta_risk <- function(cohorts, type = c("group", "continuous")) {
  type <- match.arg(type)
  sb_require_curve_plotting()
  if (!is.list(cohorts) || length(cohorts) < 2L || is.null(names(cohorts)) ||
      anyNA(names(cohorts)) || any(!nzchar(names(cohorts))) || anyDuplicated(names(cohorts))) stop("Supply at least two named independent cohorts.")
  required <- c("patient_id", "risk", "time", "event", "group")
  for (z in cohorts) if (!is.data.frame(z$patients) || !all(required %in% names(z$patients))) stop("Use surv_risk_plots results.")
  ids <- unlist(lapply(cohorts, function(z) as.character(z$patients$patient_id)), use.names = FALSE)
  if (anyNA(ids) || anyDuplicated(ids)) stop("Patient IDs overlap: do not pool combined and constituent cohorts.")
  if (type == "group" && length(unique(vapply(cohorts, function(z) z$cutoff, numeric(1)))) != 1L) stop("Use the same training cutoff in every cohort.")
  estimates <- do.call(rbind, lapply(seq_along(cohorts), function(i) {
    cbind(cohort = names(cohorts)[i], sb_risk_cox_estimate(cohorts[[i]]$patients, type))
  }))
  if (any(estimates$status != "success")) return(list(estimates = estimates,
    status = "not_evaluable", plot = NULL, reason = "Review failed cohort estimates before pooling."))
  if (!requireNamespace("metafor", quietly = TRUE)) stop("Install metafor for Cox meta-analysis.")
  fixed <- metafor::rma.uni(yi = estimates$log_hr, sei = estimates$se, method = "FE")
  random <- metafor::rma.uni(yi = estimates$log_hr, sei = estimates$se, method = "DL")
  estimates$weight_fixed <- as.numeric(stats::weights(fixed))
  estimates$weight_random <- as.numeric(stats::weights(random))
  pooled <- do.call(rbind, lapply(list(Fixed = fixed, Random = random), function(m) {
    data.frame(log_hr = as.numeric(m$b), se = m$se, hr = exp(as.numeric(m$b)),
      lower = exp(m$ci.lb), upper = exp(m$ci.ub), p = m$pval)
  }))
  pooled$cohort <- rownames(pooled)
  d <- rbind(estimates[c("cohort", "hr", "lower", "upper")], pooled[c("cohort", "hr", "lower", "upper")])
  d$cohort <- factor(d$cohort, levels = rev(d$cohort))
  cohort <- hr <- lower <- upper <- NULL
  plot <- ggplot2::ggplot(d, ggplot2::aes(hr, cohort)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper), width = .15, orientation = "y") +
    ggplot2::geom_point(size = 2.5) + ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Hazard ratio (95% CI)", y = NULL) + sb_viz_theme()
  list(estimates = estimates, pooled = pooled, plot = plot, type = type, status = "success",
    heterogeneity = data.frame(Q = random$QE, p = random$QEp, I2 = random$I2, tau2 = random$tau2))
}
