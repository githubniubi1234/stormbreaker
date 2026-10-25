sb_risk_plot_task <- function(task) {
  task <- sb_validate_task_contract(task)
  if (!inherits(task, "stormbreaker_surv_task") || !is.null(task$group_col)) {
    stop("Supply a survival task with one independent row per patient.")
  }
  task
}

sb_risk_plot_scores <- function(risk, task) {
  ids <- as.character(task$data[[task$id_col]])
  if (!is.numeric(risk) || !is.null(dim(risk)) || length(risk) != length(ids) || any(!is.finite(risk))) {
    stop("Supply one finite numeric risk score per patient, higher means worse prognosis.")
  }
  if (!is.null(names(risk))) {
    if (anyDuplicated(names(risk)) || !setequal(names(risk), ids)) stop("Named risk IDs must match task patients.")
    risk <- risk[match(ids, names(risk))]
  }
  unname(risk)
}

sb_risk_predict <- function(fit, task) {
  p <- sb_predict_final(fit, task$data, task$id_col)
  if (p$status != "success" || !"risk" %in% names(p$predictions)) stop("Successful survival risk predictions are required.")
  risk <- stats::setNames(p$predictions$risk, p$predictions$patient_id)
  sb_risk_plot_scores(risk, task)
}

#' Lock risk grouping and expression scaling from the development cohort
#' @param fit Successful survival final fit.
#' @param development Original development survival task; one row per patient.
#' @return Reference containing the saved fit, training-median risk threshold,
#'   selected genes, training expression means/SDs and time unit. No refitting.
#' @details Use the original unchanged development data. IDs are checked but
#'   cannot prove that supplied feature values are unchanged. Scores strictly
#'   above the training median are High; all other scores are Low. Constant
#'   training genes use scale 1 for display. No outcome-based cutpoint search.
#' @export
surv_risk_reference <- function(fit, development) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  task <- sb_risk_plot_task(development)
  ids <- as.character(task$data[[task$id_col]])
  if (fit$pipeline$task_type != "survival" || !setequal(ids, fit$development$patient_ids)) {
    stop("Supply the original development cohort for this survival fit.")
  }
  risk <- sb_risk_predict(fit, task)
  genes <- fit$selection$selected_features
  x <- task$data[genes]
  scale <- vapply(x, stats::sd, numeric(1))
  scale[scale == 0] <- 1
  structure(list(fit = fit, cutoff = stats::median(risk), genes = genes,
    center = vapply(x, mean, numeric(1)), scale = scale, time_unit = task$time_unit),
    class = "stormbreaker_risk_reference")
}

#' Draw linked expression, risk, outcome and Kaplan-Meier panels
#' @param reference Result of surv_risk_reference, locked on development data.
#' @param task Development or validation survival task, with the same time unit.
#' @param genes Optional subset of selected genes for the heatmap. NULL uses all.
#' @return List with plots, a combined plot, patient table in risk order,
#'   standardized expression matrix, KM estimates, at-risk table, log-rank P
#'   and grouping status. `logrank_status` is success, single_group, no_events
#'   or failed. `logrank_error` retains the original error on failure, otherwise
#'   it is NA. Failed tests warn without discarding usable plots.
#'   KM curves include pointwise 95% log confidence bands;
#'   cox_group contains the High-versus-Low Cox HR, Wald CI and P value.
#'   Patient-level results are sensitive; review sharing.
#' @details Each cohort is sorted by increasing saved-model risk, breaking ties
#'   by input row order. Expression, risk and outcome panels share this order.
#'   Grouping and gene scaling never use validation outcomes. KM and log-rank
#'   describe the supplied cohort; apparent development separation is not
#'   independent validation. Log-rank P is returned, not printed as tiny text.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("patchwork", quietly = TRUE)) {
#'   d <- sb_example_survival(60, seed = 42)
#'   task <- sb_surv_task(d, "patient_id", "time", "event",
#'     c("gene_a", "gene_b"), "month")
#'   fit <- sb_fit_final(task, surv_pipeline("Cox"))
#'   reference <- surv_risk_reference(fit, task)
#'   panels <- surv_risk_plots(reference, task)
#'   panels[c("logrank_p", "logrank_status", "logrank_error")]
#' }
#' @export
surv_risk_plots <- function(reference, task, genes = NULL) {
  sb_require_curve_plotting()
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Install patchwork.")
  if (!inherits(reference, "stormbreaker_risk_reference")) stop("Use surv_risk_reference first.")
  task <- sb_risk_plot_task(task)
  if (!identical(task$time_unit, reference$time_unit)) stop("Time units must match the development reference.")
  if (is.null(genes)) genes <- reference$genes
  if (!is.character(genes) || !length(genes) || anyNA(genes) || anyDuplicated(genes) || !all(genes %in% reference$genes)) stop("Choose unique selected genes.")
  risk <- sb_risk_predict(reference$fit, task)
  d <- data.frame(patient_id = as.character(task$data[[task$id_col]]), risk = risk,
    time = task$data[[task$time_col]], event = task$data[[task$event_col]],
    group = factor(ifelse(risk > reference$cutoff, "High", "Low"), levels = c("Low", "High")))
  ord <- order(d$risk, seq_len(nrow(d)))
  d <- d[ord, ]; d$position <- seq_len(nrow(d))
  d$outcome <- factor(ifelse(d$event == 1, "Event", "Censored"), levels = c("Censored", "Event"))
  z <- sweep(as.matrix(task$data[ord, genes, drop = FALSE]), 2, reference$center[genes], "-")
  z <- sweep(z, 2, reference$scale[genes], "/")
  rownames(z) <- d$patient_id
  heat <- data.frame(position = rep(d$position, times = length(genes)),
    gene = factor(rep(genes, each = nrow(d)), levels = rev(genes)), value = as.vector(z))
  position <- gene <- value <- group <- time <- outcome <- survival <- n_risk <- lower <- upper <- NULL
  colours <- c(Low = "#2878B5", High = "#D45064")
  limits <- c(.5, nrow(d) + .5)
  common_x <- ggplot2::scale_x_continuous(limits = limits, expand = c(0, 0))
  plots <- list(
    expression = ggplot2::ggplot(heat, ggplot2::aes(position, gene, fill = value)) +
      ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "#2878B5", mid = "white", high = "#D45064", name = "Z-score") +
      common_x + ggplot2::labs(x = NULL, y = NULL) + sb_viz_theme(),
    risk = ggplot2::ggplot(d, ggplot2::aes(position, risk, colour = group)) +
      ggplot2::geom_point(size = 1.5) + ggplot2::geom_hline(yintercept = reference$cutoff, linetype = 2) +
      common_x + ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
      ggplot2::labs(x = "Patient order", y = "Risk score", colour = "Group") + sb_viz_theme(),
    outcome = ggplot2::ggplot(d, ggplot2::aes(position, time, colour = outcome)) +
      ggplot2::geom_point(size = 1.5) + common_x +
      ggplot2::scale_colour_manual(values = c(Censored = "#2878B5", Event = "#D45064"), drop = FALSE) +
      ggplot2::labs(x = "Patient order", y = paste0("Follow-up (", task$time_unit, ")"), colour = "Outcome") + sb_viz_theme())
  km <- survival::survfit(survival::Surv(time, event) ~ group, data = d)
  sm <- summary(km, censored = TRUE)
  groups <- levels(droplevels(d$group))
  kmdata <- data.frame(time = sm$time, survival = sm$surv, lower = sm$lower, upper = sm$upper,
    group = if (is.null(sm$strata)) rep(groups, length(sm$time)) else sub("^group=", "", as.character(sm$strata)))
  kmdata <- rbind(data.frame(time = 0, survival = 1, lower = 1, upper = 1, group = groups), kmdata)
  kmdata <- kmdata[order(match(kmdata$group, groups), kmdata$time), ]
  ticks <- pretty(c(0, max(d$time)), n = 5)
  ticks <- sort(unique(c(0, ticks[ticks >= 0 & ticks <= max(d$time)])))
  at_risk <- expand.grid(time = ticks, group = groups, stringsAsFactors = FALSE)
  at_risk$n_risk <- mapply(function(t, g) sum(d$time >= t & d$group == g), at_risk$time, at_risk$group)
  band <- do.call(rbind, lapply(split(kmdata, kmdata$group), function(z) {
    n <- nrow(z)
    if (n < 2L) return(z)
    out <- z[rep(seq_len(n), each = 2)[- (2 * n)], ]
    out$time <- c(z$time[1], rep(z$time[-1], each = 2))
    out
  }))
  plots$km <- ggplot2::ggplot(kmdata, ggplot2::aes(time, survival, colour = group)) +
    ggplot2::geom_ribbon(data = band, ggplot2::aes(ymin = lower, ymax = upper, fill = group),
      alpha = .18, colour = NA, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::geom_step() + ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_x_continuous(limits = c(0, max(d$time)), breaks = ticks) +
    ggplot2::labs(x = paste0("Time (", task$time_unit, ")"), y = "Survival probability", colour = "Group") + sb_viz_theme()
  plots$at_risk <- ggplot2::ggplot(at_risk, ggplot2::aes(time, group, label = n_risk, colour = group)) +
    ggplot2::geom_text(show.legend = FALSE) + ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_x_continuous(limits = c(0, max(d$time)), breaks = ticks) +
    ggplot2::labs(x = NULL, y = "At risk") + sb_viz_theme()
  logrank <- sb_risk_logrank(d)
  list(plots = plots, plot = patchwork::wrap_plots(plots, ncol = 1, heights = c(2, 1, 1, 2, .6)),
    patients = d, expression = z, km = km, at_risk = at_risk,
    logrank_p = logrank$p, logrank_status = logrank$status, logrank_error = logrank$error,
    cox_group = sb_risk_cox_estimate(d),
    cutoff = reference$cutoff, status = if (length(groups) == 2) "two_groups" else "single_group")
}

sb_risk_logrank <- function(d) {
  if (length(unique(d$group)) != 2L) {
    return(list(p = NA_real_, status = "single_group", error = NA_character_))
  }
  if (!any(d$event == 1)) {
    return(list(p = NA_real_, status = "no_events", error = NA_character_))
  }
  tryCatch({
    test <- survival::survdiff(survival::Surv(time, event) ~ group, data = d)
    p <- stats::pchisq(test$chisq, df = 1, lower.tail = FALSE)
    if (length(p) != 1L || !is.finite(p) || p < 0 || p > 1) {
      stop("The log-rank test did not return a finite P value in [0, 1].")
    }
    list(p = p, status = "success", error = NA_character_)
  }, error = function(e) {
    message <- conditionMessage(e)
    warning(paste0("Log-rank test failed: ", message), call. = FALSE)
    list(p = NA_real_, status = "failed", error = message)
  })
}
