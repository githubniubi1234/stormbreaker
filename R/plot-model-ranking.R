#' Rank models across explicitly selected cohorts and draw a Top-N panel
#' @param performance Data frame with model, cohort and auc columns; one row per
#'   model/cohort. Failed or missing AUCs must be NA, never zero.
#' @param ranking_cohorts Explicit cohort names defining an unweighted mean AUC.
#'   Use development or designated validation cohorts, not a held-out final test.
#' @param top_n Maximum number of complete models to display; Inf shows all.
#' @param metric "auc" (default) or "c_index"; supply the matching numeric
#'   column. All C-index cells must use the same estimator (e.g. Harrell),
#'   endpoint and comparable evaluation scope; do not mix Harrell and Uno.
#' @return List with full ranking, displayed data, heatmap, mean-AUC bar plot and
#'   a combined patchwork plot. Models missing any ranking cohort are unranked;
#'   other missing cells remain missing. Ties retain input model order.
#' @export
sb_plot_model_ranking <- function(performance, ranking_cohorts, top_n = 20L,
                                   metric = c("auc", "c_index")) {
  metric <- match.arg(metric)
  if (metric == "c_index") {
    if (!is.data.frame(performance) || !"c_index" %in% names(performance)) stop("Supply a c_index column.")
    d <- performance[c("model", "cohort", "c_index")]
    names(d)[3] <- "auc"
    result <- sb_plot_model_ranking(d, ranking_cohorts, top_n, metric = "auc")
    names(result$ranking)[names(result$ranking) == "mean_auc"] <- "mean_c_index"
    names(result$data)[names(result$data) == "auc"] <- "c_index"
    names(result$heatmap$data)[names(result$heatmap$data) == "auc"] <- "c_index"
    names(result$bar$data)[names(result$bar$data) == "mean_auc"] <- "mean_c_index"
    c_index <- mean_c_index <- NULL
    result$heatmap <- result$heatmap + ggplot2::aes(fill = c_index) + ggplot2::labs(fill = "C-index")
    result$bar <- result$bar + ggplot2::aes(x = mean_c_index) + ggplot2::labs(x = "Mean C-index")
    result$bar$layers[[2]]$mapping$label <- ggplot2::aes(label = sprintf("%.3f", mean_c_index))$label
    result$plot <- patchwork::wrap_plots(result$heatmap, result$bar, widths = c(4, 1.5))
    result$metric <- metric
    return(result)
  }
  sb_require_curve_plotting()
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Install patchwork for the combined ranking plot.")
  d <- performance
  if (!is.data.frame(d) || !all(c("model", "cohort", "auc") %in% names(d)) || !nrow(d)) stop("Supply model, cohort and auc columns.")
  d$model <- as.character(d$model)
  d$cohort <- as.character(d$cohort)
  if (anyNA(d[c("model", "cohort")]) || any(!nzchar(d$model)) || any(!nzchar(d$cohort)) ||
      anyDuplicated(d[c("model", "cohort")]) || !is.numeric(d$auc) ||
      any(!is.na(d$auc) & (!is.finite(d$auc) | d$auc < 0 | d$auc > 1))) stop("Use unique model/cohort cells and AUC in [0,1] or NA.")
  if (!is.character(ranking_cohorts) || !length(ranking_cohorts) || anyNA(ranking_cohorts) ||
      anyDuplicated(ranking_cohorts) || !all(ranking_cohorts %in% d$cohort)) stop("Choose explicit unique ranking cohorts from the table.")
  if (length(top_n) != 1L || is.na(top_n) || top_n <= 0 || (is.finite(top_n) && top_n != floor(top_n))) stop("top_n must be a positive integer or Inf.")
  models <- unique(d$model)
  scores <- vapply(models, function(m) {
    z <- d[d$model == m, ]
    a <- z$auc[match(ranking_cohorts, z$cohort)]
    if (anyNA(a)) NA_real_ else mean(a)
  }, numeric(1))
  ranking <- data.frame(model = models, mean_auc = unname(scores), status = ifelse(is.na(scores), "incomplete", "ranked"))
  ranking <- ranking[order(-ranking$mean_auc, seq_along(models), na.last = TRUE), ]
  selected <- utils::head(ranking$model[ranking$status == "ranked"], if (is.infinite(top_n)) length(models) else top_n)
  if (!length(selected)) stop("No model has complete AUCs for every ranking cohort.")
  grid <- expand.grid(model = selected, cohort = unique(d$cohort), stringsAsFactors = FALSE)
  z <- merge(grid, d, by = c("model", "cohort"), all.x = TRUE, sort = FALSE)
  z$model <- factor(z$model, levels = rev(selected))
  z$cohort <- factor(z$cohort, levels = unique(d$cohort))
  z$label <- ifelse(is.na(z$auc), "", sprintf("%.3f", z$auc))
  bars <- ranking[match(selected, ranking$model), ]
  bars$model <- factor(bars$model, levels = rev(selected))
  model <- cohort <- auc <- label <- mean_auc <- NULL
  heatmap <- ggplot2::ggplot(z, ggplot2::aes(cohort, model, fill = auc)) +
    ggplot2::geom_tile(colour = "white") + ggplot2::geom_text(ggplot2::aes(label = label), size = 3) +
    ggplot2::scale_fill_gradient(low = "#eff3ff", high = "#6a51a3", limits = c(0, 1), na.value = "grey90") +
    ggplot2::labs(x = NULL, y = NULL, fill = "AUC") + sb_viz_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  bar <- ggplot2::ggplot(bars, ggplot2::aes(mean_auc, model)) +
    ggplot2::geom_col(fill = "#a63683", width = .85) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", mean_auc)), hjust = -.1, size = 3) +
    ggplot2::scale_x_continuous(limits = c(0, 1.12), breaks = c(0, .5, 1)) +
    ggplot2::labs(x = "Mean AUC", y = NULL) + sb_viz_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank())
  list(ranking = ranking, data = z, ranking_cohorts = ranking_cohorts,
    heatmap = heatmap, bar = bar, plot = patchwork::wrap_plots(heatmap, bar, widths = c(4, 1.5)))
}
