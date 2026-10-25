#' Export a completed workflow into a new directory
#' @param x Result from [sb_complete_run()].
#' @param directory A new output directory. Existing paths are never overwritten.
#' @return Invisible absolute output directory. Tables use CSV, figures use PDF.
#'   Cohort/plot indexes map numbered files to their original names. Patient
#'   predictions are identifiable data; store this directory appropriately.
#'   Original model predictions and evaluated predictions are separate files.
#'   Cohort processing records and original metrics accompany evaluations.
#' @examples
#' development <- bi_example(40, seed = 42)
#' run <- bi_run(development, outcome = "outcome", id = "patient_id",
#'   predictors = c("gene_a", "gene_b"), models = "Logistic", folds = 2,
#'   seed = 42, plots = FALSE)
#' finished <- finish_run(run, "Logistic", plots = FALSE)
#' destination <- tempfile("stormbreaker-export-")
#' export_run(finished, destination)
#' list.files(destination)
#' @export
sb_export_run <- function(x, directory) {
  if (file.exists(directory)) stop("Choose a new output directory; existing paths are not overwritten.")
  if (!dir.create(directory, recursive = TRUE)) stop("Could not create output directory.")
  write <- function(d, name) utils::write.csv(d, file.path(directory, name), row.names = FALSE)
  write(x$comparison, "model-comparison.csv")
  write(x$genes, "genes-by-stage.csv")
  index <- data.frame(
    cohort = as.character(names(x$evaluations)), file = sprintf("cohort-%03d", seq_along(x$evaluations)),
    status = vapply(x$evaluations, function(z) z$status, character(1))
  )
  write(index, "cohorts.csv")
  if (length(x$risk_groups)) {
    write(do.call(rbind, lapply(names(x$risk_groups), function(name) {
      data.frame(cohort = name, x$risk_groups[[name]], row.names = NULL)
    })), "risk-groups.csv")
  }
  for (i in seq_len(nrow(index))) {
    name <- index$cohort[i]
    write(x$predictions[[name]], paste0(index$file[i], "-predictions.csv"))
    evaluation <- x$evaluations[[name]]
    if (is.data.frame(evaluation$scored_predictions)) {
      write(evaluation$scored_predictions, paste0(index$file[i], "-evaluated-predictions.csv"))
    }
    if (is.data.frame(evaluation$direction_correction)) {
      write(evaluation$direction_correction, paste0(index$file[i], "-direction-correction.csv"))
      original <- evaluation$original_metrics
      original_tab <- do.call(rbind, lapply(names(original), function(m) {
        z <- original[[m]]
        if (!is.list(z) || is.null(z$estimate) || length(z$estimate) != 1L) return(NULL)
        data.frame(metric = m, estimate = z$estimate,
          status = if (is.null(z$status)) NA_character_ else z$status)
      }))
      if (!is.null(original_tab)) write(original_tab, paste0(index$file[i], "-original-metrics.csv"))
    }
    metrics <- x$metrics[[name]]
    tab <- do.call(rbind, lapply(names(metrics), function(m) {
      z <- metrics[[m]]
      if (!is.list(z) || is.null(z$estimate) || length(z$estimate) != 1L) {
        return(NULL)
      }
      data.frame(
        metric = m, estimate = z$estimate,
        status = if (is.null(z$status)) NA_character_ else z$status
      )
    }))
    if (!is.null(tab)) write(tab, paste0(index$file[i], "-metrics.csv"))
  }
  figures <- sb_plot_pages(sb_export_plot_list(x$plots))
  plot_index <- data.frame(plot = as.character(names(figures)), file = sprintf("figure-%03d.pdf", seq_along(figures)))
  for (i in seq_len(nrow(plot_index))) {
    size <- attr(figures[[i]], "stormbreaker_suggested_size")
    if (is.null(size)) size <- c(width = 11, height = if (inherits(figures[[i]], "patchwork")) 12 else 8.5)
    ggplot2::ggsave(file.path(directory, plot_index$file[i]),
      figures[[i]],
      device = grDevices::cairo_pdf, width = size[[1]], height = size[[2]]
    )
  }
  write(plot_index, "figures.csv")
  notes <- vapply(figures, function(p) {
    paste(unique(c(
      attr(p, "stormbreaker_notes"),
      p$labels$caption
    )), collapse = "\n")
  }, character(1))
  write(data.frame(plot = as.character(names(figures)), note = notes), "figure-notes.csv")
  write(sb_export_issue_rows(x$issues), "issues.csv")
  saveRDS(x$issues, file.path(directory, "issues.rds"))
  if (!is.null(x$explanation)) saveRDS(x$explanation, file.path(directory, "explanation.rds"))
  invisible(normalizePath(directory, winslash = "/"))
}

sb_export_plot_list <- function(plots, prefix = "") {
  if (inherits(plots, "ggplot")) {
    return(stats::setNames(list(plots), prefix))
  }
  collected <- list()
  if (is.list(plots)) {
    for (name in names(plots)) {
      collected <- c(collected, sb_export_plot_list(plots[[name]], paste(prefix, name, sep = "/")))
    }
  }
  collected
}

sb_export_issue_rows <- function(issues, scope = "") {
  empty <- data.frame(scope = character(), message = character())
  if (!length(issues)) {
    return(empty)
  }
  if (is.character(issues)) {
    return(data.frame(scope = scope, message = paste(issues, collapse = "; ")))
  }
  if (is.list(issues) && !is.null(issues$status)) {
    return(data.frame(scope = scope, message = paste(c(issues$status, issues$error_message), collapse = "; ")))
  }
  if (!is.list(issues) || is.null(names(issues))) {
    return(empty)
  }
  rows <- lapply(names(issues), function(name) {
    sb_export_issue_rows(issues[[name]], paste(scope, name, sep = "/"))
  })
  do.call(rbind, rows)
}
