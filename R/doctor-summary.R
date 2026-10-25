#' @export
print.stormbreaker_run <- function(x, ...) {
  comparison <- x$comparison$comparison
  successful <- sum(comparison$status == "success")
  cat("<stormbreaker:", x$task_type, "comparison>\n")
  cat(
    "Models:", nrow(comparison), "| CV successful:", successful,
    "| Other:", nrow(comparison) - successful, "\n"
  )
  cat("Final fits retained:", length(x$details), "| Issues:", nrow(x$issues), "\n")
  cat("Model names: names(run$pipelines)\n")
  cat("Full comparison: run$comparison$comparison | Issues: run$issues\n")
  cat("Next: choose a pipeline with sb_complete_run().\n")
  invisible(x)
}

#' @export
print.stormbreaker_completed_run <- function(x, ...) {
  successful <- sum(vapply(x$evaluations, \(evaluation) evaluation$status == "success", logical(1)))
  cat("<stormbreaker:", x$task_type, "completed run>\n")
  cat("Pipeline:", x$pipeline, "| Fit:", x$fit$status, "\n")
  cat(
    "External cohorts evaluated:", length(x$evaluations),
    "| Successful:", successful, "\n"
  )
  cat("Results: $genes, $predictions, $metrics, $plots\n")
  cat("Unavailable results and failures: $issues\n")
  cat("Export: sb_export_run(finished, directory = \"results-new\")\n")
  invisible(x)
}
