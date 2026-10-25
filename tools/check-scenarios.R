# Load the source package with pkgload::load_all() or load an isolated installation first.
library(stormbreaker)
scenario_files <- c("01-binary.R", "02-compare.R", "03-external.R",
  "04-survival.R", "05-diagnostics.R", "06-shap.R")
scenario_results <- list()
scenario_checks <- list()
for (file in scenario_files) {
  env <- new.env(parent = globalenv())
  elapsed <- system.time(sys.source(system.file("examples", "scenarios", file,
    package = "stormbreaker", mustWork = TRUE), envir = env))[["elapsed"]]
  result <- env$example_result
  stopifnot(is.list(result), length(result$figures) > 0L)
  if (!is.null(result$fit)) stopifnot(identical(result$fit$status, "success"))
  if (!is.null(result$comparison)) {
    stopifnot(all(result$comparison$comparison$status == "success"))
  }
  if (!is.null(result$cv)) stopifnot(identical(result$cv$status, "success"))
  if (!is.null(result$evaluation)) stopifnot(identical(result$evaluation$status, "success"))
  if (!is.null(result$finished)) {
    stopifnot(all(vapply(result$finished$evaluations,
      \(x) identical(x$status, "success"), logical(1))))
  }
  for (figure in result$figures) {
    stopifnot(inherits(figure, "ggplot"))
    invisible(ggplot2::ggplot_build(figure))
  }
  scenario_results[[file]] <- result
  scenario_checks[[file]] <- data.frame(scenario = file, status = "success",
    figures = length(result$figures), elapsed_seconds = elapsed)
  cat(file, "success; figures", length(result$figures), "seconds", elapsed, "\n")
}
scenario_checks <- do.call(rbind, scenario_checks)
rownames(scenario_checks) <- NULL
