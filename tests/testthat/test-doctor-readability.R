readability_binary_run <- function(plots = FALSE) {
  data <- sb_example_binary(80, seed = 11)
  bi_run(data, "outcome", c("gene_a", "gene_b"), "patient_id",
    models = "Logistic", folds = 3, plots = plots
  )
}

test_that("an optional plot failure does not discard a successful model run", {
  local_mocked_bindings(sb_plots = function(...) stop("Plot device unavailable"))

  run <- readability_binary_run(plots = TRUE)

  expect_identical(run$comparison$status, "success")
  expect_identical(run$details$Logistic$fit$status, "success")
  expect_contains(run$issues$reason, "Plot device unavailable")
  expect_contains(run$issues$stage, "plot algorithm")
})

test_that("completion reuses saved algorithm figures without calling the plotter", {
  run <- readability_binary_run(plots = TRUE)
  local_mocked_bindings(sb_plots = function(...) stop("Unexpected second plot call"))

  completed <- sb_complete_run(run, "Logistic")

  expect_identical(completed$fit, run$details$Logistic$fit)
  expect_identical(completed$plots$algorithm, run$details$Logistic$figures$plots)
  expect_identical(completed$issues$algorithm, run$details$Logistic$figures$unavailable)
})

test_that("printing a comparison shows a short summary instead of patient data", {
  run <- readability_binary_run()
  run$development_task$data$patient_id[1] <- "PRIVATE_PATIENT_731"

  printed <- capture.output(returned <- withVisible(print(run)))

  expect_s3_class(run, "stormbreaker_run")
  expect_identical(returned$value, run)
  expect_identical(returned$visible, FALSE)
  expect_lte(length(printed), 12L)
  expect_match(paste(printed, collapse = "\n"), "sb_complete_run", fixed = TRUE)
  expect_identical(any(grepl("PRIVATE_PATIENT_731", printed, fixed = TRUE)), FALSE)
})

test_that("printing a completed run leaves its detailed results untouched", {
  run <- readability_binary_run()
  external <- run$development_task$data
  external$patient_id <- paste0("PRIVATE_PATIENT_952_", seq_len(nrow(external)))
  completed <- sb_complete_run(run, "Logistic", list(validation = external), plots = FALSE)

  printed <- capture.output(returned <- withVisible(print(completed)))

  expect_s3_class(completed, "stormbreaker_completed_run")
  expect_identical(returned$value, completed)
  expect_identical(returned$visible, FALSE)
  expect_lte(length(printed), 12L)
  expect_match(paste(printed, collapse = "\n"), "Logistic", fixed = TRUE)
  expect_match(paste(printed, collapse = "\n"), "sb_export_run", fixed = TRUE)
  expect_identical(any(grepl("PRIVATE_PATIENT_952", printed, fixed = TRUE)), FALSE)
})

test_that("exported indexes retain their column names when no cohorts or plots exist", {
  completed <- sb_complete_run(readability_binary_run(), "Logistic", plots = FALSE)
  directory <- file.path(withr::local_tempdir(), "results")

  sb_export_run(completed, directory)

  cohorts <- read.csv(file.path(directory, "cohorts.csv"))
  figures <- read.csv(file.path(directory, "figures.csv"))
  notes <- read.csv(file.path(directory, "figure-notes.csv"))
  expect_named(cohorts, c("cohort", "file", "status"))
  expect_named(figures, c("plot", "file"))
  expect_named(notes, c("plot", "note"))
  expect_equal(nrow(cohorts), 0L)
  expect_equal(nrow(figures), 0L)
  expect_equal(nrow(notes), 0L)
})

test_that("skipping completion plots neither redraws nor discards a saved fit", {
  run <- readability_binary_run(plots = TRUE)
  local_mocked_bindings(
    sb_plots = function(...) stop("Unexpected plotting"),
    sb_fit_final = function(...) stop("Unexpected refit")
  )

  completed <- sb_complete_run(run, "Logistic", plots = FALSE)

  expect_length(completed$plots, 0L)
  expect_identical(completed$fit, run$details$Logistic$fit)
  expect_setequal(completed$genes$gene, c("gene_a", "gene_b"))
})

test_that("a failed final fit remains visible in the completed summary", {
  run <- readability_binary_run()
  local_mocked_bindings(sb_fit_final = function(...) {
    list(status = "fit_failed", error_message = "No usable features")
  })

  completed <- sb_complete_run(run, "Logistic", plots = FALSE)

  expect_s3_class(completed, "stormbreaker_completed_run")
  expect_identical(completed$issues$fit$error_message, "No usable features")
  expect_match(paste(capture.output(print(completed)), collapse = "\n"),
    "fit_failed",
    fixed = TRUE
  )
  expect_length(completed$predictions, 0L)
})

test_that("export preserves the original reasons for nested failures", {
  completed <- sb_complete_run(readability_binary_run(), "Logistic", plots = FALSE)
  completed$issues <- list(
    external = list(broken = "Missing outcome column"),
    algorithm = list(path = "Saved path unavailable"),
    fit = list(status = "fit_failed", error_message = "Original engine error")
  )
  directory <- file.path(withr::local_tempdir(), "results")

  sb_export_run(completed, directory)

  issues <- read.csv(file.path(directory, "issues.csv"))
  expect_identical(issues$scope, c("/external/broken", "/algorithm/path", "/fit"))
  expect_identical(issues$message, c(
    "Missing outcome column", "Saved path unavailable", "fit_failed; Original engine error"
  ))
  expect_identical(readRDS(file.path(directory, "issues.rds")), completed$issues)
})

test_that("export never replaces an existing directory or its contents", {
  directory <- withr::local_tempdir()
  marker <- file.path(directory, "keep.txt")
  writeLines("Original contents", marker)

  expect_error(sb_export_run(list(), directory), "existing paths are not overwritten", fixed = TRUE)

  expect_identical(readLines(marker), "Original contents")
  expect_identical(list.files(directory), "keep.txt")
})

test_that("survival completion preserves risk groups while reusing algorithm plots", {
  data <- sb_example_survival(80, seed = 21)
  run <- surv_run(data, "time", "event", c("gene_a", "gene_b"), "patient_id",
    models = "Cox", folds = 3, plots = TRUE
  )
  external <- data
  external$patient_id <- paste0("external-", external$patient_id)
  local_mocked_bindings(sb_plots = function(...) stop("Unexpected second plot call"))

  completed <- sb_complete_run(run, "Cox", list(validation = external))

  expect_s3_class(run, "stormbreaker_run")
  expect_s3_class(completed, "stormbreaker_completed_run")
  expect_identical(completed$plots$algorithm, run$details$Cox$figures$plots)
  expect_identical(completed$evaluations$validation$status, "success")
  expect_equal(completed$risk_groups$validation$n_low + completed$risk_groups$validation$n_high, 80L)
  expect_match(paste(capture.output(print(completed)), collapse = "\n"), "survival", fixed = TRUE)
})
