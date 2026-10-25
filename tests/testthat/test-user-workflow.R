workflow_test_environment <- function() {
  env <- new.env(parent = globalenv())
  source(system.file("examples", "user-workflow.R", package = "stormbreaker"),
    local = env)
  env
}

test_that("CSV helper preserves IDs and requires explicit model confirmation", {
  helper <- workflow_test_environment()
  directory <- withr::local_tempdir()
  path <- file.path(directory, "id.csv")
  writeLines(c("patient_id,gene_a", "00123,2.5"), path)
  expect_identical(helper$workflow_read_csv(path)$patient_id, "00123")
  expect_error(helper$workflow_finalize("does-not-exist.rds"), "Explicit")
  expect_error(helper$workflow_new_directory(directory), "new output directory")
})

test_that("both workflow templates complete outcome-blind prediction and evaluation", {
  helper <- workflow_test_environment()
  directory <- withr::local_tempdir()
  input <- file.path(directory, "input")
  helper$workflow_make_templates(input)
  expect_length(list.files(input), 6L)
  for (kind in c("binary", "survival")) {
    new_data <- helper$workflow_read_csv(file.path(input,
      paste0(kind, "_external_predictors.csv")))
    expect_false(any(c("outcome", "time", "event") %in% names(new_data)))
    development_dir <- file.path(directory, paste0(kind, "_development"))
    state <- suppressMessages(helper$workflow_develop(file.path(input,
      paste0(kind, "_development.csv")), development_dir, kind))
    expect_true(state$suggested_candidate %in% names(state$pipelines))
    expect_equal(nrow(state$comparison), 2L)
    for (stem in c("gene_selection", "gene_nonzero")) {
      ranking_file <- file.path(development_dir, paste0(stem, "_ranking.csv"))
      expect_identical(file.exists(ranking_file), TRUE)
      ranking <- read.csv(ranking_file)
      expect_setequal(ranking$feature, state$predictors)
      expect_identical(file.exists(file.path(development_dir,
        paste0(stem, "_counting_groups.csv"))), TRUE)
    }
    expect_true(all(state$comparison$eligible_for_primary_rule))
    expect_true(all(vapply(state$cv$results, function(x)
      identical(x$provenance$fold_assignment, state$cv$fold_assignment), logical(1))))
    final_dir <- file.path(directory, paste0(kind, "_final"))
    result <- helper$workflow_finalize(file.path(development_dir, "development.rds"),
      state$suggested_candidate, "Synthetic test of the declared rule.",
      file.path(input, paste0(kind, "_external_predictors.csv")),
      file.path(input, paste0(kind, "_external_outcomes.csv")), final_dir,
      render_report = FALSE)
    expect_identical(result$fit$status, "success")
    expect_identical(result$prediction$status, "success")
    expect_identical(result$evaluation$status, "success")
    expect_equal(nrow(result$prediction$predictions), 40L)
    reloaded <- readRDS(file.path(final_dir, "final_model.rds"))
    expect_identical(sb_predict_final(reloaded, new_data, "patient_id"), result$prediction)
    expect_true(all(file.exists(file.path(final_dir,
      c("selection_decision.rds", "final_model.rds", "external_predictions.csv",
        "external_predictions.rds", "external_evaluation.rds")))))
    # A missing outcome file fails only after outcome-blind predictions are saved.
    blind_dir <- file.path(directory, paste0(kind, "_blind"))
    expect_error(suppressWarnings(helper$workflow_finalize(
      file.path(development_dir, "development.rds"), state$suggested_candidate,
      "Synthetic outcome isolation test.",
      file.path(input, paste0(kind, "_external_predictors.csv")),
      file.path(directory, "absent_outcomes.csv"), blind_dir,
      render_report = FALSE)))
    expect_identical(readRDS(file.path(blind_dir, "external_predictions.rds")),
      result$prediction)
  }
})
