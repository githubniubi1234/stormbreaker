test_that("binary comparison reuses fold-local selectors with identical results", {
  d <- sb_example_binary(120, 37)
  task <- sb_binary_task(d, "patient_id", "outcome",
    c("gene_a", "gene_b", "gene_c"))
  pipelines <- list(
    logistic = sb_pipeline("univ_binary_auc_top_k", "logistic", selector_n_features = 2L),
    nb = sb_pipeline("univ_binary_auc_top_k", "gaussian_naive_bayes", selector_n_features = 2L)
  )
  original <- sb_fit_runtime_selector_uncached
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    if (identical(selector$id, "univ_binary_auc_top_k")) calls <<- calls + 1L
    original(selector, x, endpoint)
  })
  plain <- sb_compare_cv(task, unname(pipelines), folds = 3L, seed = 19L,
    reuse_selectors = FALSE, progress = FALSE)
  plain_calls <- calls
  calls <- 0L
  shared <- sb_compare_cv(task, unname(pipelines), folds = 3L, seed = 19L,
    reuse_selectors = TRUE, progress = FALSE)
  expect_identical(shared, plain)
  expect_identical(shared$comparison$status, rep("success", 2L))
  expect_lt(calls, plain_calls)
  expect_identical(calls, 3L)
})

test_that("binary completion exports usable results and preserves cohort failures", {
  d <- sb_example_binary(120, 38)
  run <- bi_run(d, "outcome", c("gene_a", "gene_b", "gene_c"),
    "patient_id", models = c("Logistic", "NaiveBayes"), folds = 3L,
    plots = FALSE)
  external <- sb_example_binary(60, 39)
  external$patient_id <- paste0("external_", external$patient_id)
  result <- sb_complete_run(run, "Logistic",
    external = list(valid = external, incomplete = external["patient_id"]),
    plots = FALSE)
  expect_identical(result$evaluations$valid$status, "success")
  expect_identical(nrow(result$predictions$valid), 60L)
  expect_identical(result$fit$status, "success")
  expect_type(result$issues$external$incomplete, "character")
  target <- tempfile("binary-alignment-export-")
  sb_export_run(result, target)
  scores <- read.csv(file.path(target, "cohort-001-predictions.csv"))
  expect_identical(nrow(scores), 60L)
  expect_setequal(read.csv(file.path(target, "genes-by-stage.csv"))$gene,
    c("gene_a", "gene_b", "gene_c"))
})
