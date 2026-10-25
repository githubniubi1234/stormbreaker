test_that("two workers preserve repeated binary comparisons and caller RNG", {
  d <- sb_example_binary(60L, 1101L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  pp <- list(sb_pipeline("none", "logistic"), sb_pipeline("none", "gaussian_naive_bayes"))
  withr::local_seed(119L)
  before <- .Random.seed
  serial <- sb_compare_repeated_cv(task, pp, 3L, c(19L, 20L, 21L))
  actual <- sb_compare_repeated_cv(task, pp, 3L, c(19L, 20L, 21L), workers = 2L)
  expect_identical(actual, serial)
  expect_identical(.Random.seed, before)
  plans <- lapply(serial$repetitions, function(x) setNames(x$fold_assignment$fold, x$fold_assignment$patient_id))
  expect_identical(sb_compare_repeated_cv(task, pp, plans, NULL, workers = 2L),
    sb_compare_repeated_cv(task, pp, plans, NULL))
})

test_that("two workers preserve survival horizon and failed routes", {
  d <- sb_example_survival(72L, 11032L)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"),
    "month", prediction_horizon = 7)
  pp <- list(sb_pipeline("none", "coxph"), sb_pipeline("univ_cox_top_k", "coxph", selector_n_features = 2L))
  expect_identical(sb_compare_repeated_cv(task, pp, 3L, c(19L, 20L), workers = 2L),
    sb_compare_repeated_cv(task, pp, 3L, c(19L, 20L)))
  bad <- task
  bad$data$gene_b <- bad$data$gene_a
  failed <- sb_compare_repeated_cv(bad, pp[1], 3L, c(19L, 20L))
  expect_true(any(failed$comparison$status != "success"))
  expect_identical(sb_compare_repeated_cv(bad, pp[1], 3L, c(19L, 20L), workers = 2L),
    failed)
})

test_that("parallel repeated CV preserves absent caller RNG", {
  d <- sb_example_binary(48L, 1101L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  withr::local_preserve_seed()
  if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  result <- sb_compare_repeated_cv(task, list(sb_pipeline("none", "logistic")),
    3L, c(19L, 20L), workers = 2L)
  expect_identical(result$n_repeats, 2L)
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
})

test_that("worker count is bounded and worker errors keep their class", {
  d <- sb_example_binary(48L, 12L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  pp <- list(sb_pipeline("none", "logistic"))
  for (workers in list(0L, 3L, NA_real_, c(1L, 2L)))
    expect_error(sb_compare_repeated_cv(task, pp, workers = workers), class = "stormbreaker_bad_workers")
  connections <- nrow(showConnections(all = TRUE))
  expect_error(sb_compare_repeated_cv(task, list("invalid"), seeds = 1:2, workers = 2L),
    class = "stormbreaker_bad_pipeline_collection")
  expect_identical(nrow(showConnections(all = TRUE)), connections)
  expect_error(sb_compare_repeated_cv(task, list(sb_pipeline("none", "adaboost_binary")),
    seeds = 1:2, workers = 2L), class = "stormbreaker_parallel_sequential_rng")
})

test_that("parallel stack repetitions keep their complete results", {
  skip_if_not_installed("glmnet")
  d <- sb_example_binary(72L, 1101L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  stack <- sb_stack(list(sb_pipeline("none", "logistic"),
    sb_pipeline("none", "gaussian_naive_bayes")), inner_folds = 3L, inner_seed = 41L)
  expect_identical(sb_compare_repeated_cv(task, stack, 3L, c(19L, 20L), workers = 2L),
    sb_compare_repeated_cv(task, stack, 3L, c(19L, 20L)))
  d <- sb_example_survival(72L, 11032L)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"),
    "month", prediction_horizon = 7)
  stack <- sb_stack(list(sb_pipeline("none", "coxph"),
    sb_pipeline("univ_cox_top_k", "coxph", selector_n_features = 2L)),
    inner_folds = 3L, inner_seed = 41L)
  expect_identical(sb_compare_repeated_cv(task, stack, 3L, c(19L, 20L), workers = 2L),
    sb_compare_repeated_cv(task, stack, 3L, c(19L, 20L)))
})
