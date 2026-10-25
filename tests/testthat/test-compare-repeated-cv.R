test_that("repeated binary comparison preserves every seeded comparison", {
  d <- make_sb_binary_fixture(n = 48L, seed = 7101L)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  )
  seeds <- c(7102L, 7103L)

  observed <- sb_compare_repeated_cv(
    task,
    pipelines,
    folds = 3L,
    seeds = seeds
  )

  expect_s3_class(observed, "stormbreaker_repeated_cv_comparison")
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expect_identical(observed$n_patients, 48L)
  expect_identical(observed$n_repeats, 2L)
  expect_identical(names(observed$repetitions), c("repeat_1", "repeat_2"))
  expect_identical(
    observed$comparison$repeat_id,
    rep(1:2, each = length(pipelines))
  )
  expect_identical(
    observed$comparison$seed,
    rep(seeds, each = length(pipelines))
  )
  expect_identical(
    observed$comparison$fold_source,
    rep("generated_stratified", 4L)
  )
  expect_identical(
    observed$fold_assignments$repeat_id,
    rep(1:2, each = nrow(d))
  )

  for (repeat_id in seq_along(seeds)) {
    direct <- sb_compare_cv(
      task,
      pipelines,
      folds = 3L,
      seed = seeds[[repeat_id]]
    )
    expect_identical(observed$repetitions[[repeat_id]], direct)
    expected_assignment <- data.frame(
      repeat_id = repeat_id,
      direct$fold_assignment,
      stringsAsFactors = FALSE
    )
    actual_assignment <- observed$fold_assignments[
      observed$fold_assignments$repeat_id == repeat_id,
      ,
      drop = FALSE
    ]
    rownames(actual_assignment) <- NULL
    expect_identical(actual_assignment, expected_assignment)
  }

  expect_length(
    grep(
      "mean|standard|interval|rank|winner|best",
      names(observed),
      ignore.case = TRUE,
      value = TRUE
    ),
    0L
  )
  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(
    printed,
    "No cross-repeat mean, confidence interval, ranking, or winner",
    fixed = TRUE
  )
})

test_that("repeated survival comparison preserves explicit fold plans", {
  d <- make_sb_fixture(n = 36L, seed = 7201L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month"
  )
  pipelines <- list(
    sb_pipeline(),
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 2L
    )
  )
  first <- sb_make_stratified_folds(task, 3L, 7202L)
  second <- sb_make_stratified_folds(task, 3L, 7203L)
  fold_plans <- list(
    setNames(first$assignment, d$patient_id),
    rev(setNames(second$assignment, d$patient_id))
  )

  observed <- sb_compare_repeated_cv(
    task,
    pipelines,
    folds = fold_plans,
    seeds = NULL
  )

  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "survival")
  expect_identical(observed$n_repeats, 2L)
  expect_identical(observed$comparison$seed, rep(NA_integer_, 4L))
  expect_identical(
    observed$comparison$fold_source,
    rep("preassigned_patient_id", 4L)
  )

  for (repeat_id in seq_along(fold_plans)) {
    direct <- sb_compare_cv(
      task,
      pipelines,
      folds = fold_plans[[repeat_id]],
      seed = 999L
    )
    expect_identical(observed$repetitions[[repeat_id]], direct)
    expected_fold <- as.integer(fold_plans[[repeat_id]][match(
      d$patient_id,
      names(fold_plans[[repeat_id]])
    )])
    actual <- observed$fold_assignments[
      observed$fold_assignments$repeat_id == repeat_id,
      ,
      drop = FALSE
    ]
    expect_identical(actual$patient_id, d$patient_id)
    expect_identical(actual$fold, expected_fold)
  }
})

test_that("repeated comparison preserves incomplete pipeline evidence", {
  d <- make_sb_fixture(n = 48L, seed = 7251L)
  d$gene_b <- d$gene_a
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month"
  )
  pipelines <- list(
    sb_pipeline(
      selector = "univ_cox_top_k",
      selector_n_features = 1L
    ),
    sb_pipeline()
  )

  observed <- sb_compare_repeated_cv(
    task,
    pipelines,
    folds = 3L,
    seeds = c(7252L, 7253L)
  )

  expect_identical(observed$status, "incomplete")
  expect_identical(
    unname(vapply(
      observed$repetitions,
      `[[`,
      character(1L),
      "status"
    )),
    rep("incomplete", 2L)
  )
  expect_identical(
    observed$comparison$status,
    rep(c("success", "failed"), 2L)
  )
  expect_identical(
    observed$comparison$complete_oof,
    rep(c(TRUE, FALSE), 2L)
  )
})

test_that("repeated comparison keeps generated and explicit modes unambiguous", {
  d <- make_sb_binary_fixture(n = 24L, seed = 7301L)
  task <- sb_binary_task(
    d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(sb_pipeline(learner = "logistic"))
  generated <- sb_make_stratified_folds(task, 3L, 7302L)
  plan <- setNames(generated$assignment, d$patient_id)

  missing_seeds <- capture_sb_error(sb_compare_repeated_cv(
    task, pipelines, folds = 3L, seeds = NULL
  ))
  conflicting_seeds <- capture_sb_error(sb_compare_repeated_cv(
    task, pipelines, folds = list(plan), seeds = 7303L
  ))
  bare_plan <- capture_sb_error(sb_compare_repeated_cv(
    task, pipelines, folds = plan, seeds = NULL
  ))
  invalid_seed <- capture_sb_error(sb_compare_repeated_cv(
    task, pipelines, folds = 3L, seeds = c(7303L, NA_real_)
  ))
  empty_plans <- capture_sb_error(sb_compare_repeated_cv(
    task, pipelines, folds = list(), seeds = NULL
  ))
  scalar_plan <- capture_sb_error(sb_compare_repeated_cv(
    task, pipelines, folds = list(3L), seeds = NULL
  ))
  mixed_plans <- capture_sb_error(sb_compare_repeated_cv(
    task, pipelines, folds = list(plan, 3L), seeds = NULL
  ))

  expect_s3_class(missing_seeds, "stormbreaker_bad_repeat_seeds")
  expect_s3_class(
    conflicting_seeds,
    "stormbreaker_repeated_fold_seed_conflict"
  )
  expect_s3_class(bare_plan, "stormbreaker_bad_repeated_folds")
  expect_s3_class(invalid_seed, "stormbreaker_bad_repeat_seeds")
  expect_s3_class(empty_plans, "stormbreaker_bad_repeated_folds")
  expect_s3_class(scalar_plan, "stormbreaker_bad_repeated_folds")
  expect_s3_class(mixed_plans, "stormbreaker_bad_repeated_folds")
})

test_that("repeated comparison preserves complete generated stack runs", {
  skip_if_not_installed("glmnet")

  data <- make_sb_binary_fixture(n = 48L, seed = 7351L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 7352L)
  seeds <- c(7353L, 7354L)

  observed <- sb_compare_repeated_cv(
    task,
    stack,
    folds = 3L,
    seeds = seeds
  )

  expect_identical(observed$status, "success")
  expect_identical(observed$strategy_type, "stack")
  expect_identical(observed$n_repeats, 2L)
  expect_identical(
    observed$comparison$pipeline_id,
    rep(c(stack$pipeline_ids, stack$id), 2L)
  )
  for (repeat_id in seq_along(seeds)) {
    direct <- sb_run_stack_cv(
      task,
      stack,
      folds = 3L,
      seed = seeds[[repeat_id]]
    )
    expect_identical(observed$repetitions[[repeat_id]], direct)
    actual_assignment <- observed$fold_assignments[
      observed$fold_assignments$repeat_id == repeat_id,
      c("patient_id", "fold")
    ]
    rownames(actual_assignment) <- NULL
    expect_identical(
      actual_assignment,
      direct$provenance$fold_assignment
    )
  }
})

test_that("repeated stack comparison keeps explicit patient fold plans", {
  skip_if_not_installed("glmnet")

  data <- make_sb_binary_fixture(n = 42L, seed = 7361L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  stack <- sb_stack(list(
    sb_pipeline(learner = "logistic"),
    sb_pipeline(learner = "gaussian_naive_bayes")
  ), inner_folds = 3L, inner_seed = 7362L)
  generated <- sb_make_stratified_folds(task, 3L, 7363L)
  plan <- rev(stats::setNames(generated$assignment, data$patient_id))

  observed <- sb_compare_repeated_cv(
    task,
    stack,
    folds = list(plan),
    seeds = NULL
  )
  expected <- sb_run_stack_cv(task, stack, folds = plan, seed = 1L)

  expect_identical(observed$status, "success")
  expect_identical(observed$repetitions[[1L]], expected)
  expect_identical(observed$comparison$seed, rep(NA_integer_, 3L))
  expect_identical(
    observed$comparison$fold_source,
    rep("preassigned_patient_id", 3L)
  )
  expect_identical(
    observed$fold_assignments$patient_id,
    data$patient_id
  )
  expect_identical(
    observed$fold_assignments$fold,
    as.integer(plan[match(data$patient_id, names(plan))])
  )
})
