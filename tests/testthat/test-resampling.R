test_that("mixed-outcome patients stay together in stratified inner folds", {
  groups <- c(rep(paste0("paired", 1:7), each = 2),
              rep("large_negative", 9), paste0("negative", 1:6),
              rep("large_positive", 6), paste0("positive", 1:6))
  outcome <- c(rep(c(0L, 1L), 7), rep(0L, 15), rep(1L, 12))
  withr::local_seed(3749L)
  before <- .Random.seed
  observed <- sb_make_stratified_fold_id(
    outcome, 5L, 3748L, group_ids = groups
  )
  expect_identical(.Random.seed, before)
  expect_identical(observed, sb_make_stratified_fold_id(
    outcome, 5L, 3748L, group_ids = groups
  ))
  expect_identical(unname(vapply(split(observed, groups), function(x) {
    length(unique(x))
  }, integer(1))), rep(1L, length(unique(groups))))
  expect_equal(sort(unique(observed)), 1:5)
  expect_gt(min(table(observed, outcome)), 0)
  few_mixed_y <- c(rep(c(0L, 1L), 3L), rep(0L, 8L), rep(1L, 8L))
  few_mixed_groups <- c(rep(paste0("mixed", 1:3), each = 2L),
    rep(paste0("negative", 1:4), each = 2L),
    rep(paste0("positive", 1:4), each = 2L))
  few_mixed <- sb_make_stratified_fold_id(few_mixed_y, 5L, 3748L,
    group_ids = few_mixed_groups)
  expect_gt(min(table(few_mixed, few_mixed_y)), 0)
  expect_identical(unname(vapply(split(few_mixed, few_mixed_groups), function(x) {
    length(unique(x))
  }, integer(1))), rep(1L, 11L))
  unique_ids <- seq_along(outcome)
  expect_identical(
    sb_make_stratified_fold_id(outcome, 5L, 3748L, group_ids = unique_ids),
    sb_make_stratified_fold_id(outcome, 5L, 3748L)
  )
  error <- capture_sb_error(sb_make_stratified_fold_id(
    rep(c(0L, 1L), each = 10), 5L, 1L,
    group_ids = rep(c("only_negative", "only_positive"), each = 10)
  ))
  expect_s3_class(error, "stormbreaker_infeasible_inner_folds")
})

test_that("task grouping survives validation and protects preassigned folds", {
  d <- sb_example_binary(60L, seed = 3745L)
  d$group <- rep(paste0("patient", 1:30), each = 2)
  d$outcome <- rep(c(0L, 1L), 30)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"),
                         group_col = "group")
  expect_identical(sb_validate_binary_task_contract(task), task)
  subset <- sb_subset_task(task, 1:40)
  expect_identical(subset$group_col, "group")
  expect_identical(subset$data$group, d$group[1:40])
  folds <- sb_make_outer_folds(task, 5L, 3746L)
  expect_identical(folds$source, "generated_group_stratified")
  expect_identical(unname(vapply(split(folds$assignment, d$group), function(x) {
    length(unique(x))
  }, integer(1))), rep(1L, 30))
  supplied <- setNames(folds$assignment, d$patient_id)
  expect_identical(sb_make_outer_folds(task, rev(supplied), 1L)$assignment,
                   folds$assignment)
  supplied[2] <- supplied[1] %% 5L + 1L
  error <- capture_sb_error(sb_make_outer_folds(task, supplied, 1L))
  expect_s3_class(error, "stormbreaker_group_fold_overlap")
  d$group <- paste0("unique", seq_len(nrow(d)))
  unique_task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"),
                                group_col = "group")
  legacy_task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  expect_identical(sb_make_outer_folds(unique_task, 5L, 3746L),
                   sb_make_outer_folds(legacy_task, 5L, 3746L))
})

test_that("outer-fold dispatch preserves the scalar stratified path", {
  d <- make_sb_fixture(n = 30L, seed = 3701L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month"
  )
  expected <- sb_make_stratified_folds(task, folds = 3L, seed = 3702L)

  observed <- sb_make_outer_folds(task, folds = 3L, seed = 3702L)
  named_scalar <- sb_make_outer_folds(
    task,
    folds = c(requested = 3L),
    seed = 3702L
  )

  expect_identical(observed$assignment, expected$assignment)
  expect_identical(observed$splits, expected$splits)
  expect_identical(named_scalar$assignment, expected$assignment)
  expect_identical(named_scalar$splits, expected$splits)
  expect_identical(observed$source, "generated_stratified")
  expect_identical(observed$seed, 3702L)
  expect_identical(observed$n_folds, 3L)
})

test_that("preassigned outer folds align by patient ID rather than position", {
  d <- make_sb_fixture(n = 12L, seed = 3711L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month"
  )
  task_ids <- as.character(task$data[[task$id_col]])
  task_assignment <- rep(1:3, length.out = task$n)
  supplied_order <- c(12L, 2L, 9L, 1L, 7L, 4L, 10L, 3L, 6L, 8L, 5L, 11L)
  supplied <- setNames(
    task_assignment[supplied_order],
    task_ids[supplied_order]
  )

  observed <- sb_make_outer_folds(task, folds = supplied, seed = 999L)
  expected_assignment <- as.integer(supplied[match(task_ids, names(supplied))])

  expect_identical(observed$assignment, expected_assignment)
  expect_identical(observed$source, "preassigned_patient_id")
  expect_identical(observed$seed, NA_integer_)
  expect_identical(observed$n_folds, 3L)
  for (fold in seq_len(3L)) {
    expected_assessment <- which(expected_assignment == fold)
    expected_analysis <- which(expected_assignment != fold)
    expect_identical(observed$splits[[fold]]$fold, fold)
    expect_identical(
      observed$splits[[fold]]$assessment_rows,
      expected_assessment
    )
    expect_identical(observed$splits[[fold]]$analysis_rows, expected_analysis)
    expect_identical(
      observed$splits[[fold]]$assessment_ids,
      task_ids[expected_assessment]
    )
    expect_identical(
      observed$splits[[fold]]$analysis_ids,
      task_ids[expected_analysis]
    )
  }
})

test_that("ambiguous preassigned fold vectors fail before modeling", {
  d <- make_sb_binary_fixture(n = 12L, seed = 3721L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  ids <- as.character(task$data[[task$id_col]])
  valid <- setNames(rep(1:3, length.out = task$n), ids)
  duplicate_name <- valid
  names(duplicate_name)[[2L]] <- names(duplicate_name)[[1L]]
  blank_name <- valid
  names(blank_name)[[2L]] <- ""
  wrong_ids <- valid
  names(wrong_ids)[[2L]] <- "not_a_task_patient"
  missing_value <- valid
  missing_value[[2L]] <- NA_real_
  nonwhole_value <- valid
  nonwhole_value[[2L]] <- 1.5
  nonpositive_value <- valid
  nonpositive_value[[2L]] <- 0
  one_fold <- setNames(rep(1L, task$n), ids)
  nonconsecutive <- valid
  nonconsecutive[nonconsecutive == 3L] <- 4L
  character_value <- setNames(as.character(valid), ids)
  cases <- list(
    unnamed = list(value = unname(valid), pattern = "named"),
    duplicate = list(value = duplicate_name, pattern = "unique"),
    blank = list(value = blank_name, pattern = "non-empty"),
    wrong_ids = list(value = wrong_ids, pattern = "patient IDs"),
    missing = list(value = missing_value, pattern = "finite positive"),
    nonwhole = list(value = nonwhole_value, pattern = "finite positive"),
    nonpositive = list(value = nonpositive_value, pattern = "finite positive"),
    one_fold = list(value = one_fold, pattern = "at least two"),
    nonconsecutive = list(value = nonconsecutive, pattern = "1 through K"),
    character = list(value = character_value, pattern = "numeric")
  )

  for (case in cases) {
    error <- capture_sb_error(
      sb_make_outer_folds(task, folds = case$value, seed = 1L)
    )
    expect_s3_class(error, "stormbreaker_bad_fold_assignment")
    expect_match(conditionMessage(error), case$pattern)
  }
})

test_that("structurally valid outcome-segregated folds are not repaired", {
  d <- make_sb_binary_fixture(n = 20L, seed = 3731L)
  task <- sb_binary_task(
    d,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b")
  )
  ids <- as.character(task$data[[task$id_col]])
  supplied <- setNames(ifelse(d$outcome == 1L, 1L, 2L), ids)

  observed <- sb_make_outer_folds(task, folds = supplied, seed = 3732L)

  expect_identical(observed$assignment, as.integer(supplied))
  expect_identical(
    task$data[[task$outcome_col]][observed$splits[[1L]]$assessment_rows],
    rep(1L, task$n_positive)
  )
  expect_identical(
    task$data[[task$outcome_col]][observed$splits[[2L]]$assessment_rows],
    rep(0L, task$n_negative)
  )
})

test_that("preassigned fold construction does not consume the random stream", {
  d <- make_sb_fixture(n = 18L, seed = 3741L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month"
  )
  ids <- as.character(task$data[[task$id_col]])
  supplied <- setNames(rep(1:3, length.out = task$n), ids)
  set.seed(3742L)
  before <- .Random.seed

  first <- sb_make_outer_folds(task, folds = supplied, seed = 1L)
  after <- .Random.seed
  second <- sb_make_outer_folds(task, folds = supplied, seed = 999L)

  expect_identical(after, before)
  expect_identical(.Random.seed, before)
  expect_identical(first, second)
  expect_identical(first$seed, NA_integer_)
})
