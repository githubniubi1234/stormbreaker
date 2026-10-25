test_that("the synthetic survival example is valid and reproducible", {
  first <- sb_example_survival(n = 50, seed = 42)
  second <- sb_example_survival(n = 50, seed = 42)

  expect_identical(first, second)
  expect_identical(
    names(first),
    c("patient_id", "time", "event", "gene_a", "gene_b", "gene_c")
  )
  expect_identical(nrow(first), 50L)
  expect_identical(sum(first$event == 1L), 25L)
  expect_identical(sum(first$event == 0L), 25L)
  expect_identical(anyNA(first), FALSE)
  expect_identical(all(first$time > 0), TRUE)

  task <- sb_surv_task(
    first, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  expect_s3_class(task, "stormbreaker_surv_task")
})

test_that("synthetic data arguments are checked", {
  n_error <- capture_sb_error(sb_example_survival(n = 9, seed = 1))
  expect_s3_class(n_error, "stormbreaker_bad_integer")

  seed_error <- capture_sb_error(sb_example_survival(n = 20, seed = 0))
  expect_s3_class(seed_error, "stormbreaker_bad_integer")
})

test_that("the synthetic binary example is valid and reproducible", {
  first <- sb_example_binary(n = 50, seed = 42)
  second <- sb_example_binary(n = 50, seed = 42)

  expect_identical(first, second)
  expect_identical(
    names(first),
    c("patient_id", "outcome", "gene_a", "gene_b", "gene_c")
  )
  expect_identical(nrow(first), 50L)
  expect_identical(sum(first$outcome == 1L), 25L)
  expect_identical(sum(first$outcome == 0L), 25L)
  expect_identical(anyNA(first), FALSE)

  task <- sb_binary_task(
    first, "patient_id", "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  expect_s3_class(task, "stormbreaker_binary_task")
})

test_that("synthetic binary data arguments are checked", {
  n_error <- capture_sb_error(sb_example_binary(n = 9, seed = 1))
  expect_s3_class(n_error, "stormbreaker_bad_integer")

  seed_error <- capture_sb_error(sb_example_binary(n = 20, seed = 0))
  expect_s3_class(seed_error, "stormbreaker_bad_integer")
})
