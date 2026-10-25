probability_distribution_fixture <- function() {
  structure(list(status = "success", task_type = "binary", pipeline_id = "distribution",
    scored_predictions = data.frame(patient_id = letters[1:6], outcome = c(1L, 0L, 1L, 0L, 1L, 0L),
      probability = c(0, .5, 1, 1, .5, 0), status = "success"), metrics = list()),
    class = "stormbreaker_external_evaluation")
}

test_that("probability bins preserve zero, one, ties and absent classes", {
  x <- probability_distribution_fixture()
  before <- x
  r <- sb_probability_distribution(x, c(0, .5, 1))
  expect_equal(r$table$count, c(1, 2, 1, 2))
  expect_equal(r$table$fraction, c(1, 2, 1, 2) / 3)
  expect_equal(sum(r$table$count), 6)
  expect_identical(x, before)
  expect_length(intersect("patient_id", names(r$table)), 0)
  expect_equal(sb_probability_distribution(x, c(0, 1))$table$count, c(3, 3))
  x$scored_predictions$outcome[] <- 1L
  q <- sb_probability_distribution(x, c(0, .5, 1))
  expect_equal(q$table$count, c(0, 0, 2, 4))
  expect_equal(q$table$fraction[1:2], c(NA_real_, NA_real_))
})

test_that("distributions preserve invalid groups and cohort separation", {
  x <- probability_distribution_fixture()
  bad <- x
  bad$scored_predictions$probability[1] <- NA_real_
  r <- sb_probability_distribution(list(A = x, B = bad), c(0, .5, 1))
  expect_equal(sum(r$table$count), 6)
  expect_equal(nrow(r$diagnostics), 2)
  expect_equal(sort(r$diagnostics$status), c("invalid_predictions", "success"))
  both <- sb_probability_distribution(list(A = x, B = x), c(0, .5, 1))
  expect_equal(length(unique(both$table$panel)), 2)
  expect_equal(sum(both$table$count), 12)
  for (b in list(c(0, 0, 1), c(.1, 1), c(0, NA, 1), c(0, .8, .5, 1))) {
    err <- tryCatch(sb_probability_distribution(x, b), error = identity)
    expect_s3_class(err, "stormbreaker_bad_visualization_argument")
  }
})

test_that("probability graphics use aggregate bins including empty bins", {
  skip_if_not_installed("ggplot2")
  r <- sb_probability_distribution(probability_distribution_fixture(), c(0, .2, .5, 1))
  p <- sb_plot_probability_distribution(r)
  expect_equal(nrow(p$data), 6)
  expect_equal(sum(p$data$value), 6)
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  expect_s3_class(ggplot2::ggplot_build(sb_plot_probability_distribution(r, "fraction")), "ggplot_built")
})
