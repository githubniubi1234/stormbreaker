test_that("GBM CV reuse preserves the original two-call native fit", {
  skip_if_not_installed("gbm")
  d <- sb_example_survival(90, seed = 61)
  x <- d[c("gene_a", "gene_b")]
  for (steps in c(1L, 2L, 100L)) {
    p <- sb_gbm_cox_cv_params(inner_folds = 3L, maximum_trees = steps)
    p$shrinkage <- 0.1
    data <- sb_additional_survival_formula_data(x, d$time, d$event)
    native <- function(nt, folds = p$inner_folds) {
      withr::with_seed(p$seed, gbm::gbm(data$formula, data = data$data,
        distribution = p$distribution, n.trees = nt,
        interaction.depth = p$interaction_depth, n.minobsinnode = p$min_obs_node,
        shrinkage = p$shrinkage, bag.fraction = p$bag_fraction,
        cv.folds = folds, n.cores = 1L, keep.data = FALSE, verbose = FALSE))
    }
    withr::local_seed(91)
    before <- .Random.seed
    invisible(capture.output({
      search <- if (steps == 1L) NULL else native(steps)
      best <- if (steps == 1L) 1L else which.min(search$cv.error)
      expected <- native(best, if (best == 1L) 0L else p$inner_folds)
      actual <- sb_fit_gbm_cox_cv(x, d$time, d$event, p)
    }))
    expect_identical(.Random.seed, before)
    expect_identical(actual$n_trees, best)
    expect_identical(actual$tuning$cv_error, as.numeric(search$cv.error))
    expect_identical(actual$booster$trees, expected$trees)
    expect_identical(actual$booster$cv.error, expected$cv.error)
    expect_identical(actual$booster$cv.fitted, expected$cv.fitted)
    expect_identical(actual$booster$train.error, expected$train.error)
    expect_identical(actual$booster$oobag.improve, expected$oobag.improve)
    expect_identical(sb_predict_gbm_cox(actual, x),
      as.numeric(predict(expected, data$data, n.trees = best, type = "link")))
    restored <- unserialize(serialize(actual, NULL))
    expect_identical(sb_predict_gbm_cox(restored, x),
      sb_predict_gbm_cox(actual, x))
  }
})

test_that("warned GBM searches are not reused", {
  skip_if_not_installed("gbm")
  calls <- 0L
  local_mocked_bindings(gbm = function(...) {
    calls <<- calls + 1L
    warning("native search warning")
    structure(list(n.trees = 2L, cv.error = c(2, 1),
      call = quote(gbm()), m = quote(gbm())), class = "gbm")
  }, .package = "gbm")
  d <- sb_example_survival(60, seed = 1)
  seen <- character()
  withCallingHandlers(sb_fit_gbm_cox_cv(d[c("gene_a", "gene_b")], d$time,
    d$event, sb_gbm_cox_cv_params(inner_folds = 3L, maximum_trees = 2L)),
    warning = function(w) { seen <<- c(seen, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_identical(calls, 2L)
  expect_identical(seen, rep("native search warning", 2L))
})

test_that("an unremarkable upper-bound GBM optimum needs only one native call", {
  skip_if_not_installed("gbm")
  calls <- 0L
  local_mocked_bindings(gbm = function(...) {
    calls <<- calls + 1L
    structure(list(n.trees = 2L, cv.error = c(2, 1),
      call = quote(gbm()), m = quote(gbm())), class = "gbm")
  }, .package = "gbm")
  d <- sb_example_survival(60, seed = 1)
  fit <- sb_fit_gbm_cox_cv(d[c("gene_a", "gene_b")], d$time, d$event,
    sb_gbm_cox_cv_params(inner_folds = 3L, maximum_trees = 2L))
  expect_identical(calls, 1L)
  expect_identical(fit$n_trees, 2L)
})
