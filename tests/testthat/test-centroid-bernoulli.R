test_that("centroid scores match Euclidean geometry including ties", {
  x <- data.frame(a = c(0, 2, 4, 6), constant = 1)
  y <- c(0L, 0L, 1L, 1L)
  id <- "nearest_centroid_binary"
  m <- sb_fit_classical_binary(x, y, id, sb_classical_binary_params(id))
  expect_equal(unname(m$engine$centroids), matrix(c(1, 5, 1, 1), 2))
  q <- data.frame(a = c(1, 3, 5), constant = 1)
  p <- sb_predict_classical_binary(m, q)
  expect_equal(p, plogis(c(-8, 0, 8)))
  saved <- sb_extract_classical_binary_learner_state(m)
  expect_equal(sb_predict_classical_binary_learner_state(saved, q)$probability, p)
})

test_that("Bernoulli NB includes absent features, priors and Laplace smoothing", {
  x <- data.frame(a = c(-2, 0, 1, 2), b = c(0, 1, 0, 1), constant = 0)
  id <- "bernoulli_naive_bayes_binary"
  m <- sb_fit_classical_binary(x, c(0, 0, 0, 1), id,
    sb_classical_binary_params(id))
  expect_equal(m$engine$prior, c(.75, .25))
  expect_equal(unname(m$engine$feature_probability), rbind(c(2/5, 2/5, 1/5), c(2/3, 2/3, 1/3)))
  q <- data.frame(a = 0, b = 1, constant = 0)
  likelihood <- c(.75 * 3/5 * 2/5 * 4/5, .25 * 1/3 * 2/3 * 2/3)
  expect_equal(sb_predict_classical_binary(m, q), likelihood[2]/sum(likelihood))
})

test_that("new learners support CV, selection, saved prediction and external AUC", {
  d <- sb_example_binary(100, 17)
  genes <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", genes)
  external <- sb_example_binary(50, 18)
  for (id in c("nearest_centroid_binary", "bernoulli_naive_bayes_binary")) {
    for (selector in c("none", "univ_binary_auc_top_k", "lasso_binary_nonzero")) {
      p <- sb_pipeline(selector, id,
        selector_n_features = if (selector == "univ_binary_auc_top_k") 2L else NULL,
        selector_options = if (selector == "lasso_binary_nonzero") list(lambda_rule = "lambda.min") else NULL)
      expect_identical(sb_lint_pipeline(p)$status, "admitted")
      fit <- sb_fit_final(task, p)
      expect_identical(fit$status, "success")
      pred <- sb_predict_final(unserialize(serialize(fit, NULL)), external[c("patient_id", genes)], "patient_id")
      ev <- sb_evaluate_external(pred, sb_binary_endpoint(external, "patient_id", "outcome"))
      expect_identical(ev$status, "success")
      if (selector != "lasso_binary_nonzero") {
        cv <- sb_run_cv(task, p, folds = 3, seed = 17)
        expect_identical(cv$status, "success")
        expect_s3_class(sb_plot_roc(cv), "ggplot")
      }
    }
  }
})

test_that("CV learns class statistics only from analysis rows", {
  d <- sb_example_binary(90, 22)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  original <- sb_fit_classical_binary
  seen <- list()
  local_mocked_bindings(sb_fit_classical_binary = function(x, outcome, learner, params) {
    seen[[length(seen) + 1L]] <<- rownames(x)
    original(x, outcome, learner, params)
  })
  for (id in c("nearest_centroid_binary", "bernoulli_naive_bayes_binary")) {
    seen <- list()
    cv <- sb_run_cv(task, sb_pipeline(learner = id), folds = 3, seed = 22)
    expect_identical(cv$status, "success")
    expect_length(seen, 3L)
    for (i in seq_along(seen)) {
      expect_setequal(seen[[i]], rownames(d)[cv$folds[[i]]$analysis_rows])
    }
  }
})
