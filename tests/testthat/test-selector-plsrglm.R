test_that("PLS logistic selection uses native coefficients in the analysis scale", {
  skip_if_not_installed("plsRglm")
  d <- sb_example_survival(180, 82); d$outcome <- d$event
  x <- d[c("gene_a", "gene_b", "gene_c")]
  names(x) <- c("gene A", "gene-B", "gene C")
  p <- sb_pipeline("plsrglm_binary_weights", "logistic")
  a <- sb_resolve_selector_adapter(p$selector, p$selector_params)
  set.seed(19); before <- .Random.seed
  selected <- a$fit(x, d$outcome)
  expect_identical(.Random.seed, before)
  native <- plsRglm::PLS_glm(dataY = d$outcome, dataX = x, nt = 2,
    modele = "pls-glm-logistic", typeVC = "none", scaleX = TRUE,
    pvals.expli = FALSE, MClassed = FALSE, sparse = FALSE, naive = FALSE, verbose = FALSE)
  score <- abs(drop(native$Coeffs)[-1] * vapply(x, sd, numeric(1)))
  expect_equal(selected$scores$score, unname(score), tolerance = 1e-12)
  expect_setequal(selected$selected_features, names(x)[order(-score)[1:2]])
  expect_null(selected$engine)
  expect_identical(a$apply(unserialize(serialize(selected, NULL)), x), x[selected$selected_features])
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  cv <- sb_run_cv(task, p, folds = 3, seed = 82)
  expect_identical(cv$status, "success")
  expect_identical(cv$provenance$plsrglm_version, as.character(packageVersion("plsRglm")))
  fit <- sb_fit_final(task, p)
  expect_identical(fit$status, "success")
  expect_identical(fit$provenance$plsrglm_version, as.character(packageVersion("plsRglm")))
  pred <- sb_predict_final(unserialize(serialize(fit, NULL)), d[c("patient_id", features)], "patient_id")
  expect_identical(pred$status, "success")
})
