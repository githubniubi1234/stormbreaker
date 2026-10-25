test_that("tree SHAP preserves prediction additivity and seeds", {
  d <- sb_example_binary(90, 42)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  for (model in c("RF", "GBM", "XGBoost")) {
    fit <- sb_fit_final(task, bi_pipeline(model))
    expect_identical(fit$status, "success")
    s <- sb_tree_shap(fit, d[1:2, ], "patient_id", d, permutations = 3, background_n = 4)
    expect_s3_class(s, "stormbreaker_tree_shap")
    expect_equal(unname(s$reference_value + rowSums(s$contributions)), s$scores$score)
    expect_equal(s, sb_tree_shap(fit, d[1:2, ], "patient_id", d, permutations = 3, background_n = 4))
    expect_identical(s$method, "approximate_permutation_interventional_shap")
    expect_s3_class(sb_plot_shap_importance(s), "ggplot")
    expect_s3_class(sb_plot_shap_beeswarm(s), "ggplot")
    expect_s3_class(sb_plot_shap_waterfall(s), "ggplot")
    expect_s3_class(sb_plot_shap_dependence(s, "gene_a", d[1:2, ], "patient_id"), "ggplot")
  }
})

test_that("one-feature SHAP equals prediction minus empirical reference", {
  d <- sb_example_binary(70, 12)
  t <- sb_binary_task(d, "patient_id", "outcome", "gene_a")
  f <- sb_fit_final(t, bi_pipeline("RF"))
  s <- sb_tree_shap(f, d[1:2, ], "patient_id", d[3:7, ], permutations = 1, background_n = 5)
  expected <- sb_predict_final(f, d[3:7, c("patient_id", "gene_a")], "patient_id")$predictions$probability
  expect_equal(s$reference_value, mean(expected))
  expect_equal(as.numeric(s$contributions), s$scores$score - mean(expected))
})

test_that("survival tree SHAP explains saved risk not event probability", {
  d <- sb_example_survival(70, 12)
  t <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month")
  f <- sb_fit_final(t, surv_pipeline("RF"))
  s <- sb_tree_shap(f, d[1:2, ], "patient_id", d, permutations = 2, background_n = 3)
  expect_identical(s$score_scale, "saved_risk_score")
  expect_equal(unname(rowSums(s$contributions) + s$reference_value), s$scores$score)
})
