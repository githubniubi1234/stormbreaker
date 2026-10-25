test_that("LogitBoost and RDA SHAP use original probability predictions", {
  render <- function(p) ggplot2::ggplotGrob(p + ggplot2::theme(text = ggplot2::element_text(family = "sans")))
  d <- sb_example_binary(100, 23)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  for (id in c("logitboost_binary", "rda_binary")) {
    fit <- sb_fit_final(task, sb_pipeline(learner = id))
    expect_identical(fit$status, "success")
    s <- sb_model_shap(fit, d[1:3, ], "patient_id", d, permutations = 3, background_n = 4)
    expect_s3_class(s, "stormbreaker_model_shap")
    expect_identical(s$score_scale, "class_1_probability")
    expect_equal(unname(s$reference_value + rowSums(s$contributions)), s$scores$score)
    expect_equal(s$scores$score, sb_predict_final(fit, d[1:3, ], "patient_id")$predictions$probability)
    expect_equal(s, sb_model_shap(fit, d[1:3, ], "patient_id", d, permutations = 3, background_n = 4))
    expect_s3_class(render(sb_plot_shap_beeswarm(s, feature_values = d[1:3, ], id_col = "patient_id")), "gtable")
    expect_s3_class(render(sb_plot_shap_importance(s, groupby = c("A", "B", "A"))), "gtable")
    expect_s3_class(render(sb_plot_shap_waterfall(s)), "gtable")
  }
})
