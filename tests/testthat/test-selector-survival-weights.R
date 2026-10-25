weight_screen_data <- function() {
  withr::local_seed(92)
  d <- as.data.frame(matrix(rnorm(160 * 8), 160, 8))
  names(d) <- paste0("gene", 1:8)
  d$id <- paste0("p", 1:160)
  d$time <- rexp(160, exp(d$gene1 - d$gene2) / 10) + .01
  d$event <- rbinom(160, 1, .8)
  d
}

test_that("default survival weight screens select original genes and render", {
  d <- weight_screen_data()
  t <- sb_surv_task(d, "id", "time", "event", paste0("gene", 1:8), "month")
  for (method in c("Ridge", "PLSR", "SVM")) {
    p <- surv_pipeline("Cox", surv_selector(method))
    f <- sb_fit_final(t, p, capture_diagnostics = TRUE)
    expect_identical(f$status, "success")
    s <- f$selection$state
    expect_length(s$selected_features, 4)
    expect_equal(s$scores$score, abs(s$scores$coefficient * s$scores$training_sd))
    expect_setequal(s$selected_features, s$scores$feature[order(-s$scores$score)[1:4]])
    expect_identical(sb_predict_final(f, d[1:5, ], "id")$status, "success")
    z <- sb_plots(f)
    expect_length(z$unavailable, 0)
    expect_s3_class(z$plots$selector_1, "ggplot")
    for (plot in z$plots) expect_s3_class(ggplot2::ggplotGrob(plot +
      ggplot2::theme(text = ggplot2::element_text(family = "sans"))), "gtable")
    cv <- sb_run_cv(t, p, folds = 3, seed = 12)
    expect_identical(cv$status, "success")
  }
  f <- sb_fit_final(t, surv_pipeline("Cox", surv_selector(c("Ridge", "PLSR"))))
  expect_identical(f$status, "success")
  expect_length(f$selection$selected_features, 2)
})

test_that("Ridge screening keeps repeated patients in the same tuning fold", {
  d <- weight_screen_data()
  groups <- rep(seq_len(80), each = 2)
  p <- sb_survival_weight_selector_params("ridge_cox_weights")
  s <- sb_survival_weight_selector_adapter("ridge_cox_weights", p)$fit(
    d[paste0("gene", 1:8)], d$time, d$event, group_ids = groups)
  expect_equal(vapply(split(s$evidence$fold_id, groups), function(x) length(unique(x)), integer(1)),
    stats::setNames(rep(1L, 80), as.character(1:80)))
})
