selector_diagnostic_fixture <- function(kind) {
  withr::with_seed(881L, {
    n <- 120L
    d <- data.frame(patient_id = paste0("s", seq_len(n)),
      gene_a = rnorm(n), gene_b = rnorm(n), gene_c = rnorm(n),
      gene_d = rnorm(n), gene_e = rnorm(n), gene_f = rnorm(n))
    if (kind == "binary") {
      d$outcome <- rbinom(n, 1, plogis(2 * d$gene_a - 1.3 * d$gene_b))
      sb_binary_task(d, "patient_id", "outcome", paste0("gene_", letters[1:6]))
    } else {
      times <- rexp(n, rate = .08 * exp(1.8 * d$gene_a - 1.3 * d$gene_b)) + .1
      d$time <- pmin(times, median(times))
      d$event <- as.integer(times <= median(times))
      sb_surv_task(d, "patient_id", "time", "event", paste0("gene_", letters[1:6]), "month")
    }
  })
}

test_that("Boruta capture preserves selection and real predictor/shadow histories", {
  skip_if_not_installed("Boruta")
  skip_if_not_installed("ranger")
  for (kind in c("binary", "survival")) {
    task <- selector_diagnostic_fixture(kind)
    pipeline <- sb_pipeline(selector = paste0("boruta_", kind, "_roughfix"),
      learner = if (kind == "binary") "logistic" else "coxph")
    plain <- sb_fit_final(task, pipeline)
    captured <- sb_fit_final(task, pipeline, capture_diagnostics = TRUE)
    expect_identical(plain$status, "success")
    expect_identical(captured$status, "success")
    expect_identical(plain$learner_state, captured$learner_state)
    expect_identical(plain$selection$scores, captured$selection$scores)
    expect_null(plain$selection$state$diagnostic_history)
    q <- sb_selector_diagnostics(captured, "boruta_summary")
    expect_identical(q$table, captured$selection$scores)
    h <- sb_selector_diagnostics(captured, "boruta_history")
    native <- captured$selection$state$diagnostic_history$importance
    expect_identical(h$table$raw_importance, as.vector(native))
    expect_equal(sum(h$table$not_evaluated_after_rejection), sum(native == -Inf, na.rm = TRUE))
    expect_equal(nrow(h$table), length(native))
    expect_identical(sort(unique(h$table$kind)), c("Predictor", "Shadow reference"))
    err <- tryCatch(sb_selector_diagnostics(plain, "boruta_history"), error = identity)
    expect_s3_class(err, "stormbreaker_missing_diagnostic_artifacts")
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      for (type in c("boruta_summary", "boruta_hits", "boruta_history", "boruta_boxplot")) {
        expect_s3_class(ggplot2::ggplot_build(sb_plot_selector_diagnostic(
          sb_selector_diagnostics(captured, type))), "ggplot_built")
      }
    }
  }
})

test_that("varSelRF exposes saved OOB sequence without relabelling it CV", {
  skip_if_not_installed("varSelRF")
  fit <- sb_fit_final(selector_diagnostic_fixture("binary"),
    sb_pipeline(selector = "varselrf_binary_oob_1se", learner = "logistic"))
  expect_identical(fit$status, "success")
  q <- sb_selector_diagnostics(fit, "selection_history")
  expect_identical(q$table, fit$selection$state$selection_history)
  expect_identical(q$selected_features, fit$required_features)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(ggplot2::ggplot_build(sb_plot_selector_diagnostic(q)), "ggplot_built")
  }
})
