test_that("stack comparison reports actual prediction variation for both tasks", {
  skip_if_not_installed("glmnet")
  for (type in c("binary", "survival")) {
    if (type == "binary") {
      d <- make_sb_binary_fixture(n = 60L, seed = 101L)
      task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
      learners <- c("logistic", "gaussian_naive_bayes")
    } else {
      d <- make_sb_fixture(n = 72L, seed = 103L)
      task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
      learners <- c("coxph", "survreg_weibull")
    }
    stack <- sb_stack(lapply(learners, function(id) sb_pipeline(learner = id)),
      inner_folds = 3L, inner_seed = 31L)
    result <- sb_run_stack_cv(task, stack, folds = 3L, seed = 37L)
    expect_identical(result$status, "success")
    expect_identical(tail(result$comparison$prediction_variation, 1L),
      sb_prediction_variation(result$predictions, type)$status)
    value <- if (type == "binary") "probability" else "risk"
    result$predictions[[value]][] <- 0.5
    expect_identical(sb_stack_comparison_row(result)$prediction_variation, "constant")
    result$predictions$status[1L] <- "failed"
    expect_identical(sb_stack_comparison_row(result)$prediction_variation, "incomplete")
  }
})
