doctor_chain_task <- function(grouped = FALSE) {
  withr::local_seed(81)
  d <- as.data.frame(matrix(rnorm(300 * 6), 300, 6))
  names(d) <- paste0("gene", 1:6)
  d$id <- paste0("p", 1:300)
  d$y <- rbinom(300, 1, plogis(2 * d$gene1 - 2 * d$gene2))
  if (grouped) {
    d$patient <- d$id
    d <- d[rep(seq_len(nrow(d)), each = 2), ]
    d$id <- paste0("r", seq_len(nrow(d)))
  }
  sb_binary_task(d, "id", "y", paste0("gene", 1:6), group_col = if (grouped) "patient" else NULL)
}

test_that("short menus cover all existing binary and survival learners", {
  for (task in c("binary", "survival")) {
    menu <- sb_models(task)
    predictors <- if (task == "binary") menu$model_id[menu$can_predict] else menu$model_id
    expect_setequal(predictors, sb_model_catalog(task)$model_id)
    expect_identical(anyDuplicated(tolower(menu$name)), 0L)
    plan <- sb_doctor_pipeline_set(task, if (task == "binary") c("RF", "GBM") else c("Cox", "CoxBoost"), c(1, 3))
    expect_length(plan$unavailable, 0L)
    expect_gte(length(plan$pipelines), 3L)
    for (p in plan$pipelines) expect_identical(sb_lint_pipeline(p)$status, "admitted")
  }
})

test_that("ordered binary chains keep stage evidence and frozen predictions", {
  task <- doctor_chain_task()
  p <- bi_pipeline("GBM", bi_selector(c("Lasso", "RF")))
  f <- sb_fit_final(task, p, capture_diagnostics = TRUE)
  expect_identical(f$status, "success")
  stages <- f$selection$state$stages
  expect_length(stages, 2L)
  expect_identical(stages[[2]]$source_features, stages[[1]]$selected_features)
  copied <- unserialize(serialize(f, NULL))
  new <- task$data[c("id", task$predictors)]
  expect_equal(sb_predict_final(copied, new, "id")$predictions,
    sb_predict_final(f, new, "id")$predictions)
  z <- sb_plots(f)
  expect_identical(all(c("selector_1", "selector_2") %in% names(z$plots)), TRUE)
  cv <- sb_run_cv(task, p, folds = 3, seed = 7)
  expect_identical(cv$status, "success")
  reverse <- bi_pipeline("GBM", bi_selector(c("RF", "Lasso")))
  expect_identical(p$id == reverse$id, FALSE)
})

test_that("chains refit within glmnet and SVM inner folds with grouped patients", {
  task <- doctor_chain_task(TRUE)
  original <- sb_fit_binary_model_selector
  sizes <- integer()
  local_mocked_bindings(sb_fit_binary_model_selector = function(x, outcome, id, params, group_ids = NULL) {
    result <- original(x, outcome, id, params, group_ids)
    sizes <<- c(sizes, nrow(x))
    if (id == "lasso_binary_nonzero") {
      expect_length(group_ids, nrow(x))
      expect_identical(all(vapply(split(result$evidence$fold_id, group_ids),
        function(v) length(unique(v)) == 1L, logical(1))), TRUE)
    }
    result
  })
  for (model in c("Lasso", "SVM")) {
    sizes <- integer()
    f <- sb_fit_final(task, bi_pipeline(model, bi_selector(c("Lasso", "RF"))))
    expect_identical(f$status, "success")
    expect_gt(length(sizes), 2L)
    expect_identical(any(sizes < task$n), TRUE)
  }
})

test_that("survival selector lists use saved stages in final and CV paths", {
  d <- sb_example_survival(120, 37)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
  p <- surv_pipeline("Cox", surv_selector(c("StepCox", "StepCox")))
  f <- sb_fit_final(task, p)
  expect_identical(f$status, "success")
  expect_length(f$selection$state$stages, 2)
  expect_identical(sb_run_cv(task, p, folds = 3)$status, "success")
})

test_that("doctor runs produce small-model figures without a browser", {
  d <- sb_example_binary(80, 11)
  r <- bi_run(d, "outcome", id = "patient_id", models = c("NearestCentroid", "BernoulliNB"), stages = 1, folds = 3, plots = TRUE)
  expect_length(r$details, 2)
  for (detail in r$details) {
    expect_identical(detail$fit$status, "success")
    expect_length(detail$figures$plots, 1)
    expect_s3_class(detail$figures$plots[[1]], "ggplot")
  }
  expect_length(sb_doctor_pipeline_set("binary", "RF", 2)$pipelines, 1L)
})

test_that("multiclass short entrypoint returns CV and algorithm figures", {
  d <- iris
  r <- multi_run(d, "Species", models = "RF", folds = 3)
  expect_identical(r$details$status$status, "success")
  expect_identical(r$comparison$comparison$status, "success")
  expect_gt(r$details$status$n_plots, 0)
})
