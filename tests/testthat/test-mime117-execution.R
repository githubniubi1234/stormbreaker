make_mime117_execution_fixture <- function(n = 150L, p = 15L,
                                           seed = 8801L) {
  set.seed(seed)
  x <- replicate(p, stats::rnorm(n))
  colnames(x) <- sprintf("gene_%02d", seq_len(p))
  linear_predictor <- 1.4 * x[, 1L] - 1.1 * x[, 2L] +
    0.8 * x[, 3L] - 0.5 * x[, 4L]
  event_time <- stats::rexp(n, rate = exp(linear_predictor) / 18)
  censor_time <- stats::rexp(n, rate = 1 / 26)
  data.frame(
    patient_id = sprintf("mime117_%03d", seq_len(n)),
    time = pmax(pmin(event_time, censor_time), 0.05),
    event = as.integer(event_time <= censor_time),
    x,
    check.names = FALSE
  )
}

test_that("representative MIME117 pipelines really fit and predict", {
  skip_if_not_installed("randomForestSRC")
  skip_if_not_installed("CoxBoost")
  skip_if_not_installed("glmnet")

  data <- make_mime117_execution_fixture()
  predictors <- grep("^gene_", names(data), value = TRUE)
  task <- sb_surv_task(
    data,
    "patient_id",
    "time",
    "event",
    predictors,
    "month"
  )
  suite <- sb_mime117_pipelines(
    inner_folds = 3L,
    inner_seed = 8802L,
    selector_inner_folds = 3L,
    selector_seed = 8803L,
    rsf_nodesize = 5L,
    learner_tuning_folds = 3L,
    plsrcox_tuning_folds = 3L,
    learner_seed = 8804L,
    gbm_cv_max_trees = 50L
  )
  labels <- c(
    "RSF + CoxBoost",
    "CoxBoost + Enet[\u03b1=0.5]",
    "Lasso + StepCox[forward]",
    "StepCox[backward] + RSF"
  )

  fits <- lapply(labels, function(label) {
    index <- match(label, suite$official_model_label)
    expect_false(is.na(index), info = label)
    fit <- sb_fit_final(task, suite$pipeline[[index]])
    expect_identical(fit$status, "success", info = label)
    prediction <- sb_predict_final(
      fit,
      data[1:5, , drop = FALSE],
      "patient_id"
    )
    expect_identical(prediction$status, "success", info = label)
    expect_true(
      all(is.finite(prediction$predictions$risk)),
      info = label
    )
    fit
  })

  nested <- fits[[match("CoxBoost + Enet[\u03b1=0.5]", labels)]]
  expect_identical(
    nested$tuning$inner_selection_method,
    "glmnet_exclude_function"
  )
  expect_length(nested$tuning$inner_selections, 3L)
  expect_true(all(vapply(
    nested$tuning$inner_selections,
    function(selection) {
      length(selection$analysis_patient_ids) < nrow(data) &&
        length(selection$selected_features) > 1L
    },
    logical(1)
  )))
})
