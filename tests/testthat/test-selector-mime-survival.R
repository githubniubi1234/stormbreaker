make_mime_selector_fixture <- function(n = 90L, seed = 8601L) {
  set.seed(seed)
  signal_a <- stats::rnorm(n)
  signal_b <- stats::rnorm(n)
  signal_c <- stats::rnorm(n)
  noise_a <- stats::rnorm(n)
  noise_b <- stats::rnorm(n)
  noise_c <- stats::rnorm(n)
  linear_predictor <- 1.2 * signal_a - 0.9 * signal_b + 0.6 * signal_c
  event_time <- stats::rexp(n, rate = exp(linear_predictor) / 15)
  censor_time <- stats::rexp(n, rate = 1 / 20)
  data.frame(
    patient_id = sprintf("mime_selector_%03d", seq_len(n)),
    time = pmax(pmin(event_time, censor_time), 0.05),
    event = as.integer(event_time <= censor_time),
    signal_a = signal_a,
    signal_b = signal_b,
    signal_c = signal_c,
    noise_a = noise_a,
    noise_b = noise_b,
    noise_c = noise_c,
    stringsAsFactors = FALSE
  )
}

test_that("MIME selection failures report how many predictors were returned", {
  source_features <- c("gene_a", "gene_b", "gene_c")
  for (selected in list(character(), NA_character_, c("gene_a", "gene_a"))) {
    error <- capture_sb_error(
      sb_validate_mime_selector_features(selected, source_features, "Lasso")
    )
    expect_s3_class(error, "stormbreaker_mime_selector_too_few_features")
    expect_match(
      conditionMessage(error),
      paste0("returned ", length(selected), " of 3 candidate predictors"),
      fixed = TRUE
    )
  }
})

test_that("Lasso can retain one predictor for a single-input survival learner", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("randomForestSRC")
  d <- sb_example_survival(n = 200L, seed = 20260905L)
  noise <- withr::with_seed(20260906L, matrix(rnorm(200L * 17L), nrow = 200L))
  for (j in seq_len(ncol(noise))) d[[paste0("noise_", j)]] <- noise[, j]
  features <- setdiff(names(d), c("patient_id", "time", "event"))
  task <- sb_surv_task(d[1:160, ], "patient_id", "time", "event", features,
                       "month", prediction_horizon = 6)
  pipeline <- sb_pipeline(
    selector = "lasso_nonzero", learner = "rsf_survival",
    selector_inner_folds = 3L, selector_seed = 20260905L
  )
  fit <- sb_fit_final(task, pipeline)
  expect_identical(fit$status, "success")
  expect_identical(fit$selection$selected_features, "gene_a")
  prediction <- sb_predict_final(fit, d[161:200, c("patient_id", features)], "patient_id")
  expect_identical(prediction$status, "success")
  expect_equal(nrow(prediction$predictions), 40L)
  expect_equal(sum(is.finite(prediction$predictions$risk)), 40L)
  expect_identical(prediction$predictions$patient_id, d$patient_id[161:200])
})

test_that("RSF minimal-depth selection uses native top variables", {
  skip_if_not_installed("randomForestSRC")
  data <- make_mime_selector_fixture()
  x <- data[4:9]
  adapter <- sb_rsf_minimal_depth_selector_adapter(
    nodesize = 5L,
    seed = 8602L
  )

  fitted <- adapter$fit(x, data$time, data$event)
  selected <- adapter$apply(fitted, x[rev(seq_len(nrow(x))), ])

  expect_s3_class(fitted, "stormbreaker_selector_fit")
  expect_identical(fitted$id, "rsf_minimal_depth")
  expect_identical(fitted$source_features, names(x))
  expect_gt(length(fitted$selected_features), 1L)
  expect_true(all(fitted$selected_features %in% names(x)))
  expect_identical(names(selected), fitted$selected_features)
  expect_identical(
    selected,
    x[rev(seq_len(nrow(x))), fitted$selected_features, drop = FALSE]
  )
})

test_that("CoxBoost nonzero selection records tuned native coefficients", {
  skip_if_not_installed("CoxBoost")
  data <- make_mime_selector_fixture(seed = 8611L)
  x <- data[4:9]
  adapter <- sb_coxboost_nonzero_selector_adapter(
    inner_folds = 3L,
    seed = 8612L
  )

  fitted <- adapter$fit(x, data$time, data$event)

  expect_s3_class(fitted, "stormbreaker_selector_fit")
  expect_identical(fitted$id, "coxboost_nonzero")
  expect_gt(length(fitted$selected_features), 1L)
  expect_true(all(fitted$selected_features %in% names(x)))
  expect_true(is.finite(fitted$penalty))
  expect_true(is.finite(fitted$optimal_step))
  expect_identical(
    fitted$scores$selected,
    fitted$scores$coefficient != 0
  )
})

test_that("Lasso nonzero selection records fold-local lambda", {
  skip_if_not_installed("glmnet")
  data <- make_mime_selector_fixture(seed = 8621L)
  x <- data[4:9]
  adapter <- sb_lasso_nonzero_selector_adapter(
    inner_folds = 3L,
    seed = 8622L
  )

  fitted <- adapter$fit(x, data$time, data$event)

  expect_s3_class(fitted, "stormbreaker_selector_fit")
  expect_identical(fitted$id, "lasso_nonzero")
  expect_gt(length(fitted$selected_features), 1L)
  expect_true(all(fitted$selected_features %in% names(x)))
  expect_true(is.finite(fitted$lambda))
  expect_identical(fitted$alpha, 1)
  expect_identical(fitted$selection_rule, "lambda.min")
  expect_identical(
    fitted$scores$selected,
    fitted$scores$coefficient != 0
  )
})

test_that("sequential selection fits StepCox only on stage-one features", {
  skip_if_not_installed("randomForestSRC")
  data <- make_mime_selector_fixture(n = 120L, seed = 8631L)
  x <- data[4:9]
  adapter <- sb_rsf_then_stepcox_selector_adapter(
    nodesize = 5L,
    seed = 8632L,
    direction = "forward"
  )

  fitted <- adapter$fit(x, data$time, data$event)
  selected <- adapter$apply(fitted, x)

  expect_s3_class(fitted, "stormbreaker_selector_fit")
  expect_identical(fitted$id, "rsf_then_stepcox")
  expect_identical(
    fitted$stages$second$source_features,
    fitted$stages$first$selected_features
  )
  expect_true(all(
    fitted$selected_features %in% fitted$stages$first$selected_features
  ))
  expect_identical(names(selected), fitted$selected_features)
  expect_identical(fitted$stages$second$direction, "forward")
})

test_that("all MIME first-stage selectors compose with StepCox", {
  adapters <- list(
    rsf_then_stepcox = sb_rsf_then_stepcox_selector_adapter(
      nodesize = 5L,
      seed = 8641L,
      direction = "both"
    ),
    coxboost_then_stepcox = sb_coxboost_then_stepcox_selector_adapter(
      inner_folds = 3L,
      seed = 8642L,
      direction = "backward"
    ),
    lasso_then_stepcox = sb_lasso_then_stepcox_selector_adapter(
      inner_folds = 3L,
      seed = 8643L,
      direction = "forward"
    )
  )

  expect_identical(
    unname(vapply(adapters, `[[`, character(1), "id")),
    names(adapters)
  )
  expect_true(all(vapply(adapters, `[[`, logical(1), "needs_outcome")))
  expect_true(all(vapply(adapters, function(adapter) {
    is.function(adapter$fit) && is.function(adapter$apply)
  }, logical(1))))
})

test_that("MIME survival selector controls reconstruct registered pipelines", {
  pipelines <- list(
    sb_pipeline(
      selector = "rsf_minimal_depth",
      learner = "coxboost_cv",
      selector_rsf_nodesize = 5L,
      selector_seed = 8651L,
      learner_tuning_folds = 3L,
      learner_seed = 8652L
    ),
    sb_pipeline(
      selector = "coxboost_nonzero",
      learner = "glmnet_cox",
      selector_inner_folds = 3L,
      selector_seed = 8652L,
      inner_folds = 3L,
      inner_seed = 8653L,
      glmnet_alpha = 0.5
    ),
    sb_pipeline(
      selector = "lasso_then_stepcox",
      learner = "coxph",
      selector_inner_folds = 3L,
      selector_seed = 8654L,
      stepcox_direction = "both"
    )
  )

  expect_length(unique(vapply(pipelines, `[[`, character(1), "id")), 3L)
  for (pipeline in pipelines) {
    runtime <- sb_resolve_pipeline_runtime(pipeline)
    expect_identical(runtime$selector$params, pipeline$selector_params)
    expect_identical(
      sb_validate_pipeline_contract(pipeline, "survival"),
      pipeline
    )
  }

  missing_seed <- capture_sb_error(sb_pipeline(
    selector = "lasso_nonzero",
    learner = "coxboost_cv",
    selector_inner_folds = 3L
  ))
  irrelevant_control <- capture_sb_error(sb_pipeline(
    learner = "coxph",
    selector_seed = 8655L
  ))
  admitted_shape <- sb_pipeline(
    selector = "rsf_minimal_depth",
    learner = "bart_survival",
    selector_rsf_nodesize = 5L,
    selector_seed = 8656L
  )

  expect_s3_class(missing_seed, "stormbreaker_missing_selector_parameter")
  expect_s3_class(
    irrelevant_control,
    "stormbreaker_irrelevant_selector_parameter"
  )
  expect_s3_class(admitted_shape, "stormbreaker_pipeline")
  expect_identical(admitted_shape$selector, "rsf_minimal_depth")
  expect_identical(admitted_shape$learner, "bart_survival")
})
