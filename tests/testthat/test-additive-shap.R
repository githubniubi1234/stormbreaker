test_that("logistic additive SHAP exactly reconstructs class-one log odds", {
  development <- make_sb_binary_fixture(n = 80L, seed = 9601L)
  predictors <- c("gene_c", "gene_a", "gene_b")
  task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    predictors
  )
  pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "logistic",
    selector_n_features = 2L
  )
  fit <- sb_fit_final(task, pipeline)

  new_data <- make_sb_binary_fixture(n = 8L, seed = 9602L)
  new_data$patient_id <- paste0("new_", new_data$patient_id)
  background <- development[1:30, c(
    "outcome", rev(fit$required_features)
  )]
  observed <- sb_additive_shap(
    fit,
    new_data,
    "patient_id",
    background
  )

  state <- fit$learner_state
  x <- as.matrix(new_data[fit$required_features])
  background_x <- as.matrix(background[fit$required_features])
  new_terms <- sweep(x, 2L, state$coefficient, "*")
  background_terms <- sweep(
    background_x, 2L, state$coefficient, "*"
  )
  reference_terms <- colMeans(background_terms)
  expected_contributions <- sweep(
    new_terms, 2L, reference_terms, "-"
  )
  expected_reference <- state$intercept + sum(reference_terms)
  expected_score <- expected_reference + rowSums(expected_contributions)
  direct_prediction <- sb_predict_final(
    fit,
    new_data,
    "patient_id"
  )

  expect_s3_class(observed, "stormbreaker_additive_shap")
  expect_identical(observed$status, "success")
  expect_identical(observed$task_type, "binary")
  expect_identical(observed$pipeline_id, pipeline$id)
  expect_identical(observed$learner, "logistic")
  expect_identical(
    observed$method,
    "exact_interventional_shap_for_additive_score"
  )
  expect_identical(observed$score_scale, "class_1_log_odds")
  expect_identical(observed$background$n, 30L)
  expect_identical(
    observed$background$feature_summary$feature,
    fit$required_features
  )
  expect_equal(
    observed$background$feature_summary$mean_feature_value,
    unname(colMeans(background_x)),
    tolerance = 0
  )
  expect_equal(
    observed$background$feature_summary$mean_score_term,
    unname(reference_terms),
    tolerance = 0
  )
  expect_equal(observed$reference_value, expected_reference, tolerance = 0)
  expect_equal(
    unname(observed$contributions),
    unname(expected_contributions),
    tolerance = 0
  )
  expect_identical(
    rownames(observed$contributions),
    as.character(new_data$patient_id)
  )
  expect_identical(colnames(observed$contributions), fit$required_features)
  expect_equal(observed$scores$score, expected_score, tolerance = 1e-15)
  expect_equal(
    stats::plogis(observed$scores$score),
    direct_prediction$predictions$probability,
    tolerance = 1e-15
  )
  expect_identical(observed$prediction, direct_prediction)
  expect_false(any(c("outcome", "time", "event") %in% names(observed)))

  changed_new <- new_data
  changed_new$outcome <- 1L - changed_new$outcome
  changed_background <- background
  changed_background$outcome <- 1L - changed_background$outcome
  expect_identical(
    sb_additive_shap(
      fit,
      changed_new,
      "patient_id",
      changed_background
    ),
    observed
  )
})

test_that("Cox additive SHAP stays on the saved risk-score scale", {
  development <- make_sb_fixture(n = 72L, seed = 9701L)
  task <- sb_surv_task(
    development,
    "patient_id",
    "time",
    "event",
    c("gene_b", "gene_a", "gene_c"),
    "month",
    prediction_horizon = 6
  )
  fit <- sb_fit_final(task, sb_pipeline())
  new_data <- make_sb_fixture(n = 7L, seed = 9702L)
  new_data$patient_id <- paste0("new_", new_data$patient_id)
  background <- development[1:24, rev(fit$required_features)]

  observed <- sb_additive_shap(
    fit,
    new_data,
    "patient_id",
    background
  )
  state <- fit$learner_state
  x <- as.matrix(new_data[fit$required_features])
  background_x <- as.matrix(background[fit$required_features])
  expected_reference <- sum(colMeans(background_x) * state$coefficient)
  expected_contributions <- sweep(
    sweep(x, 2L, state$coefficient, "*"),
    2L,
    colMeans(background_x) * state$coefficient,
    "-"
  )

  expect_identical(observed$task_type, "survival")
  expect_identical(
    observed$score_scale,
    "zero_reference_log_hazard_ratio"
  )
  expect_equal(observed$reference_value, expected_reference, tolerance = 0)
  expect_equal(
    unname(observed$contributions),
    unname(expected_contributions),
    tolerance = 1e-15
  )
  expect_equal(
    observed$scores$score,
    observed$prediction$predictions$risk,
    tolerance = 1e-15
  )
  expect_identical(
    "survival_probability" %in% names(observed$prediction$predictions),
    TRUE
  )
  printed <- paste(capture.output(print(observed)), collapse = "\n")
  expect_match(printed, "not probability-scale effects", fixed = TRUE)
  expect_match(printed, "No feature is ranked", fixed = TRUE)
})

test_that("Gaussian Naive Bayes SHAP reconstructs posterior log odds", {
  development <- make_sb_binary_fixture(n = 84L, seed = 9801L)
  task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  pipeline <- sb_pipeline(
    learner = "gaussian_naive_bayes",
    naive_bayes_var_smoothing = 1e-8
  )
  fit <- sb_fit_final(task, pipeline)
  new_data <- make_sb_binary_fixture(n = 9L, seed = 9802L)
  new_data$patient_id <- paste0("new_", new_data$patient_id)
  background <- development[1:28, fit$required_features]

  observed <- sb_additive_shap(
    fit,
    new_data,
    "patient_id",
    background
  )
  model <- structure(
    fit$learner_state$payload,
    class = "stormbreaker_gaussian_naive_bayes"
  )
  term_matrix <- function(data) {
    x <- as.matrix(data[model$features])
    class_term <- function(class) {
      centered <- sweep(x, 2L, model$means[class, ], "-")
      sweep(
        centered^2,
        2L,
        model$variances[class, ],
        "/"
      ) + matrix(
        log(2 * pi * model$variances[class, ]),
        nrow = nrow(x),
        ncol = ncol(x),
        byrow = TRUE
      )
    }
    -0.5 * class_term("1") + 0.5 * class_term("0")
  }
  new_terms <- term_matrix(new_data)
  background_terms <- term_matrix(background)
  reference_terms <- colMeans(background_terms)
  expected_reference <- log(model$priors[["1"]]) -
    log(model$priors[["0"]]) + sum(reference_terms)
  expected_contributions <- sweep(
    new_terms, 2L, reference_terms, "-"
  )

  expect_identical(
    observed$score_scale,
    "class_1_log_posterior_odds"
  )
  expect_equal(observed$reference_value, expected_reference, tolerance = 1e-15)
  expect_equal(
    unname(observed$contributions),
    unname(expected_contributions),
    tolerance = 1e-15
  )
  expect_equal(
    stats::plogis(observed$scores$score),
    observed$prediction$predictions$probability,
    tolerance = 1e-14
  )
})

test_that("glmnet logistic uses the same exact additive explanation", {
  skip_if_not_installed("glmnet")

  development <- make_sb_binary_fixture(n = 72L, seed = 9851L)
  task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 0.5,
    inner_folds = 3L,
    inner_seed = 9852L
  )
  fit <- sb_fit_final(task, pipeline)
  observed <- sb_additive_shap(
    fit,
    development[1:5, ],
    "patient_id",
    development[6:30, ]
  )

  expect_identical(observed$learner, "glmnet_logistic")
  expect_equal(
    stats::plogis(observed$scores$score),
    observed$prediction$predictions$probability,
    tolerance = 1e-14
  )
})

test_that("additive SHAP rejects invalid references and unsupported learners", {
  development <- make_sb_binary_fixture(n = 60L, seed = 9901L)
  task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  new_data <- development[1:4, ]
  missing_background <- development[5:20, c("gene_a", "gene_b")]
  duplicated_background <- development[5:20, fit$required_features]
  duplicated_background[[ncol(duplicated_background) + 1L]] <-
    duplicated_background[[fit$required_features[[1L]]]]
  names(duplicated_background)[[ncol(duplicated_background)]] <-
    fit$required_features[[1L]]

  expect_s3_class(
    capture_sb_error(sb_additive_shap(
      fit,
      new_data,
      "patient_id",
      missing_background
    )),
    "stormbreaker_missing_final_features"
  )
  expect_s3_class(
    capture_sb_error(sb_additive_shap(
      fit,
      new_data,
      "patient_id",
      duplicated_background
    )),
    "stormbreaker_duplicate_final_feature_columns"
  )
  expect_s3_class(
    capture_sb_error(sb_additive_shap(
      fit,
      new_data,
      "patient_id",
      development[FALSE, fit$required_features]
    )),
    "stormbreaker_bad_adapter_predictors"
  )

  if (requireNamespace("ranger", quietly = TRUE)) {
    ranger_fit <- sb_fit_final(
      task,
      sb_pipeline(
        learner = "ranger_binary",
        ranger_num_trees = 20L,
        ranger_seed = 9902L
      )
    )
    expect_s3_class(
      capture_sb_error(sb_additive_shap(
        ranger_fit,
        new_data,
        "patient_id",
        development[5:20, ]
      )),
      "stormbreaker_unsupported_additive_shap"
    )
  }
})
